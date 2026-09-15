/**
 * CFD Solver Application
 * ImGui + OpenGL UI – responsive fullscreen-friendly layout
 * Adds material presets, per-side thermal BCs, and auto simulation modes.
 */
#include "renderer/Window.hpp"
#include "renderer/MeshEditor2D.hpp"
#include "solvers/HeatSolver2D.h"
#include "linearAlgebra/ConjugateGradient.hpp"
#include "mesh/mesh2D.h"
#include "IO/MeshReader.hpp"

#include <imgui.h>
#include <iostream>
#include <string>
#include <vector>
#include <thread>
#include <mutex>
#include <atomic>
#include <deque>
#include <algorithm>
#include <cmath>

#include <nfd.h>

struct PendingFrame {
    std::vector<double> values;
    double time = 0.0;
    double globalMin = 0.0;
    double globalMax = 1.0;
};

static std::mutex               g_frameMutex;
static std::deque<PendingFrame> g_frameQueue;
static std::atomic<bool>        g_solverRunning{ false };
static std::atomic<float>       g_solverProgress{ 0.0f };
static std::atomic<bool>        g_solverError{ false };
static std::string              g_solverErrorMsg;
static std::atomic<float>       g_suggestedAnimSpeed{ 0.0f };

enum class SimMode {
    MaterialComparison = 0,
    Transient = 1,
    NearSteady = 2,
    Manual = 3
};

struct MaterialPreset {
    const char* name;
    float alpha;
};

static void solverThread(
    int    meshNx, int meshNy,
    double width, double height,
    double alpha,
    double T_west,
    double T_east,
    double T_south,
    double T_north,
    double tEnd,
    int    targetFrames)
{
    g_solverError = false;
    g_solverProgress = 0.0f;
    g_suggestedAnimSpeed = 0.0f;

    try {
        Mesh mesh(meshNx, meshNy, width, height);

        const double dx = width / static_cast<double>(meshNx);
        const double dy = height / static_cast<double>(meshNy);
        const double dMin = std::min(dx, dy);
        const double Fo_target = 5.0;
        const double dt = Fo_target * dMin * dMin / alpha;

        CFD::HeatSolver2D solver(mesh, alpha, dt, tEnd);

        CFD::HeatSolver2D::BoundaryConditions bcs;
        bcs.T_west = T_west;
        bcs.T_east = T_east;
        bcs.T_south = T_south;
        bcs.T_north = T_north;

        const double icMean = 0.25 * (T_west + T_east + T_south + T_north);

        solver.setBCs(bcs);
        solver.setIC(icMean);
        solver.setOutputFreq(1);
        solver.enableSteadyStop(true);
        solver.setSteadyTolerance(1.0e-6);

        const int safeFrameCount = std::max(2, targetFrames);
        const double frameDt = tEnd / static_cast<double>(safeFrameCount - 1);
        double nextFrameTime = 0.0;

        CFD::ConjugateGradient cg(1e-8, 5000);
        std::vector<int> cgIters;

        struct RawFrame { std::vector<double> values; double time = 0.0; };
        std::vector<RawFrame> rawFrames;
        rawFrames.reserve(static_cast<std::size_t>(safeFrameCount + 4));

        double globalMin = 1e30;
        double globalMax = -1e30;

        auto recordFrame = [&](double timeValue) {
            const CFD::ScalarField& T = solver.T();
            RawFrame rf;
            rf.values.resize(T.size());
            for (std::size_t k = 0; k < T.size(); ++k) {
                rf.values[k] = T[k];
                globalMin = std::min(globalMin, T[k]);
                globalMax = std::max(globalMax, T[k]);
            }
            rf.time = timeValue;
            rawFrames.push_back(std::move(rf));
            };

        recordFrame(0.0);
        nextFrameTime += frameDt;

        while (!solver.finished() && g_solverRunning) {
            solver.step(cg, cgIters);

            while (solver.time() + 1.0e-12 >= nextFrameTime &&
                static_cast<int>(rawFrames.size()) < safeFrameCount)
            {
                recordFrame(solver.time());
                nextFrameTime += frameDt;
            }

            g_solverProgress = static_cast<float>(
                std::min(solver.time() / std::max(tEnd, 1.0e-12), 1.0));
        }

        if (rawFrames.empty() || rawFrames.back().time < solver.time()) {
            recordFrame(solver.time());
        }

        if (globalMax <= globalMin) globalMax = globalMin + 1.0;

        const int   nFrames = static_cast<int>(rawFrames.size());
        const float desiredPlaySecs = 8.0f;
        const float autoSpeed = static_cast<float>(nFrames) / desiredPlaySecs;
        g_suggestedAnimSpeed = std::clamp(autoSpeed, 1.0f, 60.0f);

        {
            std::lock_guard<std::mutex> lock(g_frameMutex);
            for (auto& rf : rawFrames) {
                PendingFrame pf;
                pf.values = std::move(rf.values);
                pf.time = rf.time;
                pf.globalMin = globalMin;
                pf.globalMax = globalMax;
                g_frameQueue.push_back(std::move(pf));
            }
        }

        g_solverProgress = 1.0f;
    }
    catch (const std::exception& e) {
        g_solverErrorMsg = e.what();
        g_solverError = true;
    }

    g_solverRunning = false;
}

int main()
{
    CFD::Window window(1280, 720, "CFD Solver");
    if (!window.init()) {
        std::cerr << "Failed to init!\n";
        return -1;
    }

    if (NFD_Init() != NFD_OKAY) {
        std::cerr << "NFD init failed: " << NFD_GetError() << '\n';
        window.cleanup();
        return -1;
    }

    static const MaterialPreset materials[] = {
        { "Custom",     1.0e-4f },
        { "Air",        2.1e-5f },
        { "Water",      1.4e-7f },
        { "Steel",      1.2e-5f },
        { "Aluminium",  8.4e-5f },
        { "Copper",     1.11e-4f }
    };

    int materialIndex = 0;
    float alpha = materials[materialIndex].alpha;

    float T_west = 1.0f;
    float T_east = 0.0f;
    float T_south = 0.0f;
    float T_north = 0.0f;

    int simModeIndex = 0; // Material Comparison
    float manualEndTime = 10.0f;
    int manualFrames = 50;
    bool loopAnimation = false;

    CFD::OBJMesh loadedMesh;
    CFD::UI::MeshEditor2D meshEditor;
    meshEditor.showEditorWindow(false);
    meshEditor.showViewportWindow(true);

    std::cout << "App running!\n";

    while (!window.shouldClose())
    {
        window.beginFrame();

        {
            std::lock_guard<std::mutex> lock(g_frameMutex);
            while (!g_frameQueue.empty()) {
                auto& pf = g_frameQueue.front();
                meshEditor.pushAnimationFrame(
                    pf.values, pf.time, pf.globalMin, pf.globalMax);
                g_frameQueue.pop_front();
            }
        }

        {
            const float spd = g_suggestedAnimSpeed.exchange(0.0f);
            if (spd > 0.0f)
                meshEditor.setAnimSpeed(spd);
        }

        if (ImGui::BeginMainMenuBar()) {
            if (ImGui::BeginMenu("File")) {
                ImGui::MenuItem("New");
                ImGui::Separator();
                if (ImGui::MenuItem("Exit")) break;
                ImGui::EndMenu();
            }

            if (ImGui::BeginMenu("View")) {
                bool showEditor = meshEditor.isEditorWindowVisible();
                bool showViewport = meshEditor.isViewportWindowVisible();
                if (ImGui::MenuItem("Mesh Editor", nullptr, showEditor))
                    meshEditor.showEditorWindow(!showEditor);
                if (ImGui::MenuItem("Mesh Viewport", nullptr, showViewport))
                    meshEditor.showViewportWindow(!showViewport);
                ImGui::EndMenu();
            }

            if (ImGui::BeginMenu("Mesh")) {
                if (ImGui::MenuItem("Import OBJ...")) {
                    nfdchar_t* outPath = nullptr;
                    nfdfilteritem_t filter = { "OBJ Files", "obj" };
                    nfdresult_t result = NFD_OpenDialog(&outPath, &filter, 1, nullptr);

                    if (result == NFD_OKAY) {
                        try {
                            loadedMesh = CFD::loadOBJ(outPath);
                            meshEditor.setMesh(&loadedMesh);
                            meshEditor.showViewportWindow(true);
                        }
                        catch (const std::exception& e) {
                            std::cerr << "Failed to load OBJ: " << e.what() << '\n';
                        }
                        NFD_FreePath(outPath);
                    }
                }

                if (ImGui::MenuItem("Clear Mesh", nullptr, false, meshEditor.hasMesh())) {
                    meshEditor.setMesh(nullptr);
                }

                ImGui::EndMenu();
            }

            if (ImGui::BeginMenu("Help")) {
                if (ImGui::MenuItem("About")) window.openAbout();
                ImGui::EndMenu();
            }

            ImGui::EndMainMenuBar();
        }
        window.renderAbout();

        const float W = static_cast<float>(window.width());
        const float H = static_cast<float>(window.height());
        const float menuBarH = ImGui::GetFrameHeight();
        const float leftToolbarW = 55.0f;
        const float rightPanelW = 330.0f;
        const float bottomTimelineH = 60.0f;
        const float topY = menuBarH;
        const float mainH = H - topY;
        const float centerX = leftToolbarW;
        const float centerW = std::max(200.0f, W - leftToolbarW - rightPanelW);
        const float centerH = std::max(200.0f, mainH - bottomTimelineH);

        meshEditor.drawUI();
        meshEditor.drawViewport();

        ImGui::SetNextWindowPos({ 0.0f, topY }, ImGuiCond_Always);
        ImGui::SetNextWindowSize({ leftToolbarW, mainH }, ImGuiCond_Always);
        ImGui::Begin("##toolbar", nullptr,
            ImGuiWindowFlags_NoTitleBar |
            ImGuiWindowFlags_NoResize |
            ImGuiWindowFlags_NoMove |
            ImGuiWindowFlags_NoScrollbar |
            ImGuiWindowFlags_NoCollapse);

        ImGui::PushStyleColor(ImGuiCol_Button, { 0.3f, 0.5f, 0.9f, 1.0f });
        if (ImGui::Button("GEO\n   ", { 40, 50 })) {
            meshEditor.showEditorWindow(true);
            meshEditor.showViewportWindow(true);
        }
        ImGui::PopStyleColor();
        ImGui::SetItemTooltip("Geometry Mode");

        ImGui::Button("PHY\n   ", { 40, 50 }); ImGui::SetItemTooltip("Physics Setup");
        ImGui::Button("MSH\n   ", { 40, 50 }); ImGui::SetItemTooltip("Mesh");

        ImGui::PushStyleColor(ImGuiCol_Button, { 0.2f, 0.7f, 0.2f, 1.0f });
        ImGui::Button("RUN\n   ", { 40, 50 });
        ImGui::PopStyleColor();
        ImGui::SetItemTooltip("Run Solver");

        ImGui::Button("RES\n   ", { 40, 50 }); ImGui::SetItemTooltip("Results");
        ImGui::End();

        ImGui::SetNextWindowPos({ W - rightPanelW, topY }, ImGuiCond_Always);
        ImGui::SetNextWindowSize({ rightPanelW, mainH }, ImGuiCond_Always);
        ImGui::Begin("Properties", nullptr,
            ImGuiWindowFlags_NoMove |
            ImGuiWindowFlags_NoResize |
            ImGuiWindowFlags_NoCollapse);

        ImGui::Text("Physics");
        ImGui::Separator();
        const char* physics[] = { "Heat Diffusion 2D", "Navier-Stokes 2D", "Navier-Stokes 3D" };
        static int physType = 0;
        ImGui::Combo("##phys", &physType, physics, 3);

        ImGui::Spacing();
        ImGui::Text("Material");
        ImGui::Separator();
        const char* materialNames[] = { "Custom", "Air", "Water", "Steel", "Aluminium", "Copper" };
        static const float materialAlpha[] = { 1.0e-4f, 2.1e-5f, 1.4e-7f, 1.2e-5f, 8.4e-5f, 1.11e-4f };

        if (ImGui::Combo("##material", &materialIndex, materialNames, IM_ARRAYSIZE(materialNames))) {
            if (materialIndex != 0) {
                alpha = materialAlpha[materialIndex];
            }
        }
        ImGui::SameLine();
        ImGui::Text("Material");

        ImGui::BeginDisabled(materialIndex != 0);
        ImGui::InputFloat("Alpha [m2/s]", &alpha, 0, 0, "%.2e");
        ImGui::EndDisabled();

        ImGui::Spacing();
        ImGui::Text("Boundary Conditions");
        ImGui::Separator();
        ImGui::PushStyleColor(ImGuiCol_FrameBg, { 0.4f, 0.1f, 0.1f, 1.0f });
        ImGui::InputFloat("##west", &T_west, 0.1f, 1.0f, "%.2f");
        ImGui::PopStyleColor();
        ImGui::SameLine(); ImGui::Text("T west [K]");

        ImGui::PushStyleColor(ImGuiCol_FrameBg, { 0.1f, 0.2f, 0.4f, 1.0f });
        ImGui::InputFloat("##east", &T_east, 0.1f, 1.0f, "%.2f");
        ImGui::PopStyleColor();
        ImGui::SameLine(); ImGui::Text("T east [K]");

        ImGui::PushStyleColor(ImGuiCol_FrameBg, { 0.1f, 0.2f, 0.4f, 1.0f });
        ImGui::InputFloat("##south", &T_south, 0.1f, 1.0f, "%.2f");
        ImGui::PopStyleColor();
        ImGui::SameLine(); ImGui::Text("T south [K]");

        ImGui::PushStyleColor(ImGuiCol_FrameBg, { 0.1f, 0.2f, 0.4f, 1.0f });
        ImGui::InputFloat("##north", &T_north, 0.1f, 1.0f, "%.2f");
        ImGui::PopStyleColor();
        ImGui::SameLine(); ImGui::Text("T north [K]");

        ImGui::Spacing();
        ImGui::Text("Simulation");
        ImGui::Separator();
        const char* simModes[] = { "Material Comparison", "Transient", "Near Steady", "Manual" };
        ImGui::Combo("##sim_mode", &simModeIndex, simModes, IM_ARRAYSIZE(simModes));
        ImGui::SameLine();
        ImGui::Text("Mode");

        if (simModeIndex == static_cast<int>(SimMode::Manual)) {
            ImGui::InputFloat("End Time [s]", &manualEndTime, 0.1f, 1.0f, "%.3f");
            ImGui::SliderInt("Frames", &manualFrames, 2, 120);
            ImGui::Checkbox("Loop Animation", &loopAnimation);
        }
        else {
            ImGui::BeginDisabled();
            ImGui::InputFloat("End Time [s]", &manualEndTime, 0.1f, 1.0f, "%.3f");
            ImGui::SliderInt("Frames", &manualFrames, 2, 120);
            ImGui::Checkbox("Loop Animation", &loopAnimation);
            ImGui::EndDisabled();
        }

        ImGui::Spacing();
        ImGui::Separator();

        const bool solverRunning = g_solverRunning.load();
        if (!solverRunning) {
            ImGui::PushStyleColor(ImGuiCol_Button, { 0.2f, 0.7f, 0.2f, 1.0f });
            if (ImGui::Button("▶  RUN SOLVER", { -1, 45 })) {
                if (!g_solverRunning.exchange(true)) {
                    meshEditor.clearAnimation();
                    meshEditor.clearScalarField();
                    meshEditor.showViewportWindow(true);

                    const auto settings = meshEditor.meshSettings();
                    const double L = std::max(settings.width, settings.height);
                    double tEnd = manualEndTime;
                    int frameCount = manualFrames;

                    if (simModeIndex == static_cast<int>(SimMode::MaterialComparison)) {
                        tEnd = 0.2 * L * L / 1.0e-4;
                        frameCount = 40;
                        loopAnimation = false;
                    }
                    else if (simModeIndex == static_cast<int>(SimMode::Transient)) {
                        tEnd = 0.5 * L * L / alpha;
                        frameCount = 50;
                        loopAnimation = false;
                    }
                    else if (simModeIndex == static_cast<int>(SimMode::NearSteady)) {
                        tEnd = 2.0 * L * L / alpha;
                        frameCount = 60;
                        loopAnimation = false;
                    }

                    std::thread(solverThread,
                        settings.nx, settings.ny,
                        settings.width, settings.height,
                        static_cast<double>(alpha),
                        static_cast<double>(T_west),
                        static_cast<double>(T_east),
                        static_cast<double>(T_south),
                        static_cast<double>(T_north),
                        tEnd,
                        frameCount).detach();
                }
            }
            ImGui::PopStyleColor();
        }
        else {
            ImGui::PushStyleColor(ImGuiCol_Button, { 0.8f, 0.2f, 0.2f, 1.0f });
            if (ImGui::Button("■  STOP", { -1, 45 }))
                g_solverRunning = false;
            ImGui::PopStyleColor();
        }

        if (g_solverError) {
            ImGui::TextColored(ImVec4(1, 0.3f, 0.3f, 1), "Error: %s", g_solverErrorMsg.c_str());
        }

        ImGui::Spacing();
        ImGui::Text("Progress:");
        ImGui::ProgressBar(g_solverProgress.load(), { -1, 20 });
        ImGui::End();

        ImGui::SetNextWindowPos({ centerX, topY }, ImGuiCond_Always);
        ImGui::SetNextWindowSize({ centerW, centerH }, ImGuiCond_Always);
        ImGui::Begin("Viewport", nullptr,
            ImGuiWindowFlags_NoMove |
            ImGuiWindowFlags_NoResize |
            ImGuiWindowFlags_NoCollapse |
            ImGuiWindowFlags_NoScrollbar);

        const ImVec2 size = ImGui::GetContentRegionAvail();
        const ImVec2 centre = {
            ImGui::GetCursorPosX() + size.x * 0.5f - 150,
            ImGui::GetCursorPosY() + size.y * 0.5f - 30 };

        ImGui::SetCursorPos(centre);
        ImGui::TextDisabled("OpenGL viewport renders here");
        ImGui::SetCursorPosX(centre.x + 20);
        ImGui::TextDisabled(meshEditor.hasMesh()
            ? "OBJ mesh loaded — see Mesh Viewport"
            : meshEditor.hasScalarField()
            ? "Heat map active — see Mesh Viewport"
            : "Heat map / streamlines");
        ImGui::End();

        ImGui::SetNextWindowPos({ centerX, topY + centerH }, ImGuiCond_Always);
        ImGui::SetNextWindowSize({ centerW, bottomTimelineH }, ImGuiCond_Always);
        ImGui::Begin("##timeline", nullptr,
            ImGuiWindowFlags_NoTitleBar |
            ImGuiWindowFlags_NoResize |
            ImGuiWindowFlags_NoMove |
            ImGuiWindowFlags_NoCollapse);

        const int totalFrames = meshEditor.animFrameCount();
        if (totalFrames > 0) {
            ImGui::Text("t=%.2fs", meshEditor.animTime(0));
            ImGui::SameLine();

            int frameIdx = meshEditor.currentAnimFrame();
            ImGui::SetNextItemWidth(std::max(120.0f, centerW - 280.0f));
            if (ImGui::SliderInt("##tslider", &frameIdx, 0, totalFrames - 1, ""))
                meshEditor.setAnimFrame(frameIdx);

            ImGui::SameLine();
            ImGui::Text("t=%.2fs", meshEditor.animEndTime());
            ImGui::SameLine();
            ImGui::TextDisabled("[t=%.3fs]", meshEditor.animTime(meshEditor.currentAnimFrame()));
        }
        else {
            ImGui::TextDisabled("t = 0.0s");
            ImGui::SameLine();
            ImGui::BeginDisabled();
            float dummy = 0.0f;
            ImGui::SetNextItemWidth(std::max(120.0f, centerW - 280.0f));
            ImGui::SliderFloat("##tslider_empty", &dummy, 0.0f, 1.0f, "");
            ImGui::EndDisabled();
            ImGui::SameLine();
            ImGui::TextDisabled("t = --");
        }

        ImGui::End();
        window.endFrame();
    }

    g_solverRunning = false;
    NFD_Quit();
    window.cleanup();
    return 0;
}
