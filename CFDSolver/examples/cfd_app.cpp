/**
 * CFD Solver Application
 * ImGui + OpenGL UI  –  with real-time heat-diffusion animation
 */
#include "renderer/Window.hpp"
#include "renderer/MeshEditor2D.hpp"
#include "solvers/HeatSolver2D.h"
#include "linearAlgebra/ConjugateGradient.hpp"
#include "mesh/mesh2D.h"

#include <imgui.h>
#include <iostream>
#include <string>
#include <vector>
#include <thread>
#include <mutex>
#include <atomic>
#include <deque>

 // ── Thread-safe frame queue ───────────────────────────────────────────────────
struct PendingFrame {
    std::vector<double> values;
    double              time = 0.0;
    double              globalMin = 0.0;
    double              globalMax = 1.0;
};

static std::mutex              g_frameMutex;
static std::deque<PendingFrame> g_frameQueue;   // solver thread pushes, main pops
static std::atomic<bool>        g_solverRunning{ false };
static std::atomic<float>       g_solverProgress{ 0.0f };
static std::atomic<bool>        g_solverError{ false };
static std::string              g_solverErrorMsg;

// ── Solver thread entry point ─────────────────────────────────────────────────
static void solverThread(
    int    meshNx, int meshNy,
    double width, double height,
    double alpha,
    double T_hot, double T_cold)
{
    g_solverError = false;
    g_solverProgress = 0.0f;

    try {
        Mesh mesh(meshNx, meshNy, width, height);

        // ── Compute a stable, physically meaningful dt ────────────────────
        // Fourier number Fo = alpha*dt/dx^2 — keep it at 0.1
        // so diffusion is slow enough to animate clearly.
        const double dx = width / static_cast<double>(meshNx);
        const double dy = height / static_cast<double>(meshNy);
        const double dMin = std::min(dx, dy);

        // Target Fo = 0.1  →  dt = 0.1 * dMin^2 / alpha
        const double Fo_target = 0.1;
        const double dt = Fo_target * dMin * dMin / alpha;

        // Run for enough time to reach steady state visually.
        // Diffusion time scale: L^2 / alpha
        const double L = std::max(width, height);
        const double tEnd = 2.0 * L * L / alpha;

        CFD::HeatSolver2D solver(mesh, alpha, dt, tEnd);

        CFD::HeatSolver2D::BoundaryConditions bcs;
        bcs.T_west = T_hot;
        bcs.T_east = T_cold;
        bcs.T_south = T_cold;
        bcs.T_north = T_cold;

        solver.setBCs(bcs);
        solver.setIC(T_cold);
        solver.setOutputFreq(10);
        solver.enableSteadyStop(true);
        solver.setSteadyTolerance(1.0e-6);

        // ── Collect ALL frames first, then set global min/max ──
        // Target ~60 animation frames regardless of step count.
        const int estimatedSteps =
            static_cast<int>(std::ceil(tEnd / dt));
        const int recordEvery =
            std::max(1, estimatedSteps / 60);

        CFD::ConjugateGradient cg(1e-10, 10000);
        std::vector<int> cgIters;

        // Collect raw snapshots (no normalisation yet)
        struct RawFrame {
            std::vector<double> values;
            double              time = 0.0;
        };
        std::vector<RawFrame> rawFrames;
        rawFrames.reserve(64);

        double globalMin = 1e30;
        double globalMax = -1e30;

        while (!solver.finished() && g_solverRunning) {
            solver.step(cg, cgIters);

            if (solver.steps() % recordEvery == 0 || solver.steps() == 1) {
                const CFD::ScalarField& T = solver.T();
                RawFrame rf;
                rf.values.resize(T.size());
                for (std::size_t k = 0; k < T.size(); ++k) {
                    rf.values[k] = T[k];
                    globalMin = std::min(globalMin, T[k]);
                    globalMax = std::max(globalMax, T[k]);
                }
                rf.time = solver.time();
                rawFrames.push_back(std::move(rf));
            }

            g_solverProgress = static_cast<float>(solver.time() / tEnd);
        }

        // Clamp global range
        if (globalMax <= globalMin) globalMax = globalMin + 1.0;

        // Push all frames into the shared queue at once
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

// ─────────────────────────────────────────────────────────────────────────────
int main()
{
    CFD::Window window(1280, 720, "CFD Solver");
    if (!window.init()) {
        std::cerr << "Failed to init!\n";
        return -1;
    }

    float alpha = 1e-4f;
    float T_hot = 1.0f;
    float T_cold = 0.0f;
    int   nx = 50, ny = 50;

    CFD::UI::MeshEditor2D meshEditor;
    meshEditor.showEditorWindow(false);
    meshEditor.showViewportWindow(false);

    std::cout << "App running!\n";

    while (!window.shouldClose())
    {
        window.beginFrame();

        // ── Drain the frame queue into meshEditor (main thread only) ──────
        {
            std::lock_guard<std::mutex> lock(g_frameMutex);
            while (!g_frameQueue.empty()) {
                auto& pf = g_frameQueue.front();
                meshEditor.pushAnimationFrame(
                    pf.values, pf.time, pf.globalMin, pf.globalMax);
                g_frameQueue.pop_front();
            }
        }

        // ── Menu bar ──────────────────────────────────────────────────────
        if (ImGui::BeginMainMenuBar()) {
            if (ImGui::BeginMenu("File")) {
                ImGui::MenuItem("New");
                ImGui::MenuItem("Open");
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
            if (ImGui::BeginMenu("Help")) {
                if (ImGui::MenuItem("About")) window.openAbout();
                ImGui::EndMenu();
            }
            ImGui::EndMainMenuBar();
        }
        window.renderAbout();

        // ── Mesh editor windows ───────────────────────────────────────────
        meshEditor.drawUI();
        meshEditor.drawViewport();

        // ── Left toolbar ──────────────────────────────────────────────────
        ImGui::SetNextWindowPos({ 0, 20 });
        ImGui::SetNextWindowSize({ 55, static_cast<float>(window.height()) - 20 });
        ImGui::Begin("##toolbar", nullptr,
            ImGuiWindowFlags_NoTitleBar |
            ImGuiWindowFlags_NoResize |
            ImGuiWindowFlags_NoScrollbar);

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

        // ── Right properties panel ────────────────────────────────────────
        ImGui::SetNextWindowPos({
            static_cast<float>(window.width()) - 270, 20 });
        ImGui::SetNextWindowSize({
            270, static_cast<float>(window.height()) - 20 });
        ImGui::Begin("Properties");

        ImGui::Text("Physics");
        ImGui::Separator();
        const char* physics[] = {
            "Heat Diffusion 2D", "Navier-Stokes 2D", "Navier-Stokes 3D" };
        static int physType = 0;
        ImGui::Combo("##phys", &physType, physics, 3);

        ImGui::Spacing();
        ImGui::Text("Material");
        ImGui::Separator();
        ImGui::InputFloat("Alpha [m2/s]", &alpha, 0, 0, "%.2e");

        ImGui::Spacing();
        ImGui::Text("Boundary Conditions");
        ImGui::Separator();
        ImGui::PushStyleColor(ImGuiCol_FrameBg, { 0.4f, 0.1f, 0.1f, 1.0f });
        ImGui::InputFloat("T hot [K]", &T_hot, 0.1f, 1.0f, "%.2f");
        ImGui::PopStyleColor();
        ImGui::PushStyleColor(ImGuiCol_FrameBg, { 0.1f, 0.2f, 0.4f, 1.0f });
        ImGui::InputFloat("T cold [K]", &T_cold, 0.1f, 1.0f, "%.2f");
        ImGui::PopStyleColor();

        ImGui::Spacing();
        ImGui::Text("Mesh");
        ImGui::Separator();
        ImGui::SliderInt("Cells X", &nx, 10, 200);
        ImGui::SliderInt("Cells Y", &ny, 10, 200);
        ImGui::Text("Total: %d cells", nx * ny);

        ImGui::Spacing();
        ImGui::Separator();

        const bool solverRunning = g_solverRunning.load();

        if (!solverRunning) {
            ImGui::PushStyleColor(ImGuiCol_Button, { 0.2f, 0.7f, 0.2f, 1.0f });
            if (ImGui::Button("▶  RUN SOLVER", { -1, 45 }))
            {
                // Clear previous results
                meshEditor.clearAnimation();
                meshEditor.clearScalarField();
                {
                    std::lock_guard<std::mutex> lock(g_frameMutex);
                    g_frameQueue.clear();
                }
                g_solverProgress = 0.0f;
                g_solverError = false;

                const auto   settings = meshEditor.meshSettings();
                g_solverRunning = true;

                // Launch solver on background thread
                std::thread(solverThread,
                    settings.nx, settings.ny,
                    settings.width, settings.height,
                    static_cast<double>(alpha),
                    static_cast<double>(T_hot),
                    static_cast<double>(T_cold)
                ).detach();

                meshEditor.showViewportWindow(true);
            }
            ImGui::PopStyleColor();

            if (g_solverError) {
                ImGui::TextColored({ 1,0.3f,0.3f,1 },
                    "Error: %s", g_solverErrorMsg.c_str());
            }
        }
        else {
            // Solver is running — show stop button + spinner
            ImGui::PushStyleColor(ImGuiCol_Button, { 0.8f, 0.2f, 0.2f, 1.0f });
            if (ImGui::Button("■  STOP", { -1, 45 }))
                g_solverRunning = false;   // signals thread to exit
            ImGui::PopStyleColor();
        }

        ImGui::Spacing();
        ImGui::Text("Progress:");
        ImGui::ProgressBar(g_solverProgress.load(), { -1, 20 });

        ImGui::End();

        // ── Centre viewport placeholder ───────────────────────────────────
        ImGui::SetNextWindowPos({ 55, 20 });
        ImGui::SetNextWindowSize({
            static_cast<float>(window.width()) - 325,
            static_cast<float>(window.height()) - 80 });
        ImGui::Begin("Viewport", nullptr, ImGuiWindowFlags_NoScrollbar);

        const ImVec2 size = ImGui::GetContentRegionAvail();
        const ImVec2 centre = {
            ImGui::GetCursorPosX() + size.x * 0.5f - 150,
            ImGui::GetCursorPosY() + size.y * 0.5f - 30 };

        ImGui::SetCursorPos(centre);
        ImGui::TextDisabled("OpenGL viewport renders here");
        ImGui::SetCursorPosX(centre.x + 20);
        ImGui::TextDisabled(meshEditor.hasScalarField()
            ? "Showing real heat map in Mesh Viewport"
            : "Heat map / streamlines");

        ImGui::End();

        // ── Bottom timeline – drives animation frame ──────────────────────
        ImGui::SetNextWindowPos({
            55, static_cast<float>(window.height()) - 60 });
        ImGui::SetNextWindowSize({
            static_cast<float>(window.width()) - 325, 60 });
        ImGui::Begin("##timeline", nullptr,
            ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoResize);

        const int totalFrames = meshEditor.animFrameCount();

        if (totalFrames > 1) {
            ImGui::Text("t=%.2fs", meshEditor.animTime(0));
            ImGui::SameLine();

            int frameIdx = meshEditor.currentAnimFrame();
            ImGui::SetNextItemWidth(
                static_cast<float>(window.width()) - 500.0f);
            if (ImGui::SliderInt("##tslider", &frameIdx,
                0, totalFrames - 1, ""))
            {
                meshEditor.setAnimFrame(frameIdx);
            }
            ImGui::SameLine();
            ImGui::Text("t=%.2fs", meshEditor.animEndTime());
            ImGui::SameLine();
            ImGui::TextDisabled("[t=%.3fs]",
                meshEditor.animTime(meshEditor.currentAnimFrame()));
        }
        else {
            // No animation yet
            ImGui::TextDisabled("t = 0.0s");
            ImGui::SameLine();
            ImGui::BeginDisabled();
            float dummy = 0.0f;
            ImGui::SetNextItemWidth(
                static_cast<float>(window.width()) - 500.0f);
            ImGui::SliderFloat("##tslider_empty", &dummy, 0.0f, 1.0f, "");
            ImGui::EndDisabled();
            ImGui::SameLine();
            ImGui::TextDisabled(solverRunning
                ? "Solving..." : "t = --");
        }

        ImGui::End();

        window.endFrame();
    }

    // Signal solver thread to stop cleanly before exit
    g_solverRunning = false;
    window.cleanup();
    return 0;
}
