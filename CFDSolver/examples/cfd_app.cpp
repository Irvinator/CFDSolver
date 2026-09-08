/**
 * CFD Solver Application
 * ImGui + OpenGL UI
 */
#include "renderer/Window.hpp"
#include "renderer/MeshEditor2D.hpp"
#include "solvers/HeatSolver2D.h"
#include "mesh/mesh2D.h"

#include <imgui.h>
#include <iostream>
#include <string>

int main() {

    CFD::Window window(1280, 720,
        "CFD Solver");

    if (!window.init()) {
        std::cerr << "Failed to init!\n";
        return -1;
    }

    float alpha = 1e-4f;
    float T_hot = 1.0f;
    float T_cold = 0.0f;
    int   nx = 50, ny = 50;
    bool  running = false;
    float progress = 0.0f;

    std::cout << "App running!\n";
    CFD::UI::MeshEditor2D meshEditor;
    meshEditor.showEditorWindow(false);
    meshEditor.showViewportWindow(false);

    while (!window.shouldClose()) {
        window.beginFrame();

        if (ImGui::BeginMainMenuBar()) {
            if (ImGui::BeginMenu("File")) {
                ImGui::MenuItem("New");
                ImGui::MenuItem("Open");
                ImGui::Separator();
                if (ImGui::MenuItem("Exit"))
                    break;
                ImGui::EndMenu();
            }
            if (ImGui::BeginMenu("View")) {
                bool showEditor = meshEditor.isEditorWindowVisible();
                bool showViewport = meshEditor.isViewportWindowVisible();
                if (ImGui::MenuItem("Mesh Editor", nullptr, showEditor)) {
                    meshEditor.showEditorWindow(!showEditor);
                }
                if (ImGui::MenuItem("Mesh Viewport", nullptr, showViewport)) {
                    meshEditor.showViewportWindow(!showViewport);
                }
                ImGui::EndMenu();
            }
            if (ImGui::BeginMenu("Help")) {
                if (ImGui::MenuItem("About"))
                    window.openAbout();
                ImGui::EndMenu();
            }
            ImGui::EndMainMenuBar();
        }
        window.renderAbout();

        meshEditor.drawUI();
        meshEditor.drawViewport();

        ImGui::SetNextWindowPos({ 0, 20 });
        ImGui::SetNextWindowSize(
            { 55, (float)window.height() - 20 });
        ImGui::Begin("##toolbar", nullptr,
            ImGuiWindowFlags_NoTitleBar |
            ImGuiWindowFlags_NoResize |
            ImGuiWindowFlags_NoScrollbar);

        ImGui::PushStyleColor(
            ImGuiCol_Button,
            { 0.3f, 0.5f, 0.9f, 1.0f });
        if (ImGui::Button("GEO\n   ", { 40,50 })) {
            meshEditor.showEditorWindow(true);
            meshEditor.showViewportWindow(true);
        }
        ImGui::PopStyleColor();
        ImGui::SetItemTooltip(
            "Geometry Mode");

        ImGui::Button("PHY\n   ", { 40,50 });
        ImGui::SetItemTooltip(
            "Physics Setup");

        ImGui::Button("MSH\n   ", { 40,50 });
        ImGui::SetItemTooltip("Mesh");

        ImGui::PushStyleColor(
            ImGuiCol_Button,
            { 0.2f, 0.7f, 0.2f, 1.0f });
        ImGui::Button("RUN\n   ", { 40,50 });
        ImGui::PopStyleColor();
        ImGui::SetItemTooltip("Run Solver");

        ImGui::Button("RES\n   ", { 40,50 });
        ImGui::SetItemTooltip("Results");

        ImGui::End();

        ImGui::SetNextWindowPos(
            { (float)window.width() - 270, 20 });
        ImGui::SetNextWindowSize(
            { 270, (float)window.height() - 20 });
        ImGui::Begin("Properties");

        ImGui::Text("Physics");
        ImGui::Separator();

        const char* physics[] = {
            "Heat Diffusion 2D",
            "Navier-Stokes 2D",
            "Navier-Stokes 3D"
        };
        static int physType = 0;
        ImGui::Combo("##phys",
            &physType, physics, 3);

        ImGui::Spacing();
        ImGui::Text("Material");
        ImGui::Separator();
        ImGui::InputFloat("Alpha [m2/s]",
            &alpha, 0, 0, "%.2e");

        ImGui::Spacing();
        ImGui::Text("Boundary Conditions");
        ImGui::Separator();

        ImGui::PushStyleColor(
            ImGuiCol_FrameBg,
            { 0.4f, 0.1f, 0.1f, 1.0f });
        ImGui::InputFloat("T hot [K]",
            &T_hot, 0.1f, 1.0f, "%.2f");
        ImGui::PopStyleColor();

        ImGui::PushStyleColor(
            ImGuiCol_FrameBg,
            { 0.1f, 0.2f, 0.4f, 1.0f });
        ImGui::InputFloat("T cold [K]",
            &T_cold, 0.1f, 1.0f, "%.2f");
        ImGui::PopStyleColor();

        ImGui::Spacing();
        ImGui::Text("Mesh");
        ImGui::Separator();
        ImGui::SliderInt("Cells X",
            &nx, 10, 200);
        ImGui::SliderInt("Cells Y",
            &ny, 10, 200);
        ImGui::Text("Total: %d cells",
            nx * ny);

        ImGui::Spacing();
        ImGui::Separator();

        if (!running) {
            ImGui::PushStyleColor(
                ImGuiCol_Button,
                { 0.2f,0.7f,0.2f,1.0f });
            if (ImGui::Button(
                "▶  RUN SOLVER", { -1,45 }))
            {
                running = true;
                progress = 0.0f;

                try
                {
                    const auto settings = meshEditor.meshSettings();
                    const double width = settings.width;
                    const double height = settings.height;
                    const int meshNx = settings.nx;
                    const int meshNy = settings.ny;

                    Mesh mesh(meshNx, meshNy, width, height);
                    CFD::HeatSolver2D solver(mesh, static_cast<double>(alpha), 0.1, 10.0);

                    CFD::HeatSolver2D::BoundaryConditions bcs;
                    bcs.T_west = static_cast<double>(T_hot);
                    bcs.T_east = static_cast<double>(T_cold);
                    bcs.T_south = static_cast<double>(T_cold);
                    bcs.T_north = static_cast<double>(T_cold);

                    solver.setBCs(bcs);
                    solver.setIC(static_cast<double>(T_cold));
                    solver.setOutputFreq(1000);
                    solver.enableSteadyStop(true);
                    solver.setSteadyTolerance(1.0e-10);

                    solver.run();
                    meshEditor.setScalarField(solver.T(), meshNx, meshNy);
                    progress = 1.0f;
                    running = false;
                    meshEditor.showViewportWindow(true);
                }
                catch (const std::exception& e)
                {
                    std::cerr << "HeatSolver2D failed: " << e.what() << "\n";
                    running = false;
                    meshEditor.clearScalarField();
                }
            }
            ImGui::PopStyleColor();
        }
        else {
            ImGui::PushStyleColor(
                ImGuiCol_Button,
                { 0.8f,0.2f,0.2f,1.0f });
            if (ImGui::Button(
                "■  STOP", { -1,45 }))
                running = false;
            ImGui::PopStyleColor();
        }

        ImGui::Spacing();
        ImGui::Text("Progress:");
        ImGui::ProgressBar(progress,
            { -1, 20 });

        ImGui::End();

        ImGui::SetNextWindowPos({ 55, 20 });
        ImGui::SetNextWindowSize({
            (float)window.width() - 325,
            (float)window.height() - 80 });
        ImGui::Begin("Viewport", nullptr,
            ImGuiWindowFlags_NoScrollbar);

        ImVec2 size =
            ImGui::GetContentRegionAvail();
        ImVec2 centre = {
            ImGui::GetCursorPosX() +
            size.x * 0.5f - 150,
            ImGui::GetCursorPosY() +
            size.y * 0.5f - 30
        };

        ImGui::SetCursorPos(centre);
        ImGui::TextDisabled(
            "OpenGL viewport renders here");
        ImGui::SetCursorPosX(centre.x + 20);
        ImGui::TextDisabled(
            meshEditor.hasScalarField() ? "Showing real heat map in Mesh Viewport" : "Heat map / streamlines");

        ImGui::End();

        ImGui::SetNextWindowPos({
            55,
            (float)window.height() - 60 });
        ImGui::SetNextWindowSize({
            (float)window.width() - 325, 60 });
        ImGui::Begin("##timeline", nullptr,
            ImGuiWindowFlags_NoTitleBar |
            ImGuiWindowFlags_NoResize);

        ImGui::Text("t = 0.0s");
        ImGui::SameLine();
        ImGui::SetNextItemWidth(
            window.width() - 500.0f);
        ImGui::SliderFloat("##tslider",
            &progress, 0.0f, 1.0f, "");
        ImGui::SameLine();
        ImGui::Text("t = 100.0s");

        ImGui::End();

        window.endFrame();
    }

    window.cleanup();
    return 0;
}
