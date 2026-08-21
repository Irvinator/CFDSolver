/**
 * CFD Solver Application
 * ImGui + OpenGL UI
 */
#include "renderer/Window.hpp"
#include <imgui.h>
#include <iostream>

int main() {

    CFD::Window window(1280, 720,
        "CFD Solver");

    if (!window.init()) {
        std::cerr << "Failed to init!\n";
        return -1;
    }

    // App state
    float alpha = 1e-4f;
    float T_hot = 1.0f;
    float T_cold = 0.0f;
    int   nx = 50, ny = 50;
    bool  running = false;
    float progress = 0.0f;

    std::cout << "App running!\n";

    while (!window.shouldClose()) {
        window.beginFrame();

        // ── Menu bar ──────────────────────────
        if (ImGui::BeginMainMenuBar()) {
            if (ImGui::BeginMenu("File")) {
                ImGui::MenuItem("New");
                ImGui::MenuItem("Open");
                ImGui::MenuItem("Save");
                ImGui::Separator();
                if (ImGui::MenuItem("Exit"))
                    break;
                ImGui::EndMenu();
            }
            if (ImGui::BeginMenu("Physics")) {
                ImGui::MenuItem(
                    "Heat Diffusion");
                ImGui::MenuItem(
                    "Navier-Stokes");
                ImGui::EndMenu();
            }
            if (ImGui::BeginMenu("Help")) {
                ImGui::MenuItem("About");
                ImGui::EndMenu();
            }
            ImGui::EndMainMenuBar();
        }

        // ── Left toolbar ──────────────────────
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
        ImGui::Button("GEO\n   ", { 40,50 });
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

        // ── Properties panel ──────────────────
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

        // Run button
        if (!running) {
            ImGui::PushStyleColor(
                ImGuiCol_Button,
                { 0.2f,0.7f,0.2f,1.0f });
            if (ImGui::Button(
                "▶  RUN SOLVER", { -1,45 }))
                running = true;
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

            // Animate progress
            progress += 0.001f;
            if (progress > 1.0f) {
                progress = 0.0f;
                running = false;
            }
        }

        ImGui::Spacing();
        ImGui::Text("Progress:");
        ImGui::ProgressBar(progress,
            { -1, 20 });

        ImGui::End();

        // ── Main viewport ─────────────────────
        ImGui::SetNextWindowPos({ 55, 20 });
        ImGui::SetNextWindowSize({
            (float)window.width() - 325,
            (float)window.height() - 80 });
        ImGui::Begin("Viewport", nullptr,
            ImGuiWindowFlags_NoScrollbar);

        // Placeholder text
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
            "Heat map / streamlines");

        ImGui::End();

        // ── Bottom timeline ───────────────────
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