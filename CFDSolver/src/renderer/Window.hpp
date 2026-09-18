#pragma once
#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <imgui.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl3.h>
#include <iostream>
#include <string>

namespace CFD {

    class Window {
    public:

        Window(int width, int height,
            const std::string& title)
            : initWidth_(width)
            , initHeight_(height)
            , title_(title) {
        }

        // ── Initialise ────────────────────────────────────────────
        bool init() {
            if (!glfwInit()) {
                std::cerr << "GLFW init failed\n";
                return false;
            }

            glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
            glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
            glfwWindowHint(GLFW_OPENGL_PROFILE,
                GLFW_OPENGL_CORE_PROFILE);

            window_ = glfwCreateWindow(
                initWidth_, initHeight_,
                title_.c_str(),
                nullptr, nullptr);

            if (!window_) {
                std::cerr << "Window failed\n";
                glfwTerminate();
                return false;
            }

            glfwMakeContextCurrent(window_);
            glfwSwapInterval(1);

            // Keep OpenGL viewport in sync with OS resize / fullscreen
            glfwSetFramebufferSizeCallback(window_,
                [](GLFWwindow*, int w, int h) {
                    glViewport(0, 0, w, h);
                });

            if (!gladLoadGLLoader(
                (GLADloadproc)glfwGetProcAddress)) {
                std::cerr << "GLAD failed\n";
                return false;
            }

            IMGUI_CHECKVERSION();
            ImGui::CreateContext();
            ImGui::StyleColorsDark();
            applyStyle();

            ImGui_ImplGlfw_InitForOpenGL(window_, true);
            ImGui_ImplOpenGL3_Init("#version 330");

            std::cout << "Window ready!\n";
            std::cout << "OpenGL: "
                << glGetString(GL_VERSION) << "\n";
            return true;
        }

        // ── Frame ─────────────────────────────────────────────────
        void beginFrame() {
            glfwPollEvents();
            ImGui_ImplOpenGL3_NewFrame();
            ImGui_ImplGlfw_NewFrame();
            ImGui::NewFrame();
            glClearColor(0.08f, 0.08f, 0.10f, 1.0f);
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        }

        void endFrame() {
            ImGui::Render();
            ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
            glfwSwapBuffers(window_);
        }

        // ── State ─────────────────────────────────────────────────
        bool shouldClose() const {
            return glfwWindowShouldClose(window_);
        }

        // LIVE size query — always correct after fullscreen / OS resize
        int width() const {
            int w, h;
            glfwGetWindowSize(window_, &w, &h);
            return w;
        }
        int height() const {
            int w, h;
            glfwGetWindowSize(window_, &w, &h);
            return h;
        }

        // ── About dialog ──────────────────────────────────────────
        void openAbout() { showAbout_ = true; }

        void renderAbout() {
            if (!showAbout_) return;

            const float W = static_cast<float>(width());
            const float H = static_cast<float>(height());

            ImGui::SetNextWindowPos(
                ImVec2(W * 0.5f, H * 0.5f),
                ImGuiCond_Always,
                ImVec2(0.5f, 0.5f));

            ImGui::SetNextWindowSize(
                ImVec2(450, 320), ImGuiCond_Always);

            ImGui::Begin("About CFD Solver",
                &showAbout_,
                ImGuiWindowFlags_NoResize |
                ImGuiWindowFlags_NoMove);

            ImGui::SetCursorPosX(
                (450 - ImGui::CalcTextSize("CFD Solver").x) * 0.5f);
            ImGui::TextColored(
                { 0.4f, 0.65f, 1.0f, 1.0f }, "CFD Solver");

            ImGui::Separator();
            ImGui::Spacing();
            ImGui::Text("Version:     0.1.0");
            ImGui::Text("Build:       x64-Debug");
            ImGui::Spacing();
            ImGui::Separator();
            ImGui::Spacing();

            ImGui::TextColored(
                { 0.4f, 0.65f, 1.0f, 1.0f }, "Development Team");
            ImGui::Spacing();
            ImGui::Text("Irvin    - Core Solver & Numerics");
            ImGui::Text("Alex&Irvin - OpenGL Visualisation");
            ImGui::Text("Alex - Mesh & IO");
            ImGui::Text("Mack - Navier-Stokes");

            ImGui::Spacing();
            ImGui::Separator();
            ImGui::Spacing();

            ImGui::TextWrapped(
                "A high-performance CFD solver "
                "built from scratch in C++17. "
                "Solving heat diffusion and "
                "Navier-Stokes equations using "
                "the Finite Volume Method.");

            ImGui::Spacing();
            ImGui::Separator();
            ImGui::Spacing();

            float buttonWidth = 120.0f;
            ImGui::SetCursorPosX((450 - buttonWidth) * 0.5f);
            if (ImGui::Button("Close", ImVec2(buttonWidth, 30)))
                showAbout_ = false;

            ImGui::End();
        }

        // ── Cleanup ───────────────────────────────────────────────
        void cleanup() {
            ImGui_ImplOpenGL3_Shutdown();
            ImGui_ImplGlfw_Shutdown();
            ImGui::DestroyContext();
            glfwDestroyWindow(window_);
            glfwTerminate();
        }

    private:
        GLFWwindow* window_ = nullptr;
        int         initWidth_, initHeight_;
        std::string title_;
        bool        showAbout_ = false;

        void applyStyle() {
            ImGuiStyle& s = ImGui::GetStyle();
            s.WindowRounding = 4.0f;
            s.FrameRounding = 3.0f;
            s.GrabRounding = 3.0f;
            s.WindowBorderSize = 0.0f;
            s.FrameBorderSize = 0.0f;
            s.WindowPadding = { 8, 8 };
            s.FramePadding = { 6, 4 };
            s.ItemSpacing = { 6, 4 };

            ImVec4* c = s.Colors;
            c[ImGuiCol_WindowBg] = { 0.15f,0.15f,0.15f,1.0f };
            c[ImGuiCol_ChildBg] = { 0.13f,0.13f,0.13f,1.0f };
            c[ImGuiCol_Button] = { 0.26f,0.26f,0.26f,1.0f };
            c[ImGuiCol_ButtonHovered] = { 0.35f,0.35f,0.35f,1.0f };
            c[ImGuiCol_ButtonActive] = { 0.45f,0.45f,0.45f,1.0f };
            c[ImGuiCol_FrameBg] = { 0.20f,0.20f,0.20f,1.0f };
            c[ImGuiCol_Header] = { 0.26f,0.26f,0.26f,1.0f };
            c[ImGuiCol_TitleBg] = { 0.10f,0.10f,0.10f,1.0f };
            c[ImGuiCol_TitleBgActive] = { 0.15f,0.15f,0.15f,1.0f };
            c[ImGuiCol_Tab] = { 0.18f,0.18f,0.18f,1.0f };
            c[ImGuiCol_TabActive] = { 0.26f,0.26f,0.26f,1.0f };
            c[ImGuiCol_Separator] = { 0.30f,0.30f,0.30f,1.0f };
            c[ImGuiCol_CheckMark] = { 0.40f,0.65f,1.0f,1.0f };
            c[ImGuiCol_SliderGrab] = { 0.40f,0.65f,1.0f,1.0f };
        }
    };

} // namespace CFD
