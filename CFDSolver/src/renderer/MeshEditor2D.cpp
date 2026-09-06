#include "MeshEditor2D.hpp"

#include <algorithm>
#include <imgui.h>

namespace {
    ImU32 scalarToColor(double value)
    {
        value = std::clamp(value, 0.0, 1.0);

        const int r = static_cast<int>(255.0 * value);
        const int g = 0;
        const int b = static_cast<int>(255.0 * (1.0 - value));

        return IM_COL32(r, g, b, 255);
    }
}

namespace CFD::UI {

    namespace {
        ImVec2 toImVec2(double x, double y)
        {
            return ImVec2(static_cast<float>(x), static_cast<float>(y));
        }

        void normalizeRect(RectDomain& rect)
        {
            if (rect.x0 > rect.x1) {
                std::swap(rect.x0, rect.x1);
            }
            if (rect.y0 > rect.y1) {
                std::swap(rect.y0, rect.y1);
            }
        }
    }

    void MeshEditor2D::drawUI()
    {
        ImGui::Begin("Mesh Editor 2D");

        if (ImGui::Button("Draw Rectangle")) {
            drawMode_ = true;
            drawing_ = false;
            domainReady_ = false;
        }

        ImGui::Separator();
        ImGui::Text("Mesh Settings");
        ImGui::InputInt("Nx", &mesh_.nx);
        ImGui::InputInt("Ny", &mesh_.ny);

        if (mesh_.nx < 1) mesh_.nx = 1;
        if (mesh_.ny < 1) mesh_.ny = 1;

        ImGui::Checkbox("Show Grid", &showGrid_);
        ImGui::Checkbox("Show Colormap", &showColormap_);

        ImGui::Separator();
        ImGui::Text("Rectangle Domain");
        ImGui::Text("x0 = %.1f, y0 = %.1f", domain_.x0, domain_.y0);
        ImGui::Text("x1 = %.1f, y1 = %.1f", domain_.x1, domain_.y1);
        ImGui::Text("Width  = %.1f", domain_.width());
        ImGui::Text("Height = %.1f", domain_.height());

        if (domainReady_) {
            ImGui::TextColored(ImVec4(0.2f, 0.8f, 0.2f, 1.0f), "Domain ready");
            
        }
        else if (drawMode_) {
            ImGui::TextColored(ImVec4(0.9f, 0.8f, 0.2f, 1.0f), "Click and drag in viewport");
        }
        else {
            ImGui::Text("No domain yet");
        }

        ImGui::End();
    }

    void MeshEditor2D::drawViewport()
    {
        ImGui::Begin("Mesh Viewport");

        const ImVec2 canvasPos = ImGui::GetCursorScreenPos();
        const ImVec2 canvasSize = ImGui::GetContentRegionAvail();
        const ImVec2 canvasEnd(canvasPos.x + canvasSize.x, canvasPos.y + canvasSize.y);

        ImDrawList* drawList = ImGui::GetWindowDrawList();

        drawList->AddRectFilled(canvasPos, canvasEnd, IM_COL32(25, 25, 30, 255));
        drawList->AddRect(canvasPos, canvasEnd, IM_COL32(200, 200, 200, 255));

        ImGui::InvisibleButton("mesh_canvas", canvasSize,
            ImGuiButtonFlags_MouseButtonLeft);
        
        const bool hovered = ImGui::IsItemHovered();
        const ImVec2 mouse = ImGui::GetIO().MousePos;

        if (drawMode_ && hovered && ImGui::IsMouseClicked(ImGuiMouseButton_Left)) {
            domain_.x0 = mouse.x;
            domain_.y0 = mouse.y;
            domain_.x1 = mouse.x;
            domain_.y1 = mouse.y;
            drawing_ = true;
            domainReady_ = false;
        }

        if (drawMode_ && drawing_) {
            domain_.x1 = mouse.x;
            domain_.y1 = mouse.y;

            if (ImGui::IsMouseReleased(ImGuiMouseButton_Left)) {
                normalizeRect(domain_);
                drawing_ = false;
                domainReady_ = domain_.valid();
                drawMode_ = false;
            }
        }

        if (drawing_ || domainReady_) {
            const ImVec2 p0 = toImVec2(domain_.x0, domain_.y0);
            const ImVec2 p1 = toImVec2(domain_.x1, domain_.y1);

            drawList->AddRect(p0, p1, IM_COL32(80, 180, 255, 255), 0.0f, 0, 2.0f);
        }

        if (domainReady_ && showGrid_) {
            const double width = domain_.width();
            const double height = domain_.height();
            const double dx = width / static_cast<double>(mesh_.nx);
            const double dy = height / static_cast<double>(mesh_.ny);

            for (int i = 0; i <= mesh_.nx; ++i) {
                const double x = domain_.x0 + static_cast<double>(i) * dx;
                drawList->AddLine(
                    toImVec2(x, domain_.y0),
                    toImVec2(x, domain_.y1),
                    IM_COL32(120, 120, 120, 255)
                );
            }

            for (int j = 0; j <= mesh_.ny; ++j) {
                const double y = domain_.y0 + static_cast<double>(j) * dy;
                drawList->AddLine(
                    toImVec2(domain_.x0, y),
                    toImVec2(domain_.x1, y),
                    IM_COL32(120, 120, 120, 255)
                );
            }
        }
        if (domainReady_ && showColormap_) {
            const double width = domain_.width();
            const double height = domain_.height();

            const double dx = width / static_cast<double>(mesh_.nx);
            const double dy = height / static_cast<double>(mesh_.ny);

            for (int j = 0; j < mesh_.ny; ++j) {
                for (int i = 0; i < mesh_.nx; ++i) {
                    const double x0 = domain_.x0 + static_cast<double>(i) * dx;
                    const double y0 = domain_.y0 + static_cast<double>(j) * dy;
                    const double x1 = x0 + dx;
                    const double y1 = y0 + dy;

                    const double xc = x0 + 0.5 * dx;

                    // Fake scalar field: left-to-right gradient
                    const double phi = (xc - domain_.x0) / width;

                    drawList->AddRectFilled(
                        ImVec2(static_cast<float>(x0), static_cast<float>(y0)),
                        ImVec2(static_cast<float>(x1), static_cast<float>(y1)),
                        scalarToColor(phi)
                    );
                }
            }
        }

        ImGui::End();
    }

} // namespace CFD::UI
