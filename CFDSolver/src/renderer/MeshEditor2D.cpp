#include "MeshEditor2D.hpp"

#include <algorithm>
#include <cmath>
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
} // namespace

namespace CFD::UI {

    namespace {
        void normalizeRect(RectDomain& rect)
        {
            if (rect.x0 > rect.x1) std::swap(rect.x0, rect.x1);
            if (rect.y0 > rect.y1) std::swap(rect.y0, rect.y1);
        }
    } // namespace

    // ─────────────────────────────────────────────────────────────────────────
    //  Animation helpers
    // ─────────────────────────────────────────────────────────────────────────

    void MeshEditor2D::pushAnimationFrame(const std::vector<double>& values,
        double                      time)
    {
        AnimFrame f;
        f.values = values;
        f.time = time;
        if (!values.empty()) {
            auto [mn, mx] = std::minmax_element(values.begin(), values.end());
            f.minValue = *mn;
            f.maxValue = *mx;
        }
        animFrames_.push_back(std::move(f));
    }

    void MeshEditor2D::clearAnimation()
    {
        animFrames_.clear();
        animFrame_ = 0;
        animPlaying_ = false;
    }

    // ─────────────────────────────────────────────────────────────────────────
    //  Centre-view helper
    //  Fits the domain rectangle into the canvas with a small margin.
    // ─────────────────────────────────────────────────────────────────────────
    void MeshEditor2D::centreView(ImVec2 canvasSize)
    {
        if (!domainReady_) return;

        const float margin = 40.0f;
        const float domW = static_cast<float>(domain_.width());
        const float domH = static_cast<float>(domain_.height());

        if (domW <= 0.0f || domH <= 0.0f) return;

        const float scaleX = (canvasSize.x - 2.0f * margin) / domW;
        const float scaleY = (canvasSize.y - 2.0f * margin) / domH;
        zoom_ = std::min(scaleX, scaleY);
        zoom_ = std::clamp(zoom_, 0.2f, 50.0f);

        // Pan so that the domain centre maps to the canvas centre
        const float domCx = static_cast<float>(domain_.x0 + domain_.x1) * 0.5f;
        const float domCy = static_cast<float>(domain_.y0 + domain_.y1) * 0.5f;
        panX_ = canvasSize.x * 0.5f - domCx * zoom_;
        panY_ = canvasSize.y * 0.5f - domCy * zoom_;
    }

    // ─────────────────────────────────────────────────────────────────────────
    //  drawUI  (small floating helper window)
    // ─────────────────────────────────────────────────────────────────────────
    void MeshEditor2D::drawUI()
    {
        if (!showEditorWindow_) return;

        ImGui::Begin("Mesh Editor 2D", &showEditorWindow_);
        ImGui::Text("Use the combined docked editor in the Viewport window.");
        ImGui::Checkbox("Show Viewport", &showViewportWindow_);
        ImGui::End();
    }

    // ─────────────────────────────────────────────────────────────────────────
    //  drawViewport  – main panel
    // ─────────────────────────────────────────────────────────────────────────
    void MeshEditor2D::drawViewport()
    {
        if (!showViewportWindow_) return;

        ImGui::Begin("Mesh Viewport", &showViewportWindow_);

        // ── Measure total available area ──────────────────────────────────
        const ImVec2 fullAvail = ImGui::GetContentRegionAvail();

        const float panelWidth = 280.0f;
        const float legendWidth = 90.0f;   // bar(22) + gap(4) + label text
        const float spacing = 8.0f;

        // ── Left control panel ────────────────────────────────────────────
        ImGui::BeginChild("mesh_controls",
            ImVec2(panelWidth, fullAvail.y), true);

        if (ImGui::Button("Draw Rectangle")) {
            drawMode_ = true;
            drawing_ = false;
            domainReady_ = false;
            clearScalarField();
            clearAnimation();
        }

        // Centre-view button (only useful when domain exists)
        if (domainReady_) {
            ImGui::SameLine();
            if (ImGui::Button("Centre View")) {
                // We don't know canvasSize here; set a flag and apply next frame
                centreViewRequested_ = true;
            }
        }

        ImGui::Separator();
        ImGui::Text("Mesh Settings");
        ImGui::InputInt("Nx", &mesh_.nx);
        ImGui::InputInt("Ny", &mesh_.ny);
        ImGui::InputDouble("Width", &mesh_.width, 0.1, 1.0, "%.3f");
        ImGui::InputDouble("Height", &mesh_.height, 0.1, 1.0, "%.3f");

        if (mesh_.nx < 1)     mesh_.nx = 1;
        if (mesh_.ny < 1)     mesh_.ny = 1;
        if (mesh_.width <= 0.0) mesh_.width = 1.0;
        if (mesh_.height <= 0.0) mesh_.height = 1.0;

        ImGui::Checkbox("Show Grid", &showGrid_);
        ImGui::Checkbox("Show Colormap", &showColormap_);

        ImGui::Separator();
        ImGui::Text("Rectangle Domain");
        ImGui::Text("x0 = %.1f, y0 = %.1f", domain_.x0, domain_.y0);
        ImGui::Text("x1 = %.1f, y1 = %.1f", domain_.x1, domain_.y1);
        ImGui::Text("Width  = %.1f", domain_.width());
        ImGui::Text("Height = %.1f", domain_.height());

        if (scalar_.loaded) {
            ImGui::Separator();
            ImGui::Text("Scalar Field Loaded");
            ImGui::Text("Min = %.6f", scalar_.minValue);
            ImGui::Text("Max = %.6f", scalar_.maxValue);
        }

        // ── Animation controls ────────────────────────────────────────────
        if (!animFrames_.empty()) {
            ImGui::Separator();
            ImGui::Text("Animation  (%d frames)", static_cast<int>(animFrames_.size()));

            if (animPlaying_) {
                if (ImGui::Button("|| Pause")) animPlaying_ = false;
            }
            else {
                if (ImGui::Button("|> Play"))  animPlaying_ = true;
            }
            ImGui::SameLine();
            if (ImGui::Button("|< Reset")) {
                animFrame_ = 0;
                animPlaying_ = false;
                applyAnimFrame(animFrame_);
            }

            int fi = animFrame_;
            ImGui::SetNextItemWidth(-1.0f);
            if (ImGui::SliderInt("##anim", &fi, 0,
                static_cast<int>(animFrames_.size()) - 1))
            {
                animFrame_ = fi;
                animPlaying_ = false;
                applyAnimFrame(animFrame_);
            }

            ImGui::Text("t = %.4f", animFrames_[animFrame_].time);

            ImGui::SliderFloat("Speed", &animSpeed_, 0.5f, 30.0f, "%.1f fps");
        }

        if (domainReady_)
            ImGui::TextColored(ImVec4(0.2f, 0.8f, 0.2f, 1.0f), "Domain ready");
        else if (drawMode_)
            ImGui::TextColored(ImVec4(0.9f, 0.8f, 0.2f, 1.0f), "Click and drag in viewport");
        else
            ImGui::Text("No domain yet");

        ImGui::EndChild();
        ImGui::SameLine(0.0f, spacing);

        // ── Centre canvas ─────────────────────────────────────────────────
        const float canvasWidth = fullAvail.x
            - panelWidth - spacing
            - legendWidth - spacing;

        ImGui::BeginChild("mesh_canvas_region",
            ImVec2(canvasWidth, fullAvail.y), false,
            ImGuiWindowFlags_NoScrollbar);

        const ImVec2 canvasPos = ImGui::GetCursorScreenPos();
        const ImVec2 canvasSize = ImGui::GetContentRegionAvail();
        const ImVec2 canvasEnd(canvasPos.x + canvasSize.x,
            canvasPos.y + canvasSize.y);

        // Apply deferred centre-view now that we know canvasSize
        if (centreViewRequested_) {
            centreView(canvasSize);
            centreViewRequested_ = false;
        }

        // ── Advance animation ─────────────────────────────────────────────
        if (animPlaying_ && !animFrames_.empty()) {
            animAccum_ += ImGui::GetIO().DeltaTime;
            const float frameDuration = 1.0f / animSpeed_;
            if (animAccum_ >= frameDuration) {
                animAccum_ -= frameDuration;
                animFrame_ = (animFrame_ + 1) %
                    static_cast<int>(animFrames_.size());
                applyAnimFrame(animFrame_);
            }
        }

        ImDrawList* drawList = ImGui::GetWindowDrawList();
        drawList->AddRectFilled(canvasPos, canvasEnd, IM_COL32(25, 25, 30, 255));
        drawList->AddRect(canvasPos, canvasEnd, IM_COL32(200, 200, 200, 255));

        ImGui::InvisibleButton("mesh_canvas", canvasSize,
            ImGuiButtonFlags_MouseButtonLeft |
            ImGuiButtonFlags_MouseButtonRight |
            ImGuiButtonFlags_MouseButtonMiddle);

        const bool   hovered = ImGui::IsItemHovered();
        const ImVec2 mouse = ImGui::GetIO().MousePos;

        // ── Zoom / pan ────────────────────────────────────────────────────
        if (hovered) {
            const float wheel = ImGui::GetIO().MouseWheel;
            if (wheel != 0.0f) {
                const float oldZoom = zoom_;
                zoom_ = std::clamp(zoom_ * (1.0f + 0.1f * wheel), 0.05f, 50.0f);
                const ImVec2 ml(mouse.x - canvasPos.x, mouse.y - canvasPos.y);
                panX_ = ml.x - (ml.x - panX_) * (zoom_ / oldZoom);
                panY_ = ml.y - (ml.y - panY_) * (zoom_ / oldZoom);
            }
            if (ImGui::IsMouseDragging(ImGuiMouseButton_Middle) ||
                ImGui::IsMouseDragging(ImGuiMouseButton_Right))
            {
                const ImVec2 delta = ImGui::GetIO().MouseDelta;
                panX_ += delta.x;
                panY_ += delta.y;
            }
        }

        // ── Coordinate helpers ────────────────────────────────────────────
        auto toScreen = [&](double x, double y) -> ImVec2 {
            return ImVec2(
                canvasPos.x + panX_ + static_cast<float>(x) * zoom_,
                canvasPos.y + panY_ + static_cast<float>(y) * zoom_);
            };

        auto toWorld = [&](float sx, float sy) -> ImVec2 {
            return ImVec2(
                (sx - canvasPos.x - panX_) / zoom_,
                (sy - canvasPos.y - panY_) / zoom_);
            };

        // ── Rectangle drawing interaction ─────────────────────────────────
        if (drawMode_ && hovered && ImGui::IsMouseClicked(ImGuiMouseButton_Left)) {
            const ImVec2 world = toWorld(mouse.x, mouse.y);
            domain_.x0 = static_cast<double>(world.x);
            domain_.y0 = static_cast<double>(world.y);
            domain_.x1 = domain_.x0;
            domain_.y1 = domain_.y0;
            drawing_ = true;
            domainReady_ = false;
            clearScalarField();
            clearAnimation();
        }

        if (drawMode_ && drawing_) {
            const ImVec2 world = toWorld(mouse.x, mouse.y);
            domain_.x1 = static_cast<double>(world.x);
            domain_.y1 = static_cast<double>(world.y);

            if (ImGui::IsMouseReleased(ImGuiMouseButton_Left)) {
                normalizeRect(domain_);
                drawing_ = false;
                domainReady_ = domain_.valid();
                drawMode_ = false;

                // Auto-centre after drawing so the rect is always visible
                centreView(canvasSize);
            }
        }

        // ── Colormap fill ─────────────────────────────────────────────────
        if (domainReady_ && showColormap_ && scalar_.loaded) {
            const double dw = domain_.width();
            const double dh = domain_.height();
            const double dx = dw / static_cast<double>(mesh_.nx);
            const double dy = dh / static_cast<double>(mesh_.ny);
            const double rng = scalar_.maxValue - scalar_.minValue;

            for (int j = 0; j < mesh_.ny; ++j) {
                for (int i = 0; i < mesh_.nx; ++i) {
                    const double x0 = domain_.x0 + static_cast<double>(i) * dx;
                    const double y0 = domain_.y0 + static_cast<double>(j) * dy;

                    const std::size_t P =
                        static_cast<std::size_t>(j * mesh_.nx + i);
                    double phi = 0.0;
                    if (P < scalar_.values.size()) {
                        const double v = scalar_.values[P];
                        phi = (rng > 0.0)
                            ? (v - scalar_.minValue) / rng
                            : 0.0;
                    }

                    drawList->AddRectFilled(
                        toScreen(x0, y0),
                        toScreen(x0 + dx, y0 + dy),
                        scalarToColor(phi));
                }
            }
        }

        // ── Domain outline ────────────────────────────────────────────────
        if (drawing_ || domainReady_) {
            drawList->AddRect(
                toScreen(domain_.x0, domain_.y0),
                toScreen(domain_.x1, domain_.y1),
                IM_COL32(80, 180, 255, 255), 0.0f, 0, 2.0f);
        }

        // ── Grid ──────────────────────────────────────────────────────────
        if (domainReady_ && showGrid_) {
            const double dw = domain_.width();
            const double dh = domain_.height();
            const double dx = dw / static_cast<double>(mesh_.nx);
            const double dy = dh / static_cast<double>(mesh_.ny);

            for (int i = 0; i <= mesh_.nx; ++i) {
                const double x = domain_.x0 + static_cast<double>(i) * dx;
                drawList->AddLine(
                    toScreen(x, domain_.y0),
                    toScreen(x, domain_.y1),
                    IM_COL32(120, 120, 120, 255));
            }
            for (int j = 0; j <= mesh_.ny; ++j) {
                const double y = domain_.y0 + static_cast<double>(j) * dy;
                drawList->AddLine(
                    toScreen(domain_.x0, y),
                    toScreen(domain_.x1, y),
                    IM_COL32(120, 120, 120, 255));
            }
        }

        // ── Hover tooltip ─────────────────────────────────────────────────
        if (domainReady_ && hovered && scalar_.loaded) {
            const ImVec2 world = toWorld(mouse.x, mouse.y);
            if (world.x >= domain_.x0 && world.x <= domain_.x1 &&
                world.y >= domain_.y0 && world.y <= domain_.y1)
            {
                const double dx = domain_.width() / static_cast<double>(mesh_.nx);
                const double dy = domain_.height() / static_cast<double>(mesh_.ny);
                int ci = std::clamp(static_cast<int>((world.x - domain_.x0) / dx),
                    0, mesh_.nx - 1);
                int cj = std::clamp(static_cast<int>((world.y - domain_.y0) / dy),
                    0, mesh_.ny - 1);
                const std::size_t P =
                    static_cast<std::size_t>(cj * mesh_.nx + ci);
                if (P < scalar_.values.size())
                    ImGui::SetTooltip("T = %.6f", scalar_.values[P]);
            }
        }

        ImGui::EndChild();   // mesh_canvas_region
        ImGui::SameLine(0.0f, spacing);

        // ── Legend strip ──────────────────────────────────────────────────
        ImGui::BeginChild("mesh_legend",
            ImVec2(legendWidth, fullAvail.y), true);

        ImGui::Text("Legend");
        ImGui::Separator();

        if (scalar_.loaded) {
            ImDrawList* ld = ImGui::GetWindowDrawList();
            const ImVec2 barMin = ImGui::GetCursorScreenPos();
            const float  barW = 22.0f;
            const float  barH = std::max(120.0f, fullAvail.y - 110.0f);
            const ImVec2 barMax(barMin.x + barW, barMin.y + barH);

            // Draw gradient bar
            for (int y = 0; y < static_cast<int>(barH); ++y) {
                const double t = 1.0 - static_cast<double>(y) /
                    static_cast<double>(barH - 1);
                ld->AddLine(
                    ImVec2(barMin.x, barMin.y + static_cast<float>(y)),
                    ImVec2(barMin.x + barW, barMin.y + static_cast<float>(y)),
                    scalarToColor(t));
            }
            ld->AddRect(barMin, barMax, IM_COL32(255, 255, 255, 200));

            // Advance cursor past the bar
            ImGui::Dummy(ImVec2(barW, barH));

            // Max label — pinned to top-right of bar
            const float lineH = ImGui::GetTextLineHeight();
            ImGui::SetCursorScreenPos(
                ImVec2(barMax.x + 4.0f, barMin.y));
            ImGui::Text("%.3f", scalar_.maxValue);

            // Mid label
            ImGui::SetCursorScreenPos(
                ImVec2(barMax.x + 4.0f, barMin.y + barH * 0.5f - lineH * 0.5f));
            const double mid = scalar_.minValue +
                (scalar_.maxValue - scalar_.minValue) * 0.5;
            ImGui::Text("%.3f", mid);

            // Min label — pinned to bottom-right of bar
            ImGui::SetCursorScreenPos(
                ImVec2(barMax.x + 4.0f, barMax.y - lineH));
            ImGui::Text("%.3f", scalar_.minValue);

            // Tick marks
            ld->AddLine(ImVec2(barMax.x, barMin.y + barH * 0.0f),
                ImVec2(barMax.x + 3.0f, barMin.y + barH * 0.0f),
                IM_COL32(255, 255, 255, 200));
            ld->AddLine(ImVec2(barMax.x, barMin.y + barH * 0.5f),
                ImVec2(barMax.x + 3.0f, barMin.y + barH * 0.5f),
                IM_COL32(255, 255, 255, 200));
            ld->AddLine(ImVec2(barMax.x, barMax.y),
                ImVec2(barMax.x + 3.0f, barMax.y),
                IM_COL32(255, 255, 255, 200));
        }
        else {
            ImGui::TextDisabled("Run solver\nto load\nscalar field");
        }

        ImGui::EndChild();   // mesh_legend

        ImGui::End();        // Mesh Viewport
    }

    // ─────────────────────────────────────────────────────────────────────────
    //  applyAnimFrame  – copy frame data into the live scalar_ struct
    // ─────────────────────────────────────────────────────────────────────────
    void MeshEditor2D::applyAnimFrame(int idx)
    {
        if (idx < 0 || idx >= static_cast<int>(animFrames_.size())) return;
        const AnimFrame& f = animFrames_[idx];
        scalar_.values = f.values;
        scalar_.minValue = f.minValue;
        scalar_.maxValue = f.maxValue;
        scalar_.loaded = true;
    }

    // ─────────────────────────────────────────────────────────────────────────
    //  setScalarField overloads
    // ─────────────────────────────────────────────────────────────────────────
    void MeshEditor2D::setScalarField(const std::vector<double>& values,
        int nx, int ny)
    {
        scalar_.values = values;
        scalar_.loaded = !values.empty();
        mesh_.nx = std::max(1, nx);
        mesh_.ny = std::max(1, ny);

        if (scalar_.loaded) {
            auto [mn, mx] = std::minmax_element(
                scalar_.values.begin(), scalar_.values.end());
            scalar_.minValue = *mn;
            scalar_.maxValue = *mx;
        }
        else {
            scalar_.minValue = 0.0;
            scalar_.maxValue = 1.0;
        }
    }

    void MeshEditor2D::setScalarField(const CFD::ScalarField& field,
        int nx, int ny)
    {
        scalar_.values.resize(field.size());
        for (std::size_t i = 0; i < field.size(); ++i)
            scalar_.values[i] = field[i];

        scalar_.minValue = field.min();
        scalar_.maxValue = field.max();
        scalar_.loaded = true;
        mesh_.nx = nx;
        mesh_.ny = ny;
    }

    void MeshEditor2D::clearScalarField()
    {
        scalar_.values.clear();
        scalar_.minValue = 0.0;
        scalar_.maxValue = 1.0;
        scalar_.loaded = false;
    }

} // namespace CFD::UI
