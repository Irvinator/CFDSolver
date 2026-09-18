#include "MeshEditor2D.hpp"

#include <algorithm>
#include <cmath>
#include <imgui.h>

namespace {

    ImU32 scalarToColor(double value)
    {
        value = std::clamp(value, 0.0, 1.0);

        double r = 0.0;
        double g = 0.0;
        double b = 0.0;

        // Blue -> Cyan
        if (value < 0.25)
        {
            const double t = value / 0.25;

            r = 0.0;
            g = t;
            b = 1.0;
        }
        // Cyan -> Green
        else if (value < 0.50)
        {
            const double t = (value - 0.25) / 0.25;

            r = 0.0;
            g = 1.0;
            b = 1.0 - t;
        }
        // Green -> Yellow
        else if (value < 0.75)
        {
            const double t = (value - 0.50) / 0.25;

            r = t;
            g = 1.0;
            b = 0.0;
        }
        // Yellow -> Red
        else
        {
            const double t = (value - 0.75) / 0.25;

            r = 1.0;
            g = 1.0 - t;
            b = 0.0;
        }

        return IM_COL32(
            static_cast<int>(255.0 * r),
            static_cast<int>(255.0 * g),
            static_cast<int>(255.0 * b),
            255);
    }

}

namespace CFD::UI {

    namespace {

        void normalizeRect(RectDomain& rect)
        {
            if (rect.x0 > rect.x1)
                std::swap(rect.x0, rect.x1);

            if (rect.y0 > rect.y1)
                std::swap(rect.y0, rect.y1);
        }

        bool nearlyEqual(
            float a,
            float b,
            float eps = 0.5f)
        {
            return std::fabs(a - b) <= eps;
        }

    }

    // ================================================================
    // OBJ MESH
    // ================================================================

    void MeshEditor2D::setMesh(
        CFD::OBJMesh* mesh)
    {
        importedMesh_ = mesh;

        objBBoxDirty_ = true;
        objCentreRequested_ = true;

        autoRefitOnResize_ = true;
    }

    static void computeBBox(
        const CFD::OBJMesh& m,
        float& minX,
        float& maxX,
        float& minY,
        float& maxY)
    {
        minX = minY = 1e30f;
        maxX = maxY = -1e30f;

        for (const auto& v : m.vertices)
        {
            minX = std::min(minX, v.x);
            maxX = std::max(maxX, v.x);

            minY = std::min(minY, v.y);
            maxY = std::max(maxY, v.y);
        }
    }

    // ================================================================
    // VIEW CENTRING
    // ================================================================

    void MeshEditor2D::centreOBJView(
        ImVec2 canvasSize)
    {
        const float w =
            objMaxX_ - objMinX_;

        const float h =
            objMaxY_ - objMinY_;

        if (w <= 0.0f || h <= 0.0f)
            return;

        const float margin = 40.0f;

        const float scaleX =
            std::max(
                0.01f,
                (canvasSize.x - 2.0f * margin) / w);

        const float scaleY =
            std::max(
                0.01f,
                (canvasSize.y - 2.0f * margin) / h);

        objZoom_ =
            std::clamp(
                std::min(scaleX, scaleY),
                0.001f,
                5000.0f);

        const float cx =
            (objMinX_ + objMaxX_) * 0.5f;

        const float cy =
            (objMinY_ + objMaxY_) * 0.5f;

        objPanX_ =
            canvasSize.x * 0.5f -
            cx * objZoom_;

        objPanY_ =
            canvasSize.y * 0.5f +
            cy * objZoom_;
    }

    void MeshEditor2D::centreView(
        ImVec2 canvasSize)
    {
        if (!domainReady_)
            return;

        const float margin = 40.0f;

        const float domW =
            static_cast<float>(domain_.width());

        const float domH =
            static_cast<float>(domain_.height());

        if (domW <= 0.0f || domH <= 0.0f)
            return;

        const float scaleX =
            std::max(
                0.01f,
                (canvasSize.x - 2.0f * margin) / domW);

        const float scaleY =
            std::max(
                0.01f,
                (canvasSize.y - 2.0f * margin) / domH);

        zoom_ =
            std::clamp(
                std::min(scaleX, scaleY),
                0.05f,
                500.0f);

        const float domCx =
            static_cast<float>(
                domain_.x0 + domain_.x1) * 0.5f;

        const float domCy =
            static_cast<float>(
                domain_.y0 + domain_.y1) * 0.5f;

        panX_ =
            canvasSize.x * 0.5f -
            domCx * zoom_;

        // Structured mesh uses mathematical Y:
        // positive Y is UP on screen.
        panY_ =
            canvasSize.y * 0.5f +
            domCy * zoom_;
    }

    // ================================================================
    // HEAT ANIMATION
    // ================================================================

    void MeshEditor2D::pushAnimationFrame(
        const std::vector<double>& values,
        double time,
        double globalMin,
        double globalMax)
    {
        AnimFrame f;

        f.values = values;
        f.time = time;

        f.minValue = globalMin;
        f.maxValue = globalMax;

        f.navierStokes = false;

        animFrames_.push_back(
            std::move(f));

        if (animFrames_.size() == 1)
        {
            animFrame_ = 0;
            animPlaying_ = true;

            applyAnimFrame(0);
        }
    }

    // ================================================================
    // NAVIER-STOKES ANIMATION
    // ================================================================

    void MeshEditor2D::pushNavierStokesAnimationFrame(
        const std::vector<double>& pressure,
        const std::vector<double>& velocityU,
        const std::vector<double>& velocityV,
        double time,
        double pressureMin,
        double pressureMax,
        double velocityUMin,
        double velocityUMax,
        double velocityVMin,
        double velocityVMax)
    {
        AnimFrame f;

        f.pressure = pressure;
        f.velocityU = velocityU;
        f.velocityV = velocityV;

        // ------------------------------------------------------------
        // Calculate cell-centred velocity magnitude
        //
        // |V| = sqrt(U^2 + V^2)
        // ------------------------------------------------------------

        f.velocityMagnitude.resize(
            velocityU.size());

        double velocityMagnitudeMin =
            1.0e30;

        double velocityMagnitudeMax =
            -1.0e30;

        for (std::size_t i = 0;
            i < velocityU.size();
            ++i)
        {
            const double u =
                velocityU[i];

            const double v =
                (i < velocityV.size())
                ? velocityV[i]
                : 0.0;

            const double magnitude =
                std::sqrt(
                    u * u +
                    v * v);

            f.velocityMagnitude[i] =
                magnitude;

            velocityMagnitudeMin =
                std::min(
                    velocityMagnitudeMin,
                    magnitude);

            velocityMagnitudeMax =
                std::max(
                    velocityMagnitudeMax,
                    magnitude);
        }

        if (f.velocityMagnitude.empty())
        {
            velocityMagnitudeMin = 0.0;
            velocityMagnitudeMax = 1.0;
        }
        else if (velocityMagnitudeMax <=
            velocityMagnitudeMin)
        {
            velocityMagnitudeMax =
                velocityMagnitudeMin + 1.0;
        }

        f.velocityMagnitudeMin =
            velocityMagnitudeMin;

        f.velocityMagnitudeMax =
            velocityMagnitudeMax;

        f.time = time;

        f.pressureMin = pressureMin;
        f.pressureMax = pressureMax;

        f.velocityUMin = velocityUMin;
        f.velocityUMax = velocityUMax;

        f.velocityVMin = velocityVMin;
        f.velocityVMax = velocityVMax;

        f.navierStokes = true;

        animFrames_.push_back(
            std::move(f));

        if (animFrames_.size() == 1)
        {
            animFrame_ = 0;
            animPlaying_ = true;

            applyAnimFrame(0);
        }
    }

    // ================================================================
    // CLEAR ANIMATION
    // ================================================================

    void MeshEditor2D::clearAnimation()
    {
        animFrames_.clear();

        animFrame_ = 0;

        animPlaying_ = false;

        animAccum_ = 0.0f;
    }

    // ================================================================
    // APPLY ANIMATION FRAME
    // ================================================================

    void MeshEditor2D::applyAnimFrame(
        int idx)
    {
        if (idx < 0 ||
            idx >= static_cast<int>(animFrames_.size()))
        {
            return;
        }

        const AnimFrame& f =
            animFrames_[idx];

        // ------------------------------------------------------------
        // Heat
        // ------------------------------------------------------------

        if (!f.navierStokes)
        {
            scalar_.values = f.values;

            scalar_.minValue =
                f.minValue;

            scalar_.maxValue =
                f.maxValue;

            scalar_.loaded =
                !f.values.empty();

            return;
        }

        // ------------------------------------------------------------
        // Navier-Stokes
        // ------------------------------------------------------------

        const std::vector<double>* selectedValues =
            nullptr;

        switch (outputField_)
        {
        case OutputField::Pressure:

            selectedValues =
                &f.pressure;

            scalar_.minValue =
                f.pressureMin;

            scalar_.maxValue =
                f.pressureMax;

            break;

        case OutputField::VelocityU:

            selectedValues =
                &f.velocityU;

            scalar_.minValue =
                f.velocityUMin;

            scalar_.maxValue =
                f.velocityUMax;

            break;

        case OutputField::VelocityV:

            selectedValues =
                &f.velocityV;

            scalar_.minValue =
                f.velocityVMin;

            scalar_.maxValue =
                f.velocityVMax;

            break;

        case OutputField::VelocityMagnitude:

            selectedValues =
                &f.velocityMagnitude;

            scalar_.minValue =
                f.velocityMagnitudeMin;

            scalar_.maxValue =
                f.velocityMagnitudeMax;

            break;

        case OutputField::Temperature:

        default:

            selectedValues =
                &f.pressure;

            scalar_.minValue =
                f.pressureMin;

            scalar_.maxValue =
                f.pressureMax;

            break;
        }

        if (selectedValues)
        {
            scalar_.values =
                *selectedValues;

            scalar_.loaded =
                !selectedValues->empty();
        }
        else
        {
            scalar_.values.clear();
            scalar_.loaded = false;
        }
    }

    // ================================================================
    // EDITOR UI
    // ================================================================

    void MeshEditor2D::drawUI()
    {
        if (!showEditorWindow_)
            return;

        ImGui::Begin(
            "Mesh Editor 2D",
            &showEditorWindow_);

        ImGui::Text(
            "Use the combined docked editor in the Viewport window.");

        ImGui::Checkbox(
            "Show Viewport",
            &showViewportWindow_);

        ImGui::End();
    }

    // ================================================================
    // MAIN VIEWPORT
    // ================================================================
    void MeshEditor2D::drawVelocityVectors(
        ImDrawList* drawList,
        ImVec2 canvasPos,
        const std::function<ImVec2(double, double)>& toScreen)
    {
        if (animFrames_.empty())
            return;

        if (animFrame_ < 0 ||
            animFrame_ >=
            static_cast<int>(animFrames_.size()))
            return;

        const AnimFrame& frame =
            animFrames_[animFrame_];

        if (!frame.navierStokes)
            return;

        const int nx = mesh_.nx;
        const int ny = mesh_.ny;

        if (nx <= 0 || ny <= 0)
            return;

        const double dx =
            domain_.width() /
            static_cast<double>(nx);

        const double dy =
            domain_.height() /
            static_cast<double>(ny);

        const double minDimension =
            std::min(dx, dy);

        /*
            Scale the arrows relative to the cell size.

            This is deliberately independent of the velocity
            magnitude range so that the arrows remain visible
            even when the flow contains small velocities.
        */
        const double vectorScale =
            0.35 * minDimension;

        const double maxArrowLength =
            0.45 * minDimension;

        const double minArrowLength =
            0.02 * minDimension;

        for (int j = 0; j < ny; ++j)
        {
            for (int i = 0; i < nx; ++i)
            {
                const std::size_t index =
                    static_cast<std::size_t>(
                        j * nx + i);

                if (index >= frame.velocityU.size() ||
                    index >= frame.velocityV.size())
                {
                    continue;
                }

                const double u =
                    frame.velocityU[index];

                const double v =
                    frame.velocityV[index];

                const double magnitude =
                    std::sqrt(
                        u * u +
                        v * v);

                if (magnitude < 1.0e-12)
                    continue;

                /*
                    Cell centre.
                */
                const double x =
                    domain_.x0 +
                    (static_cast<double>(i) + 0.5) * dx;

                const double y =
                    domain_.y0 +
                    (static_cast<double>(j) + 0.5) * dy;

                double arrowLength =
                    magnitude * vectorScale;

                arrowLength =
                    std::clamp(
                        arrowLength,
                        minArrowLength,
                        maxArrowLength);

                const double invMagnitude =
                    1.0 / magnitude;

                const double dirX =
                    u * invMagnitude;

                const double dirY =
                    v * invMagnitude;

                const double endX =
                    x + dirX * arrowLength;

                const double endY =
                    y + dirY * arrowLength;

                const ImVec2 start =
                    toScreen(x, y);

                const ImVec2 end =
                    toScreen(endX, endY);

                drawList->AddLine(
                    start,
                    end,
                    IM_COL32(
                        255,
                        255,
                        255,
                        230),
                    1.5f);

                /*
                    Arrow head.
                */
                const double headLength =
                    std::min(
                        arrowLength * 0.30,
                        minDimension * 0.12);

                const double headAngle =
                    0.5;

                const double cosA =
                    std::cos(headAngle);

                const double sinA =
                    std::sin(headAngle);

                const double backX =
                    -dirX * cosA -
                    -dirY * sinA;

                const double backY =
                    -dirX * sinA +
                    -dirY * cosA;

                const double backX2 =
                    -dirX * cosA +
                    -dirY * sinA;

                const double backY2 =
                    dirX * sinA +
                    -dirY * cosA;

                const ImVec2 head1 =
                    toScreen(
                        endX + backX * headLength,
                        endY + backY * headLength);

                const ImVec2 head2 =
                    toScreen(
                        endX + backX2 * headLength,
                        endY + backY2 * headLength);

                drawList->AddLine(
                    end,
                    head1,
                    IM_COL32(
                        255,
                        255,
                        255,
                        230),
                    1.5f);

                drawList->AddLine(
                    end,
                    head2,
                    IM_COL32(
                        255,
                        255,
                        255,
                        230),
                    1.5f);
            }
        }
    }
    void MeshEditor2D::drawStreamlines(
        ImDrawList* drawList,
        ImVec2 canvasPos,
        const std::function<ImVec2(double, double)>& toScreen)
    {
        if (animFrames_.empty())
            return;

        if (animFrame_ < 0 ||
            animFrame_ >=
            static_cast<int>(animFrames_.size()))
            return;

        const AnimFrame& frame =
            animFrames_[animFrame_];

        if (!frame.navierStokes)
            return;

        const int nx = mesh_.nx;
        const int ny = mesh_.ny;

        if (nx <= 0 || ny <= 0)
            return;

        if (frame.velocityU.size() <
            static_cast<std::size_t>(nx * ny))
            return;

        if (frame.velocityV.size() <
            static_cast<std::size_t>(nx * ny))
            return;

        const double dx =
            domain_.width() /
            static_cast<double>(nx);

        const double dy =
            domain_.height() /
            static_cast<double>(ny);

        if (dx <= 0.0 || dy <= 0.0)
            return;

        /*
            Bilinear interpolation of the cell-centred velocity field.
        */
        auto sampleVelocity =
            [&](double x,
                double y,
                double& u,
                double& v) -> bool
            {
                if (x < domain_.x0 ||
                    x > domain_.x1 ||
                    y < domain_.y0 ||
                    y > domain_.y1)
                {
                    return false;
                }

                double fx =
                    (x - domain_.x0) / dx - 0.5;

                double fy =
                    (y - domain_.y0) / dy - 0.5;

                int i0 =
                    static_cast<int>(
                        std::floor(fx));

                int j0 =
                    static_cast<int>(
                        std::floor(fy));

                double tx =
                    fx - static_cast<double>(i0);

                double ty =
                    fy - static_cast<double>(j0);

                /*
                    Clamp to the valid cell range.

                    At the boundaries this effectively reduces
                    interpolation to the nearest available cells.
                */
                i0 =
                    std::clamp(
                        i0,
                        0,
                        nx - 1);

                j0 =
                    std::clamp(
                        j0,
                        0,
                        ny - 1);

                const int i1 =
                    std::min(
                        i0 + 1,
                        nx - 1);

                const int j1 =
                    std::min(
                        j0 + 1,
                        ny - 1);

                if (i0 == i1)
                    tx = 0.0;

                if (j0 == j1)
                    ty = 0.0;

                auto velocityAt =
                    [&](int i,
                        int j,
                        double& uu,
                        double& vv)
                    {
                        const std::size_t index =
                            static_cast<std::size_t>(
                                j * nx + i);

                        uu =
                            frame.velocityU[index];

                        vv =
                            frame.velocityV[index];
                    };

                double u00, v00;
                double u10, v10;
                double u01, v01;
                double u11, v11;

                velocityAt(
                    i0,
                    j0,
                    u00,
                    v00);

                velocityAt(
                    i1,
                    j0,
                    u10,
                    v10);

                velocityAt(
                    i0,
                    j1,
                    u01,
                    v01);

                velocityAt(
                    i1,
                    j1,
                    u11,
                    v11);

                const double u0 =
                    u00 * (1.0 - tx) +
                    u10 * tx;

                const double u1 =
                    u01 * (1.0 - tx) +
                    u11 * tx;

                const double v0 =
                    v00 * (1.0 - tx) +
                    v10 * tx;

                const double v1 =
                    v01 * (1.0 - tx) +
                    v11 * tx;

                u =
                    u0 * (1.0 - ty) +
                    u1 * ty;

                v =
                    v0 * (1.0 - ty) +
                    v1 * ty;

                return true;
            };

        /*
            Streamline integration.

            A relatively small step keeps the lines smooth while
            remaining inexpensive for your current 2D solver.
        */
        const double stepSize =
            0.20 *
            std::min(dx, dy);

        const int maxSteps = 500;

        /*
            Seeds are placed throughout the domain.

            More seeds horizontally than vertically gives good
            coverage without making the viewport too cluttered.
        */
        const int seedNX =
            std::clamp(
                nx / 2,
                8,
                20);

        const int seedNY =
            std::clamp(
                ny / 2,
                8,
                20);

        auto integrate =
            [&](double startX,
                double startY,
                double direction)
            {
                double x = startX;
                double y = startY;

                ImVec2 previous =
                    toScreen(x, y);

                for (int step = 0;
                    step < maxSteps;
                    ++step)
                {
                    double u1;
                    double v1;

                    if (!sampleVelocity(
                        x,
                        y,
                        u1,
                        v1))
                    {
                        break;
                    }

                    u1 *= direction;
                    v1 *= direction;

                    const double speed1 =
                        std::sqrt(
                            u1 * u1 +
                            v1 * v1);

                    if (speed1 < 1.0e-10)
                        break;

                    /*
                        RK2 midpoint step.
                    */
                    const double midX =
                        x +
                        0.5 *
                        stepSize *
                        u1 /
                        speed1;

                    const double midY =
                        y +
                        0.5 *
                        stepSize *
                        v1 /
                        speed1;

                    double u2;
                    double v2;

                    if (!sampleVelocity(
                        midX,
                        midY,
                        u2,
                        v2))
                    {
                        break;
                    }

                    u2 *= direction;
                    v2 *= direction;

                    const double speed2 =
                        std::sqrt(
                            u2 * u2 +
                            v2 * v2);

                    if (speed2 < 1.0e-10)
                        break;

                    x +=
                        stepSize *
                        u2 /
                        speed2;

                    y +=
                        stepSize *
                        v2 /
                        speed2;

                    if (x < domain_.x0 ||
                        x > domain_.x1 ||
                        y < domain_.y0 ||
                        y > domain_.y1)
                    {
                        break;
                    }

                    const ImVec2 current =
                        toScreen(x, y);

                    drawList->AddLine(
                        previous,
                        current,
                        IM_COL32(
                            255,
                            220,
                            80,
                            210),
                        1.5f);

                    previous = current;
                }
            };

        /*
            Seed streamlines throughout the domain.

            We integrate both forward and backward from each seed.
            This makes the seed distribution much more useful for
            recirculating flows such as the lid-driven cavity.
        */
        for (int sy = 0;
            sy < seedNY;
            ++sy)
        {
            const double y =
                domain_.y0 +
                (static_cast<double>(sy) + 0.5) /
                static_cast<double>(seedNY) *
                domain_.height();

            for (int sx = 0;
                sx < seedNX;
                ++sx)
            {
                const double x =
                    domain_.x0 +
                    (static_cast<double>(sx) + 0.5) /
                    static_cast<double>(seedNX) *
                    domain_.width();

                integrate(
                    x,
                    y,
                    1.0);

                integrate(
                    x,
                    y,
                    -1.0);
            }
        }
    }
    void MeshEditor2D::drawViewport()
    {
        if (!showViewportWindow_)
            return;

        ImGui::SetNextWindowSize(
            ImVec2(900, 600),
            ImGuiCond_FirstUseEver);

        ImGui::Begin(
            "Mesh Viewport",
            &showViewportWindow_);

        const ImVec2 fullAvail =
            ImGui::GetContentRegionAvail();

        const float panelWidth = 280.0f;
        const float legendWidth = 100.0f;
        const float spacing = 8.0f;

        // ============================================================
        // LEFT CONTROL PANEL
        // ============================================================

        ImGui::BeginChild(
            "mesh_controls",
            ImVec2(panelWidth, fullAvail.y),
            true);

        if (hasMesh())
        {
            ImGui::Text(
                "Imported OBJ Mesh");

            ImGui::Separator();

            ImGui::Text(
                "Vertices : %d",
                static_cast<int>(
                    importedMesh_->vertices.size()));

            ImGui::Text(
                "Faces    : %d",
                static_cast<int>(
                    importedMesh_->faces.size()));

            ImGui::Text(
                "Normals  : %d",
                static_cast<int>(
                    importedMesh_->normals.size()));

            ImGui::Separator();

            ImGui::Text(
                "Bounds X : %.3f -> %.3f",
                objMinX_,
                objMaxX_);

            ImGui::Text(
                "Bounds Y : %.3f -> %.3f",
                objMinY_,
                objMaxY_);

            ImGui::Text(
                "Size     : %.3f x %.3f",
                objMaxX_ - objMinX_,
                objMaxY_ - objMinY_);

            if (ImGui::Button("Centre Mesh"))
            {
                objCentreRequested_ = true;
                autoRefitOnResize_ = true;
            }

            ImGui::SameLine();

            if (ImGui::Button("Clear Mesh"))
            {
                importedMesh_ = nullptr;

                objBBoxDirty_ = true;
                objCentreRequested_ = true;

                autoRefitOnResize_ = true;
            }

            ImGui::TextDisabled(
                "Scroll to zoom");

            ImGui::TextDisabled(
                "MMB / RMB drag to pan");
        }
        else
        {
            if (ImGui::Button("Draw Rectangle"))
            {
                drawMode_ = true;
                drawing_ = false;
                domainReady_ = false;

                clearScalarField();
                clearAnimation();

                autoRefitOnResize_ = true;
            }

            if (domainReady_)
            {
                ImGui::SameLine();

                if (ImGui::Button("Centre View"))
                {
                    centreViewRequested_ = true;
                    autoRefitOnResize_ = true;
                }
            }

            ImGui::Separator();

            ImGui::Text(
                "Mesh Settings");

            ImGui::InputInt(
                "Nx",
                &mesh_.nx);

            ImGui::InputInt(
                "Ny",
                &mesh_.ny);

            ImGui::InputDouble(
                "Width",
                &mesh_.width,
                0.1,
                1.0,
                "%.3f");

            ImGui::InputDouble(
                "Height",
                &mesh_.height,
                0.1,
                1.0,
                "%.3f");

            if (mesh_.nx < 1)
                mesh_.nx = 1;

            if (mesh_.ny < 1)
                mesh_.ny = 1;

            if (mesh_.width <= 0.0)
                mesh_.width = 1.0;

            if (mesh_.height <= 0.0)
                mesh_.height = 1.0;

            ImGui::Checkbox(
                "Show Grid",
                &showGrid_);

            ImGui::Checkbox(
                "Show Colormap",
                &showColormap_);

            // ------------------------------------------------------------
            // Navier-Stokes overlays
            // ------------------------------------------------------------

            if (!hasMesh() && !animFrames_.empty())
            {
                const bool hasNavierStokes =
                    animFrame_ >= 0 &&
                    animFrame_ <
                    static_cast<int>(animFrames_.size()) &&
                    animFrames_[animFrame_].navierStokes;

                if (hasNavierStokes)
                {
                    ImGui::Separator();

                    ImGui::Text("Flow Visualisation");

                    ImGui::Checkbox(
                        "Show Velocity Vectors",
                        &showVelocityVectors_);

                    ImGui::Checkbox(
                        "Show Streamlines",
                        &showStreamlines_);
                }
            }

            ImGui::Separator();

            ImGui::Text(
                "Rectangle Domain");

            ImGui::Text(
                "x0 = %.1f, y0 = %.1f",
                domain_.x0,
                domain_.y0);

            ImGui::Text(
                "x1 = %.1f, y1 = %.1f",
                domain_.x1,
                domain_.y1);

            ImGui::Text(
                "Width  = %.1f",
                domain_.width());

            ImGui::Text(
                "Height = %.1f",
                domain_.height());

            if (scalar_.loaded)
            {
                ImGui::Separator();

                ImGui::Text(
                    "Field Loaded");

                ImGui::Text(
                    "Min = %.6f",
                    scalar_.minValue);

                ImGui::Text(
                    "Max = %.6f",
                    scalar_.maxValue);
            }

            if (!animFrames_.empty())
            {
                ImGui::Separator();

                ImGui::Text(
                    "Animation (%d frames)",
                    static_cast<int>(
                        animFrames_.size()));

                if (animPlaying_)
                {
                    if (ImGui::Button("|| Pause"))
                        animPlaying_ = false;
                }
                else
                {
                    if (ImGui::Button("|> Play"))
                        animPlaying_ = true;
                }

                ImGui::SameLine();

                if (ImGui::Button("|< Reset"))
                {
                    animFrame_ = 0;
                    animPlaying_ = false;
                    animAccum_ = 0.0f;

                    applyAnimFrame(0);
                }

                int fi = animFrame_;

                ImGui::SetNextItemWidth(-1.0f);

                if (ImGui::SliderInt(
                    "##anim_inner",
                    &fi,
                    0,
                    static_cast<int>(
                        animFrames_.size()) - 1))
                {
                    animFrame_ = fi;

                    animPlaying_ = false;
                    animAccum_ = 0.0f;

                    applyAnimFrame(
                        animFrame_);
                }

                if (!animFrames_.empty())
                {
                    ImGui::Text(
                        "t = %.0f",
                        animFrames_[animFrame_].time);
                }

                ImGui::SliderFloat(
                    "Speed",
                    &animSpeed_,
                    0.5f,
                    60.0f,
                    "%.1f fps");
            }

            if (domainReady_)
            {
                ImGui::TextColored(
                    ImVec4(
                        0.2f,
                        0.8f,
                        0.2f,
                        1.0f),
                    "Domain ready");
            }
            else if (drawMode_)
            {
                ImGui::TextColored(
                    ImVec4(
                        0.9f,
                        0.8f,
                        0.2f,
                        1.0f),
                    "Click and drag in viewport");
            }
            else
            {
                ImGui::Text(
                    "No domain yet");
            }
        }

        ImGui::EndChild();

        ImGui::SameLine(
            0.0f,
            spacing);

        // ============================================================
        // CANVAS
        // ============================================================

        const float canvasWidth =
            std::max(
                100.0f,
                fullAvail.x -
                panelWidth -
                legendWidth -
                2.0f * spacing);

        ImGui::BeginChild(
            "mesh_canvas_region",
            ImVec2(canvasWidth, fullAvail.y),
            false,
            ImGuiWindowFlags_NoScrollbar);

        const ImVec2 canvasPos =
            ImGui::GetCursorScreenPos();

        const ImVec2 canvasSize =
            ImGui::GetContentRegionAvail();

        const ImVec2 canvasEnd(
            canvasPos.x + canvasSize.x,
            canvasPos.y + canvasSize.y);

        const bool canvasChanged =
            !nearlyEqual(
                canvasSize.x,
                lastCanvasSize_.x) ||
            !nearlyEqual(
                canvasSize.y,
                lastCanvasSize_.y);

        if (canvasChanged &&
            canvasSize.x > 10.0f &&
            canvasSize.y > 10.0f)
        {
            if (autoRefitOnResize_)
            {
                if (hasMesh())
                    objCentreRequested_ = true;
                else if (domainReady_)
                    centreViewRequested_ = true;
            }

            lastCanvasSize_ = canvasSize;
        }

        ImDrawList* drawList =
            ImGui::GetWindowDrawList();

        drawList->AddRectFilled(
            canvasPos,
            canvasEnd,
            IM_COL32(25, 25, 30, 255));

        drawList->AddRect(
            canvasPos,
            canvasEnd,
            IM_COL32(200, 200, 200, 255));

        ImGui::InvisibleButton(
            "mesh_canvas",
            canvasSize,
            ImGuiButtonFlags_MouseButtonLeft |
            ImGuiButtonFlags_MouseButtonRight |
            ImGuiButtonFlags_MouseButtonMiddle);

        const bool hovered =
            ImGui::IsItemHovered();

        const ImVec2 mouse =
            ImGui::GetIO().MousePos;

        // ============================================================
        // OBJ VIEW
        // ============================================================

        if (hasMesh())
        {
            if (objBBoxDirty_)
            {
                computeBBox(
                    *importedMesh_,
                    objMinX_,
                    objMaxX_,
                    objMinY_,
                    objMaxY_);

                objBBoxDirty_ = false;
            }

            if (objCentreRequested_)
            {
                centreOBJView(canvasSize);
                objCentreRequested_ = false;
            }

            if (hovered)
            {
                const float wheel =
                    ImGui::GetIO().MouseWheel;

                if (wheel != 0.0f)
                {
                    const float oldZoom =
                        objZoom_;

                    objZoom_ =
                        std::clamp(
                            objZoom_ *
                            (1.0f + 0.12f * wheel),
                            0.001f,
                            5000.0f);

                    const ImVec2 ml(
                        mouse.x - canvasPos.x,
                        mouse.y - canvasPos.y);

                    objPanX_ =
                        ml.x -
                        (ml.x - objPanX_) *
                        (objZoom_ / oldZoom);

                    objPanY_ =
                        ml.y -
                        (ml.y - objPanY_) *
                        (objZoom_ / oldZoom);

                    autoRefitOnResize_ = false;
                }

                if (ImGui::IsMouseDragging(
                    ImGuiMouseButton_Middle) ||
                    ImGui::IsMouseDragging(
                        ImGuiMouseButton_Right))
                {
                    const ImVec2 d =
                        ImGui::GetIO().MouseDelta;

                    objPanX_ += d.x;
                    objPanY_ += d.y;

                    autoRefitOnResize_ = false;
                }
            }

            auto objToScreen =
                [&](float x, float y) -> ImVec2
                {
                    return ImVec2(
                        canvasPos.x +
                        objPanX_ +
                        x * objZoom_,

                        canvasPos.y +
                        objPanY_ -
                        y * objZoom_);
                };

            for (const auto& face :
                importedMesh_->faces)
            {
                const std::size_t nv =
                    face.verticesFace.size();

                if (nv < 2)
                    continue;

                for (std::size_t i = 0;
                    i < nv;
                    ++i)
                {
                    const std::size_t j =
                        (i + 1) % nv;

                    const int ia =
                        face.verticesFace[i]
                        .vertexIndex;

                    const int ib =
                        face.verticesFace[j]
                        .vertexIndex;

                    if (ia < 0 ||
                        ib < 0 ||
                        ia >= static_cast<int>(
                            importedMesh_->vertices.size()) ||
                        ib >= static_cast<int>(
                            importedMesh_->vertices.size()))
                    {
                        continue;
                    }

                    const auto& a =
                        importedMesh_->vertices[ia];

                    const auto& b =
                        importedMesh_->vertices[ib];

                    drawList->AddLine(
                        objToScreen(a.x, a.y),
                        objToScreen(b.x, b.y),
                        IM_COL32(
                            80,
                            200,
                            255,
                            220),
                        1.0f);
                }
            }
        }

        // ============================================================
        // STRUCTURED DOMAIN
        // ============================================================

        else
        {
            if (centreViewRequested_)
            {
                centreView(canvasSize);
                centreViewRequested_ = false;
            }

            // --------------------------------------------------------
            // Animation
            // --------------------------------------------------------

            if (animPlaying_ &&
                !animFrames_.empty())
            {
                animAccum_ +=
                    ImGui::GetIO().DeltaTime;

                const float frameDur =
                    1.0f / animSpeed_;

                while (animAccum_ >= frameDur)
                {
                    animAccum_ -= frameDur;

                    const int last =
                        static_cast<int>(
                            animFrames_.size()) - 1;

                    if (animFrame_ >= last)
                    {
                        if (animLoop_)
                        {
                            animFrame_ = 0;
                        }
                        else
                        {
                            animFrame_ = last;
                            animPlaying_ = false;
                            animAccum_ = 0.0f;
                            break;
                        }
                    }
                    else
                    {
                        ++animFrame_;
                    }

                    applyAnimFrame(
                        animFrame_);
                }
            }

            // --------------------------------------------------------
            // Mouse
            // --------------------------------------------------------

            if (hovered)
            {
                const float wheel =
                    ImGui::GetIO().MouseWheel;

                if (wheel != 0.0f)
                {
                    const float oldZoom =
                        zoom_;

                    zoom_ =
                        std::clamp(
                            zoom_ *
                            (1.0f + 0.1f * wheel),
                            0.05f,
                            500.0f);

                    const ImVec2 ml(
                        mouse.x - canvasPos.x,
                        mouse.y - canvasPos.y);

                    panX_ =
                        ml.x -
                        (ml.x - panX_) *
                        (zoom_ / oldZoom);

                    panY_ =
                        ml.y -
                        (ml.y - panY_) *
                        (zoom_ / oldZoom);

                    autoRefitOnResize_ = false;
                }

                if (ImGui::IsMouseDragging(
                    ImGuiMouseButton_Middle) ||
                    ImGui::IsMouseDragging(
                        ImGuiMouseButton_Right))
                {
                    const ImVec2 d =
                        ImGui::GetIO().MouseDelta;

                    panX_ += d.x;
                    panY_ += d.y;

                    autoRefitOnResize_ = false;
                }
            }

            // --------------------------------------------------------
            // IMPORTANT:
            // Structured domain uses mathematical coordinates:
            //
            // +X -> right
            // +Y -> up
            //
            // ImGui screen coordinates have +Y downward, so Y
            // must be inverted here.
            // --------------------------------------------------------

            auto toScreen =
                [&](double x, double y) -> ImVec2
                {
                    return ImVec2(
                        canvasPos.x +
                        panX_ +
                        static_cast<float>(x) * zoom_,

                        canvasPos.y +
                        panY_ -
                        static_cast<float>(y) * zoom_);
                };

            // Inverse transform for mouse -> world coordinates.
            auto toWorld =
                [&](float sx, float sy) -> ImVec2
                {
                    return ImVec2(
                        (sx -
                            canvasPos.x -
                            panX_) / zoom_,

                        (canvasPos.y +
                            panY_ -
                            sy) / zoom_);
                };

            // --------------------------------------------------------
            // Draw domain
            // --------------------------------------------------------

            if (drawMode_ &&
                hovered &&
                ImGui::IsMouseClicked(
                    ImGuiMouseButton_Left))
            {
                const ImVec2 w =
                    toWorld(
                        mouse.x,
                        mouse.y);

                domain_.x0 =
                    static_cast<double>(w.x);

                domain_.y0 =
                    static_cast<double>(w.y);

                domain_.x1 =
                    domain_.x0;

                domain_.y1 =
                    domain_.y0;

                drawing_ = true;
                domainReady_ = false;

                clearScalarField();
                clearAnimation();

                autoRefitOnResize_ = true;
            }

            if (drawMode_ &&
                drawing_)
            {
                const ImVec2 w =
                    toWorld(
                        mouse.x,
                        mouse.y);

                domain_.x1 =
                    static_cast<double>(w.x);

                domain_.y1 =
                    static_cast<double>(w.y);

                if (ImGui::IsMouseReleased(
                    ImGuiMouseButton_Left))
                {
                    normalizeRect(domain_);

                    drawing_ = false;

                    domainReady_ =
                        domain_.valid();

                    drawMode_ = false;

                    centreView(
                        canvasSize);

                    lastCanvasSize_ =
                        canvasSize;
                }
            }

            // --------------------------------------------------------
            // Colour map
            // --------------------------------------------------------

            if (domainReady_ &&
                showColormap_ &&
                scalar_.loaded)
            {
                const double dw =
                    domain_.width();

                const double dh =
                    domain_.height();

                const double dx =
                    dw /
                    static_cast<double>(
                        mesh_.nx);

                const double dy =
                    dh /
                    static_cast<double>(
                        mesh_.ny);

                const double rng =
                    scalar_.maxValue -
                    scalar_.minValue;

                for (int j = 0;
                    j < mesh_.ny;
                    ++j)
                {
                    for (int i = 0;
                        i < mesh_.nx;
                        ++i)
                    {
                        const double x0 =
                            domain_.x0 +
                            static_cast<double>(i) *
                            dx;

                        const double y0 =
                            domain_.y0 +
                            static_cast<double>(j) *
                            dy;

                        const std::size_t P =
                            static_cast<std::size_t>(
                                j * mesh_.nx + i);

                        double phi = 0.0;

                        if (P <
                            scalar_.values.size())
                        {
                            const double v =
                                scalar_.values[P];

                            phi =
                                (rng > 0.0)
                                ? (v -
                                    scalar_.minValue) /
                                rng
                                : 0.0;
                        }

                        drawList->AddRectFilled(
                            toScreen(x0, y0),
                            toScreen(
                                x0 + dx,
                                y0 + dy),
                            scalarToColor(phi));
                    }
                }
            }

            // --------------------------------------------------------
            // Domain outline
            // --------------------------------------------------------

            if (drawing_ ||
                domainReady_)
            {
                drawList->AddRect(
                    toScreen(
                        domain_.x0,
                        domain_.y0),

                    toScreen(
                        domain_.x1,
                        domain_.y1),

                    IM_COL32(
                        80,
                        180,
                        255,
                        255),

                    0.0f,
                    0,
                    2.0f);
            }

            // --------------------------------------------------------
            // Grid
            // --------------------------------------------------------

            if (domainReady_ &&
                showGrid_)
            {
                const double dw =
                    domain_.width();

                const double dh =
                    domain_.height();

                const double dx =
                    dw /
                    static_cast<double>(
                        mesh_.nx);

                const double dy =
                    dh /
                    static_cast<double>(
                        mesh_.ny);

                for (int i = 0;
                    i <= mesh_.nx;
                    ++i)
                {
                    const double x =
                        domain_.x0 +
                        static_cast<double>(i) *
                        dx;

                    drawList->AddLine(
                        toScreen(
                            x,
                            domain_.y0),

                        toScreen(
                            x,
                            domain_.y1),

                        IM_COL32(
                            120,
                            120,
                            120,
                            255));
                }

                for (int j = 0;
                    j <= mesh_.ny;
                    ++j)
                {
                    const double y =
                        domain_.y0 +
                        static_cast<double>(j) *
                        dy;

                    drawList->AddLine(
                        toScreen(
                            domain_.x0,
                            y),

                        toScreen(
                            domain_.x1,
                            y),

                        IM_COL32(
                            120,
                            120,
                            120,
                            255));
                }
            }
            // --------------------------------------------------------
            // Navier-Stokes velocity overlays
            // --------------------------------------------------------

            if (domainReady_ &&
                !animFrames_.empty() &&
                animFrame_ >= 0 &&
                animFrame_ <
                static_cast<int>(animFrames_.size()) &&
                animFrames_[animFrame_].navierStokes)
            {
                if (showStreamlines_)
                {
                    drawStreamlines(
                        drawList,
                        canvasPos,
                        toScreen);
                }

                if (showVelocityVectors_)
                {
                    drawVelocityVectors(
                        drawList,
                        canvasPos,
                        toScreen);
                }
            }
            // --------------------------------------------------------
            // Hover information
            // --------------------------------------------------------

            if (domainReady_ &&
                hovered &&
                scalar_.loaded)
            {
                const ImVec2 w =
                    toWorld(
                        mouse.x,
                        mouse.y);

                if (w.x >= domain_.x0 &&
                    w.x <= domain_.x1 &&
                    w.y >= domain_.y0 &&
                    w.y <= domain_.y1)
                {
                    const double dx =
                        domain_.width() /
                        static_cast<double>(
                            mesh_.nx);

                    const double dy =
                        domain_.height() /
                        static_cast<double>(
                            mesh_.ny);

                    int ci =
                        std::clamp(
                            static_cast<int>(
                                (w.x -
                                    domain_.x0) /
                                dx),
                            0,
                            mesh_.nx - 1);

                    int cj =
                        std::clamp(
                            static_cast<int>(
                                (w.y -
                                    domain_.y0) /
                                dy),
                            0,
                            mesh_.ny - 1);

                    const std::size_t P =
                        static_cast<std::size_t>(
                            cj * mesh_.nx + ci);

                    if (P <
                        scalar_.values.size())
                    {
                        const double value =
                            scalar_.values[P];

                        switch (outputField_)
                        {
                        case OutputField::Pressure:

                            ImGui::SetTooltip(
                                "Cell (%d, %d)\nP = %.6f Pa",
                                ci,
                                cj,
                                value);

                            break;

                        case OutputField::VelocityU:

                            ImGui::SetTooltip(
                                "Cell (%d, %d)\nU = %.6f m/s",
                                ci,
                                cj,
                                value);

                            break;

                        case OutputField::VelocityV:

                            ImGui::SetTooltip(
                                "Cell (%d, %d)\nV = %.6f m/s",
                                ci,
                                cj,
                                value);

                            break;

                        case OutputField::VelocityMagnitude:

                            ImGui::SetTooltip(
                                "Cell (%d, %d)\n|V| = %.6f m/s",
                                ci,
                                cj,
                                value);

                            break;

                        case OutputField::Temperature:

                        default:

                            ImGui::SetTooltip(
                                "Cell (%d, %d)\nT = %.6f",
                                ci,
                                cj,
                                value);

                            break;
                        }
                    }
                }
            }
        }

        ImGui::EndChild();

        ImGui::SameLine(
            0.0f,
            spacing);

        // ============================================================
        // LEGEND
        // ============================================================

        ImGui::BeginChild(
            "mesh_legend",
            ImVec2(
                legendWidth,
                fullAvail.y),
            true);

        ImGui::Text("Legend");

        ImGui::Separator();

        if (scalar_.loaded &&
            !hasMesh())
        {
            ImDrawList* ld =
                ImGui::GetWindowDrawList();

            const ImVec2 barMin =
                ImGui::GetCursorScreenPos();

            const float barW = 22.0f;

            const float barH =
                std::max(
                    120.0f,
                    fullAvail.y - 110.0f);

            const ImVec2 barMax(
                barMin.x + barW,
                barMin.y + barH);

            for (int y = 0;
                y < static_cast<int>(barH);
                ++y)
            {
                const double t =
                    1.0 -
                    static_cast<double>(y) /
                    static_cast<double>(
                        barH - 1);

                ld->AddLine(
                    ImVec2(
                        barMin.x,
                        barMin.y +
                        static_cast<float>(y)),

                    ImVec2(
                        barMin.x + barW,
                        barMin.y +
                        static_cast<float>(y)),

                    scalarToColor(t));
            }

            ld->AddRect(
                barMin,
                barMax,
                IM_COL32(
                    255,
                    255,
                    255,
                    200));

            ld->AddLine(
                ImVec2(
                    barMax.x,
                    barMin.y),

                ImVec2(
                    barMax.x + 3.0f,
                    barMin.y),

                IM_COL32(
                    255,
                    255,
                    255,
                    200));

            ld->AddLine(
                ImVec2(
                    barMax.x,
                    barMin.y +
                    barH * 0.5f),

                ImVec2(
                    barMax.x + 3.0f,
                    barMin.y +
                    barH * 0.5f),

                IM_COL32(
                    255,
                    255,
                    255,
                    200));

            ld->AddLine(
                ImVec2(
                    barMax.x,
                    barMax.y),

                ImVec2(
                    barMax.x + 3.0f,
                    barMax.y),

                IM_COL32(
                    255,
                    255,
                    255,
                    200));

            ImGui::Dummy(
                ImVec2(
                    barW,
                    barH));

            const float lineH =
                ImGui::GetTextLineHeight();

            ImGui::SetCursorScreenPos(
                ImVec2(
                    barMax.x + 4.0f,
                    barMin.y));

            ImGui::Text(
                "%.3f",
                scalar_.maxValue);

            const double mid =
                scalar_.minValue +
                (scalar_.maxValue -
                    scalar_.minValue) *
                0.5;

            ImGui::SetCursorScreenPos(
                ImVec2(
                    barMax.x + 4.0f,
                    barMin.y +
                    barH * 0.5f -
                    lineH * 0.5f));

            ImGui::Text(
                "%.3f",
                mid);

            ImGui::SetCursorScreenPos(
                ImVec2(
                    barMax.x + 4.0f,
                    barMax.y - lineH));

            ImGui::Text(
                "%.3f",
                scalar_.minValue);
        }
        else if (hasMesh())
        {
            ImGui::TextDisabled(
                "Wireframe\nview");

            ImGui::Spacing();

            ImGui::TextColored(
                ImVec4(
                    0.3f,
                    0.78f,
                    1.0f,
                    1.0f),
                "Edges");

            ImGui::TextColored(
                ImVec4(
                    1.0f,
                    0.7f,
                    0.2f,
                    1.0f),
                "Vertices");
        }
        else
        {
            ImGui::TextDisabled(
                "Run solver\nto load\nfield");
        }

        ImGui::EndChild();

        ImGui::End();
    }

    // ================================================================
    // STATIC SCALAR FIELD
    // ================================================================

    void MeshEditor2D::setScalarField(
        const std::vector<double>& values,
        int nx,
        int ny)
    {
        scalar_.values = values;

        scalar_.loaded =
            !values.empty();

        mesh_.nx =
            std::max(1, nx);

        mesh_.ny =
            std::max(1, ny);

        if (scalar_.loaded)
        {
            auto [mn, mx] =
                std::minmax_element(
                    scalar_.values.begin(),
                    scalar_.values.end());

            scalar_.minValue = *mn;
            scalar_.maxValue = *mx;
        }
        else
        {
            scalar_.minValue = 0.0;
            scalar_.maxValue = 1.0;
        }
    }

    void MeshEditor2D::setScalarField(
        const CFD::ScalarField& field,
        int nx,
        int ny)
    {
        scalar_.values.resize(
            field.size());

        for (std::size_t i = 0;
            i < field.size();
            ++i)
        {
            scalar_.values[i] =
                field[i];
        }

        scalar_.minValue =
            field.min();

        scalar_.maxValue =
            field.max();

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