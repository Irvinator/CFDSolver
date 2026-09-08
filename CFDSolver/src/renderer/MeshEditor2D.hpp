#pragma once

#include "fields/ScalarField.hpp"   // CFD::ScalarField
#include <imgui.h>
#include <vector>
#include <string>

namespace CFD::UI {

    // ── Domain rectangle ──────────────────────────────────────────────────────
    struct RectDomain {
        double x0 = 0.0, y0 = 0.0;
        double x1 = 0.0, y1 = 0.0;

        double width()  const { return x1 - x0; }
        double height() const { return y1 - y0; }
        bool   valid()  const { return width() > 0.0 && height() > 0.0; }
    };

    // ── Mesh settings ─────────────────────────────────────────────────────────
    struct MeshSettings {
        int    nx = 20;
        int    ny = 20;
        double width = 1.0;
        double height = 1.0;
    };

    // ─────────────────────────────────────────────────────────────────────────
    class MeshEditor2D {
    public:
        // ── UI entry points ───────────────────────────────────────────────
        void drawUI();
        void drawViewport();

        // ── Visibility ────────────────────────────────────────────────────
        void showEditorWindow(bool v) { showEditorWindow_ = v; }
        void showViewportWindow(bool v) { showViewportWindow_ = v; }
        bool isEditorWindowVisible()  const { return showEditorWindow_; }
        bool isViewportWindowVisible() const { return showViewportWindow_; }

        // ── Scalar field (single frame / final result) ────────────────────
        void setScalarField(const std::vector<double>& values, int nx, int ny);
        void setScalarField(const CFD::ScalarField& field, int nx, int ny);
        void clearScalarField();
        bool hasScalarField() const { return scalar_.loaded; }

        // ── Animation API ─────────────────────────────────────────────────
        /// Call once per solver step to record a frame for playback.
        void pushAnimationFrame(const std::vector<double>& values, double time);
        void clearAnimation();

        // ── Accessors ─────────────────────────────────────────────────────
        const RectDomain& domain()       const { return domain_; }
        bool               domainReady()  const { return domainReady_; }
        MeshSettings       meshSettings() const { return mesh_; }

    private:
        // ── Animation frame storage ───────────────────────────────────────
        struct AnimFrame {
            std::vector<double> values;
            double              time = 0.0;
            double              minValue = 0.0;
            double              maxValue = 1.0;
        };

        void applyAnimFrame(int idx);
        void centreView(ImVec2 canvasSize);

        // ── State ─────────────────────────────────────────────────────────
        bool showEditorWindow_ = true;
        bool showViewportWindow_ = true;

        // Drawing
        bool       drawMode_ = false;
        bool       drawing_ = false;
        bool       domainReady_ = false;
        RectDomain domain_;
        MeshSettings mesh_;

        // Viewport transform
        float panX_ = 0.0f, panY_ = 0.0f, zoom_ = 1.0f;
        bool  centreViewRequested_ = false;

        // Display flags
        bool showGrid_ = true;
        bool showColormap_ = true;

        // Live scalar field (current display frame)
        struct ScalarData {
            std::vector<double> values;
            double minValue = 0.0;
            double maxValue = 1.0;
            bool   loaded = false;
        } scalar_;

        // Animation
        std::vector<AnimFrame> animFrames_;
        int   animFrame_ = 0;
        bool  animPlaying_ = false;
        float animSpeed_ = 10.0f;   // frames per second
        float animAccum_ = 0.0f;
    };

} // namespace CFD::UI
