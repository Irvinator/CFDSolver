#pragma once

#include "solvers/HeatSolver2D.h"
#include "IO/MeshReader.hpp"

#include <imgui.h>
#include <vector>
#include <string>
#include <algorithm>

namespace CFD::UI {

    struct RectDomain {
        double x0 = 0.0, y0 = 0.0;
        double x1 = 0.0, y1 = 0.0;

        double width()  const { return x1 - x0; }
        double height() const { return y1 - y0; }
        bool   valid()  const { return width() > 0.0 && height() > 0.0; }
    };

    struct MeshSettings {
        int    nx = 20;
        int    ny = 20;
        double width = 1.0;
        double height = 1.0;
    };

    class MeshEditor2D {
    public:
        void drawUI();
        void drawViewport();

        void showEditorWindow(bool v) { showEditorWindow_ = v; }
        void showViewportWindow(bool v) { showViewportWindow_ = v; }
        bool isEditorWindowVisible()   const { return showEditorWindow_; }
        bool isViewportWindowVisible() const { return showViewportWindow_; }

        void setScalarField(const std::vector<double>& values, int nx, int ny);
        void setScalarField(const CFD::ScalarField& field, int nx, int ny);
        void clearScalarField();
        bool hasScalarField() const { return scalar_.loaded; }

        void pushAnimationFrame(const std::vector<double>& values,
            double time,
            double globalMin,
            double globalMax);
        void clearAnimation();

        int animFrameCount() const { return static_cast<int>(animFrames_.size()); }

        double animTime(int idx) const {
            if (animFrames_.empty()) return 0.0;
            idx = std::clamp(idx, 0, static_cast<int>(animFrames_.size()) - 1);
            return animFrames_[idx].time;
        }

        double animEndTime() const {
            return animFrames_.empty() ? 0.0 : animFrames_.back().time;
        }

        void setAnimFrame(int idx) {
            if (animFrames_.empty()) return;
            animFrame_ = std::clamp(idx, 0,
                static_cast<int>(animFrames_.size()) - 1);
            applyAnimFrame(animFrame_);
        }

        int  currentAnimFrame() const { return animFrame_; }
        bool animPlaying()      const { return animPlaying_; }

        void setAnimSpeed(float fps) {
            animSpeed_ = std::clamp(fps, 0.5f, 60.0f);
        }

        void setMesh(CFD::OBJMesh* mesh);
        bool hasMesh() const {
            return importedMesh_ != nullptr &&
                !importedMesh_->vertices.empty();
        }

        const RectDomain& domain()       const { return domain_; }
        bool               domainReady()  const { return domainReady_; }
        MeshSettings       meshSettings() const { return mesh_; }

    private:
        struct AnimFrame {
            std::vector<double> values;
            double              time = 0.0;
            double              minValue = 0.0;
            double              maxValue = 1.0;
        };

        void applyAnimFrame(int idx);
        void centreView(ImVec2 canvasSize);
        void centreOBJView(ImVec2 canvasSize);

        bool showEditorWindow_ = true;
        bool showViewportWindow_ = true;

        bool         drawMode_ = false;
        bool         drawing_ = false;
        bool         domainReady_ = false;
        RectDomain   domain_;
        MeshSettings mesh_;

        float panX_ = 0.0f, panY_ = 0.0f, zoom_ = 1.0f;
        bool  centreViewRequested_ = false;
        ImVec2 lastCanvasSize_ = ImVec2(0.0f, 0.0f);
        bool autoRefitOnResize_ = true;

        bool showGrid_ = true;
        bool showColormap_ = true;

        struct ScalarData {
            std::vector<double> values;
            double minValue = 0.0;
            double maxValue = 1.0;
            bool   loaded = false;
        } scalar_;

        std::vector<AnimFrame> animFrames_;
        int   animFrame_ = 0;
        bool  animPlaying_ = false;
        bool  animLoop_ = true;
        float animSpeed_ = 10.0f;
        float animAccum_ = 0.0f;

        CFD::OBJMesh* importedMesh_ = nullptr;

        float objPanX_ = 0.0f;
        float objPanY_ = 0.0f;
        float objZoom_ = 1.0f;
        bool  objCentreRequested_ = true;

        float objMinX_ = 0.0f, objMaxX_ = 1.0f;
        float objMinY_ = 0.0f, objMaxY_ = 1.0f;
        bool  objBBoxDirty_ = true;
    };

}