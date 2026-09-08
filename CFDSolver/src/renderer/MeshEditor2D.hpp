#pragma once

#include <vector>
#include "IO/MeshReader.hpp"

namespace CFD::UI {

    struct RectDomain
    {
        double x0 = 0.0;
        double y0 = 0.0;
        double x1 = 1.0;
        double y1 = 1.0;

        bool valid() const { return x1 > x0 && y1 > y0; }
        double width() const { return x1 - x0; }
        double height() const { return y1 - y0; }
    };

    struct MeshSettings {
        int nx = 20;
        int ny = 20;
    };

    struct ScalarDisplayField {
        std::vector<double> values;
        double minValue = 0.0;
        double maxValue = 1.0;
        bool loaded = false;
    };

    class MeshEditor2D {
    public:

        void setMesh(CFD::OBJMesh* mesh);
        void drawUI();
        void drawViewport();

    private:
        bool drawMode_ = false;
        bool drawing_ = false;
        bool domainReady_ = false;
        bool showGrid_ = true;
        bool showColormap_ = true;

        RectDomain domain_;
        MeshSettings mesh_;
        ScalarDisplayField scalar_;

        CFD::OBJMesh* importedMesh_ = nullptr;
    };
}