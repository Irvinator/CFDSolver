#pragma once

#include <string>
#include <vector>

namespace CFD {

    struct ReadVertex
    {
        float x;
        float y;
        float z;

        ReadVertex(float x, float y, float z)
            : x(x), y(y), z(z) {
        }

        ReadVertex(const ReadVertex& vertex)
            : x(vertex.x), y(vertex.y), z(vertex.z) {
        }
    };

    struct ReadNormal
    {
        float x;
        float y;
        float z;

        ReadNormal(float x, float y, float z)
            : x(x), y(y), z(z) {
        }

        ReadNormal(const ReadNormal& normal)
            : x(normal.x), y(normal.y), z(normal.z) {
        }
    };

    struct FaceVertex
    {
        int vertexIndex;
        int normalIndex;
    };

    struct Face
    {
        std::vector<FaceVertex> verticesFace;
    };

    struct OBJMesh
    {
        std::vector<ReadVertex> vertices;
        std::vector<ReadNormal> normals;
        std::vector<Face> faces;
    };

    OBJMesh loadOBJ(const std::string& filepath);

} // namespace CFD