#pragma once
#include <vector>
#include <array>
#include "mesh/mesh2D.h"


class MeshReader2D
{
public:
    struct Vertex
    {
        float x;
        float y;
        float z;
    };

    struct Face
    {
        std::array<std::size_t, 4> vertexIndices;
        std::size_t cellIndex;
    };

    void readMesh(const Mesh& inputMesh) {

        const int nx = inputMesh.getNx();
        const int ny = inputMesh.getNy();

        const float dx = static_cast<float>(inputMesh.getDx());
        const float dy = static_cast<float>(inputMesh.getDy());
//Vertices
        for (int j = 0; j <= ny; ++j)
        {
            for (int i = 0; i <= nx; ++i)
            {
                vertices.push_back({
                    i * dx,
                    j * dy,
                    0.0f
                    });
            }
        }
//Faces (Quads)
        for (int j = 0; j < ny; ++j)
        {
            for (int i = 0; i < nx; ++i)
            {
                const std::size_t bottomLeft =
                    j * (nx + 1) + i;

                const std::size_t bottomRight =
                    bottomLeft + 1;

                const std::size_t topLeft =
                    (j + 1) * (nx + 1) + i;

                const std::size_t topRight =
                    topLeft + 1;

                Face face;

                face.vertexIndices = {
                    bottomLeft,
                    bottomRight,
                    topRight,
                    topLeft
                };

                face.cellIndex =
                    static_cast<std::size_t>(
                        inputMesh.cellIndex(i, j)
                        );

                faces.push_back(face);
            }
        }
    }

    const std::vector<Vertex>& getVertices() const
    {
        return vertices;
    }

    const std::vector<Face>& getFaces() const
    {
        return faces;
    }

private:
    std::vector<Vertex> vertices;
    std::vector<Face> faces;
};