#pragma once
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace CDF {

    struct Vertex
    {
        float x;
        float y;
        float z;
    };

    struct Normal
    {
        float x;
        float y;
        float z;
    };

    struct FaceVertex
    {
        int vertexIndex;
        int normalIndex;
    };

    struct Face
    {
        std::vector<FaceVertex> vertices;
    };

    struct OBJMesh
    {
        std::vector<Vertex> vertices;
        std::vector<Normal> normals;
        std::vector<Face> faces;
    };

    OBJMesh loadONJ(const std::string& filepath){
        OBJMesh mesh;

        std << ifstream file(filepath);

        std::string line;

        while (std::getline(file, line)) {
            std::istringstream stream(line);

            std::string type;
            stream >> type;

            if (type.empty())
                continue;

            if (type == "v") {
                float x, y, z;
                stream >> x >> y >> z;
                mesh.vertices.emplace_back(x, y, z);
            }

            if (type == "vn") {
                float x, y, z;
                stream >> x >> y >> z;

                mesh.normals.emplace_back(x, y, z);
            }

            if (type == "f") {
                Face face;
                std::string token;

                while (stream >> token) {
                    std::stringstream tokenStream(token);

                    std::string viString;
                    std::string textureIndexString;
                    std::string vniString;

                    std::getline(tokenStream, viString, '/');
                    std::getline(tokenStream, textureIndexString, '/');
                    std::getline(tokenStream, vniString, '/');

                    int vi =
                        std::stoi(viString) - 1;

                    int vni = -1;

                    if (!vniString.empty())
                    {
                        vni =
                            std::stoi(vniString) - 1;
                    }

                    face.vertices.push_back({
                        vi,
                        vni
                        });
                }
                mesh.faces.push_back(face);
            }
        }
        return mesh;
    }

} // namespace CDF
