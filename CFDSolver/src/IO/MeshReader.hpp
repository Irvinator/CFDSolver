#pragma once
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <filesystem>

namespace CFD {

    struct ReadVertex
    {
        float x;
        float y;
        float z;

        ReadVertex(float x, float y, float z)
            :x(x), y(y), z(z) {
        }

        ReadVertex(const ReadVertex& vertex)
            :x(vertex.x), y(vertex.y), z(vertex.z) {
        }
       
    };

    struct ReadNormal
    {
        float x;
        float y;
        float z;
        ReadNormal(float x, float y, float z)
            :x(x), y(y), z(z) {
        }

        ReadNormal(const ReadNormal& normal)
            :x(normal.x), y(normal.y), z(normal.z) {
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

    OBJMesh loadOBJ(const std::string& filepath){
        OBJMesh mesh;
        std::cout << "mesh created!!!";

        std::ifstream file(filepath);
        if (std::filesystem::exists(filepath)) {
            std::cout << "File exists!\n";
        }
        else {
            std::cout << "File NOT found!\n";
            std::cout << "Current dir: "
                << std::filesystem::current_path()
                << "\n";
        }
        std::string line;

        while (std::getline(file, line)) {
            std::istringstream stream(line);
            std::string type;
            stream >> type;

            if (type.empty())
                continue;

            else if (type == "v") {
                float x, y, z;
                stream >> x >> y >> z;
                mesh.vertices.emplace_back(x, y, z);
            }

            else if (type == "vn") {
                float x, y, z;
                stream >> x >> y >> z;

                mesh.normals.emplace_back(x, y, z);
            }

            else if (type == "f") {
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

                    face.verticesFace.push_back({
                        vi,
                        vni
                        });
                }
                mesh.faces.push_back(face);
            }
        }
        return mesh;
    }

} // namespace CFD