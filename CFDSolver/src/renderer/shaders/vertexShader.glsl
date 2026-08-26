#version 330 core

layout(location = 0) in vec3 position;
layout(location = 1) in vec3 normal;

uniform mat4 MVP;
uniform mat3 normalMatrix;

out vec3 fragmentNormal;

void main() {
    fragmentNormal = normalize(normalMatrix * normal);
    gl_Position = MVP * vec4(position, 1.0);
}