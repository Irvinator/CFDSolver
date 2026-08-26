#version 330 core

in vec3 fragmentNormal;

uniform vec3 lightDir;

out vec4 fragmentColor;

void main() {
    vec3 n = normalize(fragmentNormal);
    vec3 l = normalize(lightDir);
    float intensity = dot(n, l);

    vec3 ambient =  0.15 * vec3(1.0, 1.0, 1.0);
    vec3 color = intensity * vec3(1.0, 1.0, 1.0);

    fragmentColor = vec4(ambient + color, 1.0);
}