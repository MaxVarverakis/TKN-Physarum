#version 410 core

out vec4 FragColor;
uniform vec4 uColor;  // fade color, e.g., (0.0, 0.0, 0.0, 0.05)

void main() {
    FragColor = uColor;
}
