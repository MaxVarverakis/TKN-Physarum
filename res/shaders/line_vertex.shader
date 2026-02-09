#version 410 core

layout (location = 0) in vec2 aPos;

out vec4 v_Color;

uniform mat4 u_MVP;
uniform vec4 u_Color;

void main()
{
    gl_Position = u_MVP * vec4(aPos, 0.0, 1.0);
    v_Color = u_Color;
}
