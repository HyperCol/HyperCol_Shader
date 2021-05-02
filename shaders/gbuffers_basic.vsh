#version 120

varying vec4 texcoord;

varying vec4 color;

uniform vec2 texelSize;
uniform int frameMod8;
uniform int frameMod16;
uniform vec2 taaOffset16;

void main() {
    gl_Position = ftransform();
    texcoord = gl_MultiTexCoord0;

    color = gl_Color;
}