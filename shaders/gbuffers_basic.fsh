#version 130

varying vec4 texcoord;

varying vec4 color;

#include "libs/uniform.glsl"
#include "libs/composite.glsl"
#include "libs/noise.glsl"

/* DRAWBUFFERS:0 */

void main() {
    gl_FragData[0] = color;
}