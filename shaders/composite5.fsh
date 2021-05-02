#version 130

varying vec4 texcoord;

varying vec3 sunVector;
varying vec3 wSunVector;
varying vec3 moonVector;
varying vec3 wMoonVector;
varying vec3 upVector;
varying vec3 lightVector;
varying vec3 wLightVector;

varying vec3 baseSunColor;
varying vec3 sunColor;
varying vec3 baseMoonColor;
varying vec3 moonColor;
varying vec3 skyColor;

varying float transitionFading;

#include "libs/uniform.glsl"
#include "libs/composite.glsl"
#include "libs/noise.glsl"

#include "/libs/sky.glsl"

/* DRAWBUFFERS:3 */

void main() {
    surfaceStruct surface = GetSurfaceStruct();
	maskStruct mask = GetMaskStruct(texcoord.st);
    vec3 color = texture2D(colortex0, texcoord.st).rgb;
    gl_FragData[0] = vec4(color, 1.0f);
}