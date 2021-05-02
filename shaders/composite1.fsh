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

/* DRAWBUFFERS:7 */

void main() {
    surfaceStruct surface = GetSurfaceStruct();
	maskStruct mask = GetMaskStruct(texcoord.st);

    vec3 wavesNormal = texture2D(colortex7, texcoord.st).xyz;
	vec3 wavesNormalr = texture2D(colortex7, mod(texcoord.st + vec2(0.5, 0.0), vec2(1.0))).xyz;
	vec3 wavesNormalu = texture2D(colortex7, mod(texcoord.st + vec2(0.0, 0.5), vec2(1.0))).xyz;
	vec3 wavesNormalur = texture2D(colortex7, mod(texcoord.st + vec2(0.5, 0.5), vec2(1.0))).xyz;

	float lerpx = saturate((abs(texcoord.x - 0.5) * 2.0) * 3.0 - 2.0);
	float lerpy = saturate((abs(texcoord.y - 0.5) * 2.0) * 3.0 - 2.0);

	vec3 x0 = mix(wavesNormal, wavesNormalr, vec3(lerpx));
	vec3 x1 = mix(wavesNormalu, wavesNormalur, vec3(lerpx));
	vec3 seamlessWavesNormal = normalize(mix(x0, x1, vec3(lerpy)));

    gl_FragData[0] = encodeRGBE8(seamlessWavesNormal);
}