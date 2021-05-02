#version 120
#pragma optimize (on)

// GPU Shader 5
#ifdef MC_GL_ARB_gpu_shader5
#extension GL_ARB_gpu_shader5 : require
#else
#define fma(a,b,c) ((a)*(b)+c)
#endif

#define plus(m, n) ((m + n) - m * n)
#define select(x, def) float(x == def)
#define Cselect(x, edge0, edge1) float(x == clamp(x, edge0, edge1))

float pow2(float a) { return a*a; }
float pow3(float a) { return (a*a)*a; }

vec2 pow2(vec2 a) { return a*a; }
vec2 pow3(vec2 a) { return (a*a)*a; }

vec3 pow2(vec3 a) { return a*a; }
vec3 pow3(vec3 a) { return (a*a)*a; }

vec4 pow2(vec4 a) { return a*a; }
vec4 pow3(vec4 a) { return (a*a)*a; }

//#define TONE_DEBUG

varying vec2 tex;
#ifdef TONE_DEBUG
vec2 texcoord = vec2(0.25 + fract(tex.x * 2.0) * 0.5, tex.y);
#else
vec2 texcoord = tex;
#endif

#define ANIMATION 0 	//[0 1 2 3]
// 0 off || 1 HyperCol Logo || 2 eyes open || 3 simple animation

//#define VIBRANCE    0.0 //[-1.0 -0.9 -0.8 -0.7 -0.6 -0.5 -0.4 -0.3 -0.2 -0.1 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0]
//#define SATURATION  0.0 //[-1.0 -0.9 -0.8 -0.7 -0.6 -0.5 -0.4 -0.3 -0.2 -0.1 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0]
/*
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

varying float avgSkyBrightness;*/

const float PI = 3.14159265359;

#include "libs/uniform.glsl"
//#include "libs/composite.glsl"
//#include "libs/noise.glsl"
#include "/libs/Tone.glsl"
#if ANIMATION > 0
#include "/libs/Animation.glsl.frag"
#endif

/*vec3 vibranceSaturation(vec3 color){
	const float amountVibrance = VIBRANCE;
	const float amountSaturation = SATURATION;

	float lum = dot(color, vec3(0.2125, 0.7154, 0.0721));
	float mn = min3(color);
	float mx = max3(color);
	float sat = (1.0 - (mx - mn)) * (1.0 - mx) * lum * 5.0;
	vec3 lig = vec3((mn + mx) * 0.5);

	// Vibrance
	color = mix(color, mix(color, lig, -amountVibrance), sat);

	// Inverse Vibrance
	color = mix(color, lig, (1.0 - lig) * (1.0 - amountVibrance) * 0.5 * abs(amountVibrance));

	// saturation
	color = mix(color, vec3(lum), -amountSaturation);

	return color;
}*/

uniform float centerDepthSmooth;

#define FringeOffset 0.5 //[0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0]

vec2 dispersion(in vec2 uv, float e) {
	e = 1.0 - e * 0.04;
	return (e * uv + vec2(1.0 - e) * 0.5);
}

void simple_DOF(inout Tone t) {
	float depth = (isEyeInWater > 0) ? texture2D(depthtex1, texcoord.st).x : texture2D(depthtex0, texcoord.st).x;
	
	float naive = 0.0;
	
	naive += abs(depth - clamp(centerDepthSmooth, 0.2, 0.995)) * 12;

	naive += pow(dot(texcoord.st - vec2(0.5), texcoord.st - vec2(0.5)), 1.5) * 1.5;
	
	t.blurIndex = naive;
	//t.color = vec3(naive);
	
	t.blur.r = texture2D(gaux4, dispersion(texcoord, FringeOffset)).r;
	t.blur.g = texture2D(gaux4, texcoord).g;
	t.blur.b = texture2D(gaux4, dispersion(texcoord, -FringeOffset)).b;
}

Tone tone;

void main() {
	// build up tone
	init_tone(tone, texcoord);
    //vec3 color = texture2D(colortex0, texcoord.st).rgb;
	
	#if ANIMATION > 0
	animation(tone, texcoord);
	#endif
	
	//simple_DOF(tone);
	
	//TONED
	Hue_Adjustment(tone);

    gl_FragColor = vec4(tone.color, 1.0f);
}