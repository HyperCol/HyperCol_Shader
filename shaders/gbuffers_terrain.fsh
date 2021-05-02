#version 130
#extension GL_ARB_shader_texture_lod : enable

// Half float support
#ifdef MC_GL_AMD_shader_half_float
#extension GL_AMD_shader_half_float : require
#else

	#define float16_t float
	#define f16vec2 vec2
	#define f16vec3 vec3
	#define f16vec4 vec4
	#define f16mat2 mat2
	#define f16mat3 mat3
	#define f16mat4 mat4
	#define HF f

#endif

#define plus(m, n) ((m + n) - m * n)
#define select(x, def) float(x == def)
#define Cselect(x, edge0, edge1) float(x == clamp(x, edge0, edge1))

#define POM
#define NORMALS

#define TEXTURE_RESOLUTION 1024      //[8 16 32 64 128 256 512 1024 2048 4096]
#define PARALLAX_DEPTH 1.0          //[0.5 0.75 1.0 1.25 1.5 1.75 2.0 2.5 3.0 4.0 5.0]

varying vec4 texcoord;
varying vec4 lmcoord;

varying vec4 vtexcoord;
varying vec4 vtexcoordam; // .st for add, .pq for mul.

varying vec4 color;
varying vec3 normal;
varying vec3 tangent;
varying vec3 binormal;
varying vec3 viewVector;

varying float terrainMats;
varying float distances;;
varying vec3 wpos;
varying float top;

varying float minecraftMats;

#include "libs/uniform.glsl"
#include "libs/composite.glsl"
#include "libs/noise.glsl"

uniform int terrainIconSize;
uniform vec3 BiomeType;

float rain0 = rainStrength * smoothstep(0.0, 0.5, BiomeType.y);

#ifndef NORMALS
varying float n2;
#endif

const float normalMapAngle = 2.0; //[0.5 1.0 1.5 2.0 3.0]

vec2 dcdx = dFdx(vtexcoord.st * vtexcoordam.pq);
vec2 dcdy = dFdy(vtexcoord.st * vtexcoordam.pq);

float absoluteTexGrad = dot(abs(dcdx) + abs(dcdy), vec2(1.0));
	
vec2 OffsetCoord(in vec2 coord, in vec2 offset, in int level) {

	int tileResolution    = terrainIconSize;
	ivec2 atlasTiles      = textureSize(texture, 0) / TEXTURE_RESOLUTION;
	ivec2 atlasResolution = tileResolution * atlasTiles;

	coord *= atlasResolution;

	vec2 offsetCoord = coord + mod(offset.xy * atlasResolution, vec2(tileResolution));

	vec2 minCoord = vec2(coord.x - mod(coord.x, tileResolution), coord.y - mod(coord.y, tileResolution));
	vec2 maxCoord = minCoord + tileResolution;

	if (offsetCoord.x > maxCoord.x) {
		offsetCoord.x -= tileResolution;
	} else if (offsetCoord.x < minCoord.x) {
		offsetCoord.x += tileResolution;
	}

	if (offsetCoord.y > maxCoord.y) {
		offsetCoord.y -= tileResolution;
	} else if (offsetCoord.y < minCoord.y) {
		offsetCoord.y += tileResolution;
	}

	offsetCoord /= atlasResolution;
	return offsetCoord;
	
}

vec2 CalculateParallaxCoord(in vec2 coord, in vec3 viewVector) {

	const int maxSteps = 112;
	const float gradThreshold = 0.004;
	const float parallaxStepSize = 0.5;
	
	vec2 parallaxCoord = coord.st;
	vec3 stepSize 	   = vec3(0.001f, 0.001f, 0.15f);
	
	float parallaxDepth  = PARALLAX_DEPTH;
		  parallaxDepth *= clamp(1.0 - clamp(absoluteTexGrad / gradThreshold, 0.0f, 1.0f), 0.0, 1.0);
	
	if (absoluteTexGrad > gradThreshold) return texcoord.st;

	stepSize.xy  *= parallaxDepth;
	stepSize.xyz *= parallaxStepSize;
	
	float heightmap = texture2D(normals, coord.st, 0).a;
	vec3  pCoord    = vec3(0.0, 0.0, 1.0);

	if(heightmap < 1.0 && heightmap != 0.0) {
	
		float distAngleWeight = ((distances * 0.6) * (2.1 - viewVector.z)) / 16.0;
		vec3 step = viewVector * stepSize;
			 step *= distAngleWeight * 2.0;
		float sampleHeight = heightmap;

		for (int i = 0; sampleHeight < pCoord.z && i < maxSteps; ++i) {

			pCoord.xy = mix(pCoord.xy, pCoord.xy + step.xy, clamp((pCoord.z - sampleHeight) / (stepSize.z * 1.0 * distAngleWeight / (-viewVector.z + 0.05)), 0.0, 1.0));
			pCoord.z += step.z;
			sampleHeight = texture2DGradARB(normals, OffsetCoord(coord.st, pCoord.st, 0), dcdx, dcdy).a;
			
		}
		
		parallaxCoord.xy = OffsetCoord(coord.st, pCoord.st, 0);
		
	}
	
	return parallaxCoord;
	
}

//#define SPECULAR_TO_PBR_CONVERSION
#define CONTINUUM2_TEXTURE_FORMAT

#define NEW_RAIN_SPLASHES
#define RAIN_SPLASH_LEVEL 1.5 //[0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0]

vec2 getWetness(inout vec4 sp, float height) {
	#ifdef NEW_RAIN_SPLASHES
	float wet0 = BiomeType.x * min(wetness * 6, 1.0) * 3.8 * smoothstep(0.5, 1.0, BiomeType.y) * RAIN_SPLASH_LEVEL;
	float rH = wet0 * plus(sp.r, sp.g) * smoothstep(0.92, 1.0, lmcoord.y) * top;
	float isWater = step(height, rH * 1.2);
	wet0 *= smoothstep(0.8, 1.0, lmcoord.y);
	sp.r = mix(plus(sp.r, wet0 * sp.g * 0.7), 0.95, isWater);
	sp.g = mix(sp.g, 0.03, isWater);
	sp.b = max(sp.b - wet0 * 0.1 - isWater * 0.4, 0.0);
	return vec2(isWater, (wet0 * 0.5 + isWater) * (1.0- plus(sp.g, sp.r)));
	#else
	return vec2(0.0);
	#endif
}

#define RAIN_SPLASH_WAVE

#ifdef RAIN_SPLASH_WAVE
#define hash_fast(p) fract(mod(p.x, 1.0) * 73758.23f - p.y)

float16_t hash(f16vec2 p) {
	vec3 p3 = fract(vec3(p.xyx) * 0.2031);
	p3 += dot(p3, p3.yzx + 19.19);
	return fract((p3.x + p3.y) * p3.z);
}

float16_t noise(f16vec2 p) {
	f16vec2 i = floor(p);
	f16vec2 f = fract(p);
	f16vec2 u = (f * f) * fma(f16vec2(-2.0f), f, f16vec2(3.0f));
	return fma(2.0f, mix(
		mix(hash(i),                      hash(i + f16vec2(1.0f,0.0f)), u.x),
		mix(hash(i + f16vec2(0.0f,1.0f)), hash(i + f16vec2(1.0f,1.0f)), u.x),
	u.y), -1.0f);
}
#endif

/* DRAWBUFFERS:024 */

void main() {
	const float pomRenderDistance = 32.0;
	
	vec2 adjustedTexCoord  = vtexcoord.st * vtexcoordam.pq + vtexcoordam.st;

	#ifdef POM
		if (distances < pomRenderDistance) adjustedTexCoord = CalculateParallaxCoord(texcoord.st, viewVector);
	#endif

    vec4 fragColor = texture2DGradARB(texture, adjustedTexCoord.st, dcdx, dcdy) * color;
	//vec4 specularData = texture2DGradARB(specular, adjustedTexCoord.st, dcdx, dcdy);
	
	vec4 sp; float nor2;
	#ifdef SPECULAR_TO_PBR_CONVERSION
	vec3 spec = texture2DGradARB(specular, adjustedTexCoord.st, dcdx, dcdy);
	float spec_strength = dot(spec, mix(vec3(0.4, 0.4, 0.2), vec3(0.3, 0.6, 0.1), wetness));
	sp = vec4(spec_strength, spec_strength, 0.0, 1.0);
	#else
	#ifdef CONTINUUM2_TEXTURE_FORMAT
	sp = texture2DGradARB(specular, adjustedTexCoord.st, dcdx, dcdy).brga;
	#else
	sp = texture2DGradARB(specular, adjustedTexCoord.st, dcdx, dcdy);
	#endif
	#endif
	
	#ifdef NORMALS
		vec4 norMap = texture2D(normals, adjustedTexCoord.st);
		vec2 wet = getWetness(sp, norMap.w);
		fragColor *= 1.0 - wet.y * 0.4;
		
		vec3 N = normal;
		f16vec3 normal2 = normal;
		if (distances < 64.0) {
			#ifdef RAIN_SPLASH_WAVE
			N.x += noise(wpos.xz * 5.0 - vec2(frameTimeCounter * 2.0, 0.0)) * 0.04 * wet.x * (1.0 + rain0);
			N.y -= noise(wpos.xz * 6.0 - vec2(frameTimeCounter * 2.0, 0.0)) * 0.04 * wet.x * (1.0 + rain0);
			#endif
			normal2 = norMap.xyz * 2.0 - 1.0;
			//rainNormal = rainNormal * 2.0 + 1.0;
			const float16_t bumpmult = 0.5;
			normal2 = normal2 * bumpmult + vec3(0.0f, 0.0f, 1.0f - bumpmult);
			//rainNormal = rainNormal * bumpmult + vec3(0.0f, 0.0f, 1.0f - bumpmult);
			
			f16mat3 tbnMatrix = mat3(tangent, binormal, normal);
			normal2 = tbnMatrix * normal2;
		}
		N = normalize(N);//  + vec4(rainNormal
		float n2 = encodeNormal(N);
		float d = encodeNormal(normal2);
		//if (!(d.x > 0.0 && d.y > 0.0)) d = n2;
		nor2 = mix(n2, d, min((1.0 - wet.x), Cselect(d, 0.0, 1.0)));
	#else
		nor2 = n2;
	#endif
	
	float roughness = 1.0 - sp.x;
	float f0 = sp.g;
	//encodeNormal(normal)
	vec2 lightmap = lmcoord.st * (1.0 / 255.0);

    gl_FragData[0] = fragColor * (1.0 + sp.b * 0.5);
	gl_FragData[1] = vec4(nor2, lightmap, minecraftMats + 0.1);
    gl_FragData[2] = vec4(roughness, f0, terrainMats / 255.0f, 1.0f);
}