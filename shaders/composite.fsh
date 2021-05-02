#version 130

//#define GI
	#define GI_RENDER_RESOLUTION 9.0 // Render resolution of GI. 1.0 = Original. Set higher for faster but blurrier GI. [1.0 2.0 4.0 9.0 16.0 25.0 36.0 49.0 64.0 81.0 100.0 225.0]
	#define GI_ARTIFACT_REDUCTION
	#define GI_QUALITY 1.0	//[0.1 0.25 0.3 0.4 0.5 0.75 1.0 1.5 2.0 2.5 3.0 3.5 4.0 5.0 6.0 7.0 8.0]
	#define GI_RADIUS 20.0	//[0.25 0.5 0.75 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0 16.0 20.0 25.0 30.0 32.0 40.0 50.0 60.0 64.0 70.0 80.0 90.0 100.0 125.0 150.0 200.0 250.0 300.0 400.0 500.0]
	#define GI_STEPS 2.0	//[1.0 1.25 1.5 1.75 2.0 2.25 2.5 2.75 3.0 4.0 8.0]
	#define FAST_GI

#define SSAO
	#define AO_RENDER_RESOLUTION 1.0 // Render resolution of AO. 1.0 = Original. Set higher for faster but blurrier AO. [1.0 4.0 9.0 16.0 25.0 36.0 49.0 64.0 81.0 100.0 225.0]
	//#define AO_DEBUG

const int 		R8 						    = 0;
const int 		RG8 					    = 0;
const int 		RGB8 					    = 0;
const int 		RGBA8 					    = 0;
const int 		RGB16 					    = 0;
const int 		RGBA16 					    = 0;
const int 		RGBA16F 					= 0;
const int 		R11F_G11F_B10F 			    = 0;
const int 		colortex0Format 			= RGBA16F;
const int 		colortex1Format 			= RGBA8;
const int 		colortex2Format 			= RGBA16;
const int 		colortex3Format 			= RGBA16F;
const int 		colortex4Format 			= RGBA16;
const int 		colortex5Format 			= RGBA16F;
const int 		colortex7Format 			= RGBA8;

const bool colortex0Clear = false;
const bool colortex3Clear = true;
const bool colortex5Clear = false;

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

/* DRAWBUFFERS:037 */

float 	ExpToLinearDepth(in float depth)
{
	return 2.0f * near * far / (far + near - (2.0f * depth - 1.0f) * (far - near));
}

vec3 ProjectBack(vec3 cameraSpace) 
{
    vec4 clipSpace = gbufferProjection * vec4(cameraSpace, 1.0);
    vec3 NDCSpace = clipSpace.xyz / clipSpace.w;
    vec3 screenSpace = 0.5 * NDCSpace + 0.5;
		 //screenSpace.z = 0.1f;
    return screenSpace;
}

float GetAO(vec2 coord, float dither)
{
	const int numRays = 24;

	float depth = GetDepth2(coord);
	float linDepth = ExpToLinearDepth(depth);
	vec3 origin = GetScreenPosition(coord, depth).xyz;
	vec3 normal = GetNormal(coord);

	float aoAccum = 0.0;

	float radius = 0.15 * -origin.z;
		  radius = mix(radius, 0.8, 0.5);
	float zThickness = 0.15 * -origin.z;
		  zThickness = mix(zThickness, 1.0, 0.5);

	float aoMul = 1.0;
	
	for (int i = 0; i < numRays; i++)
	{
		//dither = fract(dither + 0.618);
		float fi = float(i) + dither * 0.5;
		float fiN = fi / float(numRays);
		float lon = goldenAngle * fi * 6.0;
		float lat = asin(fiN * 2.0 - 1.0) * 1.0;

		vec3 kernel;
		kernel.x = cos(lat) * cos(lon);
		kernel.z = cos(lat) * sin(lon);
		kernel.y = sin(lat);

		kernel.xyz = normalize(kernel.xyz + normal.xyz);

		float sampleLength = radius * mod(fiN, 0.07) / 0.07;

		vec3 samplePos = origin + kernel * sampleLength;

		vec3 samplePosProj = ProjectBack(samplePos);

		vec3 actualSamplePos = GetScreenPosition(samplePosProj.xy, GetDepth2(samplePosProj.xy)).xyz;

		vec3 sampleVector = normalize(samplePos - origin);

		float depthDiff = actualSamplePos.z - samplePos.z;

		if (depthDiff > 0.0 && depthDiff < zThickness)
		{
			//aoAccum += 1.0 * saturate(depthDiff * 100.0) * saturate(1.0 - depthDiff * 0.25 / (sampleLength + 0.001));
			float aow = 1.25 * saturate(dot(sampleVector, normal));
			//aow *= saturate(dot(sampleVector, upVector) * 0.5 + 0.5) * 1.5 + 0.0;
			aoAccum += aow;
			//aoAccum += 1.0 * saturate(dot(kernel, upVector));
			//aoMul *= mix(1.0, 0.0, saturate(dot(kernel, upVector)));
			//aoMul *= 0.45;
		}
	}

	aoAccum /= numRays;

	float ao = 1.0 - aoAccum;
	//ao = aoMul;
	ao = pow(ao, 2.5);
	//ao *= exp2(ao);

	return ao;
}

vec2 DistortShadowSpace(in vec2 pos)
{
	vec2 signedPos = pos * 2.0f - 1.0f;

	float dist = sqrt(signedPos.x * signedPos.x + signedPos.y * signedPos.y);
	float distortFactor = (1.0f - SHADOW_MAP_BIAS) + dist * SHADOW_MAP_BIAS;
	signedPos.xy *= 0.95f / distortFactor;

	pos = signedPos * 0.5f + 0.5f;

	return pos;
}

#ifdef GI
	vec3 GetLight(in vec2 screenCoord, in float range, in float quality, surfaceStruct surface, float dither)
	{
			vec3 normal 				= GetNormal(screenCoord.st);						//Gets the screen-space normals

			vec4 gn = gbufferModelViewInverse * vec4(normal.xyz, 0.0f);
				gn = shadowModelView * gn;
				gn.xyz = normalize(gn.xyz);

			vec3 shadowSpaceNormal = gn.xyz;

			vec4 screenSpacePosition 	= GetScreenPosition(screenCoord.st); 			//Gets the screen-space position
			vec3 viewVector 			= normalize(screenSpacePosition.xyz);


			float distances = length(screenSpacePosition.xyz);

			vec4 upVectorShadowSpace = shadowModelView * vec4(0.0f, 1.0, 0.0, 0.0);

			vec4 worldposition = gbufferModelViewInverse * screenSpacePosition;		//Transform from screen space to world space
				worldposition = shadowModelView * worldposition;							//Transform from world space to shadow space
			float comparedepth = -worldposition.z;											//Surface distance from sun to be compared to the shadow map
			
			worldposition = shadowProjection * worldposition;								//Transform from shadow space to shadow projection space					
			worldposition /= worldposition.w;

			float d = sqrt(worldposition.x * worldposition.x + worldposition.y * worldposition.y);
			float distortFactor = (1.0f - SHADOW_MAP_BIAS) + d * SHADOW_MAP_BIAS;
			//worldposition.xy *= 1.0f / distortFactor;
			//worldposition.z = mix(worldposition.z, 0.5, 0.8);
			worldposition = worldposition * 0.5f + 0.5f;		//Transform from shadow projection space to shadow map coordinates

			const vec2[8] offsets = vec2[8](vec2(1./8.,-3./8.),
										vec2(-1.,3.)/8.,
										vec2(5.0,1.)/8.,
										vec2(-3,-5.)/8.,
										vec2(-5.,5.)/8.,
										vec2(-7.,-1.)/8.,
										vec2(3,7.)/8.,
										vec2(7.,-7.)/8.);

			//worldposition.xy += offsets[frameMod8] / 1024.0;

			float shadowMult = 0.0f;														//Multiplier used to fade out shadows at distance
			float shad = 0.0f;
			vec3 fakeIndirect = vec3(0.0f);

			float fakeLargeAO = 0.0;


			float mcSkylight = GetLightmap(screenCoord.st).g * 0.8 + 0.2;

			float fademult = 0.15f;

			shadowMult = clamp((shadowDistance * 41.4f * fademult) - (distances * fademult), 0.0f, 1.0f);	//Calculate shadowMult to fade shadows out

			worldposition.z -= 0.0025;

			float compare = sin(frameTimeCounter) > -0.2 ? 1.0 : 0.0;

			if (shadowMult > 0.0) 
			{
				//big shadow
				float rad = range;

				int c = 0;
				float s = 2.0f * rad / 2048;

				//dither *= 5.0;
				//vec2 dither = vec2(0.0f);
				//vec2 dither = rand(screenCoord.st + sin(frameTimeCounter)).xy / 2.0f;
				//float dither = bayer128(gl_FragCoord.xy + vec2(sin(frameTimeCounter), cos(frameTimeCounter)) * 1024.0) * 0.5;

				float step = 1.0f / quality;

				for (float i = -GI_STEPS; i <= GI_STEPS; i += step) {
					for (float j = -GI_STEPS; j <= GI_STEPS; j += step) {
						//dither = fract(dither + 0.618);
						vec2 offset = (vec2(i, j) + dither * step) / 100;
						// offset = length(offset) * normalize(offset + shadowSpaceNormal.xy * 0.12);

						//float cmpt=max(tmp.x,tmp.y);
						//if(tmp<=0.05) break;
						//float GI_RADIUS = 100.0;

						offset *= pow2(length(offset) * 3.0) * 15.0 + 0.02;
						offset *= GI_RADIUS;

						// offset *= gbufferSkylight + 0.2;

						// offset += shadowSpaceNormal.xz * 0.01;
						// offset = normalize(normalize(offset) + shadowSpaceNormal.xz * vec2(1.0, 1.0)) * length(offset);

						vec2 coord =  worldposition.st + offset;
						vec2 lookupCoord = DistortShadowSpace(coord);
						// lookupCoord = lookupCoord + (vec2(i, j) * mat2(sin(dither.x), cos(dither.x), -cos(dither.x), sin(dither.x))) / vec2(viewWidth, viewHeight);

						#ifdef GI_ARTIFACT_REDUCTION
						float depthSample = texture2DLod(shadowtex1, lookupCoord, 0).x;
						#else
						float depthSample = texture2DLod(shadowtex1, lookupCoord, 3).x;
						#endif

						/*
						depthSample = depthSample * 2.0 - 1.0;
						depthSample -= 0.4;
						depthSample /= 0.2;
						depthSample = depthSample * 0.5 + 0.5;
						*/

						depthSample = -3 + 5.0 * depthSample;
						vec3 samplePos = vec3(coord.x, coord.y, depthSample);

						vec3 lightVector = normalize(samplePos.xyz - worldposition.xyz);

						#ifdef FAST_GI
						vec4 normalSample = texture2DLod(shadowcolor1, lookupCoord, 6);
						#else
						vec4 normalSample = texture2DLod(shadowcolor1, lookupCoord, 0);
						#endif
						vec3 surfaceNormal = normalSample.rgb * 2.0f - 1.0f;
							surfaceNormal.x = -surfaceNormal.x;
							surfaceNormal.y = -surfaceNormal.y;

						float surfaceSkylight = normalSample.a;

						if (surfaceSkylight < 0.2)
						{
							surfaceSkylight = mcSkylight;
						}

						float NdotL = max(0.0f, dot(shadowSpaceNormal.xyz, lightVector * vec3(1.0, 1.0, -1.0)));
							NdotL = NdotL * 0.9f + 0.15f;
							// NdotL = 1.0;

						if (NdotL > 0.0)
						{
							bool isTranslucent = length(surfaceNormal) < 0.5f;

							if (isTranslucent)
							{
								surfaceNormal.xyz = vec3(0.0f, 0.0f, 1.0f);
							}

							//float leafMix = clamp(-surfaceNormal.b * 10.0f, 0.0f, 1.0f);


							float weight = dot(lightVector, surfaceNormal);
							float rawdot = weight;
							// float aoWeight = abs(weight);
								// aoWeight *= clamp(dot(lightVector, upVectorShadowSpace.xyz), 0.0, 1.0);
							//weight = mix(weight, 1.0f, leafMix);
							if (isTranslucent)
							{
								weight = abs(weight) * 0.85f;
							}

							if (normalSample.a < 0.2)
							{
								weight = 0.5;
							}

							weight = max(weight, 0.0f);

							float dist = length(samplePos.xyz - worldposition.xyz - vec3(0.0f, 0.0f, 0.0f));
							if (dist < 0.0015f)
							{
								dist = 10000000.0f;
							}

							const float falloffPower = 1.9f;
							float distanceWeight = (1.0f / (pow(dist * (62260.0f / rad), falloffPower) + 100.1f));
								//distanceWeight = max(0.0f, distanceWeight - 0.000009f);
								distanceWeight *= pow(length(offset), 2.0) * 50000.0 + 1.01;
							

							//Leaves self-occlusion
							if (rawdot < 0.0f)
							{
								distanceWeight = max(distanceWeight * 30.0f - 0.43f, 0.0f);
								distanceWeight *= 0.04f;
							}
								

							//float skylightWeight = clamp(1.0 - abs(surfaceSkylight - mcSkylight) * 10.0, 0.0, 1.0);
							float skylightWeight = 1.0 / (max(0.0, surfaceSkylight - mcSkylight) * 15.0 + 1.0);

							#ifdef FAST_GI
							vec3 colorSample = srgbToLinear(texture2DLod(shadowcolor, lookupCoord, 6).rgb);
							#else
							vec3 colorSample = srgbToLinear(texture2DLod(shadowcolor, lookupCoord, 0).rgb);
							#endif

							colorSample /= surfaceSkylight;
							//colorSample 				= pow(colorSample, vec3(1.4f));

							//colorSample = Contrast(colorSample, 0.8f);

							// colorSample = normalize(colorSample) * pow(length(colorSample), 1.0 / 0.8);



							fakeIndirect += colorSample * weight * distanceWeight * NdotL;
							//fakeIndirect += skylightWeight * weight * distanceWeight * NdotL;
							//fakeLargeAO += aoDistanceWeight * NdotL;
						}
						c += 1;
					}
				}

				fakeIndirect /= c;
			}

			fakeIndirect = mix(vec3(0.0f), fakeIndirect, vec3(shadowMult));

			return fakeIndirect.rgb * 100.0;
	}
#endif

#define DRAG_MULT 0.048
const int Iterations2 = 24;

vec2 wavedx(vec2 position, vec2 direction, float speed, float frequency, float timeshift) {
    float x = dot(direction, position) * frequency + timeshift * speed;
    float wave = exp(sin(x) - 1.0);
    float dx = wave * cos(x);
    return vec2(wave, -dx);
}

float getwaves(vec2 position, int iterations) {
	float iter = 0.0;
    float phase = 6.0;
    float speed = 2.0;
    float weight = 1.0;
    float w = 0.0;
    float ws = 0.0;
    for(int i=0;i<iterations;i++){
        vec2 p = vec2(sin(iter), cos(iter));
        vec2 res = wavedx(position, p, speed, phase, frameTimeCounter * 1.5);
        position += normalize(p) * res.y * weight * DRAG_MULT;
        w += res.x * weight;
        iter += 12.0;
        ws += weight;
        weight = mix(weight, 0.0, 0.2);
        phase *= 1.18;
        speed *= 1.07;
    }
    return w / ws;
}

float GetWaves2(vec3 p)
{
	return -(getwaves(p.xz * 0.1, Iterations2));
}

vec3 GetWavesNormal(vec3 position) {
	float q = 0.05;
	float w1 = GetWaves2(position + vec3(q, 0.0, 0.0));
	float w2 = GetWaves2(position + vec3(0.0, 0.0, q));
	float w0 = GetWaves2(position);

	float tangent = w1 - w0;
	float bitangent = w2 - w0;

	vec3 normal;
	normal.r = tangent * 5.0;
	normal.g = bitangent * 5.0;
	normal.b = 1.0;

	return normalize(normal);
}

void main() {
    surfaceStruct surface = GetSurfaceStruct();
    maskStruct mask = GetMaskStruct(texcoord.st);
    vec3 color = texture2D(colortex0, texcoord.st).rgb;
	float tempOffsets = HaltonSeq2(frameCounter%10000);

	//体积光
	/*
	const float shadowResolution = 1.0f / shadowMapResolution;

	vec4 worldposition = vec4(0.0f);
		 worldposition = gbufferModelViewInverse * surface.viewPos2;		//Transform from screen space to world space

	float totalDist = 0.0;
	float rayDist = 0.0;

	vec3 rayPos = worldposition.xyz - vec3(0.0,1.62,0.0);
	float rayDepth = length(rayPos);
	float stepLength = min(shadowDistance, rayDepth) / ATOMSPHERIC_RAYS_STEP;
	vec3 rayDir = normalize(rayPos) * stepLength;

	float dither = calculateBlueNoise(texcoord.st, 1.0).x * 5.0;

	float dist = 1.0;

	//float sunSpot = calculateSunGlow(surface.viewDir, surface.lightVector);
	float anisoFactor = MiePhase(0.8, worldDir, worldLightVector) + MiePhase(0.5, worldDir, -worldLightVector) * 0.1;

	for(int i = 0; i < ATOMSPHERIC_RAYS_STEP; i++)
	{
		worldposition.xyz -= rayDir;
		//dither = fract(dither + 0.618);

		vec3 pos = clacWorldposToShadowpos(worldposition.xyz + rayDir * dither);

		float sample = texture2D(shadowtex0, pos.st).x;
		rayDist = float(pos.z + 0.0006 < sample);
		//rayDist = sample;
		totalDist += mix(rayDist * stepLength, rayDist * stepLength * 0.05, mask.sky);
	}
	if(isEyeInWater > 0.5 || mask.water > 0.5) totalDist /= 4.0;
	else totalDist /= 8.0;
	*/

	// vec3 rayColorMix = sunColor * transitionFading * totalDist * 0.0000005;
	// if(isEyeInWater > 0.5 || mask.water > 0.5) rayColorMix *= vec3(0.0824, 0.3333, 0.604);
	// rayColorMix *= ATOMSPHERIC_RAYS_DESENTIYE;

	// #ifdef ATOMSPHERIC_RAYS
	// 	if(isEyeInWater > 0.5 || mask.water > 0.5) color += rayColorMix;
	// 	else color += rayColorMix * sunSpot;
	// #endif
	

	#ifdef AO_DEBUG
		gl_FragData[0] = vec3(1.0);
	#else
		gl_FragData[0] = vec4(color, 1.0);
	#endif
//--------------------------------------------------------------------------------------------------
	float ao = 1.0;

	vec4 lightGI = vec4(vec3(0.0), 1.0);

	//UV
	float zoomGI = 1.0f / GI_RENDER_RESOLUTION;
	if (zoomGI < 0.0f || zoomGI > 1.0f)
	{
		zoomGI = 1.0f;
	}
	float zoomAO = 1.0f / AO_RENDER_RESOLUTION;
	if (zoomAO < 0.0f || zoomAO > 1.0f)
	{
		zoomAO = 1.0f;
	}

	vec2 coordGI = texcoord.st / sqrt(zoomGI);
	vec2 coordAO = texcoord.st / sqrt(zoomAO);

	if (coordGI.s <= 1.0f && coordGI.s >= 0.0f
	 && coordGI.t <= 1.0f && coordGI.t >= 0.0f)
	{
		#ifdef GI
			//float dither_GI = bayer128(gl_FragCoord.xy + vec2(sin(frameTimeCounter), cos(frameTimeCounter)) * 1024.0) * 0.5;
			//float dither_GI = rand(texcoord.st + vec2(sin(frameTimeCounter), cos(frameTimeCounter))).x * 0.5;
			float dither_GI = calculateBlueNoise(texcoord.st, 1.0).x * 0.5;
			//float dither_GI = 0.25f;
			lightGI.rgb  = GetLight(coordGI, 16.0f, GI_QUALITY, surface, dither_GI);
		#endif
	}
	if (coordAO.s <= 1.0f && coordAO.s >= 0.0f
	 && coordAO.t <= 1.0f && coordAO.t >= 0.0f)
	{
		#ifdef SSAO
			//ao = GetAO(coordAO.st, rand(texcoord.st + sin(frameTimeCounter)).x);
			ao = GetAO(coordAO.st, calculateBlueNoise(texcoord.st, 1.0).x * 0.5);
		#endif
	}

	vec3 normal = vec3(0.0);
		 normal = GetWavesNormal(vec3(texcoord.s * 50.0, 1.0, texcoord.t * 50.0));

	gl_FragData[1] = vec4(lightGI.rgb, ao);
	gl_FragData[2] = encodeRGBE8(normal);
}
