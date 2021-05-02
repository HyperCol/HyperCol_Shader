#version 130

#define TAA

#define TORCHLIGHT_COLOR_TEMPERATURE 3000	//[2000 2300 2500 3000]

#define COLORED_SHADOW

//#define GI
	#define GI_RENDER_RESOLUTION 9.0 // Render resolution of GI. 1.0 = Original. Set higher for faster but blurrier GI. [1.0 4.0 9.0 16.0 25.0 36.0 49.0 64.0 81.0 100.0 225.0]
	#define GI_BLUR_QUALITY 0.25 //[0.25 0.5 1.0 1.25 1.5 1.75 2.0 2.25 2.5 3.0 4.0 5.0 6.0 7.0 8.0 16.0]
	#define GI_BRIGHTNESS 1.0	//[0.25 0.5 0.75 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0 20.0 30.0 40.0 50.0 60.0 70.0 80.0 90.0 100.0]

#define SSAO
	#define AO_RENDER_RESOLUTION 1.0 // Render resolution of AO. 1.0 = Original. Set higher for faster but blurrier AO. [1.0 4.0 9.0 16.0 25.0 36.0 49.0 64.0 81.0 100.0 225.0]
	#define AO_BLUR_QUALITY 0.5 //[0.25 0.5 1.0 1.25 1.5 1.75 2.0 2.25 2.5 3.0 4.0 5.0 6.0 7.0 8.0 16.0]

//#define PLANE_CLOUD

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
varying mat3x4 skySH;

varying float transitionFading;

uniform vec3 shadowLightPosition;

#include "libs/uniform.glsl"
#include "libs/composite.glsl"
#include "libs/noise.glsl"

#include "/libs/sky.glsl"

vec3 calculateShadow(vec4 screenPosition, vec3 normal) 
{
	float distance = length(screenPosition.xyz);

	const float iStep = 2.0f;
	const float jStep = 2.0f;
	vec4 ssp = screenPosition;

	const float shadowResolution = 1.0f / shadowMapResolution;

	vec4 worldposition = vec4(0.0f);
		 worldposition = gbufferModelViewInverse * ssp;		//Transform from screen space to world space

	float worldDepth = worldposition.z;

	worldposition = shadowModelView * worldposition;	//Transform from world space to shadow space

	worldposition = shadowProjection * worldposition;
	worldposition /= worldposition.w;

	float dist = sqrt(worldposition.x * worldposition.x + worldposition.y * worldposition.y);
	float distortFactor = (1.0f - SHADOW_MAP_BIAS) + dist * SHADOW_MAP_BIAS;
	worldposition.xy *= 0.95f / distortFactor;
	worldposition.z = mix(worldposition.z, 0.5, 0.8);
	worldposition = worldposition * 0.5f + 0.5f;		//Transform from shadow space to shadow map coordinates

	const vec2[8] offsets = vec2[8](vec2(1./8.,-3./8.),
									vec2(-1.,3.)/8.,
									vec2(5.0,1.)/8.,
									vec2(-3,-5.)/8.,
									vec2(-5.,5.)/8.,
									vec2(-7.,-1.)/8.,
									vec2(3,7.)/8.,
									vec2(7.,-7.)/8.);

	float shadowMult = 0.0f;																			//Multiplier used to fade out shadows at distance

	float fademult = 0.15f;
		shadowMult = clamp((shadowDistance * 1.4f * fademult) - (distance * fademult), 0.0f, 1.0f);	//Calculate shadowMult to fade shadows out

	float shadowDepth;

	float shading;

	vec3 result;
	if(shadowMult > 0.0)
	{
		#ifdef TAA
			worldposition.z -= 0.00005;
		#endif

		float diffthresh = dist * 1.0f + 0.10f;
			  diffthresh *= 1.5f / (shadowMapResolution / 2048.0f);

		float vpsSpread = 0.145 / distortFactor;

		float avgDepth = 0.0;
		float minDepth = 11.0;
		int c;

		for (int i = -1; i <= 1; i++)
		{
			for (int j = -1; j <= 1; j++)
			{
				vec2 lookupCoord = worldposition.xy + (vec2(i, j) / shadowMapResolution) * 8.0 * vpsSpread;
				//avgDepth += pow(texture2DLod(shadowtex1, lookupCoord, 2).x, 4.1);
				float depthSample = texture2DLod(shadowtex1, lookupCoord, 2).x;
				minDepth = min(minDepth, depthSample);
				avgDepth += pow(min(max(0.0, worldposition.z - depthSample) * 1.0, 0.025), 2.0);
				c++;
			}
		}

		avgDepth /= c;
		avgDepth = pow(avgDepth, 0.5);

		float penumbraSize = avgDepth;

		int count = 0;
		float spread = penumbraSize * 0.125 * vpsSpread + 0.25 / shadowMapResolution;
		spread += dist * 2.0 / shadowMapResolution;

		diffthresh *= 0.45 + avgDepth * 80.0;

		int stepTotal = 0;
		float dither = calculateBlueNoise(texcoord.st, 1.0).x;
		for(float i = -iStep; i < iStep; i++) {
			for(float j = -jStep; j < jStep; j++) {
				dither *= TAU;
				//dither = fract(dither + 0.618);
				vec2 offset = vec2(i, j) * mat2(cos(dither), sin(dither), -sin(dither), cos(dither));
				shading += shadow2D(shadow, vec3(worldposition.xy + offset * spread, worldposition.z - 0.0012f * diffthresh)).x;
				//worldposition.z -= 0.0012f * diffthresh;
				//shading += 1.0 - float((texture2D(shadowtex1, worldposition.xy + offset * spread).x - 0.0012f * diffthresh) < worldposition.z);
				stepTotal++;
			}
		}
		shading /= stepTotal;

		float clampFactor = max(0.0, dist - 0.1) * 5.0 + 1.0;
		shading = saturate(((shading * 2.0 - 1.0) * clampFactor) * 0.5 + 0.5);


		result = vec3(shading);

		#ifdef COLORED_SHADOW
			float shadowNormalAlpha = texture2DLod(shadowcolor1, worldposition.st, 0).a;

			if (shadowNormalAlpha < 0.1)
			{
				#ifdef TAA
					vec2 noise = calculateBlueNoise(texcoord.st, 1.0).xx;
				#else
					vec2 noise = bayer128(gl_FragCoord.xy) * vec2(1.0) * 0.5;
				#endif

				float angle = noise.x * 3.14159 * 2.0;
				mat2 rot = mat2(cos(angle), -sin(angle), sin(angle), cos(angle));

				float solidShadowSum = 0.0;
				vec3 shadowColorSampleSum = vec3(0.0);

				int c = 0;

				for (float i = -1.5f; i <= 1.5f; i += 1.0f) 
				{
					for (float j = -1.5f; j <= 1.5f; j += 1.0f) 
					{
						worldposition.st += (vec2(i + noise.y - 0.5, j + noise.y - 0.5) * rot) * (0.5 / shadowMapResolution);
						vec4 shadowColorSample = texture2DLod(shadowcolor, worldposition.st, 0);
						float opacityCheck = 1.0 - saturate(pow(shadowColorSample.a * 1.1, 4.0));
						shadowColorSampleSum += pow(shadowColorSample.rgb, vec3(1.6)) * (opacityCheck);
						float solidDepth = texture2DLod(shadowtex1, worldposition.st, 0).x;
						float solidShadow = 1.0 - clamp((worldposition.z - solidDepth) * 5200.0, 0.0, 1.0); 
						solidShadowSum += solidShadow;
						c++;
					}
				}

				solidShadowSum /= c;
				shadowColorSampleSum /= c;

				result = mix(vec3(1.0), shadowColorSampleSum.rgb, vec3(1.0 - shading));
				result *= solidShadowSum;
			}
		#endif
	}
	
	result = mix(vec3(1.0), result, shadowMult);

	return result;
}

float 	ExpToLinearDepth(in float depth)
{
	return 2.0f * near * far / (far + near - (2.0f * depth - 1.0f) * (far - near));
}

void RayPlaneIntersectionWorld(vec3 dir, vec3 normal, vec3 origin, out vec3 pos, out float distance, out float angle)
{
	float rayPlaneAngle = dot(dir, normal);

	float planeRayDist = 100000000.0f;
	vec3 intersectionPos = dir * planeRayDist;

	if (rayPlaneAngle > 0.0001f || rayPlaneAngle < -0.0001f)
	{
		planeRayDist = dot((origin), normal) / rayPlaneAngle;
		intersectionPos = dir * planeRayDist;
		intersectionPos = -intersectionPos;

		intersectionPos += cameraPosition.xyz;
	}

	pos = intersectionPos;
	distance = planeRayDist;
	angle = rayPlaneAngle;
}

float GetCoverage(in float coverage, in float density, in float clouds)
{
	clouds = clamp(clouds - (1.0f - coverage), 0.0f, 1.0f -density) / (1.0f - density);
		clouds = max(0.0f, clouds * 1.1f - 0.1f);
	 clouds = clouds = clouds * clouds * (3.0f - 2.0f * clouds);
	 // clouds = pow(clouds, 1.0f);
	return clouds;
}

vec4 GetCloudPlaneColor(vec3 worldPos) {
    vec3 pos = (worldPos.xyz + worldPos.xyz * 2.5) / 5000.0;

        vec2 coord = pos.xz;
		float wind = frameTimeCounter * 0.005;
    	mat2 octave_m = mat2(1.4, 1.2, -1.1, 1.4);

		coord -= wind * 0.1;

		float noise = texture2D(noisetex, coord * vec2(0.65, 1.0) * 0.2).x * 0.45;	coord *= 3.0;
		    coord += vec2(noise * 0.6, 0.0) * octave_m;
		noise += texture2D(noisetex, coord * 0.4).x * 0.25; coord *= 3.01;
		    coord += vec2(noise * 0.4, 0.0) * octave_m + vec2(wind * 0.1, 0.2); coord += wind;
		noise += texture2D(noisetex, coord * 0.45).x * 0.105; coord *= 4.02;
		    coord += vec2(noise * 0.3, 0.0) * octave_m + vec2(wind * 0.03, 0.1);
		noise += texture2D(noisetex, coord * 0.25).x * 0.0625;
		noise /= 1.45;

	float coverage = 0.82f;
		coverage = mix(coverage, 0.97f, rainStrength);

	float dist = length(worldPos.xz * 0.5);
		coverage *= max(0.0f, 1.0f - dist / 14000.0f);
	float density = 0.1f + rainStrength * 0.3;

	noise = GetCoverage(coverage, density, noise);

	//noise = smoothstep(0.0, 1.0, noise + wetness);

	vec3 color = mix(sunColor, skyColor, 0.5);

    return vec4(color * 1.75f, noise);
}

vec3 GetCloudPlane(vec3 color, vec4 cameraPos, float sky, float depth) {
	
	vec4 worldVector = gbufferModelViewInverse * (-cameraPos);

	vec3 worldDir = normalize(worldVector.xyz);

	float cloudsAltitude = 1000.0f;
	float cloudsThickness = 150.0f;

	float cloudsUpperLimit = cloudsAltitude + cloudsThickness * 0.5f;
	float cloudsLowerLimit = cloudsAltitude - cloudsThickness * 0.5f;

	float density = 1.0f;

	float planeHeight = cloudsUpperLimit;
	float stepSize = 25.5f;
	planeHeight -= cloudsThickness * 0.85f;

	vec3 origin = vec3(0.0f, cameraPosition.y - planeHeight, 0.0f);
	vec3 normal = vec3(0.0f, 1.0f, 0.0f);
	float angle = 0.0f;
	float distance = 0.0f;
	vec3 ipos = vec3(0.0);
	RayPlaneIntersectionWorld(worldDir, normal, origin, ipos, distance, angle);

	float linearDepth = ExpToLinearDepth(depth);

	if (angle < 0.0f)
	{
		if (distance < linearDepth || sky > 0.5 || linearDepth >= far - 0.1)
		{
			
			vec4 cloudSample = GetCloudPlaneColor(vec3(ipos.xyz * 0.5f + vec3(30.0f)));
			 	 cloudSample.a = min(1.0f, cloudSample.a * density);

			float cloudDist = length(ipos.xyz - cameraPosition.xyz);

			const vec3 absorption = vec3(0.2, 0.4, 1.0);

			cloudSample.rgb *= exp(-cloudDist * absorption * 0.0001 * (1.0 - rainStrength));

			cloudSample.a *= exp(-cloudDist * (0.0002 + rainStrength * 0.0029));

			//cloudSample.rgb *= sin(cloudDist * 0.3) * 0.5 + 0.5;

			color.rgb = mix(color.rgb, cloudSample.rgb, cloudSample.a);

		}
	}
	
    return color;
}

vec3 GetWavesNormal(vec3 worldpos)
{
	vec2 coord = worldpos.xz / 50.0;
	coord.xy -= worldpos.y / 50.0;

	coord = mod(coord, vec2(1.0));

	float texelScale = 4.0;

	//to fix color error with GL_CLAMP
	coord.x = coord.x * ((viewWidth - 1 * texelScale) / viewWidth) + ((0.5 * texelScale) / viewWidth);
	coord.y = coord.y * ((viewHeight - 1 * texelScale) / viewHeight) + ((0.5 * texelScale) / viewHeight);

	vec3 normal = decodeRGBE8(texture2DLod(colortex7, coord, 2));
	return normal;
}

float CalculateWaterCaustics(vec3 worldpos, vec3 worldpos0) 
{
	vec2 dither = calculateBlueNoise(texcoord.st,1.0).xx;
	
	float waterPlaneHeight = 65.0 + distance(worldpos, worldpos0);
		
	vec4 ViewlightVector = gbufferModelViewInverse * vec4(lightVector.xyz, 0.0);
	vec3 worldLightVector = -normalize(ViewlightVector.xyz);

	float pointToWaterVerticalLength = min(abs(worldpos.y - waterPlaneHeight), 2.0);
	vec3 flatRefractVector = refract(worldLightVector, vec3(0.0, 1.0, 0.0), 1.0 / 1.3333);
	float PTWLength = pointToWaterVerticalLength / -flatRefractVector.y;
	vec3 lookupCenter = worldpos.xyz - flatRefractVector * PTWLength;

	const float distanceThreshold = 0.135;

	const int numSamples = 1;
	int c = 0;

	float caustics = 0.0;

	for (int i = -numSamples; i <= numSamples; i++) {
		for (int j = -numSamples; j <= numSamples; j++) {
			vec2 offset = vec2(i + dither.x, j + dither.y) * 0.15;
			vec3 lookupPoint = lookupCenter + vec3(offset.x, 0.0, offset.y);
			vec3 wavesNormal = GetWavesNormal(lookupPoint).xzy;
			vec3 refractVector = refract(worldLightVector.xyz, wavesNormal.xyz, 1.0 / 1.3333);
			float rayLength = pointToWaterVerticalLength / refractVector.y;
			vec3 collisionPoint = lookupPoint - refractVector * rayLength;

			float dist = distance(collisionPoint, worldpos.xyz);

			caustics += 1.0 - abs(saturate(dist / distanceThreshold));

			c++;
		}
	}

	caustics /= c;

	caustics /= distanceThreshold;

	//caustics = 1.0 - saturate(caustics);
	return smoothstep(0,1,pow(caustics,2));
	
} 

//vec3 calculateAtmosphereRays()

float calculateTorchLightAttenuation(float lightmap){
	lightmap *= 1.135;
	lightmap = clamp01(lightmap);

	float dist = (1.0 - lightmap) * 16.0 + 1.0;
	return lightmap * pow(dist, -2.0);
}

vec4 GetBlurTexture(in sampler2D tex, in vec2 coord) {
	vec2 res = vec2(viewWidth, viewHeight);
	coord = coord * res + 0.5f;
	vec2 i = floor(coord);
	vec2 f = fract(coord);
	f = f * f * (3.0f - 2.0f * f);
	coord = i + f;
	coord = (coord - 0.5f) / res;
	return BicubicTexture(tex, coord);
}

vec4 calcFilterColor(in sampler2D tex, in vec2 coord, in float quality)
{
	float steps = 1.0f / quality;
	int totalsteps = 0;

	vec4 total = vec4(0.0f);
	for(float x = -2.0f; x < 2.0f; x += steps)
		for(float y = -2.0f; y < 2.0f; y += steps)
		{
			float angle = exp(-sin(totalsteps + 1)) * TAU;
			vec2 offset = rotate(vec2(x, y), angle) * texelSize;
			total += GetBlurTexture(tex, coord + offset);
			totalsteps++;
		}

	return total / totalsteps;
}

vec3 GetLightingEffect(vec3 color, maskStruct mask)
{
	if(mask.land < 0.5) return color;

	float zoom = 1.0f / GI_RENDER_RESOLUTION;
	if (zoom < 0.0f || zoom > 1.0f)
	{
		zoom = 1.0f;
	}

	//GI
	// const float iStep = 0.5f;
	// const float jStep = 0.5f;
	// float quality = 1./GI_BLUR_QUALITY;

	// int total = 1;

	vec3 gi = vec3(0.0);

	// float dither = rand(texcoord.st + sin(frameTimeCounter)).x / 2.0f;
	// for(float i = -iStep; i < iStep; i += quality)
	// for(float j = -jStep; j < jStep; j += quality)
	// {
	// 	float angle = dither * TAU;
	// 	vec2 offset = mat2(cos(angle), -sin(angle), sin(angle), cos(angle)) * vec2(i, j) / 512.0;
	// 	gi += texture2D(colortex3, texcoord.st * sqrt(zoom) + offset * sqrt(zoom)).rgb;
	// 	total++;
	// }
	// gi /= total;

	//gi += texture2DLod(colortex3, texcoord.st * sqrt(zoom), 0).rgb;
	gi += calcFilterColor(colortex3, texcoord.st * sqrt(zoom), GI_BLUR_QUALITY).rgb;
	gi = gi * (
	mask.torch * 0.0015 + mask.lava * 0.075 + mask.glowstone * 0.03 + mask.fire * 0.045 + 
	mask.sky * 0.015 + 
	mask.ice * 0.0075 + 
	mask.goldBlock * 0.25 + mask.ironBlock * 0.25 + mask.diamondBlock * 0.45 + mask.emeraldBlock * 0.35 +
	mask.wool * 2.5 + mask.grass * 0.15 + 
	mask.land * 0.1);

	gi *= 3.0;
	gi *= GI_BRIGHTNESS;
	gi *= transitionFading;
	//gi /= pow(eyeBrightnessSmooth.y / 240.0f, 6.0f);

	#ifdef GI
		color += gi;
	#endif

	return color;
}

float GetAOTexture()
{
	float zoom = 1.0f / AO_RENDER_RESOLUTION;
	if (zoom < 0.0f || zoom > 1.0f)
	{
		zoom = 1.0f;
	}

	//SSAO
	const float iStep = 2.0f;
	const float jStep = 2.0f;
	float quality = 1./AO_BLUR_QUALITY;

	int total = 1;

	float ao = 1.0;

	ao = calcFilterColor(colortex3, texcoord.st * sqrt(zoom), AO_BLUR_QUALITY).a;
	//ao /= total;

	return ao;
}

// float calculateHypercolAO(vec2 coord, vec3 pos, vec3 normal)
// {
// 	const int steps = 6;
// 	float dither = calculateBlueNoise().x;

// 	float radius = radians(pos.y);

// 	float PdotN = dot(pos, normal);

// 	float total = 0.0;
// 	for(int i = 0; i <= steps; i++)
// 	{
// 		vec2 offset = vec2(sin(i), cos(i)) * vec2(256.0) + dither;

// 		vec3 aopos = GetWorldPosition(coord + offset).xyz;

// 		total += dot(normal, normalize(aopos - pos));
// 		dither = fract(dither + 0.618);
// 	}

// 	float ao = smoothstep(0.0, 1.0, radius / abs(pos.z - total));
// 	return total * 16000.0;
// }

bool ismmats(float mmats, int id)
{
	if(mmats >= id-0.1 && mmats < id+1)
	return true;
	else return false;
}

void main() {
    surfaceStruct surface = GetSurfaceStruct();
	maskStruct mask = GetMaskStruct(texcoord.st);
    vec3 color = vec3(0.0);
	vec3 albedo = srgbToLinear(texture2DLod(colortex0, texcoord.st, 0).rgb);

	vec4 data4 = texture2D(colortex4, texcoord.st);
	float roughness = data4.x, f0 = data4.y;
	//getRoughnessF0(texture2D(colortex4, texcoord.st).r, roughness, f0);

	vec2 planetSphere = vec2(0.0);
	vec3 sky = vec3(0.0);
	vec3 skyAbsorb = vec3(0.0);

	vec3 worldLightVector = normalize((shadowModelViewInverse * vec4(0.0, 0.0, 1.0, 0.0)).xyz);

    float vDotL = dot(surface.viewDir, sunVector);

	sky = calculateAtmosphere(vec3(0.0), surface.viewDir, upVector, sunVector, moonVector, planetSphere, skyAbsorb, 25);
	// sky = AtmosphericScattering(surface.worldDir, wSunVector, 1.0);
	// sky = mix(sky, vec3(0.6) * Luminance(skyColor), vec3(0.3 * wetness));
	// sky *= 1000.0;
	if (mask.sky > 0.5) {
		color = sky;

		color += calculateSunSpot(vDotL) * sunColorBase * skyAbsorb;
		color += calculateMoonSpot(-vDotL) * moonColorBase * skyAbsorb;
		color += calculateStars(surface.worldDir, wMoonVector) * skyAbsorb;

		#ifdef PLANE_CLOUD
			color = GetCloudPlane(color, surface.viewPos2, mask.sky, surface.depth2);
		#endif

		color = encodeColor(color);
		if(mask.translucent > 0.5) color += albedo * 10.0;
	} else {
		//Lighting
		vec3 shadows = calculateShadow(surface.viewPos2, surface.worldNormal);

		//color += surface.lightmap.r * albedo * 0.5 * 0.7;

		vec3 colorTorchlight;

		if (TORCHLIGHT_COLOR_TEMPERATURE == 2000)
			//2000k
			colorTorchlight = vec3(255, 141, 11) / 255.0;
		else if (TORCHLIGHT_COLOR_TEMPERATURE == 2300)
			//2300k
			colorTorchlight = vec3(255, 152, 54) / 255.0;
		else if (TORCHLIGHT_COLOR_TEMPERATURE == 2500)
			//2500k
			colorTorchlight = vec3(255, 166, 69) / 255.0;
		else
			//3000k
			colorTorchlight = vec3(255, 180, 107) / 255.0;

		//colorTorchlight = srgbToLinear(colorTorchlight);

		float ao = 1.0f;
		#ifdef SSAO
			if(mask.land > 0.5) ao = GetAOTexture();
		#endif

		#define lumCoeff vec3(0.2125, 0.7154, 0.0721)
		float moonlum = dot(moonColor, lumCoeff);
		float sunlum = dot(sunColor, lumCoeff);
		float albedolum = dot(albedo, lumCoeff);
		vec3 unsaturatedAlbedo = mix(albedo, vec3(albedolum), clamp01(moonlum * 2.0 - sunlum));

		vec3 viewVector = -normalize(surface.viewPos.xyz);
		vec3 shadowLightVector = normalize(shadowLightPosition);
		vec3 diffuse = GeometrySmithGGX(albedo, surface.normal, viewVector, shadowLightVector, roughness);
		//vec3 diffuse = vec3(rPI);

		vec3 skylight = FromSH(skySH[0], skySH[1], skySH[2], surface.worldNormal) * PI;
		//skylight *= 1.0 - rainStrength * 0.55;
		skylight *= surface.lightmap.y;
		//skylight *= surface.lightmap.g;

		if(mask.water > 0.5 || isEyeInWater > 0)
		{
			diffuse *= CalculateWaterCaustics(surface.worldPos2.xyz + cameraPosition, surface.worldPos.xyz + cameraPosition);
		}

		//color += (sunColor + moonColor * 100.0) * albedo * unsaturatedAlbedo * 0.25;

		//color += skylight * unsaturatedAlbedo;
		//color += diffuse * albedo * shadows * (sunColor + moonColor * 200.0);

		//color += albedo * surface.lightmap.r * colorTorchlight * 1000.0;
		//color += skylight * albedo * 2000.0;

		//float calculateHypercolAO(vec2 coord, vec3 pos, vec3 pos2, vec3 normal)
		//float ao = calculateHypercolAO(texcoord.st, surface.viewPos.xyz, surface.normal);
		vec3 torchLightColor = vec3(0.0);

		float mMats = texture2D(colortex2, texcoord.st).b * 255.0;
		vec3 lightColor = vec3(0.0);
		if(ismmats(mMats,1))
		{
			lightColor = vec3(255, 48, 48) / 255.0;
		}
		else if(ismmats(mMats,2))
		{
			lightColor = vec3(255, 251, 0) / 255.0;
		}
		else if(ismmats(mMats,3))
		{
			lightColor = vec3(30, 144, 255) / 255.0;
		}
		else if(ismmats(mMats,4))
		{
			lightColor = vec3(1);
		}
		else if(ismmats(mMats,5))
		{
			lightColor = vec3(186, 85, 211) / 255.0;
		}

		lightColor = srgbToLinear(lightColor);

		torchLightColor += colorTorchlight * 0.0075 * mask.torch;
		torchLightColor += colorTorchlight * 0.0150 * mask.lava;
		torchLightColor += colorTorchlight * 0.0085 * mask.glowstone;
		torchLightColor += colorTorchlight * 0.0100 * mask.fire;

		//torchLightColor *= 10.0;
		color += albedo * torchLightColor;

		vec3 directionalLighting = shadows * diffuse * transitionFading;
		vec3 lighting = vec3(0.0);

		lighting += calculateTorchLightAttenuation(surface.lightmap.x) * colorTorchlight * 2.0 * ao;
		lighting += 0.000001 * (-surface.lightmap.y + 1.0);
		lighting += directionalLighting * sunColor * 0.000005;

		//albedo = lightmaps.x > 0.99 ? vec3(1.0) : albedo;

		lighting *= albedo;

		lighting = GetLightingEffect(lighting, mask);

		lighting += skylight * unsaturatedAlbedo * ao;
		lighting += directionalLighting * moonColor * albedolum * 0.00005; //Fake Purkinje effect

		color += lighting * 10.0;
		//color = linearToSRGB(color);
		//color += ao * vec3(1.0);

		//color *= 0.01;
		//color = pow(color, vec3(2.2));

		//color += (sky * 0.000005) * far * 0.01;
	}

	//color = texture2D(shadowtex1, texcoord.st).rgb;

	//if(mask.water > 0.5) color *= .0;

/* DRAWBUFFERS:0 */
    gl_FragData[0] = encodeRGBE8(max0(color));
}