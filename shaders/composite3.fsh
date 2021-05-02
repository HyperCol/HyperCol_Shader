#version 130

#define WATER_REFRACT_IOR 1.2

#define WATER_PARALLAX
	#define WATER_ENHANCE

#define WATER_REFLECT
	#define REFLECT_ISTEP 16.0	//[16.0 20.0 24.0 40.0 48.0 64.0 128.0 218.0]
	#define REFLECT_QUALITY 1.0	//[0.25 0.5 1.0 1.25 1.5 1.75 2.0 2.25 2.5 3.0 4.0 5.0 6.0 7.0 8.0 16.0]

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

varying float transitionFading;

#include "libs/uniform.glsl"
#include "libs/composite.glsl"
#include "libs/noise.glsl"

#include "/libs/sky.glsl"

float curve(float x)
{
	return x * x * (3.0 - 2.0 * x);
}

vec3 GetWaterFogColor()
{
	// return vec3(0.2, 0.8, 1.0);
	// return vec3(0.8, 0.9, 1.0) * 0.7;
	return vec3(0.8, 0.9, 1.0) * 0.5;
	// return vec3(1.0, 1.9, 0.4);
}

vec3 GetWaterAbsorption()
{
	// return vec3(1.0, 0.22, 0.005);
	// return vec3(1.0, 0.32, 0.05);

	// return vec3(0.25, 0.03, 0.01);// * 0.0;
	return vec3(0.25, 0.04, 0.01);// * 0.0;
}

vec3 TintUnderwaterDepth(vec3 color)
{
	// return color;
	if (isEyeInWater > 0)
	{
		float underwaterDepth = 1.0 - (eyeBrightnessSmooth.y / 240.0);
		color *= exp((-GetWaterAbsorption()) * (underwaterDepth * 8.0 + 0.0));
	}

	return color;
}

vec3 TintUnderwaterMCSkylight(vec3 color, float mcSkylight)
{
	return color;
	if (isEyeInWater > 0)
	{
		float underwaterDepth = 1.0 - mcSkylight;
		color *= exp(-GetWaterAbsorption() * underwaterDepth * 1.0);
	}

	return color;
}

void UnderwaterFog(inout vec3 color, float eyeLength, vec3 eyeDir, vec3 skyColor, vec3 sunColor)
{
	if (isEyeInWater == 0)
	return;

	vec3 fogColor = GetWaterFogColor() * 0.05;
	// fogColor *= exp(-GetWaterAbsorption() * ((1.0 - (eyeBrightnessSmooth.y / 240.0)) * 4.0 + (-eyeDir.y * 0.5 + 0.5) * 3.0));
	fogColor *= exp(-GetWaterAbsorption() * ((-eyeDir.y * 0.5 + 0.5) * 3.0));
	fogColor *= pow(curve(saturate(mix(eyeDir.y, 1.0, 0.5))), 2.0) * 0.6 + 0.4;

	// fogColor *= (skyColor * 0.1 + sunColor * 2.0);
	fogColor *= (skyColor * 1.0 + sunColor * 2.0) * 1.4;
	fogColor = TintUnderwaterDepth(fogColor);
	// fogColor *= eyeBrightnessSmooth.y / 240.0;

	float fogFactor = (exp(-eyeLength * 0.03));
	fogFactor = saturate(fogFactor * 1.1);
	fogFactor = 1.0 - pow(1.0 - fogFactor, 2.0);

	color *= exp(-GetWaterAbsorption() * eyeLength);

	// fogColor *= exp(-GetWaterAbsorption() * min(eyeLength * 1.0, 30.2));
	// fogColor *= exp(-GetWaterAbsorption() * 30.2);
	fogColor *= exp(-GetWaterAbsorption() * (1.0 - fogFactor) * 20.0);
	color = mix(fogColor, color, vec3(fogFactor));
}

vec3 calculateWaterFog(vec3 color, vec3 viewPos, vec3 viewPos2, float water, float ice) 
{
	if(water > 0.5 || isEyeInWater > 0 || ice > 0.5) 
	{
		float waterDepth = distance(viewPos.xyz, viewPos2.xyz);

		vec3 viewVector = normalize(viewPos.xyz);

		vec3 waterNormal = normalize(GetNormal(texcoord.st));

		if (isEyeInWater > 0)
		{
			waterDepth = length(viewPos.xyz) * 0.5;		
			if (water > 0.5 || ice > 0.5)
			{
				waterDepth = length(viewPos.xyz) * 0.5;		
			}	
		}

		float fogDensity = 0.05;

		vec3 waterFogColor = vec3(0.05, 0.8, 1.0);

		if (ice > 0.5)
		{
			waterFogColor = vec3(0.2, 0.6, 1.0) * 2.0;
			fogDensity = 0.7;
		}

		waterFogColor *= 0.01 * dot(vec3(0.33333), sunColor);
		waterFogColor *= isEyeInWater * 2.0 + 1.0;

		if(isEyeInWater > 0)
		{
			waterFogColor *= 0.5;
			//waterFogColor *= pow(eyeBrightnessSmooth.y / 240.0f, 6.0f);


			vec3 waterSunlightVector = refract(-lightVector, upVector, 1.0 / WATER_REFRACT_IOR);

			//waterFogColor *= (dot(lightVector, viewVector) * 0.5 + 0.5) * 2.0 + 1.0;
			float scatter = 1.0 / (pow(saturate(dot(waterSunlightVector, viewVector) * 0.5 + 0.5) * 30.0, 1.0) + 0.1);
			vec3 waterSunlightScatter = sunColor * scatter * 1.0 * waterFogColor * 0.000625;

			float eyeWaterDepth = eyeBrightnessSmooth.y / 240.0;


			waterFogColor *= dot(viewVector, upVector) * 0.5 + 0.5;
			waterFogColor = waterFogColor * pow(eyeWaterDepth, 1.0f) + waterSunlightScatter * pow(eyeWaterDepth, 1.0);

			waterFogColor *= pow(vec3(0.4, 0.72, 1.0) * 0.99, vec3(0.2 + (1.0 - eyeWaterDepth)));

			fogDensity *= 0.5;
		}

		float visibility = 1.0f / (pow(exp(waterDepth * fogDensity), 1.0f));

		vec3 viewVectorRefracted = refract(viewVector, waterNormal, 1.0 / 1.3333);
		float scatter = 1.0 / (pow(saturate(dot(-lightVector, viewVectorRefracted) * 0.5 + 0.5) * 20.0, 2.0) + 0.1);

		if (isEyeInWater < 1)
		{
			waterFogColor = mix(waterFogColor, sunColor * 0.000025 * waterFogColor, vec3(scatter));
		}

		color *= pow(vec3(0.4, 0.75, 1.0) * 0.99, vec3(waterDepth * 0.25 + 0.25));
		// color *= pow(vec3(0.7, 1.0, 0.2) * 0.8, vec3(waterDepth * 0.15 + 0.1));
		return mix(waterFogColor * 0.00005, color, saturate(visibility));
	} else {
		return color;
	}
}

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

float GetWaves(vec3 p, float lod) { return getwaves(p.xz * 0.1, Iterations2) * lod; }

vec3 GetWavesNormal(vec3 position, float lod) {
	float q = 0.05;
	vec3 w1 = vec3(q, GetWaves(position + vec3(q, 0.0, 0.0), lod), 0.0);
	vec3 w2 = vec3(0.0, GetWaves(position + vec3(0.0, 0.0, q), lod), q);
	vec3 w0 = GetWaves(position, lod) * vec3(0.0, 1.0, 0.0);

	vec3 tangent = w1 - w0;
	vec3 bitangent = w2 - w0;
	return normalize(cross(bitangent, tangent));
}

float fov = atan(1./gbufferProjection[1][1]);
float mulfov = bool(isEyeInWater) ? gbufferProjection[1][1]*tan(fov * 0.85):1.0;

vec4 fetch_vpos (vec2 uv, float z) {
	vec4 v = gbufferProjectionInverse * vec4(fma(vec3(uv, z), vec3(2.0f), vec3(-1.0)), 1.0);
	v /= v.w;
	v.xy *= mulfov;

	return v;
}

vec4 fetch_vpos (vec2 uv, sampler2D sam) {
	return fetch_vpos(uv, texture2D(sam, uv).x);
}

vec2 getScreenCoordByViewCoord(vec3 viewCoord) {
	vec4 p = vec4(viewCoord, 1.0);
	p = gbufferProjection * p;
	p /= p.w;
	if(p.z < -1 || p.z > 1)
		return vec2(-1.0);
	p = p * 0.5f + 0.5f;
	return p.st;
}

vec2 screen_project (vec3 vpos) {
	vec4 p = mat4(gbufferProjection) * vec4(vpos, 1.0f);
	p /= p.w;
	if(abs(p.z) > 1)
		return vec2(-1.0);
	return fma(p.st, vec2(0.5f), vec2(0.5f));
}

float linearizeDepth(float depth) {
    return (2.0 * near) / (far + near - depth * (far - near));
}

float getLinearDepthOfViewCoord(vec3 viewCoord) {
	vec4 p = vec4(viewCoord, 1.0);
	p = gbufferProjection * p;
	p /= p.w;
	return linearizeDepth(p.z * 0.5 + 0.5);
}

#define BISEARCH(SEARCHPOINT, DIRVEC, SIGN) DIRVEC *= 0.5; \
					SEARCHPOINT+= DIRVEC * SIGN; \
					uv = getScreenCoordByViewCoord(SEARCHPOINT); \
					sampleDepth = linearizeDepth(texture2DLod(depthtex0, uv, 0.0).x); \
					testDepth = getLinearDepthOfViewCoord(SEARCHPOINT); \
					SIGN = sign(sampleDepth - testDepth);

vec3 waterRayTarcing(vec3 startPoint, vec3 direction, vec3 color, float jitter, vec3 skyColor) {
	const float stepBase = 0.1;
	vec3 testPoint = startPoint;
	vec3 lastPoint = testPoint;
	vec3 reflection = direction;
	direction *= stepBase;
	bool hit = false;
	vec4 hitColor = vec4(0.0);
	bool land = false;
	float quality = 1./REFLECT_QUALITY;
	for(float i = 0.0; i < REFLECT_ISTEP; i += quality)
	{
		testPoint += direction * (i+jitter);
		vec2 uv = getScreenCoordByViewCoord(testPoint);
		if(uv.x < 0.0 || uv.x > 1.0 || uv.y < 0.0 || uv.y > 1.0)
		{
			hit = true;
			break;
		}
		float sampleDepth = GetDepth2(uv);
		sampleDepth = linearizeDepth(sampleDepth);
		float testDepth = getLinearDepthOfViewCoord(testPoint);
		if(sampleDepth < testDepth && testDepth - sampleDepth < (1.0 / 2048.0) * (1.0 + testDepth * 200.0 + float(i)))
		{
			uv = getScreenCoordByViewCoord(lastPoint);
			hitColor = vec4(decodeRGBE8(texture2D(colortex0, uv)), 1.0);

			float pixeldepth = texture2D(depthtex1, uv.st).x;

			//land = pixeldepth > 0.999;
			hit = true;
			break;
		}
		lastPoint = testPoint;
	}
	if(!hit)
	{
		vec2 uv = getScreenCoordByViewCoord(lastPoint);
		float testDepth = getLinearDepthOfViewCoord(lastPoint);
		float sampleDepth = texture2DLod(depthtex0, uv, 0.0).x;
		sampleDepth = linearizeDepth(sampleDepth);
		if(testDepth - sampleDepth < 0.5)
		{
			hitColor = vec4(decodeRGBE8(texture2D(colortex0, uv)), 1.0);
		}
	}

	//if(!hit) hitColor.rgb = skyColor;

	return mix(skyColor, hitColor.rgb, hitColor.a);
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

		coord += wind * 0.1;

		float noise = texture2D(noisetex, coord * vec2(0.65, 1.0) * 0.2).x * 0.45;	coord *= 3.0;
		    coord += vec2(noise * 0.6, 0.0) * octave_m;
		noise += texture2D(noisetex, coord * 0.4).x * 0.25; coord *= 3.01;
		    coord += vec2(noise * 0.4, 0.0) * octave_m + vec2(wind * 0.1, 0.2); coord += wind;
		noise += texture2D(noisetex, coord * 0.45).x * 0.105; coord *= 4.02;
		    coord += vec2(noise * 0.3, 0.0) * octave_m + vec2(wind * 0.03, 0.1);
		noise += texture2D(noisetex, coord * 0.25).x * 0.0625;
		noise /= 1.45;

	float coverage = 0.8f;
		coverage = mix(coverage, 0.97f, rainStrength);

	float dist = length(worldPos.xz * 0.5);
		coverage *= max(0.0f, 1.0f - dist / 14000.0f);
	float density = 0.1f + rainStrength * 0.3;

	noise = GetCoverage(coverage, density, noise);

	//noise = smoothstep(0.0, 1.0, noise + wetness);

	vec3 color = mix(sunColor, skyColor, 0.5);

    return vec4(color * 0.75f, noise);
}

vec3 GetCloudPlane(vec3 color, vec3 worldDir, float sky, float depth) 
{
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

vec3 calculateSpecularBRDF(vec3 normal, vec3 lightVector, vec3 viewVector, float f0, float alpha2, const float tanSunRadius){
	vec3 H = normalize(lightVector - viewVector);
	
	float VoH = clamp01(dot(H, lightVector));
	float NoL = clamp01(dot(normal, lightVector));
	float NoV = clamp01(dot(normal, -viewVector));
	float VoL = (dot(lightVector, -viewVector));
	//float NoH = calculateNoH(tanSunRadius, NoL, NoV, VoL);
	float NoH = clamp01(dot(normal, H));

	float D = GGXDistribution(alpha2, NoH);
	float G = GSpecular(alpha2, NoV, NoL);
	vec3 F = Fresnel(f0, 1.0, VoH);

	return max0(F * D * G / (4.0 * NoL * NoV)) * NoL;
}

vec3 GetWaterEffect(vec3 color, vec3 fragpos, vec3 wave_normal, vec3 normal, float water, vec3 worldPos, float roughness, float f0) 
{
	if (f0 < 0.005 || roughness >= 0.03) return color;

	if(roughness < 0.03) {
	//if(water > 0.5) {
		vec3 reflection = normalize(reflect(fragpos, wave_normal));
		vec3 reflectionwp = mat3(gbufferModelViewInverse) * (-reflection);
		float jitter = bayer128(gl_FragCoord.xy + vec2(sin(frameTimeCounter), cos(frameTimeCounter)) * 1024.0) * 1.0;

		float VoN = clamp01(dot(wave_normal, -normalize(fragpos)));

		//vec3 fresnel = vec3(0.0);
		float fresnel = pow(clamp01(dot(-reflection, wave_normal) + 1.0), 5.0) * 0.98 + 0.02;
		//float fresnel = SchlickFresnel(f0, 1.0, VoN);

		vec2 planetSphere = vec2(0.0);
        vec3 skyAbsorb = vec3(0.0);

		float vDotL = dot(reflection, sunVector);

		vec3 reflectSky = calculateAtmosphere(vec3(0.0), reflection, upVector, sunVector, moonVector, planetSphere, skyAbsorb, 5);

		reflectSky += calculateSunSpot(vDotL) * sunColorBase * skyAbsorb;
		reflectSky += calculateMoonSpot(-vDotL) * moonColorBase * skyAbsorb;
		reflectSky += calculateStars(reflectionwp, wMoonVector) * skyAbsorb;

		reflectSky += GetCloudPlane(reflectSky, reflectionwp, 1.0, 1.0);
		reflectSky *= 0.00005;
		//reflectSky = vec3(0.0);

		vec3 reflectColor = waterRayTarcing(fragpos, reflection.rgb, color, jitter, reflectSky);
		reflectColor.rgb = mix(color, reflectColor.rgb, fresnel);

		//reflectColor *= fresnel;

		return reflectColor;
	}
	return color;
}

void WaterParallax(inout vec3 wpos, float lod) {
	const int maxLayers = 4;
	
	#ifdef WATER_ENHANCE
	wpos.y -= 0.15;
	#else
	wpos.y -= 1.62;
	#endif

	vec3 stepin = vec3(0.0);
	vec3 nwpos = normalize(wpos);
	nwpos /= max(0.01, abs(nwpos.y));

	for (int i = 0; i < maxLayers; i++) {
		#ifdef WATER_ENHANCE
		float h = GetWaves(wpos + stepin + cameraPosition, lod) * 1.5;
		#else
		float h = GetWaves(wpos + stepin + cameraPosition, lod);
		#endif

		float diff = stepin.y - h;
		if (bool(isEyeInWater)) diff = -diff;
		stepin += nwpos * diff * 0.5;
	}
	wpos += stepin;
	#ifdef WATER_ENHANCE
	wpos.y += 0.15;
	#else
	wpos.y += 1.62;
	#endif
}

void getRoughnessF0(float data, out float roughness, out float f0){
	vec2 decodedData = decodeVec2(data);

	roughness = decodedData.x;
	f0 = decodedData.y;
}

/*
vec3 GetGILightling()
{
	const float iStep = 4.0f;
	const float jStep = 4.0f;
	float quality = 1./GI_BLUR_QUALITY;

	int total = 1;

	vec3 sample = vec3(0.0);

	for(float i = -iStep; i < iStep; i += quality)
	for(float j = -jStep; j < jStep; j += quality)
	{
		vec2 offset = mat2(cos(i), -sin(i), sin(i), cos(i)) * vec2(1.0) / 1024;
		sample += texture2D(colortex3, texcoord.st + offset).rgb;
		total++;
	}
	sample /= total;
	sample *= 35.0;

	return sample;
}
*/

void main() {
    surfaceStruct surface = GetSurfaceStruct();
	maskStruct mask = GetMaskStruct(texcoord.st);
    vec3 color = decodeRGBE8(texture2D(colortex0, texcoord.st));
	//color = decodeColor(color);

	vec4 data4 = texture2D(colortex4, texcoord.st);
	float roughness = data4.x, f0 = data4.y;
	//getRoughnessF0(texture2D(colortex4, texcoord.st).r, roughness, f0);

    //水面
	vec4 lvpos = fetch_vpos(texcoord.st, depthtex1);
	vec4 vpos = gbufferModelView * vec4(surface.worldPos.xyz, 1.0);
	vec3 nvpos = normalize(vpos.xyz);
	vec4 lnvpos = vec4(0.0);

	vec3 ref_col;

	vec3 water_plain_normal = mat3(gbufferModelViewInverse) * surface.normal;
	if (bool(isEyeInWater)) water_plain_normal = -water_plain_normal;

	float lod = pow(max(water_plain_normal.y, 0.0), 4.0);

	//水体视察
	#ifdef WATER_PARALLAX
		if (lod > 0.0) WaterParallax(surface.worldPos.xyz, lod);
		float wave = GetWaves(surface.worldPos.xyz + cameraPosition, lod);
	#else
		float wave = GetWaves(surface.worldPos.xyz + cameraPosition, lod);
		vec2 p = vpos.xy / vpos.z * wave;
		vec2 wp = length(p) * normalize(surface.worldPos.xyz + cameraPosition).xz * 0.1;
		surface.worldPos.xyz -= vec3(wp.x, 0.0, wp.y);
	#endif
	
	vec3 wnormal = GetWavesNormal(surface.worldPos.xyz + cameraPosition, lod);
	vec3 vsnormal = mat3(gbufferModelView) * wnormal;

	//折射
	if(mask.water > 0.5 && isEyeInWater != 1) {
		vec3 refract_vpos = refract(nvpos, vsnormal, 1.00029 / 1.33);
		//vec3 refract_vpos = vsnormal;
		vec2 uv = screen_project(refract_vpos + surface.viewPos.xyz);
		uv = mix(uv, texcoord.st, pow(abs(uv - vec2(0.5)) * 2.0, vec2(2.0)));
		//uv = vsnormal.st;

		lnvpos.xyz = fetch_vpos(uv, depthtex1).xyz;
		float cdepth = length(lnvpos.xyz);
		lnvpos.xyz = lnvpos.xyz / cdepth;
		float cdepthN = cdepth / far;

		ref_col = decodeRGBE8(texture2D(colortex0, uv, 0)) * 0.2;
		ref_col += decodeRGBE8(texture2D(colortex0, uv, 1)) * 0.3;
		ref_col += decodeRGBE8(texture2D(colortex0, uv, 2)) * 0.5;

		color = ref_col;
    }

	//水雾
    color = calculateWaterFog(color, surface.viewPos.xyz, surface.viewPos2.xyz, mask.water, mask.ice);
	//void UnderwaterFog(inout vec3 color, float eyeLength, vec3 eyeDir, vec3 skyColor, vec3 sunColor)
	//UnderwaterFog(color, length(surface.viewPos.xyz) * 0.5, normalize(surface.viewPos.xyz), skyColor, sunColor);

	//水反射
	color = GetWaterEffect(color, surface.viewPos.xyz, vsnormal, surface.normal, mask.water, surface.worldPos.xyz, roughness, f0);

/* DRAWBUFFERS:0 */
    gl_FragData[0] = vec4(color, 1.0);
}