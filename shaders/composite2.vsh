#version 120

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

uniform mat4 gbufferModelView;
uniform mat4 gbufferModelViewInverse;

uniform mat4 shadowModelView;
uniform mat4 shadowProjection;

uniform vec3 sunPosition;
uniform vec3 upPosition;

uniform int worldTime;
uniform float eyeAltitude;

uniform float wetness;
uniform float rainStrength;

const float PI 		= acos(-1.0);
const float TAU 	= PI * 2.0;
const float hPI 	= PI * 0.5;
const float rPI 	= 1.0 / PI;
const float rTAU 	= 1.0 / TAU;

const float PHI		= sqrt(5.0) * 0.5 + 0.5;
const float rLOG2	= 1.0 / log(2.0);

const float goldenAngle = TAU / PHI / PHI;

#define clamp01(x) clamp(x, 0.0, 1.0)
#define max0(x) max(x, 0.0)
#define min0(x) min(x, 0.0)
#define max3(a) max(max(a.x, a.y), a.z)
#define min3(a) min(min(a.x, a.y), a.z)
#define max4(a, b, c, d) max(max(a, b), max(c, d))
#define min4(a, b, c, d) min(min(a, b), min(c, d))

#define fsign(x) (clamp01(x * 1e35) * 2.0 - 1.0)
#define fstep(x,y) clamp01((y - x) * 1e35)

#define HASHSCALE1 443.8975
#define HASHSCALE3 vec3(443.897, 441.423, 437.195)
#define HASHSCALE4 vec3(443.897, 441.423, 437.195, 444.129)

//  1 out, 3 in...
float hash13(vec3 p3)
{
	p3  = fract(p3 * HASHSCALE1);
    p3 += dot(p3, p3.yzx + 19.19);
    return fract((p3.x + p3.y) * p3.z);
}

vec3 hash33(vec3 p){
    p = fract(p * HASHSCALE3);
    p += dot(p.zxy, p.yxz + 19.19);
    return fract(vec3(p.x * p.y, p.z * p.x, p.y * p.z));
}

vec3 srgbToLinear(vec3 srgb){
    return mix(
        srgb / 12.92,
        pow(.947867 * srgb + .0521327, vec3(2.4) ),
        step( .04045, srgb )
    );
}

vec3 linearToSRGB(vec3 linear){
    return mix(
        linear * 12.92,
        pow(linear, vec3(1./2.4) ) * 1.055 - .055,
        step( .0031308, linear )
    );
}

vec3 blackbody(float Temp)
{
    float t = pow(Temp, -1.5);
    float lt = log(Temp);

    vec3 col = vec3(0.0);
         col.x = 220000.0 * t + 0.58039215686;
         col.y = 0.39231372549 * lt - 2.44549019608;
         col.y = Temp > 6500. ? 138039.215686 * t + 0.72156862745 : col.y;
         col.z = 0.76078431372 * lt - 5.68078431373;
         col = clamp01(col);
         col = Temp < 1000. ? col * Temp * 0.001 : col;

    return srgbToLinear(col);
}

mat3 getRotMat(vec3 x,vec3 y){
    float d = dot(x,y);
    vec3 cr = cross(y,x);
    
    float s = length(cr);
    
    float id = 1.-d;
    
    vec3 m = cr/s;
    
    vec3 m2 = m*m*id+d;
    vec3 sm = s*m;
    
    vec3 w = (m.xy*id).xxy*m.yzz;
    
    return mat3(
        m2.x,     w.x-sm.z, w.y+sm.y,
        w.x+sm.z, m2.y,     w.z-sm.x, 
        w.y-sm.y, w.z+sm.x, m2.z
    );
}

// No intersection if returned y component is < 0.0
vec2 rsi(vec3 position, vec3 direction, float radius) {
	float PoD = dot(position, direction);
	float radiusSquared = radius * radius;

	float delta = PoD * PoD + radiusSquared - dot(position, position);
	if (delta < 0.0) return vec2(-1.0);
	      delta = sqrt(delta);

	return -PoD + vec2(-delta, delta);
}

#include "/libs/sky.glsl"

vec4 ToSH(float value, vec3 dir)
{
	const float PI = 3.14159265359;
	const float N1 = sqrt(4 * PI / 3);
	const float transferl1 = (sqrt(PI) / 3.0) * N1;
	const float transferl0 = PI;

	const float sqrt1OverPI = sqrt(1.0 / PI);
	const float sqrt3OverPI = sqrt(3.0 / PI);

	vec4 coeffs;

	coeffs.x = 0.5 * sqrt1OverPI * value * transferl0;
	coeffs.y = -0.5 * sqrt3OverPI * dir.y * value * transferl1;
	coeffs.z = 0.5 * sqrt3OverPI * dir.z * value * transferl1;
	coeffs.w = -0.5 * sqrt3OverPI * dir.x * value * transferl1; //TODO: Vectorize the math so it's faster

	return coeffs;
}


vec3 FromSH(vec4 cR, vec4 cG, vec4 cB, vec3 lightDir)
{
	const float PI = 3.14159265;

	const float N1 = sqrt(4 * PI / 3);
	const float transferl1 = (sqrt(PI) / 3.0) * N1;
	const float transferl0 = PI;

	const float sqrt1OverPI = sqrt(1.0 / PI);
	const float sqrt3OverPI = sqrt(3.0 / PI);

	vec4 sh;

	sh.x = 0.5 * sqrt1OverPI;
	sh.y = -0.5 * sqrt3OverPI * lightDir.y;
	sh.z = 0.5 * sqrt3OverPI * lightDir.z;
	sh.w = -0.5 * sqrt3OverPI * lightDir.x;

	vec3 result;
	result.r = sh.x * cR.x;
	result.r += sh.y * cR.y;
	result.r += sh.z * cR.z;
	result.r += sh.w * cR.w;

	result.g = sh.x * cG.x;
	result.g += sh.y * cG.y;
	result.g += sh.z * cG.z;
	result.g += sh.w * cG.w;

	result.b = sh.x * cB.x;
	result.b += sh.y * cB.y;
	result.b += sh.z * cB.z;
	result.b += sh.w * cB.w;

	return result.rgb;
}

void CalculateSkySH(vec3 sunVector, vec3 moonVector, vec3 upVector, vec2 planetSphere, vec3 transmittance) {
	const int latSamples = 16;
	const int lonSamples = 5;
	const float rLatSamples = 1.0 / latSamples;
	const float rLonSamples = 1.0 / lonSamples;
	const float sampleCount = rLatSamples * rLonSamples;

	const float latitudeSize = rLatSamples * PI;
	const float longitudeSize = rLonSamples * TAU;

	vec4 shR = vec4(0.0), shG = vec4(0.0), shB = vec4(0.0);

	for (int i = 0; i < latSamples; ++i) {
		float latitude = float(i) * latitudeSize;

		for (int j = 0; j < lonSamples; ++j) {
			float longitude = float(j) * longitudeSize;

			float c = cos(latitude);
			vec3 kernel = vec3(c * cos(longitude), sin(latitude), c * sin(longitude));

			vec3 skyCol = calculateAtmosphere(vec3(0.0), normalize(kernel + vec3(0.0, 0.1, 0.0)), upVector, sunVector, moonVector, planetSphere, transmittance, 10);

			shR += ToSH(skyCol.r, kernel);
			shG += ToSH(skyCol.g, kernel);
			shB += ToSH(skyCol.b, kernel);
		}
	}

	skySH = mat3x4(shR, shG, shB) * sampleCount;
}

void main() {
    gl_Position = ftransform();
    texcoord = gl_MultiTexCoord0;

	const float tTime = (1.0 / 50.0);
	float wTime = float(worldTime);
	transitionFading = clamp01(clamp01((wTime - 23215.0) * tTime) + (1.0 - clamp01((wTime - 12735.0) * tTime)) + clamp01((wTime - 12925.0) * tTime) * (1.0 - clamp01((wTime - 22975.0) * tTime)));

	upVector = upPosition * 0.01;
	sunVector = sunPosition * 0.01;
	moonVector = -sunVector;

	wSunVector = mat3(gbufferModelViewInverse) * sunVector;
	wMoonVector = mat3(gbufferModelViewInverse) * moonVector;

	lightVector = (worldTime > 22975 || worldTime < 12925 ? sunVector : moonVector);
	//float timeFactor = smoothstep(0.0, 1.0, float(worldTime > 22975 || worldTime < 12925));
    //lightVector = mix(moonVector, sunVector, timeFactor);
	wLightVector = mat3(gbufferModelViewInverse) * lightVector;

	baseSunColor = sunColorBase;
	baseMoonColor = moonColorBase;

	sunColor = sky_transmittance(vec3(0.0, sky_planetRadius, 0.0), wSunVector, 8) * baseSunColor;
	moonColor = sky_transmittance(vec3(0.0, sky_planetRadius, 0.0), wMoonVector, 8) * baseMoonColor;
	// vec2 planetSphere = vec2(0.0);
	// vec3 transmittance = vec3(0.0);

	// skyColor = vec3(0.0);
	// skyColor = calculateAtmosphere(vec3(0.0), vec3(0.0, 1.0, 0.0), vec3(0.0, 1.0, 0.0), wSunVector, wMoonVector, planetSphere, transmittance, 10);

    //     const int latSamples = 5;
    //     const int lonSamples = 5;
    //     const float rLatSamples = 1.0 / latSamples;
    //     const float rLonSamples = 1.0 / lonSamples;
    //     const float sampleCount = rLatSamples * rLonSamples;

    //     const float latitudeSize = rLatSamples * PI;
    //     const float longitudeSize = rLonSamples * TAU;

    //     vec4 shR = vec4(0.0), shG = vec4(0.0), shB = vec4(0.0);
    //     const float offset = 0.1;

    //     for (int i = 0; i < latSamples; ++i) {
    //         float latitude = float(i) * latitudeSize;

    //         for (int j = 0; j < lonSamples; ++j) {
    //             float longitude = float(j) * longitudeSize;

    //             float c = cos(latitude);
    //             vec3 kernel = vec3(c * cos(longitude), sin(latitude), c * sin(longitude));

    //             vec3 skyCol = calculateAtmosphere(vec3(0.0), mat3(gbufferModelView) * normalize(kernel + vec3(0.0, offset, 0.0)), gbufferModelView[1].xyz, sunPosition * 0.01, -sunPosition * 0.01, planetSphere, transmittance, 4);

    //             shR += ToSH(skyCol.r, kernel);
    //             shG += ToSH(skyCol.g, kernel);
    //             shB += ToSH(skyCol.b, kernel);
    //         }
    //     }

    //     // fun fact, when integrating a point function over a sphere, you need to weight your sample to the size of the area (Solid Angle)
    //     // Each sample covers 2 * PI / NUM_RANDOM_SPHERE_SAMPLES solid angle  
    //     skySH = mat3x4(shR, shG, shB) * TAU * sampleCount; 

    vec2 planetSphere = vec2(0.0);
	vec3 transmittance = vec3(0.0);

	skyColor = vec3(0.0);
	skyColor = calculateAtmosphere(vec3(0.0), vec3(0.0, 1.0, 0.0), vec3(0.0, 1.0, 0.0), wSunVector, wMoonVector, planetSphere, transmittance, 10);

	CalculateSkySH(wSunVector, wMoonVector, vec3(0.0, 1.0, 0.0), planetSphere, transmittance);

}