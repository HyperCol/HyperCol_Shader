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

varying float transitionFading;

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
	wLightVector = mat3(gbufferModelViewInverse) * lightVector;

	baseSunColor = sunColorBase;
	baseMoonColor = moonColorBase;

	sunColor = sky_transmittance(vec3(0.0, sky_planetRadius, 0.0), wSunVector, 8) * baseSunColor;
	moonColor = sky_transmittance(vec3(0.0, sky_planetRadius, 0.0), wMoonVector, 8) * baseMoonColor;
	vec2 planetSphere = vec2(0.0);
	vec3 transmittance = vec3(0.0);

	skyColor = vec3(0.0);
	skyColor = calculateAtmosphere(vec3(0.0), vec3(0.0, 1.0, 0.0), vec3(0.0, 1.0, 0.0), wSunVector, wMoonVector, planetSphere, transmittance, 10);

}