#define SHADOW_MAP_BIAS 0.9

#define clamp01(x) clamp(x, 0.0, 1.0)
#define max0(x) max(x, 0.0)
#define min0(x) min(x, 0.0)
#define max3(a) max(max(a.x, a.y), a.z)
#define min3(a) min(min(a.x, a.y), a.z)
#define max4(a, b, c, d) max(max(a, b), max(c, d))
#define min4(a, b, c, d) min(min(a, b), min(c, d))

#define fsign(x) (clamp01(x * 1e35) * 2.0 - 1.0)
#define fstep(x,y) clamp01((y - x) * 1e35)

#define pow2(a) a*a
#define pow3(a) pow2(a)*a
#define pow4(a) pow3(a)*a
#define pow5(a) pow4(a)*a
#define pow6(a) pow5(a)*a
#define pow7(a) pow6(a)*a
#define pow8(a) pow7(a)*a

#define encodeColor(x) (x * 0.00005)
#define decodeColor(x) (x * 20000.0)

const int 		superSamplingLevel 		        = 0;

const float 	wetnessHalflife 			    = 1.0;
const float 	drynessHalflife 			    = 60.0;

const float		sunPathRotation					= -40.0; //[-50.0 -40.0 -30.0 -20.0 -10.0 0.0 10.0 20.0 30.0 40.0 50.0]

const float 	ambientOcclusionLevel 	        = 0.0f;

const float		eyeBrightnessHalflife		    = 10.0; //[1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0 11.0 12.0 13.0 14.0]

const int 		shadowMapResolution 	        = 1516;     //[512 1024 1516 2048 4096 16384]
const float 	shadowDistance 			        = 200.0f;   //[80.0f 100.0f 120.0f 150.0f 180.0f 200.0f 240.0f 280.0f 300.0f 400.0f 500.0f 600.0f 700.0f 800.0f 1000.0f 1200.0f 1500.0f 1600.0f]
const float 	shadowIntervalSize 		        = 4.0f;
const bool 		shadowHardwareFiltering 	    = true;
const bool 		shadowHardwareFiltering0        = true;
const bool 		shadowHardwareFiltering1        = true;

const bool 		shadowtex1Mipmap                = true;
const bool 		shadowtex1Nearest               = false;
const bool 		shadowcolor0Mipmap              = true;
const bool 		shadowcolor0Nearest             = false;
const bool 		shadowcolor1Mipmap              = true;
const bool 		shadowcolor1Nearest             = false;

const float PI 		= acos(-1.0f);
const float TAU 	= PI * 2.0;
const float hPI 	= PI * 0.5;
const float rPI 	= 1.0 / PI;
const float rTAU 	= 1.0 / TAU;

const float PHI		= sqrt(5.0) * 0.5 + 0.5;
const float rLOG2	= 1.0 / log(2.0);

const float goldenAngle = TAU / PHI / PHI;

// Faster alternative to acos when exact results are not needed.
// Intersects acos at: -1, -0.57985216912, 0, 0.57985216912, 1
// Max absolute error: 0.00426165254821 (~0.244 degrees)
// Max absolute error can be reduced if you are fine with it being discontinuous at 0 (i'm not)
float FastAcos(float x) {
	// abs() is free on input for GCN
	// I'm assuming that this is true for most modern GPUs.
	float r = sqrt(1.0 - abs(x)) * (-0.175394 * abs(x) + hPI);
	return x < 0.0 ? PI - r : r;
}
// One sub slower than facos but aside from that has the same properties.
float FastAsin(float x) {
	return hPI - FastAcos(x);
}

vec4 Cubic(float x) {
    float x2 = x * x;
    float x3 = x2 * x;

    vec4 w   = vec4(0.0);
         w.x =       -x3 + 3.0 * x2 - 3.0 * x + 1.0;
         w.y =  3.0 * x3 - 6.0 * x2           + 4.0;
         w.z = -3.0 * x3 + 3.0 * x2 + 3.0 * x + 1.0;
         w.w = x3;

    return w * 0.166666666667;
}

vec4 BicubicTexture(sampler2D tex, vec2 coord) {
    vec2 resolution = vec2(viewWidth, viewHeight);

    coord *= resolution;

    vec2 f = fract(coord);

    resolution = 1.0 / resolution;

    coord -= f;

    vec4 xCubic = Cubic(f.x);
    vec4 yCubic = Cubic(f.y);

    vec4 s = vec4(xCubic.xz + xCubic.yw, yCubic.xz + yCubic.yw);
    vec4 offset = coord.xxyy + vec4(-.5, 1.5, -.5, 1.5) + vec4(xCubic.yw, yCubic.yw) / s;

    vec4 sample0 = texture2D(tex, offset.xz * resolution);
    vec4 sample1 = texture2D(tex, offset.yz * resolution);
    vec4 sample2 = texture2D(tex, offset.xw * resolution);
    vec4 sample3 = texture2D(tex, offset.yw * resolution);

    float sx = s.x / (s.x + s.y);
    float sy = s.z / (s.z + s.w);

    return mix(mix(sample3, sample2, sx), mix(sample1, sample0, sx), sy);
}

float Luminance(in vec3 color)
{
	return dot(color.rgb, vec3(0.2125f, 0.7154f, 0.0721f));
}

float saturate(float x)
{
	return clamp(x, 0.0, 1.0);
}

vec3 saturate(vec3 x)
{
	return clamp(x, vec3(0.0), vec3(1.0));
}

vec2 saturate(vec2 x)
{
	return clamp(x, vec2(0.0), vec2(1.0));
}

float remap(float e0, float e1, float x)
{
	return saturate((x - e0) / (e1 - e0));
}

float facos(const float sx){
    float x = clamp(abs( sx ),0.,1.);
    float a = sqrt( 1. - x ) * ( -0.16882 * x + 1.56734 );
    return sx > 0. ? a : PI - a;
    //float c = clamp(-sx * 1e35, 0., 1.);
    //return c * pi + a * -(c * 2. - 1.); //no conditional version
}

vec2 sincos(float x){
    return vec2(sin(x), cos(x));
}

vec3 toGamma(vec3 color) {
    float gamma = 2.2;
    return pow(color, vec3(gamma));
}

vec3 formGamma(vec3 color) {
    float agamma = 1.0 / 2.2;
    return pow(color, vec3(agamma));
}

vec2 toGamma(vec2 color) {
    float gamma = 2.2;
    return pow(color, vec2(gamma));
}

vec2 formGamma(vec2 color) {
    float agamma = 1.0 / 2.2;
    return pow(color, vec2(agamma));
}

float encodeVec2(vec2 a) {
	const vec2 constant1 = vec2(1.0, 256.0) / 65535.0; //2^16-1
	return dot(floor(a * 255.0), constant1);
}

float encodeVec2(float x,float y) {
	return encodeVec2(vec2(x,y));
}

vec2 decodeVec2(float a) {
	const vec2 constant1 = 65535.0 / vec2(256.0, 65536.0);
	const float constant2 = 256.0 / 255.0;
	return fract(a * constant1) * constant2;
}

float encodeNormal(vec3 a) {
    vec3 b = abs(a);
    vec2 p = a.xy / dot(b, vec3(1.0));
    vec2 encoded = a.z <= 0. ? (1. - abs(p.yx)) * fsign(p) : p;
    encoded = encoded * .5 + .5;
	return encodeVec2(encoded);
}

vec3 decodeNormal(float encoded) {
	vec2 a = decodeVec2(encoded);
	     a = a * 2.0 - 1.0;
	vec2 b = abs(a);
	float z = 1.0 - b.x - b.y;

	return normalize(vec3(z < 0.0 ? (1.0 - b.yx) * fsign(a) : a, z));
}

vec3 decodeNormal(float encoded, mat4 gbufferModelView) {
	return mat3(gbufferModelView) * decodeNormal(encoded);
}

// vec2 encodeNormal(vec3 a) {
//     vec3 b = abs(a);
//     vec2 p = a.xy / dot(b, vec3(1.0));
//     vec2 encoded = a.z <= 0. ? (1. - abs(p.yx)) * fsign(p) : p;
// 	return encoded * .5 + .5;
// }

// vec3 decodeNormal(vec2 encoded) {
// 	vec2 a = encoded * 2.0 - 1.0;
// 	vec2 b = abs(a);
// 	float z = 1.0 - b.x - b.y;

// 	return normalize(vec3(z < 0.0 ? (1.0 - b.yx) * fsign(a) : a, z));
// }

// vec3 decodeNormal(vec2 encoded, mat4 gbufferModelView) {
// 	return mat3(gbufferModelView) * decodeNormal(encoded);
// }

vec4 encodeRGBE8(vec3 rgb) {
    float exponentPart = floor(log2(max3(rgb)));
    vec3  mantissaPart = clamp01((128.0 / 255.0) * rgb * exp2(-exponentPart));
          exponentPart = clamp01((exponentPart + 127.0) * (1.0 / 255.0));

    return vec4(mantissaPart, exponentPart);
}

vec3 decodeRGBE8(vec4 rgbe) {
    float exponentPart = exp2(rgbe.a * 255.0 - 127.0);
    vec3  mantissaPart = (510.0 / 256.0) * rgbe.rgb;

    return exponentPart * mantissaPart;
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

vec2 rotate(vec2 x, float r){
    vec2 sc = sincos(r);
    return mat2(sc.x, -sc.y, sc.y, sc.x) * x;
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

mat2 rotation(float x) {
    vec2 sc = sincos(x);
    return mat2(sc.y, sc.x, -sc.x, sc.y);
}

float GeometrySchlickGGX(float NdotV, float k){
    float denom = NdotV * (1.0 - k) + k;

    return NdotV / denom;
}

vec3 GeometrySmithGGX(vec3 diffuseColor, vec3 N, vec3 V, vec3 L, float r){
	float k = pow2(r + 1.0) * 0.125;
    float NdotL = clamp01(dot(N, L));

    float multiScattering = 0.1159 * r;

    return (diffuseColor * multiScattering * NdotL + GeometrySchlickGGX(NdotL, k)) * rPI;
}

float GSpecular(float alpha2, float NoV, float NoL) {
    float x = 2.0 * NoL * NoV;
    float y = (1.0 - alpha2);

    return x / (NoV * sqrt(alpha2 + y * (NoL * NoL)) + NoL * sqrt(alpha2 + y * (NoV * NoV)));
}

float GGXDistribution(const float alpha2, const float NoH) {
	float d = (NoH * alpha2 - NoH) * NoH + 1.0;

	return alpha2 / (PI * d * d);
}

float SchlickFresnel(float f0, float f90, float LoH) {
    return (f90 - f0) * pow5(1. - LoH) + f0;
}

vec3 ExactFresnel(const vec3 n, const vec3 k, float c) {
    const vec3 k2= k * k;
	const vec3 n2k2 = n * n + k2;

    vec3 c2n = (c * 2.0) * n;
    vec3 c2 = vec3(c * c);

    vec3 rs_num = n2k2 - c2n + c2;
    vec3 rs_den = n2k2 + c2n + c2;

    vec3 rs = rs_num / rs_den;

    vec3 rp_num = n2k2 * c2 - c2n + 1.0;
    vec3 rp_den = n2k2 * c2 + c2n + 1.0;

    vec3 rp = rp_num / rp_den;

    return clamp01(0.5 * (rs + rp));
}

vec3 Fresnel(float f0, float f90, float LoH) {
    /*
    if(f0 > 0.985) {
		const vec3 chromeIOR = vec3(3.1800, 3.1812, 2.3230);
        const vec3 chromeK = vec3(3.3000, 3.3291, 3.1350);
 		return ExactFresnel(chromeIOR, chromeK, LoH);
	} else if(f0 > 0.965) {
        const vec3 goldIOR = vec3(0.18299, 0.42108, 1.3734);
        const vec3 goldK = vec3(3.4242, 2.3459, 1.7704);
         return ExactFresnel(goldIOR, goldK, LoH);
    } else if(f0 > 0.45) {
        const vec3 ironIOR = vec3(2.9114, 2.9497, 2.5845);
        const vec3 ironK = vec3(3.0893, 2.9318, 2.7670);
         return ExactFresnel(ironIOR, ironK, LoH);
    } else { */
        return vec3(SchlickFresnel(f0, f90, LoH));
    //}
}

float GetDepth(vec2 coord) {
    return texture2D(depthtex0, coord).x;
}

float GetDepth2(vec2 coord) {
    return texture2D(depthtex1, coord).x;
}

vec3 GetNormal(vec2 coord) {
    return decodeNormal(texture2D(colortex2, coord).x, gbufferModelView);
}

vec2 GetLightmap(vec2 coord) {
    vec2 lightmaps = texture2D(colortex2, coord).yz;
	return vec2(lightmaps.x, lightmaps.y * lightmaps.y);
}

vec4 GetScreenPosition(vec2 coord) {
    vec4 uv = vec4(vec3(coord, GetDepth(coord)) * 2.0 - 1.0, 1.0);
    vec4 screenPosition = gbufferProjectionInverse * uv;
    return screenPosition / screenPosition.w;
}

vec4 GetScreenPosition(vec2 coord, float depth) {
    vec4 uv = vec4(vec3(coord, depth) * 2.0 - 1.0, 1.0);
    vec4 screenPosition = gbufferProjectionInverse * uv;
    return screenPosition / screenPosition.w;
}

vec4 GetWorldPosition(vec2 coord) {
    vec4 screenPosition = GetScreenPosition(coord);
    return gbufferModelViewInverse * screenPosition;
}

vec4 GetWorldPosition(vec2 coord, float depth) {
    vec4 screenPosition = GetScreenPosition(coord, depth);
    return gbufferModelViewInverse * screenPosition;
}

vec3 GetShadowPosition(vec3 worldPostion) {
	vec4 shadowPosition = shadowModelView * vec4(worldPostion, 1.0);
	shadowPosition = shadowProjection * shadowPosition;
    shadowPosition /= shadowPosition.w;
	float dist = sqrt(shadowPosition.x * shadowPosition.x + shadowPosition.y * shadowPosition.y);
	float distortFactor = (1.0f - SHADOW_MAP_BIAS) + dist * SHADOW_MAP_BIAS;
	shadowPosition.xy *= 1.0f / distortFactor;
    shadowPosition.z = mix(shadowPosition.z, 0.5, 0.8);
    return shadowPosition.xyz * 0.5 + 0.5;
}

vec4 ToSH(float value, vec3 dir) {
    const float transferl1 = 0.3849 * PI;
    const float sqrt1OverPI = sqrt(rPI);
    const float sqrt3OverPI = sqrt(rPI * 3.0);

    const vec2 halfnhalf = vec2(0.5, -0.5);
    const vec2 transfer = vec2(PI * sqrt1OverPI, transferl1 * sqrt3OverPI);

    const vec4 foo = halfnhalf.xyxy * transfer.xyyy;

    return foo * vec4(1.0, dir.yzx) * value;
}

vec3 FromSH(vec4 cR, vec4 cG, vec4 cB, vec3 lightDir) {
    const float sqrt1OverPI = sqrt(rPI);
    const float sqrt3OverPI = sqrt(3.0 * rPI);
    const vec2 halfnhalf = vec2(0.5, -0.5);
    const vec2 sqrtOverPI = vec2(sqrt1OverPI, sqrt3OverPI);
    const vec4 foo = halfnhalf.xyxy * sqrtOverPI.xyyy;

    vec4 sh = foo * vec4(1.0, lightDir.yzx);

    // know to work
    return vec3(
        dot(sh,cR),
        dot(sh,cG),
        dot(sh,cB)
    );
}

struct surfaceStruct
{
    float depth;
    float depth2;
    vec3 normal;

    vec4 viewPos;
    vec4 viewPos2;
    vec3 viewDir;
    vec4 worldPos;
    vec4 worldPos2;
    vec3 worldDir;
    vec3 worldDir2;
    vec3 shadowPos;
	vec3 lightVector;

    vec3 worldNormal;

    vec2 lightmap;
};

surfaceStruct GetSurfaceStruct() {
    surfaceStruct a;
    a.depth = GetDepth(texcoord.st);
    a.depth2 = GetDepth2(texcoord.st);
    a.normal = GetNormal(texcoord.st);

    a.viewPos = GetScreenPosition(texcoord.st);
    a.viewPos2 = GetScreenPosition(texcoord.st, a.depth2);
    a.viewDir = normalize(a.viewPos.xyz);
    a.worldPos = GetWorldPosition(texcoord.st);
    a.worldPos2 = GetWorldPosition(texcoord.st, a.depth2);
    a.worldDir = normalize(a.worldPos.xyz);
    a.worldDir2 = normalize(a.worldPos2.xyz);
    a.shadowPos = GetShadowPosition(a.worldPos2.xyz);

	a.lightVector = normalize((gbufferModelView * vec4(normalize((shadowModelViewInverse * vec4(0.0, 0.0, 1.0, 0.0)).xyz), 0.0)).xyz);

    a.worldNormal = normalize((gbufferModelViewInverse * vec4(a.normal, 0.0)).xyz);

    a.lightmap = GetLightmap(texcoord.st);

    return a;
}

struct maskStruct
{
	float matIDs;

	float sky;
	float land;
	float grass;
	float leaves;
	float ice;
	float hand;
	float translucent;
	float glow;
	float sunspot;
	float goldBlock;
	float ironBlock;
	float diamondBlock;
	float emeraldBlock;
	float sand;
	float sandstone;
	float stone;
	float cobblestone;
	float wool;

	float torch;
	float lava;
	float glowstone;
	float fire;

	float water;

	float stainedGlass;
};

float GetMaskIDs(vec2 coord) {
    return texture2D(colortex4, coord).z;
}

float 	GetMaterialMask(in vec2 coord ,const in int ID, in float matID) {
	matID = (matID * 255.0f);

	//Catch last part of sky
	if (matID > 254.0f) {
		matID = 0.0f;
	}

	if (matID == ID) {
		return 1.0f;
	} else {
		return 0.0f;
	}
}

float  	GetWaterMask(in vec2 coord, in float matID) {					//Function that returns "true" if a pixel is water, and "false" if a pixel is not water.
	matID = (matID * 255.0f);

	if (matID >= 35.0f && matID <= 51) {
		return 1.0f;
	} else {
		return 0.0f;
	}
}

float  	GetStainedGlassMask(in vec2 coord, in float matID) {					//Function that returns "true" if a pixel is water, and "false" if a pixel is not water.
	matID = (matID * 255.0f);

	if (matID >= 55.0f && matID <= 70.0f) {
		return 1.0f;
	} else {
		return 0.0f;
	}
}

float  	GetIceMask(in vec2 coord, in float matID) {					//Function that returns "true" if a pixel is water, and "false" if a pixel is not water.
	matID = (matID * 255.0f);

	if (matID >= 5.0f && matID <= 5.5f) {
		return 1.0f;
	} else {
		return 0.0f;
	}
}

float  	GetTranslucentMask(in vec2 coord, in float matID) {
	matID = (matID * 255.0f);

	if (matID >= 4.0f && matID <= 4.5f) {
		return 1.0f;
	} else {
		return 0.0f;
	}
}

maskStruct GetMaskStruct(vec2 coord) {
    maskStruct mask;
    mask.matIDs = GetMaskIDs(coord.st);

		mask.sky 			= GetMaterialMask(coord.st, 0, mask.matIDs);
        //mask.sky 			= float(GetDepth2(coord.st) >= 1);

		mask.land	 		= 1.0 - mask.sky;
		mask.grass 			= GetMaterialMask(coord.st, 2, mask.matIDs);
		mask.leaves	 		= GetMaterialMask(coord.st, 3, mask.matIDs);
		mask.hand	 		= GetMaterialMask(coord.st, 5, mask.matIDs);
//		mask.translucent	= GetMaterialMask(coord.st, 6, mask.matIDs);

		mask.glow	 		= GetMaterialMask(coord.st, 10, mask.matIDs);
		mask.sunspot 		= GetMaterialMask(coord.st, 11, mask.matIDs);

		mask.goldBlock 		= GetMaterialMask(coord.st, 20, mask.matIDs);
		mask.ironBlock 		= GetMaterialMask(coord.st, 21, mask.matIDs);
		mask.diamondBlock	= GetMaterialMask(coord.st, 22, mask.matIDs);
		mask.emeraldBlock	= GetMaterialMask(coord.st, 23, mask.matIDs);
		mask.sand	 		= GetMaterialMask(coord.st, 24, mask.matIDs);
		mask.sandstone 		= GetMaterialMask(coord.st, 25, mask.matIDs);
		mask.stone	 		= GetMaterialMask(coord.st, 26, mask.matIDs);
		mask.cobblestone	= GetMaterialMask(coord.st, 27, mask.matIDs);
		mask.wool			= GetMaterialMask(coord.st, 28, mask.matIDs);

		mask.torch 			= GetMaterialMask(coord.st, 30, mask.matIDs);
		mask.lava 			= GetMaterialMask(coord.st, 31, mask.matIDs);
		mask.glowstone 		= GetMaterialMask(coord.st, 32, mask.matIDs);
		mask.fire 			= GetMaterialMask(coord.st, 33, mask.matIDs);

		mask.water 			= GetWaterMask(coord.st, mask.matIDs);
		mask.stainedGlass 	= GetStainedGlassMask(coord.st, mask.matIDs);
		mask.ice		 	= GetIceMask(coord.st, mask.matIDs);
        mask.water          = (isEyeInWater > 0) ? 1 : mask.water;
		mask.translucent	= GetTranslucentMask(coord.st, mask.water);
        mask.translucent    = (mask.ice > 0.5 || mask.stainedGlass > 0.5 || mask.water > 0.5) ? 1 : mask.translucent;

    return mask;
}