#version 120
#define program_composite2
#define FRAGMENT

#define BLOOM
#define EXPOSURE

varying vec2 texcoord;

uniform sampler2D colortex7;
uniform sampler2D colortex1;
uniform sampler2D colortex5;

uniform float viewWidth;
uniform float viewHeight;
uniform ivec2 eyeBrightnessSmooth;

#include "/lib/utilities.glsl"
#include "/lib/fragment/camera.glsl"

const float overlap = 0.2;

const float rgOverlap = 0.1 * overlap;
const float rbOverlap = 0.01 * overlap;
const float gbOverlap = 0.04 * overlap;

const mat3 coneOverlap = mat3(1.0, 			rgOverlap, 	rbOverlap,
							  rgOverlap, 	1.0, 		gbOverlap,
							  rbOverlap, 	rgOverlap, 	1.0);

const mat3 coneOverlapInverse = mat3(	1.0 + (rgOverlap + rbOverlap), 			-rgOverlap, 	-rbOverlap,
									  	-rgOverlap, 		1.0 + (rgOverlap + gbOverlap), 		-gbOverlap,
									  	-rbOverlap, 		-rgOverlap, 	1.0 + (rbOverlap + rgOverlap));

vec3 jodieReinhardTonemap(vec3 c){
    float l = dot(c, vec3(0.2126, 0.7152, 0.0722));
    vec3 tc=c/(c+1.);
    return linearToSRGB(mix(c/(l+1.),tc,tc));
}

vec3 reinhardTonemap(vec3 c){
	return linearToSRGB(c / (c + 1.0));
}

vec3 burgressTonemap(vec3 c){
	vec3 x = c;
	return (x * (6.2 * x + 0.5)) / (x * (6.2 * x + 1.7) + 0.06);
}

//	ACES Fitting by Stephen Hill
vec3 RRTAndODTFit(vec3 v)
{
    vec3 a = v * (v + 0.0245786f) - 0.000090537f;
    vec3 b = v * (1.0f * v + 0.4329510f) + 0.238081f;
    return a / b;
}

vec3 ACESTonemap2(vec3 color)
{
	color *= 1.5;
	color = color * coneOverlap;
	//color = linearToSRGB(color);

    // Apply RRT and ODT
    color = RRTAndODTFit(color);


    // Clamp to [0, 1]
	color = color * coneOverlapInverse;
    color = clamp01(color);

    return color;
}

vec3 ACESFilmTonemap(vec3 x )
{
	x = linearToSRGB(x);
    const float a = 2.51;
    const float b = 0.03;
    const float c = 2.43;
    const float d = 0.59;
    const float e = 0.14;
    return clamp01((x*(a*x+b))/(x*(c*x+d)+e));
}

vec3 hableTonemapInjector(vec3 x){
	const float A = 0.15;
	const float B = 0.50;
	const float C = 0.10;
	const float D = 0.20;
	const float E = 0.02;
	const float F = 0.30;

	return ((x*(A*x+C*B)+D*E)/(x*(A*x+B)+D*F))-E/F;
}

vec3 hableTonemap(vec3 x){
	const float W = 11.2;

	const float exposureBias = 2.0;
	vec3 curr = hableTonemapInjector(x);

	vec3 whiteScale = 1.0 / hableTonemapInjector(vec3(W));
	vec3 color = curr * whiteScale;

	return linearToSRGB(color);
}

vec3 roboboTonemap(vec3 x) {

	const vec3 tonemapCurve = vec3(5.5);
	x *= 2.0;

	x = linearToSRGB(x);
	x = pow(x, tonemapCurve);
	x = x / (x + 1.0);
	x = pow(x, 1.0 / tonemapCurve);

	x = mix(x, x * x * (3.0 - 2.0 * x), vec3(0.2));

	return x;
}

#include "/lib/utilities/bicubic.glsl"

vec3 calculateBloomTile(sampler2D textureSample, vec2 coord, const float lod){
	float lodScale = exp2(-lod);
	float offset = lodScale * 1.5;

	return decodeRGBE8(BicubicTexture(textureSample, coord * lodScale + offset));
}

vec3 calculateBloom(vec2 coord, float EV, vec2 pixelSize){
	vec3 bloom = vec3(0.0);

	// const float lods[7] = float[7](
	// 	2.0,
	// 	3.0,
	// 	4.0,
	// 	5.0,
	// 	6.0,
	// 	7.0,
	// 	8.0
	// );

	// for (int i = 0; i < 7; ++i){
	// 	bloom += calculateBloomTile(colortex1, coord, lods[i]);
	// }
	bloom += calculateBloomTile(colortex1, coord, 2.0);
	bloom += calculateBloomTile(colortex1, coord, 3.0);
	bloom += calculateBloomTile(colortex1, coord, 4.0);
	bloom += calculateBloomTile(colortex1, coord, 5.0);
	bloom += calculateBloomTile(colortex1, coord, 6.0);
	bloom += calculateBloomTile(colortex1, coord, 7.0);
	bloom += calculateBloomTile(colortex1, coord, 8.0);

	return decodeColor(bloom) * (1.0 / 7.) * exp2(EV - 3.0);
}

vec3 calculateLowLightDesaturation(vec3 color) {
	const vec3 preTint = vec3(0.25, 0.75, 2.0);
	const vec3 saturatedPreTint = mix(preTint, vec3(dot(preTint, lumCoeff)), -0.5);

	float avg = dot(color, lumCoeff);
	float range = exp2(-avg * 17.0);

	return mix(color, avg * saturatedPreTint, range);
}

void calculateAverageExposure(vec3 color)
{
	float avglod = int(log2(min(viewWidth, viewHeight))) - 0;
	//color /= pow(Luminance(texture2DLod(gdepth, vec2(0.5, 0.5), avglod).rgb), 1.18) * 1.2 + 0.00000005;

	float avgLumPow = 1.1;
	float exposureMax = 0.9;
	float exposureMin = 0.00005;

	color /= pow(dot(texture2DLod(colortex7, vec2(0.65, 0.65), avglod).rgb, lumCoeff), avgLumPow) * exposureMax + exposureMin;
}

void calculateExposureEyeBrightness(vec3 color) 
{
	float exposureMax = 1.55f;
		  //exposureMax *= mix(1.0f, 0.25f, timeSunriseSunset);
		  //exposureMax *= mix(1.0f, 0.0f, timeMidnight);
		  //exposureMax *= mix(1.0f, 0.25f, rainStrength);
	float exposureMin = 0.07f;
	float exposure = pow(eyeBrightnessSmooth.y / 240.0f, 6.0f) * exposureMax + exposureMin;

	//exposure = 1.0f;

	color.rgb /= vec3(exposure);
	color.rgb *= 350.0;
}

/* DRAWBUFFERS:0 */
void main() {
	vec2 pixelSize = 1.0 / vec2(viewWidth, viewHeight);

	vec4 colorSample = texture2D(colortex7, texcoord);
	vec3 color = decodeColor(decodeRGBE8(colorSample));

	float avgLum = decodeColor(texture2D(colortex5, texcoord).a);
	float exposureValue = calculateExposure(avgLum);

	vec3 bloom = calculateBloom(texcoord, exposureValue, pixelSize);

	//color = calculateLowLightDesaturation(color);
	#ifdef BLOOM
		color += bloom;
	#endif

	//曝光处理
	color *= exposureValue;
	color *= 1.5;
	//bloom *= exposureValue * 3;
	//calculateExposureEyeBrightness(color);

	//color = encodeColor(color);
	color = roboboTonemap(color);
	//bloom = roboboTonemap(bloom);

	gl_FragData[0] = vec4(color, 1.0);
}
