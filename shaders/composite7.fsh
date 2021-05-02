#version 120

#extension GL_EXT_gpu_shader4 : enable

#define TAA

//#define FAST_TAA //disables bicubic resampling and closest velocity, improves fps especially at high resolutions

//TAA OPTIONS
//#define NO_CLIP	//introduces a lot of ghosting but the image will be sharper and without flickering (good for screenshots)
#define BLEND_FACTOR 0.1 //[0.01 0.02 0.03 0.04 0.05 0.06 0.08 0.1 0.12 0.14 0.16] higher values = more flickering but sharper image, lower values = less flickering but the image will be blurrier
#define BICUBIC_RESAMPLING //Can cause artifacts on high contrast edges, but looks way better in motion. Performance cost is a bit higher.
//#define FAST_BICUBIC
#define MOTION_REJECTION 0.5 //[0.0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.5] //Higher values=sharper image in motion at the cost of flickering
#define FLICKER_REDUCTION 0.75 //[0.0 0.75 0.85 0.9 0.92 0.94 0.96 0.98 1.0]  //Higher values = Blurrier image but greatly reduces flickering (0-1 range)
#define CLOSEST_VELOCITY //improves edge quality in motion at the cost of performance
//#define SMOOTHESTSTEP_INTERPOLATION //Only if not using bicubic resampling, reduces blurring in motion but might cause some "wobbling" on the image

#ifdef FAST_TAA
	#undef BICUBIC_RESAMPLING
	#undef CLOSEST_VELOCITY
#endif

varying vec4 texcoord;

#include "libs/uniform.glsl"
#include "libs/composite.glsl"

#define diagonal3(m) vec3((m)[0].x, (m)[1].y, m[2].z)
#define  projMAD(m, v) (diagonal3(m) * (v) + (m)[3].xyz)

const int noiseTextureResolution = 32;

flat varying float tempOffsets;

float luma(vec3 color) {
	return dot(color,vec3(0.21, 0.72, 0.07));
}

vec3 toClipSpace3Prev(vec3 viewSpacePosition) {
    return projMAD(gbufferPreviousProjection, viewSpacePosition) / -viewSpacePosition.z * 0.5 + 0.5;
}

vec3 toScreenSpace(vec3 p) {
	vec4 iProjDiag = vec4(gbufferProjectionInverse[0].x, gbufferProjectionInverse[1].y, gbufferProjectionInverse[2].zw);
    vec3 p3 = p * 2. - 1.;
    vec4 fragposition = iProjDiag * p3.xyzz + gbufferProjectionInverse[3];
    return fragposition.xyz / fragposition.w;
}

//returns the projected coordinates of the closest point to the camera in the 3x3 neighborhood
vec3 closestToCamera3x3()
{
	vec2 du = vec2(texelSize.x, 0.0);
	vec2 dv = vec2(0.0, texelSize.y);

	vec3 dtl = vec3(texcoord.st,0.) + vec3(-texelSize, texture2D(depthtex0, texcoord.st - dv - du).x);
	vec3 dtc = vec3(texcoord.st,0.) + vec3( 0.0, -texelSize.y, texture2D(depthtex0, texcoord.st - dv).x);
	vec3 dtr = vec3(texcoord.st,0.) +  vec3( texelSize.x, -texelSize.y, texture2D(depthtex0, texcoord.st - dv + du).x);

	vec3 dml = vec3(texcoord.st,0.) +  vec3(-texelSize.x, 0.0, texture2D(depthtex0, texcoord.st - du).x);
	vec3 dmc = vec3(texcoord.st,0.) + vec3( 0.0, 0.0, texture2D(depthtex0, texcoord.st).x);
	vec3 dmr = vec3(texcoord.st,0.) + vec3( texelSize.x, 0.0, texture2D(depthtex0, texcoord.st + du).x);

	vec3 dbl = vec3(texcoord.st,0.) + vec3(-texelSize.x, texelSize.y, texture2D(depthtex0, texcoord.st + dv - du).x);
	vec3 dbc = vec3(texcoord.st,0.) + vec3( 0.0, texelSize.y, texture2D(depthtex0, texcoord.st + dv).x);
	vec3 dbr = vec3(texcoord.st,0.) + vec3( texelSize.x, texelSize.y, texture2D(depthtex0, texcoord.st + dv + du).x);

	vec3 dmin = dmc;

	dmin = dmin.z > dtc.z? dtc : dmin;
	dmin = dmin.z > dtr.z? dtr : dmin;

	dmin = dmin.z > dml.z? dml : dmin;
	dmin = dmin.z > dtl.z? dtl : dmin;
	dmin = dmin.z > dmr.z? dmr : dmin;

	dmin = dmin.z > dbl.z? dbl : dmin;
	dmin = dmin.z > dbc.z? dbc : dmin;
	dmin = dmin.z > dbr.z? dbr : dmin;

	return dmin;
}

//Modified texture interpolation from inigo quilez
vec4 smoothfilter(in sampler2D tex, in vec2 uv)
{
	vec2 textureResolution = vec2(viewWidth,viewHeight);
	uv = uv*textureResolution + 0.5;
	vec2 iuv = floor( uv );
	vec2 fuv = fract( uv );
	#ifndef SMOOTHESTSTEP_INTERPOLATION
	uv = iuv + (fuv*fuv)*(3.0-2.0*fuv);
	#endif
	#ifdef SMOOTHESTSTEP_INTERPOLATION
	uv = iuv + fuv*fuv*fuv*(fuv*(fuv*6.0-15.0)+10.0);
	#endif
	uv = (uv - 0.5)/textureResolution;
	return texture2D( tex, uv);
}

//from : https://gist.github.com/TheRealMJP/c83b8c0f46b63f3a88a5986f4fa982b1
vec4 SampleTextureCatmullRom(sampler2D tex, vec2 uv, vec2 texSize )
{
    // We're going to sample a a 4x4 grid of texels surrounding the target UV coordinate. We'll do this by rounding
    // down the sample location to get the exact center of our "starting" texel. The starting texel will be at
    // location [1, 1] in the grid, where [0, 0] is the top left corner.
    vec2 samplePos = uv * texSize;
    vec2 texPos1 = floor(samplePos - 0.5) + 0.5;

    // Compute the fractional offset from our starting texel to our original sample location, which we'll
    // feed into the Catmull-Rom spline function to get our filter weights.
    vec2 f = samplePos - texPos1;

    // Compute the Catmull-Rom weights using the fractional offset that we calculated earlier.
    // These equations are pre-expanded based on our knowledge of where the texels will be located,
    // which lets us avoid having to evaluate a piece-wise function.
    vec2 w0 = f * ( -0.5 + f * (1.0 - 0.5*f));
    vec2 w1 = 1.0 + f * f * (-2.5 + 1.5*f);
    vec2 w2 = f * ( 0.5 + f * (2.0 - 1.5*f) );
    vec2 w3 = f * f * (-0.5 + 0.5 * f);

    // Work out weighting factors and sampling offsets that will let us use bilinear filtering to
    // simultaneously evaluate the middle 2 samples from the 4x4 grid.
    vec2 w12 = w1 + w2;
    vec2 offset12 = w2 / (w1 + w2);

    // Compute the final UV coordinates we'll use for sampling the texture
    vec2 texPos0 = texPos1 - vec2(1.0);
    vec2 texPos3 = texPos1 + vec2(2.0);
    vec2 texPos12 = texPos1 + offset12;

    texPos0 *= texelSize;
    texPos3 *= texelSize;
    texPos12 *= texelSize;

    vec4 result = vec4(0.0);
    result += texture2D(tex, vec2(texPos0.x,  texPos0.y)) * w0.x * w0.y;
    result += texture2D(tex, vec2(texPos12.x, texPos0.y)) * w12.x * w0.y;
    result += texture2D(tex, vec2(texPos3.x,  texPos0.y)) * w3.x * w0.y;

    result += texture2D(tex, vec2(texPos0.x,  texPos12.y)) * w0.x * w12.y;
    result += texture2D(tex, vec2(texPos12.x, texPos12.y)) * w12.x * w12.y;
    result += texture2D(tex, vec2(texPos3.x,  texPos12.y)) * w3.x * w12.y;

    result += texture2D(tex, vec2(texPos0.x,  texPos3.y)) * w0.x * w3.y;
    result += texture2D(tex, vec2(texPos12.x, texPos3.y)) * w12.x * w3.y;
    result += texture2D(tex, vec2(texPos3.x,  texPos3.y)) * w3.x * w3.y;

    return result;
}
//approximation from SMAA presentation from siggraph 2016
vec3 FastCatmulRom(sampler2D colorTex, vec2 texcoord, vec4 rtMetrics)
{
    vec2 position = rtMetrics.zw * texcoord;
    vec2 centerPosition = floor(position - 0.5) + 0.5;
    vec2 f = position - centerPosition;
    vec2 f2 = f * f;
    vec2 f3 = f * f2;

    float c = 0.5;
    vec2 w0 =        -c  * f3 +  2.0 * c         * f2 - c * f;
    vec2 w1 =  (2.0 - c) * f3 - (3.0 - c)        * f2         + 1.0;
    vec2 w2 = -(2.0 - c) * f3 + (3.0 -  2.0 * c) * f2 + c * f;
    vec2 w3 =         c  * f3 -                c * f2;

    vec2 w12 = w1 + w2;
    vec2 tc12 = rtMetrics.xy * (centerPosition + w2 / w12);
    vec3 centerColor = texture2D(colorTex, vec2(tc12.x, tc12.y)).rgb;

    vec2 tc0 = rtMetrics.xy * (centerPosition - 1.0);
    vec2 tc3 = rtMetrics.xy * (centerPosition + 2.0);
    vec4 color = vec4(texture2D(colorTex, vec2(tc12.x, tc0.y )).rgb, 1.0) * (w12.x * w0.y ) +
                   vec4(texture2D(colorTex, vec2(tc0.x,  tc12.y)).rgb, 1.0) * (w0.x  * w12.y) +
                   vec4(centerColor,                                      1.0) * (w12.x * w12.y) +
                   vec4(texture2D(colorTex, vec2(tc3.x,  tc12.y)).rgb, 1.0) * (w3.x  * w12.y) +
                   vec4(texture2D(colorTex, vec2(tc12.x, tc3.y )).rgb, 1.0) * (w12.x * w3.y );
	return color.rgb/color.a;

}

	vec3 clip_aabb(vec3 q,vec3 aabb_min, vec3 aabb_max)
	{
		vec3 p_clip = 0.5 * (aabb_max + aabb_min);
		vec3 e_clip = 0.5 * (aabb_max - aabb_min) + 0.00000001;

		vec3 v_clip = q - vec3(p_clip);
		vec3 v_unit = v_clip.xyz / e_clip;
		vec3 a_unit = abs(v_unit);
		float ma_unit = max(a_unit.x, max(a_unit.y, a_unit.z));

		if (ma_unit > 1.0)
			return vec3(p_clip) + v_clip / ma_unit;
		else
			return q;
	}

vec3 TAA_hq(){
	//Samples current frame 3x3 neighboorhood
	vec3 albedoCurrent0 = texture2D(colortex3, texcoord.st).rgb;
  vec3 albedoCurrent1 = texture2D(colortex3, texcoord.st + vec2(texelSize.x,texelSize.y)).rgb;
	vec3 albedoCurrent2 = texture2D(colortex3, texcoord.st + vec2(texelSize.x,-texelSize.y)).rgb;
	vec3 albedoCurrent3 = texture2D(colortex3, texcoord.st + vec2(-texelSize.x,-texelSize.y)).rgb;
	vec3 albedoCurrent4 = texture2D(colortex3, texcoord.st + vec2(-texelSize.x,texelSize.y)).rgb;
	vec3 albedoCurrent5 = texture2D(colortex3, texcoord.st + vec2(0.0,texelSize.y)).rgb;
	vec3 albedoCurrent6 = texture2D(colortex3, texcoord.st + vec2(0.0,-texelSize.y)).rgb;
	vec3 albedoCurrent7 = texture2D(colortex3, texcoord.st + vec2(-texelSize.x,0.0)).rgb;
	vec3 albedoCurrent8 = texture2D(colortex3, texcoord.st + vec2(texelSize.x,0.0)).rgb;

	//use velocity from the nearest texel from camera in a 3x3 box in order to improve edge quality in motion
	#ifdef CLOSEST_VELOCITY
	vec3 closestToCamera = closestToCamera3x3();
	#endif

	#ifndef CLOSEST_VELOCITY
	vec3 closestToCamera = vec3(texcoord.st,texture2D(depthtex0,texcoord.st).x);
	#endif

	//reproject previous frame
	vec3 fragposition = toScreenSpace(closestToCamera);
	fragposition = mat3(gbufferModelViewInverse) * fragposition + gbufferModelViewInverse[3].xyz + (cameraPosition - previousCameraPosition);
	vec3 previousPosition = mat3(gbufferPreviousModelView) * fragposition + gbufferPreviousModelView[3].xyz;
	previousPosition = toClipSpace3Prev(previousPosition);
	vec2 velocity = previousPosition.xy - closestToCamera.xy;
	previousPosition.xy = texcoord.st + velocity;

	//to reduce error propagation caused by interpolation during history resampling, we will introduce back some aliasing in motion
	vec2 d = 0.5-abs(fract(previousPosition.xy*vec2(viewWidth,viewHeight)-texcoord.st*vec2(viewWidth,viewHeight))-0.5);
	#ifdef SMOOTHESTSTEP_INTERPOLATION
	d = d*d*d*(d*(d*6.0-15.0)+10.0);
	#endif
	#ifdef BICUBIC_RESAMPLING
	d = d*d*d*(d*(d*6.0-15.0)+10.0);
	#endif
	#ifndef SMOOTHESTSTEP_INTERPOLATION
	#ifndef BICUBIC_RESAMPLING
	d = d*d*(3.0-2.0*d);
	#endif
	#endif
	float mixFactor = (d.x+d.y)*MOTION_REJECTION;

	//reject history if off-screen
	if (previousPosition.x < 0.0 || previousPosition.y < 0.0 || previousPosition.x > 1.0 || previousPosition.y > 1.0) mixFactor = 1.0;

	//Sample history
	#ifndef BICUBIC_RESAMPLING
	vec3 albedoPrev = smoothfilter(colortex5, previousPosition.xy).xyz;
	#endif
	#ifdef BICUBIC_RESAMPLING
	#ifdef FAST_BICUBIC
	vec3 albedoPrev = FastCatmulRom(colortex5, previousPosition.xy,vec4(texelSize.x,texelSize.y,viewWidth,viewHeight)).xyz;
	#else
	vec3 albedoPrev = SampleTextureCatmullRom(colortex5, previousPosition.xy,vec2(viewWidth,viewHeight)).xyz;
	#endif
	#endif


	#ifndef NO_CLIP
	//Assuming the history color is a blend of the 3x3 neighborhood, we clamp the history to the min and max of each channel in the 3x3 neighborhood
	vec3 cMax = max(max(max(albedoCurrent0,albedoCurrent1),albedoCurrent2),max(albedoCurrent3,max(albedoCurrent4,max(albedoCurrent5,max(albedoCurrent6,max(albedoCurrent7,albedoCurrent8))))));
	vec3 cMin = min(min(min(albedoCurrent0,albedoCurrent1),albedoCurrent2),min(albedoCurrent3,min(albedoCurrent4,min(albedoCurrent5,min(albedoCurrent6,min(albedoCurrent7,albedoCurrent8))))));


	vec3 finalcAcc = clip_aabb(albedoPrev,cMin,cMax);


	//increases blending factor if history is far away from aabb, reduces ghosting at the cost of some flickering
	float invCmax = 1.0/luma(cMax);
	float isclamped = abs(luma(albedoPrev)-luma(finalcAcc))*invCmax;

	//reduces blending factor if current texel is far from history, reduces flickering
	float lumDiff2 = abs(luma(albedoPrev)-luma(albedoCurrent0))*invCmax;
	lumDiff2 = 1.0-clamp(lumDiff2*lumDiff2,0.,1.)*FLICKER_REDUCTION;

	//Blend current pixel with clamped history
	vec3 supersampled =  mix(finalcAcc,albedoCurrent0,clamp(BLEND_FACTOR*lumDiff2+mixFactor+0.01,0.,1.));
	//vec3 supersampled = finalcAcc;
	#endif


	#ifdef NO_CLIP
	vec3 supersampled =  mix(albedoPrev,albedoCurrent0,clamp(BLEND_FACTOR,0.,1.));
	#endif

	//De-tonemap
	return supersampled;
}

#ifdef FAST_TAA

#endif

/* DRAWBUFFERS:5 */

void main() {
    gl_FragData[0].a = 1.0;
	#ifdef TAA
	vec3 color = TAA_hq();
	gl_FragData[0].rgb = color;
	#endif
	#ifndef TAA
	vec3 color = texture2D(colortex5,texcoord.st).rgb;
	gl_FragData[0].rgb = color;
	#endif
}
