#version 120

#define SHADOW_MAP_BIAS 0.9


varying vec4 texcoord;
varying vec4 vPosition;
varying vec4 color;
varying vec4 lmcoord;

varying vec3 normal;
varying vec3 rawNormal;

attribute vec4 mc_Entity;

varying float materialIDs;
varying float iswater;
varying float isStainedGlass;
varying vec4 viewPos;

#include "libs/uniform.glsl"

		const vec2[8] offsets = vec2[8](vec2(1./8.,-3./8.),
									vec2(-1.,3.)/8.,
									vec2(5.0,1.)/8.,
									vec2(-3,-5.)/8.,
									vec2(-5.,5.)/8.,
									vec2(-7.,-1.)/8.,
									vec2(3,7.)/8.,
									vec2(7.,-7.)/8.);

vec4 cubic(float x)
{
    float x2 = x * x;
    float x3 = x2 * x;
    vec4 w;
    w.x =   -x3 + 3*x2 - 3*x + 1;
    w.y =  3*x3 - 6*x2       + 4;
    w.z = -3*x3 + 3*x2 + 3*x + 1;
    w.w =  x3;
    return w / 6.f;
}

vec2 saturate(vec2 x)
{
	return clamp(x, 0.0, 1.0);   
}

vec4 BicubicTexture(in sampler2D tex, in vec2 coord)
{
	int resolution = 64;

	coord *= resolution;

	float fx = fract(coord.x);
    float fy = fract(coord.y);
    coord.x -= fx;
    coord.y -= fy;

    vec4 xcubic = cubic(fx);
    vec4 ycubic = cubic(fy);

    vec4 c = vec4(coord.x - 0.5, coord.x + 1.5, coord.y - 0.5, coord.y + 1.5);
    vec4 s = vec4(xcubic.x + xcubic.y, xcubic.z + xcubic.w, ycubic.x + ycubic.y, ycubic.z + ycubic.w);
    vec4 offset = c + vec4(xcubic.y, xcubic.w, ycubic.y, ycubic.w) / s;

    vec4 sample0 = texture2D(tex, vec2(offset.x, offset.z) / resolution);
    vec4 sample1 = texture2D(tex, vec2(offset.y, offset.z) / resolution);
    vec4 sample2 = texture2D(tex, vec2(offset.x, offset.w) / resolution);
    vec4 sample3 = texture2D(tex, vec2(offset.y, offset.w) / resolution);

    float sx = s.x / (s.x + s.y);
    float sy = s.z / (s.z + s.w);

    return mix( mix(sample3, sample2, sx), mix(sample1, sample0, sx), sy);
}


vec4 TextureSmooth(in sampler2D tex, in vec2 coord)
{
	int resolution = 64;

	coord *= resolution;
	vec2 i = floor(coord);
	vec2 f = fract(coord);
	     f = f * f * (3.0f - 2.0f * f);

	coord = (i + f) / resolution;

	vec4 result = texture2D(tex, coord);

	return result;
}

void main() {
	gl_Position = ftransform();

	lmcoord = gl_TextureMatrix[1] * gl_MultiTexCoord1;
	texcoord = gl_MultiTexCoord0;
	

	viewPos = gl_ModelViewMatrix * gl_Vertex;

	vec4 position = gl_Position;


		 //position *= position.w;

		 position = shadowProjectionInverse * position;
		 position = shadowModelViewInverse * position;
		 position.xyz += cameraPosition.xyz;
		 //position = gbufferModelView * position;

	//convert to world-space position

	materialIDs = 0.0f;


	iswater = 0.0;

	if (mc_Entity.x == 1971.0f)
	{
		iswater = 1.0f;
	}

	if (mc_Entity.x == 8 || mc_Entity.x == 9) {
		iswater = 1.0f;
	}

	float isice = 0.0f;


	
	if (mc_Entity.x == 79) {
		isice = 1.0f;
	}

	isStainedGlass = 0.0f;

	if (mc_Entity.x == 95 || mc_Entity.x == 160)
	{
		isStainedGlass = 1.0f;
	}


	//Grass
	if  (  mc_Entity.x == 31.0

		|| mc_Entity.x == 38.0f 	//Rose
		|| mc_Entity.x == 37.0f 	//Flower

		/*
		|| mc_Entity.x == 1925.0f 	//Biomes O Plenty: Medium Grass
		|| mc_Entity.x == 1920.0f 	//Biomes O Plenty: Thorns, barley
		|| mc_Entity.x == 1921.0f 	//Biomes O Plenty: Sunflower
		|| mc_Entity.x == 188.0f 	//Biomes O Plenty: Medium Grass
		|| mc_Entity.x == 176.0f 	//Biomes O Plenty: Desert Grass
		|| mc_Entity.x == 177.0f 	//Biomes O Plenty: Desert Grass
		|| mc_Entity.x == 178.0f 	//Lavender

		*/
		)
	{
			materialIDs = max(materialIDs, 2.0f);
	}

	//Wheat
	if (mc_Entity.x == 59.0) {
		materialIDs = max(materialIDs, 2.0f);
	}	
	
	//Leaves
	if   ( mc_Entity.x == 18.0 

		|| mc_Entity.x == 161.0f

		/*
		|| mc_Entity.x == 1962.0f //Biomes O Plenty: Leaves
		|| mc_Entity.x == 1924.0f //Biomes O Plenty: Leaves
		|| mc_Entity.x == 1923.0f //Biomes O Plenty: Leaves
		|| mc_Entity.x == 1926.0f //Biomes O Plenty: Leaves
		|| mc_Entity.x == 1936.0f //Biomes O Plenty: Giant Flower Leaves
		|| mc_Entity.x == 184.0f  //Yellow autumn leaves
		|| mc_Entity.x == 185.0f  //Dying leaves
		|| mc_Entity.x == 186.0f  //maple leaves
		|| mc_Entity.x == 187.0f  //maple leaves
		|| mc_Entity.x == 192.0f  //maple leaves
		|| mc_Entity.x == 249.0f  //Willow leaves
		|| mc_Entity.x == 248.0f  //Sacred Oak Leaves
		*/

		 ) {
		materialIDs = max(materialIDs, 3.0f);
	}
		
	//Ice
	if (  mc_Entity.x == 79.0f
	   || mc_Entity.x == 174.0f)
	{
		materialIDs = max(materialIDs, 4.0f);
	}

	//Cobweb
	if ( mc_Entity.x == 30.0f)
	{
		materialIDs = max(materialIDs, 11.0f);
	}

	if (iswater > 0.5 || isice > 0.5)
	{
		position.xyz += 10000.0;
	}

	//position = gbufferModelViewInverse * position;
	position.xyz -= cameraPosition.xyz;
	position = shadowModelView * position;
	position = shadowProjection * position;
	//position.xy += offsets[frameMod8] * position.w / 2048.0;

	vec3 worldNormal = gl_Normal;

	if (abs(materialIDs - 2.0) < 0.1)
	{
		worldNormal = vec3(0.0, 1.0, 0.0);
	}

	normal = normalize(gl_NormalMatrix * worldNormal);

	color = gl_Color;

	gl_Position = position;



	float dist = length(gl_Position.xy);
	float distortFactor = (1.0f - SHADOW_MAP_BIAS) + dist * SHADOW_MAP_BIAS;
	gl_Position.xy *= 0.95f / distortFactor;
	gl_Position.z = mix(gl_Position.z, 0.5, 0.8);

	//vec2 warp = abs(gl_Position.xy);
	

	//gl_Position.x /= warp.x * 0.90 + 0.1;
	//gl_Position.y /= warp.y * 0.9 + 0.1;




	vPosition = gl_Position;

	gl_FrontColor = gl_Color;
	
}
