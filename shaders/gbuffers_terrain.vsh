#version 120

#define TAA
#define NORMALS

uniform sampler2D noisetex;

varying vec4 texcoord;
varying vec4 lmcoord;

varying vec4 vtexcoord;
varying vec4 vtexcoordam; // .st for add, .pq for mul.

varying vec4 color;
varying vec3 normal;
varying vec3 tangent;
varying vec3 binormal;
varying vec3 viewVector;

uniform vec2 texelSize;
uniform int frameMod8;
uniform int frameMod16;
uniform vec2 taaOffset16;
uniform int frameCounter;
uniform float frameTimeCounter;
uniform float rainStrength;

attribute vec4 mc_Entity;
attribute vec4 mc_midTexCoord;
attribute vec4 at_tangent;

varying float terrainMats;
varying float distances;;
varying vec3 wpos;
varying float top;

uniform vec3 cameraPosition;
uniform vec3 upPosition;
uniform mat4 gbufferModelView;
uniform mat4 gbufferModelViewInverse;
uniform mat4 shadowModelView;
uniform mat4 shadowModelViewInverse;
uniform mat4 gbufferProjection;
uniform mat4 gbufferProjectionInverse;

varying float minecraftMats;

#ifndef NORMALS
varying float n2;
#define clamp01(x) clamp(x, 0.0, 1.0)
#define fsign(x) (clamp01(x * 1e35) * 2.0 - 1.0)

float encodeVec2(vec2 a) {
	const vec2 constant1 = vec2(1.0, 256.0) / 65535.0; //2^16-1
	return dot(floor(a * 255.0), constant1);
}

float encodeNormal(vec3 a) {
    vec3 b = abs(a);
    vec2 p = a.xy / dot(b, vec3(1.0));
    vec2 encoded = a.z <= 0. ? (1. - abs(p.yx)) * fsign(p) : p;
    encoded = encoded * .5 + .5;
	return encodeVec2(encoded);
}
#endif

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

float clacBasicWave(float pos, float offset)
{
	return exp(-sin(pos * 3.1415926)) * offset;
}

void main() {

	vec4 viewpos = gbufferModelViewInverse * gl_ModelViewMatrix * gl_Vertex;
	vec4 position = viewpos;

    //gl_Position = gl_ProjectionMatrix * gbufferModelView * position;
    texcoord = gl_MultiTexCoord0;
	lmcoord = gl_TextureMatrix[1] * gl_MultiTexCoord1;
	
	//tangent = normalize(gl_NormalMatrix * at_tangent.xyz);
    //binormal = cross(tangent, normal);


	vec2 midcoord 		  = (gl_TextureMatrix[0] *  mc_midTexCoord).st;
	vec2 texcoordminusmid = texcoord.st - midcoord;
	
	vtexcoordam.pq  = abs(texcoordminusmid) * 2;
	vtexcoordam.st  = min(texcoord.st, midcoord - texcoordminusmid);
	vtexcoord.xy    = sign(texcoordminusmid) * 0.5 + 0.5;

    color = gl_Color;

	//Gather materials
	terrainMats = 1.0f;

	float facingEast = abs(normalize(gl_Normal.xz).x);
	float facingUp = abs(gl_Normal.y);

	float waveCoeff = 0.0f;

	minecraftMats = 0.0f;
	if(mc_Entity.x==10||mc_Entity.x==11||mc_Entity.x==62||mc_Entity.x==76||mc_Entity.x==213)
	{
		minecraftMats=1.0;
	}
	if(mc_Entity.x==50||mc_Entity.x==51||mc_Entity.x==89||mc_Entity.x==91||mc_Entity.x==348)
	{
		minecraftMats=2.0;
	}
	if(mc_Entity.x==138||mc_Entity.x==169||mc_Entity.x==411)
	{
		minecraftMats=3.0;
	}
	if(mc_Entity.x==119||mc_Entity.x==124||mc_Entity.x==198||mc_Entity.x==412)
	{
		minecraftMats=4.0;
	}
	if(mc_Entity.x==90||mc_Entity.x==426)
	{
		minecraftMats=5.0;
	}

	//Grass
	if  (  mc_Entity.x == 31.0

		|| mc_Entity.x == 38.0f 	//Rose
		|| mc_Entity.x == 37.0f 	//Flower
		|| mc_Entity.x == 1925.0f 	//Biomes O Plenty: Medium Grass
		|| mc_Entity.x == 1920.0f 	//Biomes O Plenty: Thorns, barley
		|| mc_Entity.x == 1921.0f 	//Biomes O Plenty: Sunflower
		|| mc_Entity.x == 2.0 && gl_Normal.y < 0.5 && facingEast > 0.01 && facingEast < 0.99 && facingUp < 0.9

		)
	{
		terrainMats = max(terrainMats, 2.0f);
		waveCoeff = 1.0f;
	}

	if (  mc_Entity.x == 175.0f)
	{
		terrainMats = max(terrainMats, 2.0f);
	}
	
	//Wheat
	if (mc_Entity.x == 59.0) {
		terrainMats = max(terrainMats, 2.0f);
		waveCoeff = 1.0f;
	}	
	
	//Leaves
	if   ( mc_Entity.x == 18.0 

		|| mc_Entity.x == 1962.0f //Biomes O Plenty: Leaves
		|| mc_Entity.x == 1924.0f //Biomes O Plenty: Leaves
		|| mc_Entity.x == 1923.0f //Biomes O Plenty: Leaves
		|| mc_Entity.x == 1926.0f //Biomes O Plenty: Leaves
		|| mc_Entity.x == 1936.0f //Biomes O Plenty: Giant Flower Leaves
		|| mc_Entity.x == 161.0f //Biomes O Plenty: Giant Flower Leaves

		 ) {
		terrainMats = max(terrainMats, 3.0f);
	}	

	
	//Gold block
	if (mc_Entity.x == 41) {
		terrainMats = max(terrainMats, 20.0f);
	}
	
	//Iron block
	if (mc_Entity.x == 42) {
		terrainMats = max(terrainMats, 21.0f);
	}
	
	//Diamond Block
	if (mc_Entity.x == 57) {
		terrainMats = max(terrainMats, 22.0f);
	}
	
	//Emerald Block
	if (mc_Entity.x == -123) {
		terrainMats = max(terrainMats, 23.0f);
	}
	
	//sand
	if (mc_Entity.x == 12) {
		terrainMats = max(terrainMats, 24.0f);
	}

	//sandstone
	if (mc_Entity.x == 24 || mc_Entity.x == -128) {
		terrainMats = max(terrainMats, 25.0f);
	}
	
	//stone
	if (mc_Entity.x == 1) {
		terrainMats = max(terrainMats, 26.0f);
	}
	
	//cobblestone
	if (mc_Entity.x == 4) {
		terrainMats = max(terrainMats, 27.0f);
	}
	
	//wool
	if (mc_Entity.x == 35) {
		terrainMats = max(terrainMats, 28.0f);
	}


	//torch	
	if (mc_Entity.x == 50) {
		terrainMats = max(terrainMats, 30.0f);
	}

	//lava
	if (mc_Entity.x == 10 || mc_Entity.x == 11) {
		terrainMats = max(terrainMats, 31.0f);
	}

	//glowstone and lamp
	if (mc_Entity.x == 89 || mc_Entity.x == 124) {
		terrainMats = max(terrainMats, 32.0f);
	}

	//fire
	if (mc_Entity.x == 51) {
		terrainMats = max(terrainMats, 33.0f);
	}

	position.xyz += cameraPosition.xyz;
	if(waveCoeff > 0.5f || (mc_Entity.x == 59.0 && texcoord.t < 0.35) || (terrainMats == 3.0f && texcoord.t < 1.90 && texcoord.t > -1.0))
	{
		float windStrength = mix(0.85f, 1.0f, rainStrength);
		float wind = frameTimeCounter * windStrength;

		vec3 pn    = position.xyz;
			 pn.x *= 2.0f;
			 pn.x -= frameTimeCounter * 15.0f;
			 pn.z *= 8.0f;

		vec3 randomwind = BicubicTexture(noisetex, pn.xz / 512.0f).xyz;

		vec3 rp = position.xyz;
		rp += exp(-sin(rp.z / 2.0) * 4.0) + randomwind * 2.0f;

		float freq = 8.0f, offsetamp = 0.03f;
		vec3 amp = vec3(0.05, 0.0001, 0.05);
		float iSteps = 2.0f, jSteps = 2.0f;
		int total = 1;
		for(float i = -iSteps; i < iSteps; i++)
		for(float j = -jSteps; j < jSteps; j++)
		{
			vec3 offset = vec3(i, randomwind.z + i, j) * offsetamp;
			position.x += clacBasicWave(rp.x, offset.x) * amp.x;
			position.y += clacBasicWave(rp.y, offset.y) * amp.y;
			position.z += clacBasicWave(rp.z, offset.z) * amp.z;
		}
		//position.xyz *= amp;
	}
	position.xyz -= cameraPosition.xyz;

	distances = length(position.xyz);

	gl_Position = gl_ProjectionMatrix * gbufferModelView * position;

	#ifdef TAA
		gl_Position.xy += offsets[frameMod8] * gl_Position.w * texelSize;
	#endif

	normal = gl_Normal;
	normal = (gl_NormalMatrix * normal) * mat3(gbufferModelView);
	
	wpos = position.xyz;
	top = dot(normalize(upPosition), normal);

	// normal			= normalize(gl_NormalMatrix * gl_Normal);

	if (gl_Normal.y > 0.5) {
		//  0.0,  1.0,  0.0
		tangent.xyz  = normalize(gl_NormalMatrix * vec3( 1.0,  0.0,  0.0));
		binormal.xyz = normalize(gl_NormalMatrix * vec3( 0.0,  0.0,  1.0));
	} else if (gl_Normal.x >  0.5) {
		//  1.0,  0.0,  0.0
		tangent.xyz  = normalize(gl_NormalMatrix * vec3( 0.0,  0.0, -1.0));
		binormal.xyz = normalize(gl_NormalMatrix * vec3( 0.0, -1.0,  0.0));
	} else if (gl_Normal.x < -0.5) {
		// -1.0,  0.0,  0.0
		tangent.xyz  = normalize(gl_NormalMatrix * vec3( 0.0,  0.0,  1.0));
		binormal.xyz = normalize(gl_NormalMatrix * vec3( 0.0, -1.0,  0.0));
	} else if (gl_Normal.z >  0.5) {
		//  0.0,  0.0,  1.0
		tangent.xyz  = normalize(gl_NormalMatrix * vec3( 1.0,  0.0,  0.0));
		binormal.xyz = normalize(gl_NormalMatrix * vec3( 0.0, -1.0,  0.0));
	} else if (gl_Normal.z < -0.5) {
		//  0.0,  0.0, -1.0
		tangent.xyz  = normalize(gl_NormalMatrix * vec3(-1.0,  0.0,  0.0));
		binormal.xyz = normalize(gl_NormalMatrix * vec3( 0.0, -1.0,  0.0));
	} else if (gl_Normal.y < -0.5) {
		//  0.0, -1.0,  0.0
		tangent.xyz  = normalize(gl_NormalMatrix * vec3( 1.0,  0.0,  0.0));
		binormal.xyz = normalize(gl_NormalMatrix * vec3( 0.0,  0.0,  1.0));
	}
	
	mat3 tbnMatrix = mat3(tangent.x, binormal.x, normal.x,
						  tangent.y, binormal.y, normal.y,
						  tangent.z, binormal.z, normal.z);	

	viewVector = (gl_ModelViewMatrix * gl_Vertex).xyz;
	viewVector = normalize(tbnMatrix * viewVector);
	
	#ifndef NORMALS
	n2 = encodeNormal(normal);
	#endif
}