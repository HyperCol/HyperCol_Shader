#version 120

#define TAA
#define NORMALS

varying vec4 texcoord;
varying vec4 lmcoord;

varying vec4 color;
varying vec3 normal;

uniform vec2 texelSize;
uniform int frameMod8;
uniform int frameMod16;
uniform vec2 taaOffset16;

attribute vec4 mc_Entity;

varying float terrainMats;

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

		const vec2[8] offsets = vec2[8](vec2(1./8.,-3./8.),
									vec2(-1.,3.)/8.,
									vec2(5.0,1.)/8.,
									vec2(-3,-5.)/8.,
									vec2(-5.,5.)/8.,
									vec2(-7.,-1.)/8.,
									vec2(3,7.)/8.,
									vec2(7.,-7.)/8.);

void main() {
    gl_Position = ftransform();
    texcoord = gl_MultiTexCoord0;

	lmcoord = gl_TextureMatrix[1] * gl_MultiTexCoord1;

    color = gl_Color;

	normal = normalize(gl_NormalMatrix * gl_Normal);

	#ifdef TAA
	gl_Position.xy += offsets[frameMod8] * gl_Position.w * texelSize;
	#endif

	//Gather materials
	terrainMats = 1.0f;

	float facingEast = abs(normalize(gl_Normal.xz).x);
	float facingUp = abs(gl_Normal.y);

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
	}

	if (  mc_Entity.x == 175.0f)
	{
		terrainMats = max(terrainMats, 2.0f);
	}
	
	//Wheat
	if (mc_Entity.x == 59.0) {
		terrainMats = max(terrainMats, 2.0f);
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
	

	n2 = encodeNormal(normal);
}