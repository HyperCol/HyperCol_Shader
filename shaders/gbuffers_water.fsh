#version 130

varying vec4 texcoord;
varying vec4 lmcoord;

varying vec3 color;
varying vec3 normal;

varying float iswater;
varying float isice;
varying float isStainedGlass;

varying vec3 viewVector;
varying vec3 binormal;
varying vec3 tangent;
varying mat3 tbnMatrix;

varying float minecraftMats;

#include "libs/uniform.glsl"
#include "libs/composite.glsl"
#include "libs/noise.glsl"

/* DRAWBUFFERS:024 */

float GetWaves(vec3 p, float lod)
{
  vec2 pos = p.xz*4.0;
  float moving = 1.;
	vec2 movement = vec2(-0.02*frameTimeCounter*moving);
	float caustic = 0.0;
	float weightSum = 0.0;
	float radiance =  2.39996;
	mat2 rotationMatrix  = mat2(vec2(cos(radiance),  -sin(radiance)),  vec2(sin(radiance),  cos(radiance)));
	for (int i = 0; i < 4; i++){
		vec2 displ = texture2D(noisetex, pos/32.0/1.74/1.74 + movement).bb*2.0-1.0;
		pos = rotationMatrix * pos;
		caustic += sin(dot((pos+vec2(moving*frameTimeCounter))/1.74/1.74 * exp2(0.8*i) + displ*2.0,vec2(0.5)))*exp2(-0.8*i);
		weightSum += exp2(-i);
	}
	return caustic * weightSum / 300.;
}

vec3 GetWavesNormal(vec3 position, float lod) {
	float deltaPos = 0.25;

	float h0 = GetWaves(position, lod);
	float h1 = GetWaves(position + vec3(deltaPos,0.0,0.0), lod);
	float h3 = GetWaves(position + vec3(0.0,0.0,deltaPos), lod);


	float xDelta = ((h1-h0))/deltaPos*2.;
	float yDelta = ((h3-h0))/deltaPos*2.;

	vec3 wave = normalize(vec3(xDelta,yDelta,1.0-pow(abs(xDelta+yDelta),2.0)));

	return wave;
}

vec3 toScreenSpace(vec3 p) {
	vec4 iProjDiag = vec4(gbufferProjectionInverse[0].x, gbufferProjectionInverse[1].y, gbufferProjectionInverse[2].zw);
    vec3 p3 = p * 2. - 1.;
    vec4 fragposition = iProjDiag * p3.xyzz + gbufferProjectionInverse[3];
    return fragposition.xyz / fragposition.w;
}

void main() {
    vec4 fragColor = texture2D(texture, texcoord.st);

	vec4 specularData = texture2D(specular, texcoord.st);

	float roughness = 1.0 - specularData.z;
	float f0 = specularData.x;

	if(iswater > 0.5 || isice > 0.5) {
		fragColor = vec4(0.0, 0.0, 0.0, 0.2);
	}

    //store lightmap in auxilliary texture. r = torch light. g = lightning. b = sky light.
	vec2 lightmap = lmcoord.st * (1.0 / 255.0);

   	float matID = 4.0f;

	if (iswater > 0.5f)
	{
		matID = 35.0f;
	}

	if (isice > 0.5)
	{
		matID = 5.0;
	}

	if (isStainedGlass > 0.5)
	{
		matID = 55.0;
	}

	vec3 fragpos = toScreenSpace(gl_FragCoord.xyz*vec3(texelSize,1.0));

	mat3 tbnMatrix = mat3(tangent.x, binormal.x, normal.x,
						tangent.y, binormal.y, normal.y,
						tangent.z, binormal.z, normal.z);

	vec3 position = mat3(gbufferModelViewInverse) * fragpos + gbufferModelViewInverse[3].xyz;

	vec3 wavesnormal = GetWavesNormal(position + cameraPosition, 1.0);

	float bumpmult = 1.;
	vec3 bump;

	bump = normalize(GetWavesNormal(position + cameraPosition, 1.0));

	bump = bump * vec3(bumpmult, bumpmult, bumpmult) + vec3(0.0f, 0.0f, 1.0f - bumpmult);

	vec3 waternormal = normalize(bump * tbnMatrix);

	matID += 0.1f; 

	roughness = iswater > 0.5 ? 0.0 : roughness;
	f0 = iswater > 0.5 ? 0.021 : f0;

	float mMats = texture2D(colortex2, texcoord.st).b;

    gl_FragData[0] = fragColor;
	gl_FragData[1] = vec4(encodeNormal(waternormal), lightmap, 1.0);
    gl_FragData[2] = vec4(roughness, f0, matID / 255.0f, 1.0f);
}