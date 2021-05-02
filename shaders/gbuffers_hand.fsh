#version 120

varying vec4 texcoord;
varying vec4 lmcoord;

varying vec4 color;
varying vec3 normal;

varying float terrainMats;

varying float n2;

#include "libs/uniform.glsl"
#include "libs/composite.glsl"

//#define SPECULAR_TO_PBR_CONVERSION
#define CONTINUUM2_TEXTURE_FORMAT

/* DRAWBUFFERS:024 */

void main() {
    vec4 fragColor = texture2D(texture, texcoord.st) * color;
	
	vec4 sp;
	
	#ifdef SPECULAR_TO_PBR_CONVERSION
	vec3 spec = texture2D(specular, texcoord.st);
	float spec_strength = dot(spec, mix(vec3(0.4, 0.4, 0.2), vec3(0.3, 0.6, 0.1), wetness));
	sp = vec4(spec_strength, spec_strength, 0.0, 1.0);
	#else
	#ifdef CONTINUUM2_TEXTURE_FORMAT
	sp = texture2D(specular, texcoord.st).brga;
	#else
	sp = texture2D(specular, texcoord.st);
	#endif
	#endif
	
	float roughness = 1.0 - sp.x;
	float f0 = sp.g;
	//vec4 specularData = texture2D(specular, texcoord.st);

	//float roughness = 1.0 - specularData.z;
	//float f0 = specularData.x;
    
    //store lightmap in auxilliary texture. r = torch light. g = lightning. b = sky light.
	vec2 lightmap = lmcoord.st * (1.0 / 255.0);
	
    gl_FragData[0] = fragColor;
    gl_FragData[1] = vec4(n2, lightmap, 1.0);
    gl_FragData[2] = vec4(roughness, f0, terrainMats / 255.0f, 1.0f);
}