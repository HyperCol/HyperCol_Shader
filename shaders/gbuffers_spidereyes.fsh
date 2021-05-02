#version 120

varying vec4 texcoord;
varying vec4 lmcoord;

varying vec4 color;
varying vec3 normal;

varying float terrainMats;

#include "libs/uniform.glsl"
#include "libs/composite.glsl"

/* DRAWBUFFERS:024 */

void main() {
    vec4 fragColor = texture2D(texture, texcoord.st) * color;

	vec4 specularData = texture2D(specular, texcoord.st);

	float roughness = 1.0 - specularData.z;
	float f0 = specularData.x;
    
    //store lightmap in auxilliary texture. r = torch light. g = lightning. b = sky light.
	vec2 lightmap = lmcoord.st * (1.0 / 255.0);
	
    gl_FragData[0] = fragColor * 0.00000001;
    gl_FragData[1] = vec4(encodeNormal(normal), lightmap, 1.0);
    gl_FragData[2] = vec4(roughness, f0, terrainMats / 255.0f, 1.0f);
}