#version 120

#define TAA

varying vec4 texcoord;
varying vec4 lmcoord;

varying vec3 color;
varying vec3 normal;

uniform vec2 texelSize;
uniform int frameMod8;
uniform int frameMod16;
uniform vec2 taaOffset16;

attribute vec4 at_tangent;
attribute vec4 mc_Entity;

varying float iswater;
varying float isice;
varying float isStainedGlass;

varying vec3 viewVector;
varying vec3 binormal;
varying vec3 tangent;
varying mat3 tbnMatrix;

uniform mat4 gbufferModelView;
uniform mat4 gbufferModelViewInverse;

varying float minecraftMats;

varying vec4 position;

		const vec2[8] offsets = vec2[8](vec2(1./8.,-3./8.),
									vec2(-1.,3.)/8.,
									vec2(5.0,1.)/8.,
									vec2(-3,-5.)/8.,
									vec2(-5.,5.)/8.,
									vec2(-7.,-1.)/8.,
									vec2(3,7.)/8.,
									vec2(7.,-7.)/8.);

void main() {

	iswater = 0.0f;
	isice = 0.0f;
	isStainedGlass = 0.0f;

    gl_Position = ftransform();

    color = gl_Color.rgb;

	vec4 viewPos = gbufferModelViewInverse * gl_ModelViewMatrix * gl_Vertex;
	position = viewPos;

	#ifdef TAA
		gl_Position.xy += offsets[frameMod8] * gl_Position.w * texelSize;
	#endif

	if (mc_Entity.x == 79) {
		isice = 1.0f;
	}

	if (mc_Entity.x == 1971.0f)
	{
		iswater = 1.0f;
	}
	
	if (mc_Entity.x == 8 || mc_Entity.x == 9) {
		iswater = 1.0f;
	}

	if(mc_Entity.x == 90) {
		iswater = 0.0f;
	}

	if (mc_Entity.x == 95 || mc_Entity.x == 160)
	{
		isStainedGlass = 1.0f;
	}

	texcoord = gl_TextureMatrix[0] * gl_MultiTexCoord0;

	lmcoord = gl_TextureMatrix[1] * gl_MultiTexCoord1;

	gl_FogFragCoord = gl_Position.z;

	normal = gl_Normal;
	normal = (gl_NormalMatrix * normal) * mat3(gbufferModelView);

	vec3 normalMat = normalize( gl_NormalMatrix*gl_Normal);

	vec3 tangent = normalize( gl_NormalMatrix *at_tangent.rgb);
	vec3 binormal = normalize(cross(tangent.rgb,normalMat.xyz)*at_tangent.w);

	tbnMatrix = mat3(tangent.x, binormal.x, normalMat.x,
								  tangent.y, binormal.y, normalMat.y,
						     	  tangent.z, binormal.z, normalMat.z);

	viewVector = (gl_ModelViewMatrix * gl_Vertex).xyz;
	viewVector = normalize(tbnMatrix * viewVector);
}