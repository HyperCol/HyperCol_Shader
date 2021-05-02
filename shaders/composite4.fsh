#version 130

#define ATOMSPHERIC_FOG
	#define ATOMSPHERIC_FOG_BRIGHTNESS 0.5	//[0.25 0.5 0.75 1.0 1.25 1.5 1.75 2.0 2.25 2.5 2.75 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0 20.0 30.0 40.0 50.0 60.0 70.0 80.0 90.0 100.0]

#define CCVA 2000 // [60 80 100 150 200 300 356 395 400 500 600 800 1000 2000]

#define CCVT 1000 // [0 100 200 300 365 400 500 550 600 800 1000 1500 2000 3000 5000]

#define CCVC -0.25 // [-0.5 -0.4 -0.3 -0.25 -0.2 -0.15 -0.1 0.0 0.1 0.2 0.3 0.4 0.5 1.0]

#define CCVD 0.5 // [0.04 0.06 0.08 0.1 0.15 0.2 0.3 0.4 0.5 0.6 0.8 1.0]

#define CCVLD 0.7 // [0.1 0.3 0.35 0.4 0.45 0.5 0.7 0.9 1.1 1.3 1.5 1.7 1.9 2.1 2.3 2.5]

#define CCVR 0.75 // [0.0 0.2 0.4 0.6 0.75 0.9 1.0 1.1 1.2 1.3 1.4 1.6 1.8 2.0]



#define CRVA 2000 // [60 80 100 150 200 300 356 400 500 600 800 1000]

#define CRVT 1000 // [0 100 200 300 400 465 500 600 800 1000 1500 2000 3000 5000]

#define CRVC 0.4 // [-0.5 -0.4 -0.3 -0.2 -0.1 0.0 0.1 0.2 0.3 0.4 0.5 1.0]

#define CRVD 0.4 // [0.04 0.06 0.08 0.1 0.15 0.2 0.3 0.4 0.5 0.6 0.8 1.0]

#define CRVLD 1.5 // [0.1 0.3 0.5 0.7 0.9 1.1 1.3 1.5 1.7 1.9 2.1 2.3 2.5]

#define CRVR 0.8 // [0.0 0.2 0.4 0.6 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.6 1.8 2.0]

//#define VOLUMETRIC_CLOUD_DEBUG

varying vec4 texcoord;

varying vec3 sunVector;
varying vec3 wSunVector;
varying vec3 moonVector;
varying vec3 wMoonVector;
varying vec3 upVector;
varying vec3 lightVector;
varying vec3 wLightVector;

varying vec3 baseSunColor;
varying vec3 sunColor;
varying vec3 baseMoonColor;
varying vec3 moonColor;
varying vec3 skyColor;

varying float transitionFading;

varying vec3 upperCloudSunlightColor;

#include "libs/uniform.glsl"
#include "libs/composite.glsl"
#include "libs/noise.glsl"

#include "libs/sky.glsl"

#include "libs/Clouds.glsl"

vec3 clacWorldposToShadowpos(vec3 worldPos) 
{
	vec4 shadowPos = vec4(worldPos, 1.0f);

	shadowPos = shadowModelView * shadowPos;	//Transform from world space to shadow space

	shadowPos = shadowProjection * shadowPos;
	shadowPos /= shadowPos.w;

	float dist = length(shadowPos.xy);
	float distortFactor = (1.0f - SHADOW_MAP_BIAS) + dist * SHADOW_MAP_BIAS;
	shadowPos.xy *= 0.95f / distortFactor;
	shadowPos.z = mix(shadowPos.z, 0.5, 0.8);

	return shadowPos.xyz * 0.5f + 0.5f;
}

//#define GODRAYS
//#define GODRAYS_STAINED_GLASS_TINT
#define GR_VL_PERCENTAGE 1.3 // [0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.7 2.0 3.0 4.0 5.0 7.0 10.0 100.0]
#define GODRAYS_STRENGTH 1.0 // [0.5 1.0 1.5 2.0 3.0 5.0 7.5 10.0 15.0 20.0 30.0 50.0 75.0 100.0 500.0 1000.0 2000.0 5000.0 10000.0]
#define GODRAYS_LENGTH 120.0 //[120.0 160.0 200.0 300.0 500.0 700.0 1000.0]
#define GODRAYS_QUALITY 15 //[10 15 20 30 50 70 100 150 200 300 500 700 1000]
//#define NOON_GODRAY

//#define fma(a,b,c) ((a)*(b)+(c))


vec3 WorldPosToShadowProjPos(vec3 worldPos, out float dist, out float distortFactor)
{
	vec4 sp = (shadowModelView * vec4(worldPos, 1.0));
	sp = shadowProjection * sp;
	sp /= sp.w;

	dist = sqrt(sp.x * sp.x + sp.y * sp.y);
	distortFactor = (1.0f - SHADOW_MAP_BIAS) + dist * SHADOW_MAP_BIAS;
	sp.xy *= 0.95f / distortFactor;
	sp.z = mix(sp.z, 0.5, 0.8);
	sp = sp * 0.5f + 0.5f;		//Transform from shadow space to shadow map coordinates




	return sp.xyz;
}

float PhaseMie( float g, float LdotV, float LdotV2 ) {
	float gg = g * g;
	
	float a = ( 1.0 - gg ) * ( 1.0 + LdotV2 );

	float b = 1.0 + gg - 2.0 * g * LdotV;
	b *= sqrt( b );
	b *= 2.0 + gg;	
	
	return 1.5 * a / b;
}


void PeGodray(inout vec3 color, vec3 worldPos, vec3 worldDir){
	 if(isEyeInWater==0)
     {
       float Y=BlueNoise(texcoord.xy),E=GODRAYS_LENGTH;
       vec3 F=vec3(0.),Q=(gbufferModelViewInverse*vec4(0.,0.,0.,1.)).xyz;
       for(int V=0;V<GODRAYS_QUALITY;V++)
         {
           float N=float(V+Y)/float(GODRAYS_QUALITY);
           vec3 J=worldDir.xyz*E*N+Q;
           if(length(worldPos.xyz)<length(J-Q))
             {
               break;
             }
           float B,j;
           vec3 q=WorldPosToShadowProjPos(J.xyz,B,j),I=shadow2DLod(shadow,vec3(q.xy,q.z+1e-06),3).xxx;
		   
		   
		   
           #ifdef GODRAYS_STAINED_GLASS_TINT

		   
		   /*
		    float shadowNormalAlpha = texture2D(shadowcolor1, q.xy).a;

			if (shadowNormalAlpha < 0.1)
			{
				#ifdef TAA_ENABLED
					vec2 noise = rand(texcoord.st + sin(frameTimeCounter)).xy;
				#else
					vec2 noise = rand(texcoord.st).xy;
				#endif

				float angle = noise.x * 3.14159 * 2.0;
				mat2 rot = mat2(cos(angle), -sin(angle), sin(angle), cos(angle));

				float solidShadowSum = 0.0;
				vec3 shadowColorSampleSum = vec3(0.0);

				int c = 0;

				for (float i = -0.25f; i <= 0.25f; i += 0.25f) 
				{
					for (float j = -0.25f; j <= 0.25f; j += 0.25f) 
					{
						q.xy += (vec2(i, j) * rot) * (0.5 / shadowMapResolution);

						vec4 shadowColorSample = texture2D(shadowcolor, q.xy);
						float opacityCheck = 1.0 - saturate(pow(shadowColorSample.a * 1.1, 4.0));
						shadowColorSampleSum += pow(shadowColorSample.rgb, vec3(1.6)) * (opacityCheck);
						float solidDepth = texture2D(shadowtex1, q.xy).x;
						float solidShadow = 1.0 - clamp((q.z - solidDepth) * 5200.0, 0.0, 1.0); 
						solidShadowSum += solidShadow;
						c++;
					}
				}

				solidShadowSum /= c;
				shadowColorSampleSum /= c;

				I = mix(vec3(1.0), shadowColorSampleSum.rgb, vec3(1.0 - I.x));
				I *= solidShadowSum;
			}
			*/
           #endif
           
		   
		   
		   F+=I*sunColor*.2;

         }
       float I=dot(wLightVector,worldDir.xyz),q=1.;
       q=.5/(max(0.,pow(wLightVector.y,2.)*2.)+.4);
       float j=I*I,u=PhaseMie(.8,I,j);
	   u=pow(u,1.0/(GR_VL_PERCENTAGE * fma(transitionFading, 0.7, 1.0)))*3*GODRAYS_STRENGTH;
	   float mtl = 1.0;
	   #ifndef NOON_GODRAY
		   mtl *= mix(pow((1.0-transitionFading),0.2), 1.0, wetness);
	   #endif
       vec3 godRays=F*mix(sunColor*mtl,vec3(Luminance(sunColor)*10000.0),transitionFading)*2e-5*u*q*E/GODRAYS_QUALITY;

	   color += godRays;
     }
}

/* DRAWBUFFERS:0 */

void main() {
    surfaceStruct surface = GetSurfaceStruct();
	maskStruct mask = GetMaskStruct(texcoord.st);
	CloudProperties cloudProperties = GetGlobalCloudProperties();
    vec3 color = texture2D(colortex0, texcoord.st).rgb;

	vec2 planetSphere = vec2(0.0);
	vec3 skyAbsorb = vec3(0.0);

	vec3 atmosphere = vec3(0.0);
	atmosphere = calculateAtmosphere(vec3(0.0), surface.viewDir, upVector, sunVector, moonVector, planetSphere, skyAbsorb, 25) * 0.000045;

	#ifdef VOLUMETRIC_CLOUD_DEBUG
		VolumetricClouds(color, (surface.depth > 0.99999) ? surface.worldPos.xyz * 10000.0 : surface.worldPos.xyz, surface.worldDir.xyz, cloudProperties, atmosphere);
	#endif

	float totalDist = 0.0, anisoFactor = 0.0;
	#ifdef ATOMSPHERIC_RAYS
		CrepuscularRays(color, surface, mask);
	#endif

	#ifdef GODRAYS
		PeGodray(color, surface.worldPos.xyz, surface.worldDir.xyz);
	#endif

	#ifdef ATOMSPHERIC_FOG
		float cdepthN = length(surface.viewPos.xyz) * 0.075;
		float fog_density = min(cdepthN / (1024.0 - rainStrength * 512.0), 1.0);
		if(isEyeInWater == 0)
		{
			float eyeLength = length(surface.worldPos.xyz - cameraPosition);
			vec3 fogcolor = upperCloudSunlightColor * 1.215;
			fogcolor += 0.65 * skyColor * 0.00005;
			fogcolor += atmosphere * (1.0 - exp2(-eyeLength * 0.00038 * mix(1.0, 0.45, rainStrength)));
			if(mask.land > 0.5) color += fogcolor * fog_density * 4.5 * ATOMSPHERIC_FOG_BRIGHTNESS;
		}
	#endif

    gl_FragData[0] = vec4(color, texture2D(colortex5, texcoord.st).a);
}