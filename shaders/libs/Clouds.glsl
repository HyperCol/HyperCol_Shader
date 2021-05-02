#define yyqd 0.5f // [-0.5f 0f 0.3f 0.5f]

struct Ray {
	vec3 dir;
	vec3 origin;
};

struct Plane {
	vec3 normal;
	vec3 origin;
};

struct Intersection {
	vec3 pos;
	float distance;
	float angle;
};

Intersection 	RayPlaneIntersection(in Ray ray, in Plane plane)
{
	float rayPlaneAngle = dot(ray.dir, plane.normal);

	float planeRayDist = 100000000.0f;
	vec3 intersectionPos = ray.dir * planeRayDist;

	if (rayPlaneAngle > 0.0001f || rayPlaneAngle < -0.0001f)
	{
		planeRayDist = dot((plane.origin - ray.origin), plane.normal) / rayPlaneAngle;
		intersectionPos = ray.origin + ray.dir * planeRayDist;
		// intersectionPos = -intersectionPos;

		// intersectionPos += cameraPosition.xyz;
	}

	Intersection i;

	i.pos = intersectionPos;
	i.distance = planeRayDist;
	i.angle = rayPlaneAngle;

	return i;
}

//////////////////////////////////////////////////////////////
////////////////////////////Noise/////////////////////////////
//////////////////////////////////////////////////////////////

float BlueNoise(vec2 coord)
{
	vec2 noiseCoord = vec2(coord.st * vec2(viewWidth, viewHeight)) / 64.0;
	noiseCoord += vec2(sin(frameCounter * 0.75), cos(frameCounter * 0.75));

	noiseCoord = (floor(noiseCoord * 64.0) + 0.5) / 64.0;

	float blueNoise = texture2DLod(noisetex, noiseCoord.st, 0).b;

	return blueNoise;
}

float Get2DNoiseM(in vec3 pos)
{
	pos.xy = pos.xz;
	pos.xy += 0.5f;

	vec2 p = floor(pos.xy);
	vec2 f = fract(pos.xy);

	f *= f * (3.0f - 2.0f * f);

	vec2 uv =  p.xy + f.xy;

	// uv -= 0.5f;
	// uv2 -= 0.5f;

	vec2 coord =  (uv  + 0.5f) / 64.0;
	float xy1 = texture2D(noisetex, coord).x;
	return xy1;
}

float Get3DNoiseM(in vec3 pos)
{
	pos.xyz += 0.5f;

	vec3 p = floor(pos);
	vec3 f = fract(pos);

	f *= f * (3.0f - 2.0f * f);

	vec2 uv =  (p.xy + p.z * vec2(17.0f)) + f.xy;
	vec2 uv2 = (p.xy + (p.z + 1.0f) * vec2(17.0f)) + f.xy;

	vec2 coord =  (uv + 0.5f) / noiseTextureResolution;
	vec2 coord2 = (uv2 + 0.5f) / noiseTextureResolution;
	float xy1 = texture2D(noisetex, coord).x;
	float xy2 = texture2D(noisetex, coord2).x;
	return mix(xy1, xy2, f.z);
}


//Optimize by not needing to do a dot product every time
float MiePhase(float g, vec3 dir, vec3 lightDir)
{
	float VdotL = dot(dir, lightDir);

	float g2 = g * g;
	float theta = VdotL * 0.5 + 0.5;
	float anisoFactor = 1.5 * ((1.0 - g2) / (2.0 + g2)) * ((1.0 + theta * theta) / (1.0 + g2 - 2.0 * g * theta)) + g * theta;

	return anisoFactor;
}

//////////////////////////////////////////////////////////////
////////////////////////////Struct////////////////////////////
//////////////////////////////////////////////////////////////

struct CloudProperties
{
	float altitude;
	float thickness;
	float coverage;
	float density;
	float lightingDensity;
	float roughness;
};

CloudProperties CloudPropertiesLerp(CloudProperties cp1, CloudProperties cp2, float x)
{
	CloudProperties cp;

	cp.altitude = mix(cp1.altitude, cp2.altitude, x);
	cp.thickness = mix(cp1.thickness, cp2.thickness, x);
	cp.coverage = mix(cp1.coverage, cp2.coverage, x);
	cp.density = mix(cp1.density, cp2.density, x);
	cp.lightingDensity = mix(cp1.lightingDensity, cp2.lightingDensity, x);
	cp.roughness = mix(cp1.roughness, cp2.roughness, x);

	return cp;
}


CloudProperties GetGlobalCloudProperties()
{
	CloudProperties cloudPropertiesClear;
	cloudPropertiesClear.altitude = CCVA;
	cloudPropertiesClear.thickness = CCVT;
	cloudPropertiesClear.coverage = CCVC;
	// cloudPropertiesClear.coverage = mix(-0.3, 0.5, sin(frameTimeCounter * 0.5) * 0.5 + 0.5);
	cloudPropertiesClear.density = CCVD;
	cloudPropertiesClear.lightingDensity = CCVLD;
	cloudPropertiesClear.roughness = CCVR;

	CloudProperties cloudPropertiesRain;
	cloudPropertiesRain.altitude = CRVA;
	cloudPropertiesRain.thickness = CRVT;
	cloudPropertiesRain.coverage = CRVC;
	cloudPropertiesRain.density = CRVD;
	cloudPropertiesRain.lightingDensity = CRVLD;
	cloudPropertiesRain.roughness = CRVR;

	CloudProperties cp = CloudPropertiesLerp(cloudPropertiesClear, cloudPropertiesRain, rainStrength);

	return cp;
}

//////////////////////////////////////////////////////////////
////////////////////////////Clouds////////////////////////////
//////////////////////////////////////////////////////////////


float GetCloudHeightDensity(CloudProperties cp, vec3 worldPos, bool forLighting)
{
	float w = 0.5;
	float highPoint = cp.altitude + cp.thickness * 1.0;
	float lowPoint =  cp.altitude;
	float midPoint = mix(lowPoint, highPoint, 0.2);

	// float density = smoothstep(0.0, 1.0, 1.0 - abs(worldPos.y - cp.altitude) / cp.thickness);
	// density *= 1.0 - remap(lowPoint, highPoint, worldPos.y);
	// density = pow(density, 0.25);

	float density = remap(lowPoint, midPoint, worldPos.y);
	density *= remap(highPoint, midPoint, worldPos.y);
	density = smoothstep(0.0, 1.0, density);
	// density = pow(density, 0.1251);


	// density /= remap(highPoint, midPoint2, worldPos.y) + 0.1;
	// density /= remap(highPoint, midPoint2, worldPos.y) + 0.1;


	return density;
}

float GetCloudDensity(CloudProperties cp, vec3 worldPos, float roughness, vec2 clampOffsets, bool forLighting, float detailNoise)
{
	vec3 pos = worldPos * 0.00195;
	float noiseGain = 1.0;
	float noiseTotal = 0.0;

	pos.x -= frameTimeCounter * 0.01;
	pos.z -= frameTimeCounter * 0.001;
	pos.y += frameTimeCounter * 0.0025;



	roughness *= 1.0 * cp.roughness;

	float noise = 0.0;
	noise =  mix(Get2DNoiseM(pos * vec3(1.0, 1.0, 0.7)), 1.0, saturate(cp.coverage));  			 		noiseTotal += 1.0;  		pos *= 3.5; noiseGain *= 0.4 * roughness; pos.xy -= frameTimeCounter * 0.1;
	// noise =  Get2DNoiseM(pos) * (cp.coverage + 1.0);  			 		noiseTotal += 1.0;  		pos *= 4.5; noiseGain *= 0.4 * roughness; pos.xy -= frameTimeCounter * 0.01;
	noise *= GetCloudHeightDensity(cp, worldPos, forLighting);
	noise -= Get3DNoiseM(pos) * noiseGain * 1.0;  				noiseTotal-= noiseGain;    	pos *= 4.0; noiseGain *= 0.45 * roughness; 									
	noise -= Get3DNoiseM(pos) * noiseGain;  						noiseTotal-= noiseGain;    	pos *= 4.0; noiseGain *= 0.45 * roughness; 									




	noise /= noiseTotal;
	noise *= 0.9;

	noise -= saturate(-cp.coverage);


	if (!forLighting)
	{
		noise -= (detailNoise - 1.33) * 0.4 * roughness;
	}


	noise = smoothstep(0.3 - clampOffsets.x, 0.85 + clampOffsets.y, noise);



	//clouds more dense at top
	float w = forLighting ? 0.7 : 0.5;
	float highPoint = cp.altitude + cp.thickness * 1.0;
	float lowPoint =  cp.altitude;

	if (!forLighting)
	noise = saturate(noise * mix(1.0, 25.0, pow(remap(lowPoint, highPoint, worldPos.y), 2.0)));



	return noise;
}

vec4 SampleVolumetricClouds(CloudProperties cp, vec3 worldPos, vec3 worldDir, float detailNoise, vec3 atmosphere, float rayIncrementLength)
{
	float eyeLength = length(worldPos - cameraPosition);


	float detailNoiseScaled = mix(detailNoise, 2.7, 1.0 / (eyeLength * 0.00025 + 1.0)) * 0.5;
	float detailNoiseScaled2 = mix(detailNoise, 0.0, 1.0 / (eyeLength * 0.00025 + 1.0)) * 0.5;


	float 	density = GetCloudDensity(cp, worldPos, 1.0, vec2(0.0), false, detailNoiseScaled);
	float softDensity = density;
	 		density = smoothstep(0.0, mix(1.0, 0.01, cp.density), density); 	//final cloud edge resolve


	// vec3 worldPosDetailed = worldPos + worldDir * (detailNoiseScaled2) * 0.0;
	 vec3 worldPosDetailed = worldPos;



	//Early out if no cloud here
	if (density < 0.0001)
	{
		return vec4(0.0, 0.0, 0.0, 0.0);
	}

	float approxHighPoint = cp.altitude + cp.thickness * (cp.coverage + 4.0) * 0.2;
	float approxLowPoint = cp.altitude;
	float vgrad = remap(approxLowPoint, approxHighPoint, worldPosDetailed.y + detailNoise * 0.4);

	float sunAngle = dot(wLightVector, worldDir) * 0.5 + 0.5;







	//lighting..
	float sunlightExtinction = 0.0;
	float lightingRoughness = 1.0;
	for (int i = 1; i <= 3; i++)
	{
		float fi = float(i) / 2;
		fi = pow(fi, 1.5);
		vec3 checkPos = worldPosDetailed + wLightVector * fi * 150.0 /*make clouds look denser towards sun * (1.5 - sunAngle)*/;

		float densityCheck = GetCloudDensity(cp, checkPos, lightingRoughness, vec2(0.0, 0.2), true, 0.0);

		sunlightExtinction += densityCheck * 2.0;

		lightingRoughness *= 0.5;
	}
	float sunlightEnergy = 1.0 / (sunlightExtinction * 7.0 * cp.lightingDensity + 1.0);

	//powder term
	float powderFactor = exp2(-softDensity * 10.0 * cp.lightingDensity);
	powderFactor *= saturate(sunlightEnergy * sunlightEnergy * 10.5);
	sunlightEnergy *= MiePhase(powderFactor * 0.75, worldDir, wLightVector);
	// sunlightEnergy *= 1.0 - powderFactor;

	//Extra edge highlight
	sunlightEnergy *= MiePhase(0.65, worldDir, wLightVector) * sunlightEnergy * 0.0045 + 1.0;








	//fast approximate skylight energy
	// float skylightEnergy = 1.0 / (pow(approxHighPoint - worldPosDetailed.y, 2.0) * 0.0004 * cp.lightingDensity + 1.0);
	// skylightEnergy /= n2 + 0.1;
	// skylightEnergy /= n2 + 0.1;

	//Cone trace skylight energy
	float skylightExtinction = 0.0;
	lightingRoughness = 0.75;
	for (int i = 1; i < 2; i++)
	{
		float fi = float(i) / 2;
		vec3 checkPos = worldPosDetailed + vec3(0.0, 1.0, 0.0) * fi * 100.0;

		float densityCheck = GetCloudDensity(cp, checkPos, lightingRoughness, vec2(0.1, 0.0), true, 0.0);

		skylightExtinction += densityCheck;
	}
	float skylightEnergy = 1.0 / (skylightExtinction * 4.0 * cp.lightingDensity + 1.0);







	float bouncelightEnergy = pow(1.0 - skylightEnergy, 2.0);



	// vec3 sunlightColor = exp2(-vec3(0.3, 0.55, 1.0) * 0.25 * mix(1.5, 1.0, vgrad) / (worldLightVector.y + 0.01));
	vec3 sunlightColor = upperCloudSunlightColor;


	vec3 cloudColor = sunlightEnergy * sunlightColor * 1.215;
		 cloudColor += skylightEnergy * 0.65 * skyColor * 0.00005;
	// cloudColor += vec3(0.7, 0.9, 0.4) * (bouncelightEnergy) * 1.7 * (1.0 - rainStrength) * colorSunlight * worldLightVector.y;
	// cloudColor += skylightEnergy * upperCloudSunlightColor * 0.5; 	//Bounced light from upper cloud layer

	//fake scattered light
	// cloudColor += upperCloudSunlightColor * saturate(1.0 - softDensity * 0.8);


	//fake bounce lighting
	// cloudColor *= vec3(1.0) * (1.8 - softDensity);


	// rayleigh absorption and scattering
	cloudColor *= exp2(-eyeLength * vec3(0.3, 0.55, 1.0) * 0.00038 * mix(1.25, 1.0, vgrad));
	// cloudColor += 1.0 - exp2(-eyeLength * vec3(0.3, 0.55, 1.0) * 0.00048 + 0.0);
	cloudColor += atmosphere * (1.0 - exp2(-eyeLength * 0.00038 * mix(1.0, 0.45, rainStrength)));

	// cloudColor = vec3(powderFactor);


	// cloudColor = vec3(skylightEnergy * 10.0);

	// density = saturate(density * rayIncrementLength * 0.05);

	return vec4(cloudColor, density);
}

float GetDirectionalDetailNoise(vec3 worldDir)
{
	vec3 pos = worldDir * 1.0;

	pos.x -= frameTimeCounter * 0.005;
	pos.y += frameTimeCounter * 0.005;

	/*float pnoise = 	Get3DNoiseL(pos * 30.0) * 1.5; 			pos.xy += pnoise * 0.0013;
	pnoise +=  		Get3DNoiseL(pos * 60.0) * 1.5;			pos.xz -= pnoise * 0.0013;
	pnoise +=  		Get3DNoiseL(pos * 120.0) * 0.75;			pos.zy -= pnoise * 0.0013;
	pnoise +=  		Get3DNoiseL(pos * 240.0) * 0.375;			pos.xy += pnoise * 0.0013;
	pnoise +=  		Get3DNoiseL(pos * 480.0) * 0.1875;		pos.xy += pnoise * 0.0013;

	// pnoise *= 1.1;

	return pnoise;*/

	const int octaves = 10;
	const float octAlpha = 0.75; // The ratio of visibility between successive octaves
    const float octScale = 3.0; // The downscaling factor between successive octaves
    const float octShift = (octAlpha / octScale) / octaves; // Shift the FBM brightness based on how many octaves are active

    float accum = 0.0;
    float alpha = 2.0;
    
    for (int i = 0; i < octaves; ++i) {
        accum += alpha * Get3DNoiseL(pos * 20.0);
        pos = pos * octScale;
        alpha *= octAlpha;
    }
    
    return accum;
}

vec3 IntersectXZPlane(vec3 rayDir, float planeHeight)
{
	return rayDir * (planeHeight - cameraPosition.y) / rayDir.y;
}

void VolumetricClouds(inout vec3 color, in vec3 worldPos, in vec3 worldDir, CloudProperties cloudProperties, in vec3 atmosphere)
{




	float detailNoise = GetDirectionalDetailNoise(worldDir);

	//no clouds in front of sun
	// cloudProperties.coverage *= mix(0.0, 1.0, 1.0 - pow(saturate(dot(worldDir, worldLightVector) * 0.5 + 0.5), 48.0));
	CloudProperties cp = cloudProperties;

	cp.coverage += pow(dot(-worldDir, wLightVector) * 0.5 + 0.5, 2.0) * 0.3;

#if VC_MARCHING


	const int raySteps = 10;
	const float rayStepSize = 1.0 / raySteps;
	const float actualCloudHeightGuess = 0.6;

	vec3 rayStartPos = IntersectXZPlane(worldDir, cloudProperties.altitude + cloudProperties.thickness * 0.05) + cameraPosition;

	vec3 rayIncrement = (worldDir * cloudProperties.thickness * actualCloudHeightGuess * rayStepSize) / abs(worldDir.y);
	// vec3 rayIncrement = worldDir * 40.1;
	float rayIncrementLength = length(rayIncrement);

	vec4 cloudAccum = vec4(0.0);
	float visibilityAccum = 1.0;

	vec3 rayPos = rayStartPos;
	rayPos += rayIncrement * BlueNoise(texcoord.st);


	// cp.coverage -= length(rayStartPos.xz - cameraPosition.xz) * 0.0001;

	for (int i = 0; i < raySteps; i++)
	{
		//When fully opaque, stop
		if (length(rayPos - cameraPosition) > length(worldPos) || visibilityAccum < 0.00001)
		{
			break;
		}

		vec4 cloudSample = SampleVolumetricClouds(cp, rayPos, worldDir, detailNoise, atmosphere, rayIncrementLength);
		cloudAccum.rgb += (cloudSample.rgb * cloudSample.a) * visibilityAccum;
		visibilityAccum *= 1.0 - cloudSample.a;

		rayPos += rayIncrement;
	}


#else

	const int raySteps = 60;
	const float rayExtent = 24000.0f;
	const float rayStepSize = rayExtent / raySteps;


	vec4 cloudAccum = vec4(0.0);
	float visibilityAccum = 1.0;

	// float rayDepth = 0.0;
	float rayDepth = 0.0 + BlueNoise(texcoord.st) * rayStepSize * 0.2;
	// rayDepth += abs(cameraPosition.y - cloudProperties.altitude + cloudProperties.thickness * 0.2);

	for (int i = 0; i < raySteps; i++)
	{
		vec3 rayPos = (worldDir * rayDepth) + cameraPosition;

		// cloudProperties.thickness *= 1.05;
		// cp.altitude *= 0.99;
		// cloudProperties.coverage *= 1.03;
		// cloudProperties.density = saturate(cloudProperties.density * 1.05);
		// cloudProperties.density = 0.5;

		// cloudProperties.thickness = mix(200.0, 900.0, Get2DNoise(rayPos * 0.001));

		// cloudProperties.thickness = 200.0 + length(rayPos - cameraPosition) * 0.1;
		cp.altitude = 376.0 - length(rayPos - cameraPosition) * 0.0525;
		// cloudProperties.coverage = mix(-0.1, 0.5, Get2DNoise(rayPos * 0.0015));



		//When fully opaque, stop
		if (length(rayPos - cameraPosition) > length(worldPos) || visibilityAccum < 0.00001)
		{
			break;
		}

		//early out if not in correct altitude range
		if (GetCloudHeightDensity(cp, rayPos, false) < 0.5)
		{
			rayDepth += rayStepSize;
			continue;
		}

		vec4 cloudSample = SampleVolumetricClouds(cp, rayPos, worldDir, detailNoise, atmosphere, 1.0);
		cloudAccum.rgb += (cloudSample.rgb * cloudSample.a) * visibilityAccum;
		visibilityAccum *= 1.0 - cloudSample.a;

		rayDepth += rayStepSize;
	}

#endif
	//apply clouds
	color = (color * visibilityAccum) + cloudAccum.rgb;

	// color.rgb = vec3(detailNoise * 0.1);
}

//////////////////////////////////////////////////////////////
////////////////////////////Clouds Shadow/////////////////////
//////////////////////////////////////////////////////////////

float CloudShadow(vec4 worldPos, vec3 worldLightVector, CloudProperties cloudProperties)
{
	Ray ray;
	ray.dir = worldLightVector;
	ray.origin = worldPos.xyz + cameraPosition.xyz;

	Plane plane;
	plane.normal = vec3(0.0, 1.0, 0.0);
	plane.origin = vec3(0.0, cloudProperties.altitude + cloudProperties.thickness * 0.2, 0.0);

	vec3 cloudCheckPos = RayPlaneIntersection(ray, plane).pos;

	float cloudDensity = GetCloudDensity(cloudProperties, cloudCheckPos, 1.0, vec2(0.2, -0.2), true, 0.0);

	return 1.0f - cloudDensity * 0.5f + yyqd*1.0f;
}