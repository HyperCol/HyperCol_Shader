const float sunAngularSize = 0.533333;
const float moonAngularSize = 0.516667;

//Sky coefficients and heights

#define airNumberDensity 2.5035422e25 // m^3
#define ozoneConcentrationPeak 8e-6
const float ozoneNumberDensity = airNumberDensity * ozoneConcentrationPeak;
const vec3 ozoneCrossSection = vec3(4.51103766177301E-21, 3.2854797958699E-21, 1.96774621921165E-22); // cm^2 | single-wavelength values.

#define sky_planetRadius 6731e3

#define sky_atmosphereHeight 110e3
#define sky_scaleHeights vec2(8.0e3, 1.2e3)

#define sky_mieg 0.80

#define sky_coefficientRayleighR 5.8 //[0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0 2.1 2.2 2.3 2.4 2.5 2.6 2.7 2.8 2.9 3.0 3.1 3.2 3.3 3.4 3.5 3.6 3.7 3.8 3.9 4.0 4.1 4.2 4.3 4.4 4.5 4.6 4.7 4.8 4.9 5.0 5.1 5.2 5.3 5.4 5.5 5.6 5.7 5.8 5.9 6.0 6.1 6.2 6.3 6.4 6.5 6.6 6.7 6.8 6.9 7.0 7.1 7.2 7.3 7.4 7.5 7.6 7.7 7.8 7.9 8.0 8.1 8.2 8.3 8.4 8.5 8.6 8.7 8.8 8.9 9.0 9.1 9.2 9.3 9.4 9.5 9.6 9.7 9.8 9.9 10.0 ]
#define sky_coefficientRayleighG 1.35 //[0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0 2.1 2.2 2.3 2.4 2.5 2.6 2.7 2.8 2.9 3.0 3.1 3.2 3.3 3.4 3.5 3.6 3.7 3.8 3.9 4.0 4.1 4.2 4.3 4.4 4.5 4.6 4.7 4.8 4.9 5.0 5.1 5.2 5.3 5.4 5.5 5.6 5.7 5.8 5.9 6.0 6.1 6.2 6.3 6.4 6.5 6.6 6.7 6.8 6.9 7.0 7.1 7.2 7.3 7.4 7.5 7.6 7.7 7.8 7.9 8.0 8.1 8.2 8.3 8.4 8.5 8.6 8.7 8.8 8.9 9.0 9.1 9.2 9.3 9.4 9.5 9.6 9.7 9.8 9.9 10.0 ]
#define sky_coefficientRayleighB 3.31 //[0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0 2.1 2.2 2.3 2.4 2.5 2.6 2.7 2.8 2.9 3.0 3.1 3.2 3.3 3.4 3.5 3.6 3.7 3.8 3.9 4.0 4.1 4.2 4.3 4.4 4.5 4.6 4.7 4.8 4.9 5.0 5.1 5.2 5.3 5.4 5.5 5.6 5.7 5.8 5.9 6.0 6.1 6.2 6.3 6.4 6.5 6.6 6.7 6.8 6.9 7.0 7.1 7.2 7.3 7.4 7.5 7.6 7.7 7.8 7.9 8.0 8.1 8.2 8.3 8.4 8.5 8.6 8.7 8.8 8.9 9.0 9.1 9.2 9.3 9.4 9.5 9.6 9.7 9.8 9.9 10.0 ]

#define sky_coefficientRayleigh vec3(sky_coefficientRayleighR*1e-6, sky_coefficientRayleighG*1e-5, sky_coefficientRayleighB*1e-5)


#define sky_coefficientMieR 3.0 //[0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0 2.1 2.2 2.3 2.4 2.5 2.6 2.7 2.8 2.9 3.0 3.1 3.2 3.3 3.4 3.5 3.6 3.7 3.8 3.9 4.0 4.1 4.2 4.3 4.4 4.5 4.6 4.7 4.8 4.9 5.0 5.1 5.2 5.3 5.4 5.5 5.6 5.7 5.8 5.9 6.0 6.1 6.2 6.3 6.4 6.5 6.6 6.7 6.8 6.9 7.0 7.1 7.2 7.3 7.4 7.5 7.6 7.7 7.8 7.9 8.0 8.1 8.2 8.3 8.4 8.5 8.6 8.7 8.8 8.9 9.0 9.1 9.2 9.3 9.4 9.5 9.6 9.7 9.8 9.9 10.0 ]
#define sky_coefficientMieG 3.0 //[0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0 2.1 2.2 2.3 2.4 2.5 2.6 2.7 2.8 2.9 3.0 3.1 3.2 3.3 3.4 3.5 3.6 3.7 3.8 3.9 4.0 4.1 4.2 4.3 4.4 4.5 4.6 4.7 4.8 4.9 5.0 5.1 5.2 5.3 5.4 5.5 5.6 5.7 5.8 5.9 6.0 6.1 6.2 6.3 6.4 6.5 6.6 6.7 6.8 6.9 7.0 7.1 7.2 7.3 7.4 7.5 7.6 7.7 7.8 7.9 8.0 8.1 8.2 8.3 8.4 8.5 8.6 8.7 8.8 8.9 9.0 9.1 9.2 9.3 9.4 9.5 9.6 9.7 9.8 9.9 10.0 ]
#define sky_coefficientMieB 3.0 //[0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0 2.1 2.2 2.3 2.4 2.5 2.6 2.7 2.8 2.9 3.0 3.1 3.2 3.3 3.4 3.5 3.6 3.7 3.8 3.9 4.0 4.1 4.2 4.3 4.4 4.5 4.6 4.7 4.8 4.9 5.0 5.1 5.2 5.3 5.4 5.5 5.6 5.7 5.8 5.9 6.0 6.1 6.2 6.3 6.4 6.5 6.6 6.7 6.8 6.9 7.0 7.1 7.2 7.3 7.4 7.5 7.6 7.7 7.8 7.9 8.0 8.1 8.2 8.3 8.4 8.5 8.6 8.7 8.8 8.9 9.0 9.1 9.2 9.3 9.4 9.5 9.6 9.7 9.8 9.9 10.0 ]

#define sky_coefficientMie vec3(sky_coefficientMieR*1e-6, sky_coefficientMieG*1e-6, sky_coefficientMieB*1e-6) // Should be >= 2e-6


const vec3 sky_coefficientOzone = (ozoneCrossSection * (ozoneNumberDensity * 2e-6)); // ozone cross section * (ozone number density * (cm ^ 3))

const vec2 sky_inverseScaleHeights = 1.0 / sky_scaleHeights;
const vec2 sky_scaledPlanetRadius = sky_planetRadius * sky_inverseScaleHeights;
const float sky_atmosphereRadius = sky_planetRadius + sky_atmosphereHeight;
const float sky_atmosphereRadiusSquared = sky_atmosphereRadius * sky_atmosphereRadius;

#define sky_coefficientsScattering mat2x3(sky_coefficientRayleigh, sky_coefficientMie)
const mat3   sky_coefficientsAttenuation = mat3(sky_coefficientRayleigh, sky_coefficientMie * 1.11, sky_coefficientOzone); // commonly called the extinction coefficient

#define sun_illuminance 128000.0
#define moon_illuminance 2.0

#define sunColorR 1.0 //[0.0 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1 0.11 0.12 0.13 0.14 0.15 0.16 0.17 0.18 0.19 0.2 0.21 0.22 0.23 0.24 0.25 0.26 0.27 0.28 0.29 0.3 0.31 0.32 0.33 0.34 0.35 0.36 0.37 0.38 0.39 0.4 0.41 0.42 0.43 0.44 0.45 0.46 0.47 0.48 0.49 0.5 0.51 0.52 0.53 0.54 0.55 0.56 0.57 0.58 0.59 0.6 0.61 0.62 0.63 0.64 0.65 0.66 0.67 0.68 0.69 0.7 0.71 0.72 0.73 0.74 0.75 0.76 0.77 0.78 0.79 0.8 0.81 0.82 0.83 0.84 0.85 0.86 0.87 0.88 0.89 0.9 0.91 0.92 0.93 0.94 0.95 0.96 0.97 0.98 0.99 1.0 ]
#define sunColorG 0.89 //[0.0 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1 0.11 0.12 0.13 0.14 0.15 0.16 0.17 0.18 0.19 0.2 0.21 0.22 0.23 0.24 0.25 0.26 0.27 0.28 0.29 0.3 0.31 0.32 0.33 0.34 0.35 0.36 0.37 0.38 0.39 0.4 0.41 0.42 0.43 0.44 0.45 0.46 0.47 0.48 0.49 0.5 0.51 0.52 0.53 0.54 0.55 0.56 0.57 0.58 0.59 0.6 0.61 0.62 0.63 0.64 0.65 0.66 0.67 0.68 0.69 0.7 0.71 0.72 0.73 0.74 0.75 0.76 0.77 0.78 0.79 0.8 0.81 0.82 0.83 0.84 0.85 0.86 0.87 0.88 0.89 0.9 0.91 0.92 0.93 0.94 0.95 0.96 0.97 0.98 0.99 1.0 ]
#define sunColorB 0.8 //[0.0 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1 0.11 0.12 0.13 0.14 0.15 0.16 0.17 0.18 0.19 0.2 0.21 0.22 0.23 0.24 0.25 0.26 0.27 0.28 0.29 0.3 0.31 0.32 0.33 0.34 0.35 0.36 0.37 0.38 0.39 0.4 0.41 0.42 0.43 0.44 0.45 0.46 0.47 0.48 0.49 0.5 0.51 0.52 0.53 0.54 0.55 0.56 0.57 0.58 0.59 0.6 0.61 0.62 0.63 0.64 0.65 0.66 0.67 0.68 0.69 0.7 0.71 0.72 0.73 0.74 0.75 0.76 0.77 0.78 0.79 0.8 0.81 0.82 0.83 0.84 0.85 0.86 0.87 0.88 0.89 0.9 0.91 0.92 0.93 0.94 0.95 0.96 0.97 0.98 0.99 1.0 ]

#define sunColorBase (vec3(sunColorR,sunColorG,sunColorB) * sun_illuminance)
#define moonColorBase (vec3(0.36333, 0.56333, 0.92333) * moon_illuminance )    //Fake Purkinje effect

/*******************************************************************************
- Phases
******************************************************************************/

#define phaseg0 0.25

float sky_rayleighPhase(float cosTheta) {
	const vec2 mul_add = vec2(0.1, 0.28) * rPI;
	return cosTheta * mul_add.x + mul_add.y; // optimized version from [Elek09], divided by 4 pi for energy conservation
}

float sky_miePhase(float cosTheta, const float g) {
	float gg = g * g;
	return (gg * -0.25 + 0.25) * rPI * pow(-(2.0 * g) * cosTheta + (gg + 1.0), -1.5);
}

vec2 sky_phase(float cosTheta, const float g) {
	return vec2(sky_rayleighPhase(cosTheta), sky_miePhase(cosTheta, g));
}

vec3 sky_density(float centerDistance) {
	vec2 rayleighMie = exp(centerDistance * -sky_inverseScaleHeights + sky_scaledPlanetRadius);

	float ozone = exp(-max(0.0, (35000.0 - centerDistance) - sky_planetRadius) * 0.0002)
				* exp(-max(0.0, (centerDistance - 35000.0) - sky_planetRadius) / 15000.0);

	return vec3(rayleighMie, ozone);
}

vec3 sky_airmass(vec3 position, vec3 direction, float rayLength, const float steps) {
	float stepSize  = rayLength * (1.0 / steps);
	vec3  increment = direction * stepSize;
	position += increment * 0.5;

	vec3 thickness = vec3(0.0);
	for (float i = 0.0; i < steps; ++i, position += increment) {
		thickness += sky_density(length(position));
	}

	return thickness * stepSize;
}

vec3 sky_airmass(vec3 position, vec3 direction, const float steps) {
	float rayLength = dot(position, direction);
	      rayLength = rayLength * rayLength + sky_atmosphereRadiusSquared - dot(position, position);
		  if (rayLength < 0.0) return vec3(0.0);
	      rayLength = sqrt(rayLength) - dot(position, direction);

	return sky_airmass(position, direction, rayLength, steps);
}

vec3 sky_opticalDepth(vec3 position, vec3 direction, float rayLength, const float steps) {
	return sky_coefficientsAttenuation * sky_airmass(position, direction, rayLength, steps);
}
vec3 sky_opticalDepth(vec3 position, vec3 direction, const float steps) {
	return sky_coefficientsAttenuation * sky_airmass(position, direction, steps);
}

vec3 sky_transmittance(vec3 position, vec3 direction, const float steps) {
	return exp2(-sky_opticalDepth(position, direction, steps) * rLOG2);
}

float calculateSunSpot(float VdotL) {
    const float sunRadius = radians(sunAngularSize);
    const float cosSunRadius = cos(sunRadius);
    const float sunLuminance = 1.0 / ((1.0 - cosSunRadius) * TAU);

	return fstep(cosSunRadius, VdotL) * sunLuminance;
}

float   calculateSunGlow(vec3 npos, vec3 lightVector) {
	float curve = 4.0f;

	vec3 halfVector2 = normalize(-lightVector + npos);
	float factor = 1.0f - dot(halfVector2, npos);

	return factor;
}

float calculateMoonSpot(float VdotL) {
    const float moonRadius = radians(moonAngularSize);
    const float cosMoonRadius = cos(moonRadius);
	const float moonLuminance = 1.0 / ((1.0 - cosMoonRadius) * TAU);

	return fstep(cosMoonRadius, VdotL) * moonLuminance;
}

/*
float calculateStars(vec3 worldVector, vec3 moonVector){
	const int res = 256;

	worldVector = getRotMat(moonVector, vec3(0.0, 1.0, 0.0)) * worldVector;

	vec3 p = worldVector * res;
	vec3 flr = floor(p);
	vec3 fr = (p - flr) - 0.5;

	float intensity = hash13(flr);
		  intensity = fstep(intensity, 0.0025);
	float stars = smoothstep(0.5, 0.0, length(fr)) * intensity;

	return stars * 5.0;
}
*/

float calculateStars(vec3 worldVector, vec3 moonVector){
	const int res = 256;

	vec3 position = getRotMat(moonVector, upVector) * worldVector;

	vec3 p = position * res;
	vec3 flr = floor(p);
	vec3 fr = (p - flr) - 0.5;

	float hash = hash13(flr);

	float size = hash;
	size = size * size * (3.0 - 2.0 * size) * 1.2777778e2;

	float intensity = size;
		  intensity = smoothstep(0.0073111113, 0.0071111112, intensity);
	float stars = smoothstep(0.5, 0.0, length(fr)) * intensity;

	return stars * 5.0;
}

vec3 calculateAtmosphere(vec3 background, vec3 viewVector, vec3 upVector, vec3 sunVector, vec3 moonVector, out vec2 pid, out vec3 transmittance, const int iSteps) {
	const int jSteps = 4;

	vec3 viewPosition = (sky_planetRadius + eyeAltitude) * upVector;

	vec2 aid = rsi(viewPosition, viewVector, sky_atmosphereRadius);
	if (aid.y < 0.0) {transmittance = vec3(1.0); return background;}
	
	pid = rsi(viewPosition, viewVector, sky_planetRadius * 0.998);
	bool planetIntersected = pid.y >= 0.0;

	vec2 sd = vec2((planetIntersected && pid.x < 0.0) ? pid.y : max(aid.x, 0.0), (planetIntersected && pid.x > 0.0) ? pid.x : aid.y);

	float stepSize  = (sd.y - sd.x) * (1.0 / iSteps);
	vec3  increment = viewVector * stepSize;
	vec3  position  = viewVector * sd.x + (increment * 0.3 + viewPosition);

	vec2 phaseSun  = sky_phase(dot(viewVector, sunVector ), sky_mieg);
	vec2 phaseMoon = sky_phase(dot(viewVector, moonVector), sky_mieg);

	vec3 scatteringSun     = vec3(0.0);
	vec3 scatteringMoon    = vec3(0.0);
	vec3 scatteringAmbient = vec3(0.0);
	transmittance = vec3(1.0);

	for (int i = 0; i < iSteps; ++i, position += increment) {
		vec3 density          = sky_density(length(position));
		if (density.y > 1e35) break;
		vec3 stepAirmass      = density * stepSize;
		vec3 stepOpticalDepth = sky_coefficientsAttenuation * stepAirmass;

		vec3 stepTransmittance       = exp2(-stepOpticalDepth * rLOG2);
		vec3 stepTransmittedFraction = clamp01((stepTransmittance - 1.0) / -stepOpticalDepth);
		vec3 stepScatteringVisible   = transmittance * stepTransmittedFraction;

		scatteringSun  += sky_coefficientsScattering * (stepAirmass.xy * phaseSun ) * stepScatteringVisible * sky_transmittance(position, sunVector,  jSteps);
		scatteringMoon += sky_coefficientsScattering * (stepAirmass.xy * phaseMoon) * stepScatteringVisible * sky_transmittance(position, moonVector, jSteps);

		// Nice way to fake multiple scattering.
		scatteringAmbient += sky_coefficientsScattering * stepAirmass.xy * stepScatteringVisible;

		transmittance *= stepTransmittance;
	}

	vec3 scattering = scatteringSun * sunColorBase + scatteringAmbient * background + scatteringMoon*moonColorBase;

	float rainPhase = max(sky_miePhase(dot(viewVector, sunVector ),0.4),sky_miePhase(dot(viewVector, sunVector ),0.1)*0.3);
	float L = 2000.;
	float rainDensity = 800.*rainStrength;
	vec3 rainCoef = 2e-5*vec3(0.1);
	transmittance *= exp(-(rainCoef)*rainDensity*L);
	return scattering;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////SEUS/////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

const float R = 1.0;
const float R_INNER = 0.985;
const float SCALE_H = 4.0 / (R - R_INNER);
const float SCALE_L = 1.0 / (R - R_INNER);
const float rayleighAmount = 0.675;

#define fmafma(a,b,c) ((a)*(b)+(c))

//x is distance to outer surface, y is distance to inner surface
vec2 RaySphereIntersection( vec3 p, vec3 dir, float r ) 
{
	float b = dot( p, dir );
	float c = fmafma(r, -r, dot( p, p ));
	
	float d = fmafma(b, b, -c);
	if ( d < 0.0 ) 
	{
		return vec2( 10000.0, -10000.0 );
	}

	d = sqrt( d );
	
	return vec2( -b - d, -b + d );
}

// Mie
// g : ( -0.75, -0.999 )
//      3 * ( 1 - g^2 )               1 + c^2
// F = ----------------- * -------------------------------
//      2 * ( 2 + g^2 )     ( 1 + g^2 - 2 * g * c )^(3/2)
float phase_mie( float g, float c, float cc ) {
	float gg = g * g;
	
	float a = ( 1.0 - gg ) * ( 1.0 + cc );

	float b = fmafma((g * c), -2.0f, (1.0 + gg));
	b *= sqrt( b );
	b *= 2.0 + gg;	
	
	return 1.5 * a / b;
}

// Reyleigh
// g : 0
// F = 3/4 * ( 1 + c^2 )
float phase_reyleigh( float cc ) 
{
	return 0.75 * ( 1.0 + cc );
}

float density( vec3 p )
{
	return exp( -( length( p ) - R_INNER ) * SCALE_H ) * 2.0;
}

float optic( vec3 p, vec3 q ) 
{
	const int numOutscatter = 1;

	vec3 step = ( q - p ) / float(numOutscatter);
	step *= 0.3;
	vec3 v = fmafma(step, vec3(0.5f), p);
	
	float sum = 0.0;
	for ( int i = 0; i < numOutscatter; i++ ) 
	{
		sum += density( v );
		v += step;
	}
	sum *= length( step ) * SCALE_L;


	return sum;
}

vec3 in_scatter(vec3 o, vec3 rayDir, vec2 originToSkyIntersection, vec3 lightVector, const float mieAmount)
{
	const float numInscatter = 3;
	const float K_R = 0.186 * rayleighAmount;
	const float K_M_NUM = 0.020;


	const float K_M = K_M_NUM * mieAmount;
	const float E = 30;
	const vec3 C_R = vec3(0.2, 0.45, 1.0);	//Rayleigh scattering coefficients
	const float G_M = -0.75;

	//float boosty = Boosty(lightVector.y);

	//float rayStepSize = (originToSkyIntersection.y * (1.0 + boosty * 0.0)) / float(numInscatter);
	float rayStepSize = originToSkyIntersection.y / float(numInscatter);
	vec3 step = rayDir * rayStepSize;
	step *= 2.0;
	vec3 p = o;

	//vec3 skyRayPos = p + rayDir * (rayStepSize * (0.55 + boosty * 0.0));
	vec3 skyRayPos = fmafma((rayDir * vec3(rayStepSize)), vec3(0.5), p);



	vec3 sum = vec3( 0.0 );
	for ( int i = 0; i < numInscatter; i++ ) 
	{
		vec2 atmosphereIntersection = RaySphereIntersection(skyRayPos, lightVector, R);
		vec3 outerAtmospherePos = fmafma(lightVector, vec3(atmosphereIntersection.y), skyRayPos);
		
		float n = (optic(p, skyRayPos) + optic(skyRayPos, outerAtmospherePos)) * (PI * 4.0);
		
		sum += density(skyRayPos) * exp(-n * (fmafma(vec3(K_R), C_R, vec3(K_M))));

		skyRayPos += step;
	}
	sum *= rayStepSize * SCALE_L;
	
	float c  = dot(rayDir, -lightVector);
	float cc = c * c;
	
	return sum * (K_R * C_R * phase_reyleigh(cc) + K_M * phase_mie(G_M, c, cc)) * E;
}

vec3 in_scatter2(vec3 rayOrigin, vec3 rayDir, vec2 originToSkyIntersection, vec3 lightVector) 
{
	const float numInscatter = 1;

	const float K_R = 0.186;
	const float K_M = 0.00;
	const float E = 30;
	const vec3 C_R = vec3(0.05, 0.7, 1.0);	//Rayleigh scattering coefficients
	const float G_M = -0.75;

	float len = (originToSkyIntersection.y) / float(numInscatter);
	vec3 step = rayDir * len;
	step *= 2.0;
	vec3 p = rayOrigin;

	//float boosty = Boosty(lightVector.y);

	//vec3 skyRayPos = p + rayDir * (len * (0.5 + boosty * 0.0));
	vec3 skyRayPos = fmafma((rayDir * vec3(len)), vec3(0.5), p);



	vec3 sum = vec3( 0.0 );
	for ( int i = 0; i < numInscatter; i++ ) 
	{
		vec2 atmosphereIntersection = RaySphereIntersection(skyRayPos, lightVector, R);
		vec3 outerAtmospherePos = fmafma(lightVector, vec3(atmosphereIntersection.y), skyRayPos);
		
		float n = (optic(p, skyRayPos) + optic(skyRayPos, outerAtmospherePos)) * (PI * 4.0);
		
		sum += density(skyRayPos) * exp(-n * (fmafma(vec3(K_R), C_R, vec3(K_M))));

		skyRayPos += step;
	}
	sum *= len * SCALE_L;
	
	float c  = dot(rayDir, -lightVector);
	float cc = c * c;
	
	return sum * (K_R * C_R * phase_reyleigh(cc) + K_M * phase_mie(G_M, c, cc)) * E;
}


vec3 Scattering(vec3 eye, vec3 rayDir, vec2 e, vec3 lightVector, const float mieAmount, vec3 up, vec2 eup)
{
	//vec3 atmosphere = in_scatter(eye, rayDir, e, lightVector, mieAmount, 1.0);

	vec3 atmosphere = in_scatter(eye, rayDir, e, lightVector, mieAmount);

	vec3 secondary = in_scatter2(eye, up, eup, lightVector);

	//vec3 ambient = vec3(0.4, 0.55, 1.0);
	vec3 ambient = vec3(0.10, 0.30, 1.0);


	//float boosty = saturate(lightVector.y) * 0.90 + 0.10;
	//boosty = 1.0 / sin(boosty);

	//atmosphere += dot(secondary, vec3(0.06)) * ambient * boosty;
	//atmosphere += dot(secondary, vec3(0.86)) * ambient;
	//atmosphere += dot(secondary, vec3(0.55)) * ambient;
	atmosphere += dot(secondary, vec3(0.625)) * ambient;
	//atmosphere += ambient * 0.01;

	//atmosphere *= vec3(0.8, 0.89, 1.0);


	atmosphere = pow(atmosphere, vec3(1.2))*1.3;
	//atmosphere = pow(atmosphere, vec3(1.2));

	//if (originalRayDir.y < 0.0)
	//{
		//atmosphere *= curve(saturate(originalRayDir.y + 1.0));
	//}

	return atmosphere;
}


vec3 AtmosphericScattering(vec3 rayDir, vec3 lightVector, const float mieAmount)
{
	//Scatter constants
	vec3 eye = vec3(0.0, mix(R_INNER, 1.0, 0.05), 0.0);

	if (rayDir.y < 0.0)
	{
		//rayDir.y = abs(rayDir.y);
		//rayDir.y *= rayDir.y;
		rayDir.y = 0.0;
	}

	vec3 up = vec3(0.0, 1.0, 0.0);

	vec2 e = RaySphereIntersection(eye, rayDir, R);
	vec2 eup = RaySphereIntersection(eye, up, R);

	return Scattering(eye, rayDir, e, lightVector, mieAmount, up, eup);
}


vec3 AtmosphericScattering(vec3 rayDir, vec3 lightVector, const float mieAmount, float depth)
{
	//Scatter constants
	vec3 eye = vec3(0.0, mix(R_INNER, 1.0, 0.05), 0.0);

	if (rayDir.y < 0.0)
	{
		//rayDir.y = abs(rayDir.y);
		//rayDir.y *= rayDir.y;
		rayDir.y = 0.0;
	}

	vec3 up = vec3(0.0, 1.0, 0.0);

	vec2 e = RaySphereIntersection(eye, rayDir, R);
	vec2 eup = RaySphereIntersection(eye, up, R);

	e.y = depth;
	eup.y = depth;


	return Scattering(eye, rayDir, e, lightVector, mieAmount, up, eup);
}


vec3 AtmosphericScatteringSingle(vec3 rayDir, vec3 lightVector, const float mieAmount)
{
	//Scatter constants
	vec3 eye = vec3(0.0, mix(R_INNER, 1.0, 0.05), 0.0);

	if (rayDir.y < 0.0)
	{
		//rayDir.y = abs(rayDir.y);
		//rayDir.y *= rayDir.y;
		rayDir.y = 0.0;
	}

	vec3 up = vec3(0.0, 1.0, 0.0);

	vec2 e = RaySphereIntersection(eye, rayDir, R);
	vec2 eup = RaySphereIntersection(eye, up, R);


	return Scattering(eye, rayDir, e, lightVector, mieAmount, up, eup);
}
