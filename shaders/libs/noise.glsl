const int 		noiseTextureResolution 		    = 64;

#define fma(a,b,c) a*b+c

float bayer2(vec2 a){
    a = floor(a);
    return fract( dot(a, vec2(.5f, a.y * .75f)) );
}

// #define bayer4(a)   (bayer2( .5f*(a))*.25f+bayer2(a))
// #define bayer8(a)   (bayer4( .5f*(a))*.25f+bayer2(a))
// #define bayer16(a)  (bayer8( .5f*(a))*.25f+bayer2(a))
// #define bayer32(a)  (bayer16(.5f*(a))*.25f+bayer2(a))
// #define bayer64(a)  (bayer32(.5f*(a))*.25f+bayer2(a))

// float bayer_4x4(in vec2 pos, in vec2 view) {
// 	return bayer4(pos * view);
// }

// float bayer_8x8(in vec2 pos, in vec2 view) {
// 	return bayer8(pos * view);
// }

// float bayer_16x16(in vec2 pos, in vec2 view) {
// 	return bayer16(pos * view);
// }

// float bayer_32x32(in vec2 pos, in vec2 view) {
// 	return bayer32(pos * view);
// }

// float bayer_64x64(in vec2 pos, in vec2 view) {
// 	return bayer64(pos * view);
// }

float bayer4(vec2 a)   { return bayer2( 0.5  *a) * 0.25     + bayer2(a); }
float bayer8(vec2 a)   { return bayer4( 0.5  *a) * 0.25     + bayer2(a); }
float bayer16(vec2 a)  { return bayer4( 0.25 *a) * 0.0625   + bayer4(a); }
float bayer32(vec2 a)  { return bayer8( 0.25 *a) * 0.0625   + bayer4(a); }
float bayer64(vec2 a)  { return bayer8( 0.125*a) * 0.015625 + bayer8(a); }
float bayer128(vec2 a) { return bayer16(0.125*a) * 0.015625 + bayer8(a); }

vec3 rand(vec2 coord)
{
	float noiseX = saturate(fract(sin(dot(coord, vec2(12.9898, 78.223))) * 43758.5453));
	float noiseY = saturate(fract(sin(dot(coord, vec2(12.9898, 78.223)*2.0)) * 43758.5453));
	float noiseZ = saturate(fract(sin(dot(coord, vec2(12.9898, 78.223)*3.0)) * 43758.5453));

	return vec3(noiseX, noiseY, noiseZ);
}

float  	CalculateDitherPattern(vec2 coord) {
	const int[4] ditherPattern = int[4] (0, 2, 1, 4);

	vec2 count = vec2(0.0f);
	     count.x = floor(mod(coord.s * viewWidth, 2.0f));
		 count.y = floor(mod(coord.t * viewHeight, 2.0f));

	int dither = ditherPattern[int(count.x) + int(count.y) * 2];

	return float(dither) / 4.0f;
}

float  	CalculateDitherPattern1(vec2 coord) {
	const int[16] ditherPattern = int[16] (0 , 8 , 2 , 10,
									 	   12, 4 , 14, 6 ,
									 	   3 , 11, 1,  9 ,
									 	   15, 7 , 13, 5 );

	vec2 count = vec2(0.0f);
	     count.x = floor(mod(coord.s * viewWidth, 4.0f));
		 count.y = floor(mod(coord.t * viewHeight, 4.0f));

	int dither = ditherPattern[int(count.x) + int(count.y) * 4];

	return float(dither) / 16.0f;
}

vec3 	CalculateNoisePattern1(vec2 offset, float size) {
	vec2 coord = texcoord.st;

	coord *= vec2(viewWidth, viewHeight);
	coord = mod(coord + offset, vec2(size));
	coord /= noiseTextureResolution;

	return texture2D(noisetex, coord).xyz;
}

float HaltonSeq3(int index)
    {
        float r = 0.;
        float f = 1.;
        int i = index;
        while (i > 0)
        {
            f /= 3.0;
            r += f * (i % 3);
            i = int(i / 3.0);
        }
        return r;
    }
float HaltonSeq2(int index)
    {
        float r = 0.;
        float f = 1.;
        int i = index;
        while (i > 0)
        {
            f /= 2.0;
            r += f * (i % 2);
            i = int(i / 2.0);
        }
        return r;
    }

vec2 tempOffsets = vec2(HaltonSeq2(frameCounter%10000),HaltonSeq3(frameCounter%10000));

vec2 calculateBlueNoise(vec2 coord, float size) {
	vec2 noiseCoord = vec2(coord.st * vec2(viewWidth, viewHeight)) / 64.0 / size;
	vec2 noiseCoord2 = noiseCoord + vec2(sin(frameCounter * 0.75), cos(frameCounter * 0.75));

	noiseCoord = (floor(noiseCoord * 64.0 * size) + 0.5) / 64.0 / size;
	noiseCoord2 = (floor(noiseCoord2 * 64.0 * size) + 0.5) / 64.0 / size;

	return vec2(texture2DLod(noisetex, noiseCoord2.st, 0).b, texture2DLod(noisetex, noiseCoord.st, 2).b);
}

#define HASHSCALE1 443.8975
#define HASHSCALE3 vec3(443.897, 441.423, 437.195)
#define HASHSCALE4 vec3(443.897, 441.423, 437.195, 444.129)

//  1 out, 3 in...
float hash13(vec3 p3)
{
	p3  = fract(p3 * HASHSCALE1);
    p3 += dot(p3, p3.yzx + 19.19);
    return fract((p3.x + p3.y) * p3.z);
}

vec3 hash33(vec3 p){
    p = fract(p * HASHSCALE3);
    p += dot(p.zxy, p.yxz + 19.19);
    return fract(vec3(p.x * p.y, p.z * p.x, p.y * p.z));
}

////////////////////Worley Noise//////////////////////
vec2 Hash22(vec2 p){
	p=vec2(dot(p,vec2(127.1,311.7)),
			dot(p,vec2(269.5,183.3)));
	return -1.0+2.0*fract(sin(p)*43758.5453123);
}

float WNoise_2V( in vec2 x )
{
    vec2 p = floor( x );
    vec2  f = fract( x );

    float res = 8.0;
    for( int j=-1; j<=1; j++ )
    for( int i=-1; i<=1; i++ )
    {
        vec2 b = vec2( i, j );
        vec2  r = vec2( b ) - f + Hash22( p + b )*0.33;
        float d = dot( r, r );

        res = min( res, d );
    }
    return sqrt( res );
}

float Get3DWorleyNoise( in vec3 x )
{
    vec3 p = floor( x );
    vec3 f = fract( x );

    float res = 8.0;
    for( int j=-1; j<=1; j++ )
    for( int i=-1; i<=1; i++ )
	for( int k=-1; k<=1; k++ )
    {
        vec3 b = vec3( i, j, k );
        vec3  r = vec3( b ) - f + hash33( p + b )*0.33;
        float d = dot( r, r );

        res = min( res, d );
    }
    return sqrt( res );
}

float WFBM_2V(vec2 p) 
{
	float f=0.0;
	f+=0.50000*WNoise_2V(p*1.0);
	f+=0.25000*WNoise_2V(p*2.03);
	f+=0.12500*WNoise_2V(p*4.01);
	f+=0.06250*WNoise_2V(p*8.05);
	f+=0.03125*WNoise_2V(p*16.02);
	return f/0.984375;
}
float Get3DWorleyFBM(vec3 p) 
{
	float f=0.0;
	f+=0.50000*Get3DWorleyNoise(p*1.0);
	f+=0.25000*Get3DWorleyNoise(p*2.03);
	f+=0.12500*Get3DWorleyNoise(p*4.01);
	f+=0.06250*Get3DWorleyNoise(p*8.05);
	f+=0.03125*Get3DWorleyNoise(p*16.02);
	return f/0.984375;
}
////////////////////Parlin Noise//////////////////////
vec2 Hash22P(vec2 p){
	p=vec2(dot(p,vec2(127.1,311.7)),
			dot(p,vec2(269.5,183.3)));
	return -1.0+2.0*fract(sin(p)*43758.5453123);
}

float PNoise_2V( in vec2 p )
{
    vec2 i = floor( p );
    vec2 f = fract( p );
	vec2 u = f*f*(3.0-2.0*f);
    return mix( mix( dot( Hash22P( i + vec2(0.0,0.0) ), f - vec2(0.0,0.0) ), 
                     dot( Hash22P( i + vec2(1.0,0.0) ), f - vec2(1.0,0.0) ), u.x),
                mix( dot( Hash22P( i + vec2(0.0,1.0) ), f - vec2(0.0,1.0) ), 
                     dot( Hash22P( i + vec2(1.0,1.0) ), f - vec2(1.0,1.0) ), u.x), u.y);
}

float FBM_2V(vec2 p){
	float f=0.0;
	f+=0.50000*PNoise_2V(p*1.0);
	f+=0.25000*PNoise_2V(p*2.03);
	f+=0.12500*PNoise_2V(p*4.01);
	f+=0.06250*PNoise_2V(p*8.05);
	f+=0.03125*PNoise_2V(p*16.02);
	return f/0.984375;
}
////////////////////////////////////////////////////////
float Get3DNoise(in vec3 pos)
{
	pos.xyz += 0.5f;

	vec3 p = floor(pos);
	vec3 f = fract(pos);

	f = smoothstep(vec3(0.0), vec3(1.0), f);

	vec2 uv =  (p.xy + p.z * vec2(-17.0f, -17.0f)) + f.xy;

	vec2 coord =  (uv + 0.5f) / 64.0;
	vec2 noiseSample = texture2D(noisetex, coord).xy;
	float xy1 = noiseSample.x;
	float xy2 = noiseSample.y;
	return mix(xy1, xy2, f.z);
}

float Get3DNoiseFBM(vec3 p){
	float f=0.0;
	f+=0.50000*Get3DNoise(p*1.0);
	f+=0.25000*Get3DNoise(p*2.03);
	f+=0.12500*Get3DNoise(p*4.01);
	f+=0.06250*Get3DNoise(p*8.05);
	f+=0.03125*Get3DNoise(p*16.02);
	return f/0.984375;
}

float Get3DNoiseL(in vec3 pos)
{
	pos.xyz += 0.5f;

	vec3 p = floor(pos);
	vec3 f = fract(pos);


	vec2 uv =  (p.xy + p.z * vec2(-17.0f, -17.0f)) + f.xy;

	vec2 coord =  (uv + 0.5f) / 64.0;
	vec2 noiseSample = texture2D(noisetex, coord).xy;
	float xy1 = noiseSample.x;
	float xy2 = noiseSample.y;
	return mix(xy1, xy2, f.z);
}

float Get2DNoise(in vec3 pos)
{
	pos.xy = pos.xz;
	pos.xy += 0.5f;

	vec2 p = floor(pos.xy);
	vec2 f = fract(pos.xy);

	f = smoothstep(vec2(0.0), vec2(1.0), f);

	vec2 uv =  p.xy + f.xy;

	// uv -= 0.5f;
	// uv2 -= 0.5f;

	vec2 coord =  (uv  + 0.5f) / 64.0;
	float xy1 = texture2D(noisetex, coord).x;
	return xy1;
}

float Get3DNoiseRidged(in vec3 pos)
{
	pos.xyz += 0.5f;

	vec3 p = floor(pos);
	vec3 f = fract(pos);

	vec2 uv =  (p.xy + p.z * vec2(-17.0f, -17.0f)) + f.xy;

	vec2 coord =  (uv + 0.5f) / 64.0;
	vec2 noiseSample = texture2D(noisetex, coord).xy;
	float xy1 = noiseSample.x;
	float xy2 = noiseSample.y;
	return abs(mix(xy1, xy2, f.z) * 2.0 - 1.0);
}


float Get3DNoiseBillows(in vec3 pos)
{
	pos.xyz += 0.5f;

	vec3 p = floor(pos);
	vec3 f = fract(pos);

	vec2 uv =  (p.xy + p.z * vec2(-17.0f, -17.0f)) + f.xy;

	vec2 coord =  (uv + 0.5f) / 64.0;
	vec2 noiseSample = texture2D(noisetex, coord).xy;
	float xy1 = noiseSample.x;
	float xy2 = noiseSample.y;
	return 1.0 - abs(mix(xy1, xy2, f.z) * 2.0 - 1.0);
}