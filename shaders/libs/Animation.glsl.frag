#if !(defined _INCLUDE_ANIMATION)
#define _INCLUDE_ANIMATION

/* 
 * Copyright 2019 Ovizro
 *
 * This is the first shader effect made by Ovizro.
 * It includes a quantity of animation.
 * Some of them might seem to be a little crazy.
 * By the way, a LOGO of my team, HyperCol Studio, has been included in it.
 * Wish you can enjoy it.
 */

uniform bool hideGUI;
 
//#define ANIMATION_DEBUG
#define DELAY 5.0 			//[0.0 0.5 1.0 1.5 2.0 2.5 3.0 3.5 4.0 4.5 5.0 6.0 7.0 8.0 9.0 10.0 11.0 12.0 13.0 14.0 15.0 16.0 17.0 18.0 19.0 20.0]
#define FRAME_COUNTERS 50.0	//[20.0 30.0 50.0 70.0 100.0 150.0 200.0]

//AnimationTime
float fFrameCounter = float(frameCounter);
#ifndef ANIMATION_DEBUG
float animationTimeCounter = mix(max(frameTimeCounter - DELAY, 0.0), 100.0, step(FRAME_COUNTERS * (12.0 + DELAY), fFrameCounter));
#else
float animationTimeCounter = fract(frameTimeCounter * 0.1) * 10.0;//5.0;//
#endif

#define CORNER_MARK 0 //[0 1 2 3]

#if ANIMATION == 1
varying vec2 logoPos[8];
varying float Ltheter;
#endif

#if CORNER_MARK > 0
varying mat3 cLogo;
#endif

#ifdef _VERTEX_SHADER_

/*
 *==============================================================================
 *------------------------------------------------------------------------------
 *
 * 								~Vertex stuff~
 *
 *------------------------------------------------------------------------------
 *==============================================================================
 */

#if ANIMATION == 1
vec2 track(in float animationTime, out float theter) {				//Running track of graph
    const float size = 4.854369;
    
    const float edge0 = floor(log(size) / log(1.618) * 1.3333) * PI * 0.25;	//The angle of beginning. Correct rotation angle
    //edge0 = floor(edge0);					//
    
    float frametime = (smoothstep(-6.0, 6.0, animationTimeCounter) - 0.5) * 8.0;
    float frametime1 = (smoothstep(-5.0, 5.0, animationTime) - 0.5) * 8.0;
    
    theter = edge0 * (3.0 - frametime);
    float theter1 = edge0 * (3.0 - frametime1);
    float r = pow(1.618, 2.0 * theter1 / PI);			//Fibonacci helix
    r = max(r, 1.0);
	
    return r * vec2(cos(theter), sin(theter));
}
#endif

#if CORNER_MARK > 0
mat3 cTrack() {						//Running track of corner mark  
	const float size = 6.6666;
	float theter = -frameTimeCounter * 0.6;

	float s = sin(theter);
	float c = cos(theter);
	const float l = 1.0 / 0.4142;
	mat3 r0 = mat3(	size, 	0.0, 	0.0, 
					0.0, 	size, 	0.0, 
					-c, 	-s, 	1.0);
	mat3 r1 = mat3(	l, 		0.0, 	0.0, 
					0.0, 	l, 		0.0, 
					0.0, 	0.0, 	1.0);
	mat3 r2 = mat3(	s,		c,		0.0, 
					-c, 	s,		0.0, 
					0.0, 	0.0, 	1.0);

	return r2 * r1 * r0;
}
#endif

void animationCommons() {
	//HyperCol Logo
	#if ANIMATION == 1
	for (int i = 0; i < 8; ++i) {
		logoPos[i] = track(animationTimeCounter - float(i) * 0.2, Ltheter);
	}
	#endif
	
	//corner mark
	#if CORNER_MARK > 0
	cLogo = cTrack();
	#endif
}
#else

/*
 *==============================================================================
 *------------------------------------------------------------------------------
 *
 * 								~Fragment stuff~
 *
 *------------------------------------------------------------------------------
 *==============================================================================
 */

vec2 fuv_build(in vec2 uv) {                //Establish coordinate system with screen as center
    vec2 fuv = uv * 2.0 - 1.0;
    fuv.x *= aspectRatio;
    return fuv;
}

uniform sampler2D depthtex2;

const vec4 totalTexPos[] = vec4[2] (
	vec4(0.0, 0.5, 0.25, 0.1),
	vec4(0.0, 0.8, 0.1, 0.2)
);

vec4 getTexColor(in vec2 uv, int key, vec2 start, float size) {
	vec4 tPos = totalTexPos[key];
	vec2 fuv = (uv - start) * tPos.z / size;
	fuv.y *= 2.0;
	if (min(fuv.x, fuv.y) >= 0.0 && fuv.x <= tPos.z && fuv.y <= tPos.w) {
		return texture2D(depthtex2, fuv + tPos.xy);
	} else {
		return vec4(0.0);
	}	
}

/*
 *==============================================================================
 *[																				]
 *[		----------------		Simple Animation		----------------		]
 *[																				]
 *==============================================================================
 */

#define ROTATE        0.0     //[-2.0 -1.5 -1.0 -0.5 0.0 0.5 1.0 1.5 2.0]
#define ROTATING_TIME 3.0   //[0.5 1.0 1.5 2.0 2.5 3.0 3.5 4.0 4.5 5.0 5.5 6.0]
#define ROTATING_SCALE 1.0 //[0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.2 1.5 1.7 2.0 2.5 3.0 4.0 5.0 8.0]

#define SHADE_ROTATE 0.0     //[-2.0 -1.5 -1.0 -0.5 0.0 0.5 1.0 1.5 2.0]
#define SHADE_ROTATING_TIME 3.0   //[0.5 1.0 1.5 2.0 2.5 3.0 3.5 4.0 4.5 5.0 5.5 6.0]
#define SHADE_ROTATING_SCALE 1.0 //[0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.2 1.5 1.7 2.0 2.5 3.0 4.0 5.0 8.0]

//#define WHITE_SHADE
#define TRANSLUCENT_SHADE  				//Only for black shade
//#define TRANSLUCENT_SHADE_BLUR

//#define RECTANGULAR_SHADE
#define ANAMORPHIC_EDGE 1.0				//[0.0 0.5 0.55 0.65 0.7 0.75 0.8 0.85 0.9 0.95 1.0]
#define TRANSLUCENT_ANAMORPHIC_EDGE 1.0	//[0.0 0.5 0.55 0.65 0.7 0.75 0.8 0.85 0.9 0.95 1.0]
#define ANAMORPHIC_CONTRACT_SPEED 0.5	//[0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0]
#define TRANSLUCENT_ANAMORPHIC_CONTRACT_SPEED 0.4	//[0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0]

#define LOZENGULAR_SHADE
#define LOZENGULAR_SHADE_PAUSE
#define LOZENGULAR_SHADE_MIDDLE_TIME 1.0 				//[0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0]
#define LOZENGULAR_SHADE_MIDDLE_POSITION 0.6				//[0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0]
#define TRANSLUCENT_LOZENGULAR_SHADE_MIDDLE_POSITION 0.4	//[0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0]

//#define ROUND_SHADE
#define ROUND_SHADE_PAUSE
#define ROUND_SHADE_MIDDLE_TIME 1.0 				//[0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0]
#define ROUND_SHADE_MIDDLE_POSITION 0.6				//[0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0]
#define TRANSLUCENT_ROUND_SHADE_MIDDLE_POSITION 0.5	//[0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0]

mat2 mRotate(float theter) {
	float s = sin(theter);
	float c = cos(theter);
	return mat2(c, -s,
				s,  c);
}

void rotate(inout vec2 uv, in float theter) {
	uv = mRotate(theter) * uv;
}

float lozenge(in vec2 puv, float edge) {
	float e0 = puv.x + puv.y;
	return step(edge, e0);
}

float round(in vec2 puv, float r) {
	float e0 = dot(puv, puv);
	return step(pow2(r), e0);
}

float triangle(in vec2 puv) {					//Build an equilateral triangle
	float e1 = 0.57735 * (1.0 - puv.x) - abs(puv.y);
	float e2 = puv.x + 1.0; 
	return min(smoothstep(0.0, 0.05, e1), smoothstep(0.0, 0.0443, e2));
}

float func3(in float x, float m, float mY) {
	float a = mY / pow3(m);
	float b = -a * m * 3.0;
	float c = -b * m;
	return fma(x, fma(x, a, b), c) * x;
}

void simple_animation(inout Tone t, vec2 fuv) {
	//rotation
	vec2 uv = fuv;
	vec2 rM = vec2(ROTATE * min(animationTimeCounter - ROTATING_TIME, 0.0), SHADE_ROTATE * min(animationTimeCounter - SHADE_ROTATING_TIME, 0.0));
	vec2 l = vec2(smoothstep(- ROTATING_TIME, ROTATING_TIME, animationTimeCounter), smoothstep(-SHADE_ROTATING_TIME, SHADE_ROTATING_TIME, animationTimeCounter)) * 2.0 - 1.0;
	l = mix(vec2(ROTATING_SCALE, SHADE_ROTATING_SCALE), vec2(1.0), l);
	if (SHADE_ROTATE != 0 || SHADE_ROTATING_SCALE != 1) {
		rotate(fuv, rM.y);
		fuv /= l.y;
	}
	if (ROTATE != 0 || ROTATING_SCALE != 1) {
		rotate(uv, rM.x);
		uv /= vec2(aspectRatio, 1.0) * l.x;
		vec3 color = texture2D(composite, uv * 0.5 + 0.5).rgb * Cselect(uv, -1.0, 1.0);
		float t0 = smoothstep(ROTATING_TIME, ROTATING_TIME + 1.0, animationTimeCounter);
		t.color = mix(color, t.color, t0);
		t.blurIndex *= t0;
	}
	
	//shade
	fuv = abs(fuv);
	//rectangle
	#ifdef RECTANGULAR_SHADE
	float e0 = step(fuv.y, min(animationTimeCounter * ANAMORPHIC_CONTRACT_SPEED, ANAMORPHIC_EDGE));
	float e1 = 1.0 - step(min(animationTimeCounter * TRANSLUCENT_ANAMORPHIC_CONTRACT_SPEED, TRANSLUCENT_ANAMORPHIC_EDGE), fuv.y) * 0.5;
	#else
	float e0 = 1.0;
	float e1 = 1.0;
	#endif
	
	//lozenge
	#ifdef LOZENGULAR_SHADE
	#ifdef LOZENGULAR_SHADE_PAUSE
	float l0 = 1.0 - lozenge(fuv, func3(animationTimeCounter, LOZENGULAR_SHADE_MIDDLE_TIME, LOZENGULAR_SHADE_MIDDLE_POSITION));
	float l1 = 1.0 - lozenge(fuv, func3(animationTimeCounter, LOZENGULAR_SHADE_MIDDLE_TIME, TRANSLUCENT_LOZENGULAR_SHADE_MIDDLE_POSITION)) * 0.5;
	#else
	float l0 = 1.0 - lozenge(fuv, animationTimeCounter / LOZENGULAR_SHADE_MIDDLE_TIME * LOZENGULAR_SHADE_MIDDLE_POSITION);
	float l1 = 1.0 - lozenge(fuv, animationTimeCounter / LOZENGULAR_SHADE_MIDDLE_TIME * TRANSLUCENT_LOZENGULAR_SHADE_MIDDLE_POSITION) * 0.5;
	#endif
	#else
	float l0 = 1.0;
	float l1 = 1.0;
	#endif
	
	//Round
	#ifdef ROUND_SHADE
	#ifdef ROUND_SHADE_PAUSE
	float r0 = 1.0 - round(fuv, func3(animationTimeCounter, ROUND_SHADE_MIDDLE_TIME, ROUND_SHADE_MIDDLE_POSITION));
	float r1 = 1.0 - round(fuv, func3(animationTimeCounter, ROUND_SHADE_MIDDLE_TIME, TRANSLUCENT_ROUND_SHADE_MIDDLE_POSITION)) * 0.5;
	#else
	float r0 = 1.0 - round(fuv, animationTimeCounter / ROUND_SHADE_MIDDLE_TIME * ROUND_SHADE__MIDDLE_POSITION);
	float r1 = 1.0 - round(fuv, animationTimeCounter / ROUND_SHADE_MIDDLE_TIME * TRANSLUCENT_ROUND_SHADE_MIDDLE_POSITION) * 0.5;
	#endif
	#else
	float r0 = 1.0;
	float r1 = 1.0;
	#endif
	
	float c0 = min(min(e0, l0), r0);
	float c1 = min(min(e1, l1), r1);
	#ifdef TRANSLUCENT_SHADE_BLUR
	t.blurIndex = plus(t.blurIndex, (1.0 - c1) * 1.6);
	#endif
	#ifdef TRANSLUCENT_SHADE
	float c = c0 * c1;
	#else
	float c = c0;
	#endif
	#ifndef WHITE_SHADE
	t.color *= c;
	t.blur *= c;
	#else
	//c = step(1.0, c);
	t.color = mix(vec3(1.0), t.color, c0);
	t.blur = mix(vec3(1.0), t.blur, c0);
	t.useAdjustment *= c0;
	#endif
}

/*
 *==============================================================================
 *[																				]
 *[		----------------		HyperCol Logo			----------------		]
 *[																				]
 *==============================================================================
 */

#if ANIMATION == 1
float smooth_dot(in vec2 puv) {			//Just a light dot
	puv.x += 0.3333;
    float dist = dot(puv, puv);
	dist = smoothstep(0.0, 1.6, 1.6 - dist);
    return pow3(dist);
}

void ot(inout float tri, vec2 fuv) {
	const vec4 trianglePos[] = vec4[16] (		//Center position(x,y), 1/size, angle
		vec4(-1.17, 	-0.464, 4.348, 	1.4),
		vec4(-0.3356, 	1.176, 	2.1739,	1.05),
		vec4(-1.73, 	0.3, 	5.376, 	-0.4),
		vec4(1.42578,	-0.042,	3.760, 	1.24),
		vec4(0.711,		-0.92,	4.7619,	-0.3),
		vec4(0.551, 	0.64, 	16.667,	-0.5),
		vec4(-0.117, 	-0.672,	25.0, 	1.1),
		vec4(-1.56, 	-0.17,	24.39,	1.15),
		vec4(-1.26, 	0.728,	26.3,	-0.2),
		vec4(0.86,		-0.7,	12.5,	0.0),
		vec4(1.4933,	1.14,	3.70,	-0.9),
		vec4(1.36889,	-0.362,	14.3,	1.26),
		vec4(-0.89245,	0.1,	50.0,	0.1),
		vec4(1.3511,	-0.814,	37.453,	0.01),
		vec4(-1.45778,	-1.16,	58.82,	-0.2),
		vec4(1.39,		-0.722,	100.0,	1.25)
	);

	float frametime = smoothstep(0.0, 6.0, animationTimeCounter);
	
	for (int i = 0; i < trianglePos.length(); ++i) {
		vec4 pos = trianglePos[i];
		pos.w -= (1.0 - frametime) * pos.z;
		vec2 puv = (fuv - pos.xy / frametime) * pos.z;
		rotate(puv, pos.w);
		
		tri += triangle(puv - vec2(0.3333, 0.0));
	}
}

vec2 logo_puv_build(in vec2 fuv, in vec2 basePoint, in float theter) {
	float size = 4.854369;
	
    vec2 puv = (fuv * size - basePoint) / 0.4142;
    rotate(puv, theter - PI * 0.5);
	return puv;
}

void HyperCol_Logo(inout vec3 color, vec2 fuv0) {
    vec2 fuv = fuv0;
    const vec3 logoColor = vec3(0.0, 0.62, 0.9);

    float logo = 0.0;
	mat2 rot45 = mRotate(PI * 0.25);
	
	fuv.y -= 0.039 * smoothstep(5.5, 6.0, animationTimeCounter);

    for (int i = 0; i < 8; ++i) {
        fuv = rot45 * fuv; 
        vec2 puv = logo_puv_build(fuv, logoPos[i], Ltheter);

        logo += mix(smooth_dot(puv), triangle(puv), smoothstep(1.0, 3.0, animationTimeCounter));
    }
    
	ot(logo, fuv0);
    color = mix(color, logoColor, logo);
	
	vec3 tex = getTexColor(fuv0, 0, vec2(-0.35, -0.396), 0.7).rgb * smoothstep(5.5, 6.5, animationTimeCounter);
	color = mix(color, logoColor, smoothstep(0.0, 0.5, luma(tex)));
}

vec4 HyperCol_Logo_Build(in vec2 uv) {
    vec2 fuv = fuv_build(uv);

    vec3 background = vec3(1.0) * smoothstep(1.5, 3.0, animationTimeCounter);
    HyperCol_Logo(background, fuv);

    return vec4(background, 1.0) * (1.0 - smoothstep(8.0, 10.0, animationTimeCounter));
}
#endif

/*
 *==============================================================================
 *[																				]
 *[		----------------		 Eyes Open				----------------		]
 *[																				]
 *==============================================================================
 */

#if ANIMATION == 2
float E_ellipse(in vec2 puv, float eY) {
	float e0 = pow2(puv.x) / 2.6 + pow2(puv.y) / eY;
	return smoothstep(0.3, 1.0, e0);
}

void eyes_open(inout Tone t, vec2 fuv) {
	float eY = (sin(PI * 0.2 * animationTimeCounter) * 0.714 + sin(PI * 0.6 * animationTimeCounter) * 1.42857 + pow(1.3, animationTimeCounter) - 1.0) * 0.2353;
	
	float e = 1.0 - E_ellipse(fuv, eY) * (1.0 - smoothstep(6.5, 7.5, animationTimeCounter));
	t.color *= e;
	t.blur *= e;
	t.blurIndex = plus(t.blurIndex, (1.0 - eY) * (1.0 - smoothstep(6.5, 7.5, animationTimeCounter)));
}
#endif

/*
 *==============================================================================
 *[																				]
 *[		----------------		Corner Marks			------------------		]
 *[																				]
 *==============================================================================
 */

#if CORNER_MARK > 0

void cornerMark(inout Tone t, in vec2 fuv0) {
    vec2 fuv = fuv0;
	fuv += vec2(aspectRatio - 0.3, -0.7);

	if (hideGUI) {
		//HyperCol Logo
		const vec3 logoColor = vec3(0.0, 0.62, 0.9);
		float logo = 0.0;
		mat2 r45 = mRotate(PI * 0.25);
			
		//float Ctheter;
		//vec2 basePoint0 = cTrack(Ctheter);
	
		for (int i = 0; i < 8; ++i) {
			fuv = r45 * fuv;
			vec3 puv = cLogo * vec3(fuv, 1.0);//cornermark_puv_build(fuv, basePoint0, Ctheter);//
			logo += triangle(puv.st);
		}
		
		logo *= smoothstep(8.0, 10.0, animationTimeCounter);
		
		#if CORNER_MARK == 1
		t.color = mix(t.color, logoColor, logo);
		t.blurIndex *= (1.0 - logo);
		#elif CORNER_MARK == 2
		t.brightness *= (1.0 - logo * 0.1);
		t.blurIndex = plus(t.blurIndex, 0.3 * logo);
		#elif CORNER_MARK == 3
		t.color += vec3(0.1) * logo;
		t.brightness *= (1.0 + logo * 0.1);
		t.blurIndex = plus(t.blurIndex, 0.7 * logo);
		#endif
	}
}
#endif

/*
 *==============================================================================
 *[																				]
 *[		----------------		Main Animation			------------------		]
 *[																				]
 *==============================================================================
 */

void animation(inout Tone t, vec2 uv) {
	vec2 fuv = fuv_build(uv);
	
	//corner mark
	#if CORNER_MARK >= 1
	cornerMark(t, fuv);
	#endif
	
	//HyperCol Logo
	#if ANIMATION == 1
	vec3 background = vec3(1.0) * smoothstep(1.5, 3.0, animationTimeCounter);
    HyperCol_Logo(background, fuv);
    t.color = mix(background, t.color, smoothstep(8.0, 10.0, animationTimeCounter));
	t.useAdjustment *= smoothstep(8.0, 10.0, animationTimeCounter);
	#endif
	
	#if ANIMATION == 2
	eyes_open(t, fuv);
	#endif
	
	#if ANIMATION == 3
	simple_animation(t, fuv);
	#endif
}
#endif
#endif