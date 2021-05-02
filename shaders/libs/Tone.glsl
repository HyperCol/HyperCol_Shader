#if !(defined _INCLUDE_TONE)
#define _INCLUDE_TONE

//#define VIGNETTE
#ifdef VIGNETTE
uniform vec3 vignetteColor;

vec3 vignette(vec3 color) {
    float dist = distance(texcoord, vec2(0.5f));
    dist = dist * 1.7 - 0.65;
    dist = smoothstep(0.0, 1.0, dist);
    return mix(color, vignetteColor, dist);//vec3(hurt);//
}
#endif

struct Tone {
	//float exposure;
	float brightness;
	float contrast;
	float saturation;
	float vibrance;
	float hue;
	
	vec3 s;
	vec3 m;
	vec3 h;
	bool p;
	
	vec3 color;
	vec3 blur;
	vec3 hsl;
	
	float useAdjustment;
	float blurIndex;
	//float useToneMap;
};

//uniform float nightVision;
//uniform float blindness;
uniform float valLive;

float luma(in vec3 color) { return dot(color,vec3(0.2126, 0.7152, 0.0722)); }

float getLightness(vec3 rgbColor) {
    float r = rgbColor.r, g = rgbColor.g, b = rgbColor.b;
    float minval = min(r, min(g, b));
    float maxval = max(r, max(g, b));
    return ( maxval + minval ) / 2.0;
}

vec3 rgbToHsl(vec3 rgbColor) {
    rgbColor = clamp(rgbColor, vec3(0.0), vec3(1.0));
    float h, s, l;
    float r = rgbColor.r, g = rgbColor.g, b = rgbColor.b;
    float minval = min(r, min(g, b));
    float maxval = max(r, max(g, b));

    float delta = maxval - minval;
    l = ( maxval + minval ) / 2.0;  
    s = delta / mix( maxval + minval , ( 2.0 - maxval - minval ), step(0.5, l))
		* step(0.00001, delta);
	
    float deltaR = (((maxval - r) / 6.0) + (delta / 2.0)) / delta;
    float deltaG = (((maxval - g) / 6.0) + (delta / 2.0)) / delta;
    float deltaB = (((maxval - b) / 6.0) + (delta / 2.0)) / delta;

    h = mix(
			mix(
				(( 2.0 / 3.0 ) + deltaG - deltaR) * select(b, maxval),
				(( 1.0 / 3.0 ) + deltaR - deltaB), select(g, maxval)
				),
			(deltaB - deltaG), select(r, maxval)
			);
			
    h += step(0.0, -h) - step(1.0, h);
    h *= step(0.00001, delta);

    return vec3(h, s, l);
}

float hueToRgb(float v1, float v2, float vH) {
	vH += step(0.0, -vH) - step(1.0, vH);

    if ((6.0 * vH) < 1.0)
        return (v1 + (v2 - v1) * 6.0 * vH);
    if ((2.0 * vH) < 1.0)
        return v2;
    if ((3.0 * vH) < 2.0)
        return (v1 + ( v2 - v1 ) * ( ( 2.0 / 3.0 ) - vH ) * 6.0);
    return v1;
}

vec3 hslToRgb(vec3 hslColor) {
    hslColor = clamp(hslColor, vec3(0.0), vec3(1.0));
    float r, g, b;
    float h = hslColor.r, s = hslColor.g, l = hslColor.b;
    
    float v1, v2;
    v2 = mix(l * (1.0 + s), (l + s) - (s * l), step(0.5,  l));
    v1 = 2.0 * l - v2;
     
    r = hueToRgb(v1, v2, h + (1.0 / 3.0));
    g = hueToRgb(v1, v2, h);
    b = hueToRgb(v1, v2, h - (1.0 / 3.0));
    
    return mix(vec3(r, g, b), vec3(l), step(s, 0.0));
}

/*vec3 rgbToHsl0(vec3 rgbColor) {
    rgbColor = clamp(rgbColor, vec3(0.0), vec3(1.0));
    float h, s, l;
    float r = rgbColor.r, g = rgbColor.g, b = rgbColor.b;
    float minval = min(r, min(g, b));
    float maxval = max(r, max(g, b));
    float delta = maxval - minval;
    l = ( maxval + minval ) / 2.0;  
    if (delta == 0.0) 
    {
        h = 0.0;
        s = 0.0;
    }
    else
    {
        if ( l < 0.5 )
            s = delta / ( maxval + minval );
        else 
            s = delta / ( 2.0 - maxval - minval );
             
        float deltaR = (((maxval - r) / 6.0) + (delta / 2.0)) / delta;
        float deltaG = (((maxval - g) / 6.0) + (delta / 2.0)) / delta;
        float deltaB = (((maxval - b) / 6.0) + (delta / 2.0)) / delta;
         
        if(r == maxval)
            h = deltaB - deltaG;
        else if(g == maxval)
            h = ( 1.0 / 3.0 ) + deltaR - deltaB;
        else if(b == maxval)
            h = ( 2.0 / 3.0 ) + deltaG - deltaR;
             
        if ( h < 0.0 )
            h += 1.0;
        if ( h > 1.0 )
            h -= 1.0;
    }
    return vec3(h, s, l);
}

float hueToRgb0(float v1, float v2, float vH) {
    if (vH < 0.0)
        vH += 1.0;
    if (vH > 1.0)
        vH -= 1.0;
    if ((6.0 * vH) < 1.0)
        return (v1 + (v2 - v1) * 6.0 * vH);
    if ((2.0 * vH) < 1.0)
        return v2;
    if ((3.0 * vH) < 2.0)
        return (v1 + ( v2 - v1 ) * ( ( 2.0 / 3.0 ) - vH ) * 6.0);
    return v1;
}
 
vec3 hslToRgb0(vec3 hslColor) {
    hslColor = clamp(hslColor, vec3(0.0), vec3(1.0));
    float r, g, b;
    float h = hslColor.r, s = hslColor.g, l = hslColor.b;
    if (s == 0.0)
    {
        r = l;
        g = l;
        b = l;
    }
    else
    {
        float v1, v2;
        if (l < 0.5)
            v2 = l * (1.0 + s);
        else
            v2 = (l + s) - (s * l);
     
        v1 = 2.0 * l - v2;
     
        r = hueToRgb0(v1, v2, h + (1.0 / 3.0));
        g = hueToRgb0(v1, v2, h);
        b = hueToRgb0(v1, v2, h - (1.0 / 3.0));
    }
    return vec3(r, g, b);
}*/
#define HUE_ADJUSTMENT

#define TONE 0 //[0 1 2 3 4 5 6 7]

#define BRIGHTNESS 	1.0 	//[0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5]
#define CONTRAST 	1.0   	//[0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5]
#define SATURATION 	1.0 	//[0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.2 1.5 2.0 2.5 3.0]
#define VIBRANCE 	1.0 	//[0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.2 1.5 2.0 2.5 3.0]
#define HUE 		0.0		//[-0.5 -0.45 -0.4 -0.35 -0.3 -0.25 -0.2 -0.15 -0.13 -0.1 -0.09 -0.08 -0.07 -0.06 -0.05 -0.04 -0.03 -0.02 -0.01 0.0 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1 0.13 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5]

#define COLOR_BALANCE_S_R 0.0 //[-1.0 -0.9 -0.8 -0.7 -0.6 -0.5 -0.4 -0.3 -0.2 -0.1 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0]
#define COLOR_BALANCE_S_G 0.0 //[-1.0 -0.9 -0.8 -0.7 -0.6 -0.5 -0.4 -0.3 -0.2 -0.1 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0]
#define COLOR_BALANCE_S_B 0.0 //[-1.0 -0.9 -0.8 -0.7 -0.6 -0.5 -0.4 -0.3 -0.2 -0.1 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0]
#define COLOR_BALANCE_M_R 0.0 //[-1.0 -0.9 -0.8 -0.7 -0.6 -0.5 -0.4 -0.3 -0.2 -0.1 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0]
#define COLOR_BALANCE_M_G 0.0 //[-1.0 -0.9 -0.8 -0.7 -0.6 -0.5 -0.4 -0.3 -0.2 -0.1 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0]
#define COLOR_BALANCE_M_B 0.0 //[-1.0 -0.9 -0.8 -0.7 -0.6 -0.5 -0.4 -0.3 -0.2 -0.1 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0]
#define COLOR_BALANCE_H_R 0.0 //[-1.0 -0.9 -0.8 -0.7 -0.6 -0.5 -0.4 -0.3 -0.2 -0.1 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0]
#define COLOR_BALANCE_H_G 0.0 //[-1.0 -0.9 -0.8 -0.7 -0.6 -0.5 -0.4 -0.3 -0.2 -0.1 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0]
#define COLOR_BALANCE_H_B 0.0 //[-1.0 -0.9 -0.8 -0.7 -0.6 -0.5 -0.4 -0.3 -0.2 -0.1 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0]
//#define KEEP_BROGHTNESS

void init_tone(out Tone t, vec2 texcoord) {
	//t.exposure = get_exposure();
	
	t.color = texture2D(colortex0, texcoord).rgb;
	t.blur = texture2D(gaux4, texcoord).rgb;// * (1.0 + t.exposure);
	//t.hsl = rgbToHsl(col);
	
	t.useAdjustment = 1.0;
	t.blurIndex = 0.0;
	//t.useToneMap = 1.0;
	
	#if TONE == 1
	t.brightness = 1.05;
	t.contrast = 0.8;
	t.saturation = 1.5;
	t.vibrance = 1.2;
	t.hue = 0.0;
	
	t.s = vec3(0.15, 0.1, 0.2);
	t.m = vec3(0.1, 0.0, 0.0);
	t.h = vec3(0.2, 0.1, 0.0);
	t.p = true;
	#elif TONE == 2
	t.brightness = 1.0;
	t.contrast = 0.9;
	t.saturation = 1.0;
	t.vibrance = 0.8;
	t.hue = 0.0;
	
	t.s = vec3(0.1, 0.0, 0.0);
	t.m = vec3(0.0);
	t.h = vec3(0.0, -0.1, -0.1);
	t.p = true;
	#elif TONE == 3
	t.brightness = 0.95;
	t.contrast = 1.0;
	t.saturation = 1.0;
	t.vibrance = 1.1;
	t.hue = 0.0;
	
	t.s = vec3(0.4, 0.2, 0.9);
	t.m = vec3(0.4, 0.4, 0.8);
	t.h = vec3(-0.1, 0.1, 0.6);
	t.p = true;
	#elif TONE == 4
	t.brightness = 1.0;
	t.contrast = 0.9;
	t.saturation = 1.1;
	t.vibrance = 1.0;
	t.hue = 0.0;
	
	t.s = vec3(-0.1, 0.0, 0.1);
	t.m = vec3(0.0, 0.05, 0.1);
	t.h = vec3(-0.1, 0.1, 0.3);
	t.p = true;
	#elif TONE == 5
	t.brightness = 0.9;
	t.contrast = 0.8;
	t.saturation = 0.4;
	t.vibrance = 0.8;
	t.hue = -0.01;
	
	t.s = vec3(0.0);
	t.m = vec3(0.0);
	t.h = vec3(-0.08);
	t.p = false;
	#elif TONE == 6
	t.brightness = 1.0;
	t.contrast = 1.2;
	t.saturation = 2.0;
	t.vibrance = 2.0;
	t.hue = 0.0;
	
	t.s = vec3(0.0);
	t.m = vec3(0.0);
	t.h = vec3(0.0);
	t.p = false;
	#elif TONE == 7
	t.brightness = BRIGHTNESS;
	t.contrast = CONTRAST;
	t.saturation = SATURATION;
	t.vibrance = 1.0 / VIBRANCE;
	t.hue = HUE;
	
	t.s = vec3(COLOR_BALANCE_S_R, COLOR_BALANCE_S_G, COLOR_BALANCE_S_B);
	t.m = vec3(COLOR_BALANCE_M_R, COLOR_BALANCE_M_G, COLOR_BALANCE_M_B);
	t.h = vec3(COLOR_BALANCE_H_R, COLOR_BALANCE_H_G, COLOR_BALANCE_H_B);
	#ifdef KEEP_BROGHTNESS
	t.p = true;
	#else
	t.p = false;
	#endif
	#else
	t.brightness = 1.0;
	t.contrast = 1.0;
	t.saturation = 1.0;
	t.vibrance = 1.0;
	t.hue = 0.0;
	
	t.s = vec3(0.0);
	t.m = vec3(0.0);
	t.h = vec3(0.0);
	t.p = false;
	#endif
}

void contrast(inout vec3 rgbColor, float c) {
    rgbColor = mix(vec3(0.5), rgbColor, c);
}

void saturation(inout vec3 rgbColor, float s) {
	rgbColor = mix(vec3(luma(rgbColor)), rgbColor, s);
}

void vibrance(inout vec3 hslColor, float v) {
    hslColor.g = pow(hslColor.g, v);
    //return hslColor;
}

vec3 colorBalance(vec3 rgbColor, float l, vec3 s, vec3 m, vec3 h, bool p) {
    float r = rgbColor.r, g = rgbColor.g, b = rgbColor.b;
    
    s *= clamp((l - 0.333) / -0.25 + 0.5, 0.0, 1.0) * 0.7;
    m *= clamp((l - 0.333) /  0.25 + 0.5, 0.0, 1.0) *
         clamp((l + 0.333 - 1.0) / -0.25 + 0.5, 0.0, 1.0) * 0.7;
    h *= clamp((l + 0.333 - 1.0) /  0.25 + 0.5, 0.0, 1.0) * 0.7;
    
    vec3 newColor = rgbColor;
    newColor += s;
    newColor += m;
    newColor += h;
    newColor = clamp(newColor, vec3(0.0), vec3(1.0));
    
    if(p) {
        float nl = getLightness(newColor);
        newColor *= l / nl;
    }
    return newColor;
}

vec3 colorBalance(vec3 rgbColor, vec3 s, vec3 m, vec3 h, bool p) {
    float r = rgbColor.r, g = rgbColor.g, b = rgbColor.b;
    float minval = min(r, min(g, b));
    float maxval = max(r, max(g, b));
    float l = ( maxval + minval ) / 2.0;
	return colorBalance(rgbColor, l, s, m, h, p);
}

/*vec3 tonemap(in vec3 color, float adapted_lum) {
	color *= adapted_lum;

	const float a = 2.51f;
	const float b = 0.03f;
	const float c = 2.43f;
	const float d = 0.59f;
	const float e = 0.14f;
	const vec3 f = vec3(13.134f);
	
	color = pow(color, vec3(1.4));
	color *= (2.0 - rainStrength);
	//color = clamp(color, vec3(0.0), vec3(1.0));
	//color = pow(color, vec3(1.07, 1.04, 1.0));
	
	vec3 curr = (color*(a*color+b))/(color*(c*color+d)+e);
	vec3 whiteScale = 1.0f/((f*(a*f+b))/(f*(c*f+d)+e));
	color = curr*whiteScale;

	return pow(color, agamma);
}*/

void Hue_Adjustment(inout Tone t) {
	//blur & dof
	//t.blurIndex = clamp(t.blurIndex, 0.0, 1.0);
	t.color = mix(t.color, t.blur, t.blurIndex);
	
	// This will turn it into gamma space
	#ifdef BLACK_AND_WHITE
	t.saturation = 0.0;
	#elif MC_VERSION >= 11202
	t.saturation *= valLive;
	#endif
	
	//tonemap
	vec3 color = t.color;
	//t.color = mix(t.color, tonemap(color, t.exposure), t.useToneMap);*/
	
	//hue
	#ifdef HUE_ADJUSTMENT
	if (t.useAdjustment > 0) {
		#ifdef TONE_DEBUG
		if (tex.x < 0.5) {
		#endif
			contrast(t.color, t.contrast);
			saturation(t.color, t.saturation);
			if (t.hue != 0 || t.vibrance != 0 || t.s != vec3(0.0) || t.m != vec3(0.0) || t.h != vec3(0.0)) {
				t.hsl = rgbToHsl(t.color);
			
				t.hsl.r += t.hue;
				vibrance(t.hsl, t.vibrance);
				t.hsl.b *= t.brightness;
			
				t.color = hslToRgb(t.hsl);
	
				t.color = colorBalance(t.color, t.hsl.b, t.s, t.m, t.h, t.p);
			}
		#ifdef TONE_DEBUG
		}
		#endif
	}
	#endif
	
	#ifdef VIGNETTE
	t.color = vignette(t.color);
	#endif
	
	// Apply night vision gamma and blindness
	t.color = pow(t.color, vec3(1.0 - nightVision * 0.6 + blindness));	
	
	t.color = mix(color, t.color, t.useAdjustment);
}
#endif