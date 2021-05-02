#version 120
#extension GL_EXT_gpu_shader4 : enable

#define FXAA

varying vec2 texcoord;
flat varying vec3 zMults;
uniform sampler2D depthtex0;
uniform sampler2D colortex3;
uniform sampler2D colortex2;
uniform sampler2D colortex1;
uniform sampler2D colortex0;


uniform float far;
uniform float near;
uniform int isEyeInWater;

uniform vec2 texelSize;


#ifdef FXAA
    #define FXAA_SPAN_MAX 4.0 // [1.0 2.0 3.0 4.0 8.0 16.0 32.0 64.0 128.0 1024.0]
    #define FXAA_SPAN_MIN 128.0 // [1.0 2.0 3.0 4.0 8.0 16.0 32.0 64.0 128.0 1024.0]
    #define FXAA_REDUCE_SHIFT 8.0 // [1.0 2.0 3.0 4.0 8.0 16.0 32.0 64.0 128.0 1024.0]

    #define FXAA_REDUCE_MUL   (1.0/FXAA_SPAN_MAX)
    #define FXAA_REDUCE_MIN   (1.0/FXAA_SPAN_MIN)
    #define FXAA_SUBPIX_SHIFT (1.0/FXAA_REDUCE_SHIFT)

    vec3 FxaaPixelShader( vec4 uv, sampler2D tex, vec2 rcpFrame) {
        
        vec3 rgbNW = texture2DLod(tex, uv.zw, 0.0).xyz;
        vec3 rgbNE = texture2DLod(tex, uv.zw + vec2(1,0)*rcpFrame.xy, 0.0).xyz;
        vec3 rgbSW = texture2DLod(tex, uv.zw + vec2(0,1)*rcpFrame.xy, 0.0).xyz;
        vec3 rgbSE = texture2DLod(tex, uv.zw + vec2(1,1)*rcpFrame.xy, 0.0).xyz;
        vec3 rgbM  = texture2DLod(tex, uv.xy, 0.0).xyz;

        vec3 luma = vec3(0.299, 0.587, 0.114);
        float lumaNW = dot(rgbNW, luma);
        float lumaNE = dot(rgbNE, luma);
        float lumaSW = dot(rgbSW, luma);
        float lumaSE = dot(rgbSE, luma);
        float lumaM  = dot(rgbM,  luma);

        float lumaMin = min(lumaM, min(min(lumaNW, lumaNE), min(lumaSW, lumaSE)));
        float lumaMax = max(lumaM, max(max(lumaNW, lumaNE), max(lumaSW, lumaSE)));

        vec2 dir;
        dir.x = -((lumaNW + lumaNE) - (lumaSW + lumaSE));
        dir.y =  ((lumaNW + lumaSW) - (lumaNE + lumaSE));

        float dirReduce = max(
            (lumaNW + lumaNE + lumaSW + lumaSE) * (0.25 * FXAA_REDUCE_MUL),
            FXAA_REDUCE_MIN);
        float rcpDirMin = 1.0/(min(abs(dir.x), abs(dir.y)) + dirReduce);
        
        dir = min(vec2( FXAA_SPAN_MAX,  FXAA_SPAN_MAX),
              max(vec2(-FXAA_SPAN_MAX, -FXAA_SPAN_MAX),
              dir * rcpDirMin)) * rcpFrame.xy;

        vec3 rgbA = (1.0/2.0) * (
            texture2DLod(tex, uv.xy + dir * (1.0/3.0 - 0.5), 0.0).xyz +
            texture2DLod(tex, uv.xy + dir * (2.0/3.0 - 0.5), 0.0).xyz);
        vec3 rgbB = rgbA * (1.0/2.0) + (1.0/4.0) * (
            texture2DLod(tex, uv.xy + dir * (0.0/3.0 - 0.5), 0.0).xyz +
            texture2DLod(tex, uv.xy + dir * (3.0/3.0 - 0.5), 0.0).xyz);
        
        float lumaB = dot(rgbB, luma);

        if((lumaB < lumaMin) || (lumaB > lumaMax)) return rgbA;
        
        return rgbB; 
    }
#endif

void main() {
/* DRAWBUFFERS:3 */
    vec2 uv2 = texcoord.st;

    vec3 color = texture2D(colortex3, texcoord.st).rgb;

    #ifdef FXAA
      vec4 uv = vec4( uv2, uv2 - (texelSize * (0.5 + FXAA_SUBPIX_SHIFT)));
      color = FxaaPixelShader( uv, colortex3, texelSize );
    #endif

    //vec4 shadow = texture2D(colortex2,texcoord.st);
    //color = color*(1.0-shadow.a)*shadow.rgb;

    gl_FragData[0].rgb = clamp(color,0.0,65000.);
}
