#version 120
varying vec2 texcoord;

uniform float viewWidth;
uniform float viewHeight;

uniform int frameCounter;

void main() {
	gl_Position = ftransform();
	texcoord = gl_MultiTexCoord0.xy;
}
