#version 120
#extension GL_EXT_gpu_shader4 : enable

varying vec4 texcoord;

flat varying float tempOffsets;

uniform int frameCounter;

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

void main() {
    tempOffsets = HaltonSeq2(frameCounter%10000);
    gl_Position = ftransform();
    texcoord = gl_MultiTexCoord0;
}