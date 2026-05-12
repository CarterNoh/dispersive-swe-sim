// Butterfly FFT compute shader for 2-D fluid surface simulation.

#ifndef FFT_SIZE
#define FFT_SIZE 512 
#endif
#define FFT_PI 3.14159265358979323846f

// Constant Buffer
cbuffer FFTConstants : register(b1) {
    int cb_N;          // Transform size (gridsize) (must be power of 2)
    int cb_Bits;       // log2(cb_N) 
    int cb_Inverse;    // 0 = forward DFT, 1 = inverse DFT
    int cb_IsRow;     // Row = 1, Col = 0
};

// Texture Registers
#ifdef IS_ARRAY
    Texture2DArray<float2> inputTex : register(t0);
    RWTexture2DArray<float2> outputTex : register(u0);
#else
    Texture2D<float2> inputTex : register(t0);
    RWTexture2D<float2> outputTex : register(u0);
#endif

// Groupshared memory
groupshared float2 gs_Data[FFT_SIZE];

//////////////////// HELPER FUNCTIONS /////////////////////////

float2 complex_mul(float2 a, float2 b) {
    return float2(a.x * b.x - a.y * b.y,
                  a.x * b.y + a.y * b.x);
}

float2 twiddle_factor(uint k, uint N, bool inverse) {
    // Twiddle factor W_N^k = e^{-2πi k/N}  (forward); conjugated for inverse.
    float angle = -2.0f * FFT_PI * (float)k / (float)N;
    if (inverse) angle = -angle;
    return float2(cos(angle), sin(angle));
}

uint bit_reverse(uint x, uint bits) {
    // Bit-reversal of the lower `bits` bits of x.
    uint result = 0;
    for (uint i = 0; i < bits; ++i)
    {
        result = (result << 1) | (x & 1);
        x >>= 1;
    }
    return result;
}

// ---------------------------------------------------------------------------
// Kernel entry point
// ---------------------------------------------------------------------------

[numthreads(FFT_SIZE, 1, 1)]
void RowFFT(uint3 tid : SV_GroupThreadID,
            uint3 gid : SV_GroupID) {
    uint x = tid.x;
    uint y = gid.x;
    uint rev = bit_reverse(x, cb_Bits);
    #ifdef IS_ARRAY
        uint z = gid.y;
        uint3 coord = uint3(x, y, z);
    #else
        uint2 coord = uint2(x, y);
    #endif
    gs_Data[rev] = inputTex[coord];

    GroupMemoryBarrierWithGroupSync();
    for(uint p = 0; p < (uint)cb_Bits; ++p) {
        uint span = 1 << p;
        uint step = span << 1;
        float2 outA, outB;
        uint even, odd;
        if(x < (uint)cb_N / 2) {
            uint k = x % span;
            even = (x / span) * step + k;
            odd  = even + span;
            float2 W = twiddle_factor(k, step, cb_Inverse);
            float2 a = gs_Data[even];
            float2 b = complex_mul(W, gs_Data[odd]);
            outA = a + b;
            outB = a - b;
        }
        GroupMemoryBarrierWithGroupSync();
        if(x < (uint)cb_N / 2) {
            gs_Data[even] = outA;
            gs_Data[odd]  = outB;
        }
        GroupMemoryBarrierWithGroupSync();
    }

    float scale = cb_Inverse ? (1.0f / (float)cb_N) : 1.0f;
    outputTex[coord] = gs_Data[x] * scale;
}


[numthreads(16, 16, 1)]
void Transpose(uint3 id : SV_DispatchThreadID) {
    #ifdef IS_ARRAY
        uint3 dst = uint3(id.y, id.x, id.z);
        outputTex[dst] = inputTex[id];
    #else
        uint2 dst = uint2(id.y, id.x);
        outputTex[dst] = inputTex[id.xy];
    #endif
}