// Butterfly FFT compute shader for 2-D fluid surface simulation.

#ifndef FFT_SIZE
#define FFT_SIZE 512 
#endif
#define FFT_PI 3.14159265358979323846f

// Constant Buffer
cbuffer FFTConstants : register(b1) {
    int cb_Nx;         // Transform size X (width)
    int cb_Ny;         // Transform size Y (height)
    int cb_BitsX;      // log2(cb_Nx)
    int cb_BitsY;      // log2(cb_Ny)
    int cb_Inverse;    // 0 = forward DFT, 1 = inverse DFT
    int cb_IsRow;      // Row = 1, Col = 0
};

// Texture Registers
#ifdef IS_ARRAY
    RWTexture2DArray<float2> fft : register(u0);
#else
    RWTexture2D<float2> fft : register(u0);
#endif

// Groupshared memory
groupshared float2 gs_Data[FFT_SIZE];

//////////////////// HELPER FUNCTIONS /////////////////////////

float2 complex_mul(float2 a, float2 b) {
    return float2(a.x * b.x - a.y * b.y,
                  a.x * b.y + a.y * b.x);
}

float2 complex_add(float2 a, float2 b) {
    return float2(a.x + b.x, a.y + b.y);
}

float2 complex_sub(float2 a, float2 b) {
    return float2(a.x - b.x, a.y - b.y);
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

#define THREAD_COUNT (FFT_SIZE / 2)

[numthreads(THREAD_COUNT, 1, 1)]
void FFTKernel_1D(
    uint3 GI  : SV_GroupThreadID,   // tid = GI.x  ∈ [0, FFT_MAX_N)
    uint3 GID : SV_GroupID          // GID.y = row or column index
) {
    const uint tid = GI.x;
    const uint N = (cb_IsRow == 1) ? (uint)cb_Nx : (uint)cb_Ny;
    const uint bits = (cb_IsRow == 1) ? (uint)cb_BitsX : (uint)cb_BitsY;

    // -------------------------------------------------------------------------
    // PHASE 1 — Load with bit-reversal permutation into groupshared memory
    // -------------------------------------------------------------------------

    // Global element coordinate: row pass uses (tid, row), column pass uses (column, tid).
    // Each thread processes 2 elements
    const uint idx1 = tid;
    const uint idx2 = tid + THREAD_COUNT;

    #ifdef IS_ARRAY
        const uint3 coord1 = (cb_IsRow == 1) ? uint3(idx1, GID.y, GID.z) : uint3(GID.y, idx1, GID.z);
        const uint3 coord2 = (cb_IsRow == 1) ? uint3(idx2, GID.y, GID.z) : uint3(GID.y, idx2, GID.z);
    #else
        const uint2 coord1 = (cb_IsRow == 1) ? uint2(idx1, GID.y) : uint2(GID.y, idx1);
        const uint2 coord2 = (cb_IsRow == 1) ? uint2(idx2, GID.y) : uint2(GID.y, idx2);
    #endif

    if (idx1 < N) {
        const uint rev1 = bit_reverse(idx1, bits);
        gs_Data[rev1] = fft[coord1];
    }
    if (idx2 < N) {
        const uint rev2 = bit_reverse(idx2, bits);
        gs_Data[rev2] = fft[coord2];
    }
    GroupMemoryBarrierWithGroupSync();

    // -------------------------------------------------------------------------
    // PHASE 2 — Butterfly passes
    // -------------------------------------------------------------------------
    
    for (uint passNum = 0; passNum < bits; ++passNum) {
        const uint span      = 1u << passNum;
        const uint groupSize = span << 1;
        
        // Since THREAD_COUNT is N/2, every thread handles exactly one butterfly.
        const uint k       = tid % span;
        const uint evenIdx = (tid / span) * groupSize + k;
        const uint oddIdx  = evenIdx + span;
        
        float2 tw   = twiddle_factor(k, groupSize, cb_Inverse);
        float2 even = gs_Data[evenIdx];
        float2 odd  = complex_mul(tw, gs_Data[oddIdx]);
        
        float2 newEven = complex_add(even, odd);
        float2 newOdd  = complex_sub(even, odd);
        
        GroupMemoryBarrierWithGroupSync();
        gs_Data[evenIdx] = newEven;
        gs_Data[oddIdx]  = newOdd;
        GroupMemoryBarrierWithGroupSync();
    }

    // -------------------------------------------------------------------------
    // PHASE 3 — Normalize (inverse only) and write back to global memory
    // -------------------------------------------------------------------------
    
    const float scale = cb_Inverse ? (1.0f / (float)N) : 1.0f;

    if (idx1 < N) {
        fft[coord1] = gs_Data[idx1] * scale;
    }
    if (idx2 < N) {
        fft[coord2] = gs_Data[idx2] * scale;
    }
}




