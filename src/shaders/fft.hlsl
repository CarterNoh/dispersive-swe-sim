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

[numthreads(FFT_SIZE, 1, 1)]
void FFTKernel_1D(
    uint3 GI  : SV_GroupThreadID,   // tid = GI.x  ∈ [0, FFT_MAX_N)
    uint3 GID : SV_GroupID          // GID.y = row or column index
) {
    const uint tid = GI.x;

    // -------------------------------------------------------------------------
    // PHASE 1 — Load with bit-reversal permutation into groupshared memory
    // -------------------------------------------------------------------------

    // Global element coordinate: row pass uses (tid, row), column pass uses (column, tid).
    // if (tid < cb_N) {
        const uint rev = bit_reverse(tid, cb_Bits);
        #ifdef IS_ARRAY
            uint3 coord = (cb_IsRow == 1) ? uint3(tid, GID.y, GID.z) : uint3(GID.y, tid, GID.z);
        #else
            uint2 coord = (cb_IsRow == 1) ? uint2(tid, GID.y) : uint2(GID.y, tid);
        #endif
        gs_Data[rev] = fft[coord];
    // }
    GroupMemoryBarrierWithGroupSync();

    // -------------------------------------------------------------------------
    // PHASE 2 — Butterfly passes
    // -------------------------------------------------------------------------

    for (uint passNum = 0; passNum < cb_Bits; ++passNum) {
        const uint span      = 1u << passNum;
        const uint groupSize = span << 1;
        // Declare variables outside the scope so we can write them later
        float2 newEven, newOdd;
        uint evenIdx, oddIdx;
        if (tid < cb_N / 2) {
            const uint k       = tid % span;
            evenIdx            = (tid / span) * groupSize + k;
            oddIdx             = evenIdx + span;
            float2 tw   = twiddle_factor(k, groupSize, cb_Inverse);
            float2 even = gs_Data[evenIdx];
            float2 odd  = complex_mul(tw, gs_Data[oddIdx]);
            newEven = complex_add(even, odd);
            newOdd  = complex_sub(even, odd);
        }
        // All threads sync before we overwrite the old data
        GroupMemoryBarrierWithGroupSync();
        if (tid < cb_N / 2) {
            gs_Data[evenIdx] = newEven;
            gs_Data[oddIdx]  = newOdd;
        }
        // All threads sync before the next pass loop begins
        GroupMemoryBarrierWithGroupSync();
    }

    // -------------------------------------------------------------------------
    // PHASE 3 — Normalize (inverse only) and write back to global memory
    // -------------------------------------------------------------------------

    float scale = cb_Inverse ? (1.0f / (float)cb_N) : 1.0f;

    float2 result;
    result.x = gs_Data[tid].x * scale;
    result.y = gs_Data[tid].y * scale;
    fft[coord] = result;
}




