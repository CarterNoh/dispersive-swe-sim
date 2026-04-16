// fft_gpu.hlsl
// Butterfly FFT compute shader for 2-D fluid surface simulation.
//
// CORRESPONDENCE:
//   Every function and index expression here has a line-for-line equivalent
//   in fft_gpu_cpu_emu.cpp.  Search for "HLSL:" comments in that file to
//   find the matching CPU code.
//
// DISPATCH GUIDE (from C++/D3D12):
//   // 1-D row pass (process all rows)
//   cb.N       = width;
//   cb.Inverse = 0;
//   context->Dispatch(1, height, 1);        // one group per row
//
//   // 1-D column pass (process all columns)
//   cb.N       = height;
//   context->Dispatch(1, width, 1);         // one group per column
//
// REQUIREMENTS:
//   - N (cb.N) must be a power of two, <= THREADS_PER_GROUP.
//   - The input/output UAV must be bound as RWBuffer<float2> or
//     RWStructuredBuffer<float2> with width*height elements in row-major order.

// ---------------------------------------------------------------------------
// Compile-time constants — keep in sync with fft_gpu_cpu_emu.cpp
// ---------------------------------------------------------------------------

#define THREADS_PER_GROUP 512   // Must match [numthreads] below.
#define FFT_PI 3.14159265358979323846f

// ---------------------------------------------------------------------------
// Constant buffer
// ---------------------------------------------------------------------------

cbuffer FFTConstants : register(b0)
{
    uint  cb_N;           // Transform size (power of two, <= THREADS_PER_GROUP)
    uint  cb_Inverse;     // 0 = forward DFT, 1 = inverse DFT
    uint  cb_RowOffset;   // Base element index for the current row/column group.
                          // Row pass:    cb_RowOffset = groupID.y * cb_N
                          // Column pass: cb_RowOffset = groupID.y (column index)
                          //              stride        = total_width (set via separate constant)
    uint  cb_Stride;      // Element stride between samples.
                          // Row pass:    1  (contiguous)
                          // Column pass: total_width  (row-major column stride)
};

// ---------------------------------------------------------------------------
// Buffers
// ---------------------------------------------------------------------------

RWStructuredBuffer<float2> g_Data : register(u0);  // In-place read/write

// ---------------------------------------------------------------------------
// Groupshared memory
// ---------------------------------------------------------------------------

groupshared float2 gs_Data[THREADS_PER_GROUP];

// ---------------------------------------------------------------------------
// Helper functions — identical signatures to fft_gpu_cpu_emu.cpp
// ---------------------------------------------------------------------------

float2 complex_mul(float2 a, float2 b)
{
    return float2(a.x * b.x - a.y * b.y,
                  a.x * b.y + a.y * b.x);
}

float2 complex_add(float2 a, float2 b)
{
    return float2(a.x + b.x, a.y + b.y);
}

float2 complex_sub(float2 a, float2 b)
{
    return float2(a.x - b.x, a.y - b.y);
}

// Twiddle factor W_N^k = e^{-2πi k/N}  (forward); conjugated for inverse.
float2 twiddle_factor(uint k, uint N, bool inverse)
{
    float angle = -2.0f * FFT_PI * (float)k / (float)N;
    if (inverse) angle = -angle;
    return float2(cos(angle), sin(angle));
}

// Bit-reversal of the lower `bits` bits of x.
uint bit_reverse(uint x, uint bits)
{
    uint result = 0;
    for (uint i = 0; i < bits; ++i)
    {
        result = (result << 1) | (x & 1);
        x >>= 1;
    }
    return result;
}

// Integer log2 for power-of-two values.
uint log2_pot(uint n)
{
    uint k = 0;
    while ((1u << k) < n) ++k;
    return k;
}

// ---------------------------------------------------------------------------
// Kernel entry point
// ---------------------------------------------------------------------------
// One thread group handles one row (row pass) or one column (column pass).
// Each thread handles one element of the N-point transform.
//
// Threads with tid >= cb_N are idle (padding to THREADS_PER_GROUP).

[numthreads(THREADS_PER_GROUP, 1, 1)]
void FFTKernel_1D(
    uint3 GI  : SV_GroupThreadID,   // tid = GI.x  ∈ [0, THREADS_PER_GROUP)
    uint3 GID : SV_GroupID          // GID.y = row or column index
)
{
    const uint tid = GI.x;

    // Threads beyond the transform size do nothing.
    if (tid >= cb_N) return;

    // -------------------------------------------------------------------------
    // PHASE 1 — Load with bit-reversal permutation into groupshared memory
    // -------------------------------------------------------------------------
    // CPU equiv: fft_gpu_cpu_emu.cpp  PHASE 1 block

    const uint bits        = log2_pot(cb_N);
    const uint rev         = bit_reverse(tid, bits);

    // Global element index: row pass uses stride=1, column pass uses stride=total_width.
    const uint globalIdx   = cb_RowOffset + tid * cb_Stride;

    gs_Data[rev] = g_Data[globalIdx];

    GroupMemoryBarrierWithGroupSync();

    // -------------------------------------------------------------------------
    // PHASE 2 — Butterfly passes
    // -------------------------------------------------------------------------
    // CPU equiv: fft_gpu_cpu_emu.cpp  PHASE 2 block

    const bool inverse     = (cb_Inverse != 0);
    const uint numPasses   = bits;

    for (uint pass = 0; pass < numPasses; ++pass)
    {
        const uint span      = 1u << pass;
        const uint groupSize = span << 1;

        // Only the lower N/2 threads participate in each pass.
        if (tid < cb_N / 2)
        {
            const uint k       = tid % span;
            const uint evenIdx = (tid / span) * groupSize + k;
            const uint oddIdx  = evenIdx + span;

            float2 tw   = twiddle_factor(k, groupSize, inverse);
            float2 even = gs_Data[evenIdx];
            float2 odd  = complex_mul(tw, gs_Data[oddIdx]);

            // Sync before write to avoid race with threads reading the same slots.
            GroupMemoryBarrierWithGroupSync();

            gs_Data[evenIdx] = complex_add(even, odd);
            gs_Data[oddIdx]  = complex_sub(even, odd);
        }

        GroupMemoryBarrierWithGroupSync();
    }

    // -------------------------------------------------------------------------
    // PHASE 3 — Normalize (inverse only) and write back to global memory
    // -------------------------------------------------------------------------
    // CPU equiv: fft_gpu_cpu_emu.cpp  PHASE 3 block

    float scale = inverse ? (1.0f / (float)cb_N) : 1.0f;

    float2 result;
    result.x = gs_Data[tid].x * scale;
    result.y = gs_Data[tid].y * scale;

    g_Data[globalIdx] = result;
}
