// fft_gpu_cpu_emu.cpp
// CPU emulator of fft_gpu.hlsl.
//
// PORTING GUARANTEE:
//   Every function, index expression, and constant in this file has a
//   line-for-line equivalent in fft_gpu.hlsl.  Differences are limited to:
//     - C++ types  vs  HLSL types        (std::vector vs RWBuffer, etc.)
//     - Thread dispatch loop  vs  SV_GroupThreadID / SV_GroupID hardware
//     - groupshared emulated as a plain array on the stack
//
// Validate this emulator against fft_reference.cpp.  Once they agree,
// translate fft_gpu.hlsl mechanically using the comments as a guide.

#include "fft_gpu_cpu_emu.h"
#include "fft_reference.h"   // for float2, FFTDirection, FFT_PI

#include <cassert>
#include <cmath>
#include <cstring>
#include <stdexcept>
#include <vector>

// ---------------------------------------------------------------------------
// GPU constants (mirror [numthreads] in fft_gpu.hlsl)
// ---------------------------------------------------------------------------

// Number of threads per group — must equal the FFT size N for 1-D,
// and the row/column length for 2-D.  Change both here and in the HLSL.
static constexpr uint32_t THREADS_PER_GROUP = 512;   // = numthreads(512, 1, 1) in HLSL

// ---------------------------------------------------------------------------
// Shared primitives — identical to helper functions at top of fft_gpu.hlsl
// ---------------------------------------------------------------------------

// HLSL: float2 complex_mul(float2 a, float2 b)
static inline float2 complex_mul(float2 a, float2 b)
{
    return { a.x * b.x - a.y * b.y,
             a.x * b.y + a.y * b.x };
}

// HLSL: float2 complex_add(float2 a, float2 b)
static inline float2 complex_add(float2 a, float2 b)
{
    return { a.x + b.x, a.y + b.y };
}

// HLSL: float2 complex_sub(float2 a, float2 b)
static inline float2 complex_sub(float2 a, float2 b)
{
    return { a.x - b.x, a.y - b.y };
}

// HLSL: float2 twiddle_factor(uint k, uint N, bool inverse)
static inline float2 twiddle_factor(uint32_t k, uint32_t N, bool inverse)
{
    float angle = -2.0f * FFT_PI * static_cast<float>(k) / static_cast<float>(N);
    if (inverse) angle = -angle;
    return { std::cos(angle), std::sin(angle) };
}

// HLSL: uint bit_reverse(uint x, uint bits)
static inline uint32_t bit_reverse(uint32_t x, uint32_t bits)
{
    uint32_t result = 0;
    for (uint32_t i = 0; i < bits; ++i)
    {
        result = (result << 1) | (x & 1);
        x >>= 1;
    }
    return result;
}

// log2 for power-of-two values.
// HLSL: uint log2_pot(uint n)
static inline uint32_t log2_pot(uint32_t n)
{
    uint32_t k = 0;
    while ((1u << k) < n) ++k;
    return k;
}

// ---------------------------------------------------------------------------
// 1-D FFT kernel — emulates one compute shader dispatch
// ---------------------------------------------------------------------------
// HLSL equivalent:  FFTKernel_1D in fft_gpu.hlsl
//
// In HLSL this function is the body of [numthreads(N,1,1)].
// Here we loop over thread IDs explicitly to emulate the dispatch.
//
// groupshared float2 shared_data[THREADS_PER_GROUP]  →  stack array below.

static void FFTKernel_1D(float2* data, uint32_t N, bool inverse)
{
    assert(N <= THREADS_PER_GROUP && "Increase THREADS_PER_GROUP to match N.");

    // --- groupshared memory emulation ---
    // HLSL: groupshared float2 shared_data[THREADS_PER_GROUP];
    std::vector<float2> shared_data(N);

    // -------------------------------------------------------------------------
    // PHASE 1 — Load with bit-reversal permutation
    // Each thread loads one element into the bit-reversed slot.
    // -------------------------------------------------------------------------
    // HLSL:
    //   uint tid = GI.x;                         // SV_GroupThreadID.x
    //   uint rev = bit_reverse(tid, log2_pot(N));
    //   shared_data[rev] = input[tid];
    //   GroupMemoryBarrierWithGroupSync();

    const uint32_t bits = log2_pot(N);
    for (uint32_t tid = 0; tid < N; ++tid)   // <-- emulated thread loop
    {
        uint32_t rev = bit_reverse(tid, bits);
        shared_data[rev] = data[tid];
    }
    // GroupMemoryBarrierWithGroupSync() — implicit: sequential loop boundary.

    // -------------------------------------------------------------------------
    // PHASE 2 — Butterfly passes
    // -------------------------------------------------------------------------
    // HLSL:
    //   for (uint pass = 0; pass < log2_pot(N); ++pass)
    //   {
    //       uint span      = 1u << pass;
    //       uint groupSize = span << 1;
    //       uint k         = tid % span;
    //       uint evenIdx   = (tid / span) * groupSize + k;
    //       uint oddIdx    = evenIdx + span;
    //       float2 twiddle = twiddle_factor(k, groupSize, inverse);
    //       float2 even    = shared_data[evenIdx];
    //       float2 odd     = complex_mul(twiddle, shared_data[oddIdx]);
    //       GroupMemoryBarrierWithGroupSync();
    //       shared_data[evenIdx] = complex_add(even, odd);
    //       shared_data[oddIdx]  = complex_sub(even, odd);
    //       GroupMemoryBarrierWithGroupSync();
    //   }

    const uint32_t numPasses = bits;
    for (uint32_t pass = 0; pass < numPasses; ++pass)
    {
        const uint32_t span      = 1u << pass;
        const uint32_t groupSize = span << 1;

        // Pre-compute butterfly operands for all threads before writing back.
        // This emulates the mandatory GroupMemoryBarrierWithGroupSync() between
        // read and write phases in HLSL.
        std::vector<float2> even_vals(N / 2);
        std::vector<float2> odd_vals(N / 2);

        for (uint32_t tid = 0; tid < N / 2; ++tid)   // N/2 active threads per pass
        {
            const uint32_t k       = tid % span;
            const uint32_t evenIdx = (tid / span) * groupSize + k;
            const uint32_t oddIdx  = evenIdx + span;

            float2 tw = twiddle_factor(k, groupSize, inverse);
            even_vals[tid] = shared_data[evenIdx];
            odd_vals[tid]  = complex_mul(tw, shared_data[oddIdx]);
        }
        // GroupMemoryBarrierWithGroupSync()

        for (uint32_t tid = 0; tid < N / 2; ++tid)
        {
            const uint32_t k       = tid % span;
            const uint32_t evenIdx = (tid / span) * groupSize + k;
            const uint32_t oddIdx  = evenIdx + span;

            shared_data[evenIdx] = complex_add(even_vals[tid], odd_vals[tid]);
            shared_data[oddIdx]  = complex_sub(even_vals[tid], odd_vals[tid]);
        }
        // GroupMemoryBarrierWithGroupSync()
    }

    // -------------------------------------------------------------------------
    // PHASE 3 — Normalize (inverse only) and write back to global memory
    // -------------------------------------------------------------------------
    // HLSL:
    //   float scale = inverse ? (1.0f / (float)N) : 1.0f;
    //   output[tid] = shared_data[tid] * scale;   // float2 * scalar

    const float scale = inverse ? (1.0f / static_cast<float>(N)) : 1.0f;
    for (uint32_t tid = 0; tid < N; ++tid)
    {
        data[tid].x = shared_data[tid].x * scale;
        data[tid].y = shared_data[tid].y * scale;
    }
}

// ---------------------------------------------------------------------------
// 2-D FFT — row pass then column pass (two separate dispatches in HLSL)
// ---------------------------------------------------------------------------
// HLSL equivalent: two Dispatch() calls in C++, one for rows, one for columns.
//   Dispatch(width / THREADS_PER_GROUP, height, 1)   -- rows
//   Dispatch(height / THREADS_PER_GROUP, width, 1)   -- columns

static void FFT2D_GPU(float2* data, uint32_t width, uint32_t height, bool inverse)
{
    // --- Row pass ---
    // HLSL: each group processes one row.  GI.y = row index.
    for (uint32_t row = 0; row < height; ++row)
        FFTKernel_1D(data + row * width, width, inverse);

    // --- Column pass ---
    // HLSL: each group processes one column.
    // Columns are non-contiguous in row-major layout, so we copy to a temp buffer,
    // transform, and write back — same cost as a column-major UAV access on GPU.
    std::vector<float2> col_buf(height);
    for (uint32_t col = 0; col < width; ++col)
    {
        for (uint32_t row = 0; row < height; ++row)
            col_buf[row] = data[row * width + col];

        FFTKernel_1D(col_buf.data(), height, inverse);

        for (uint32_t row = 0; row < height; ++row)
            data[row * width + col] = col_buf[row];
    }
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

void FFT_GPU_Emu_1D(float2* data, uint32_t N, FFTDirection dir)
{
    if (N == 0 || (N & (N - 1)) != 0)
        throw std::invalid_argument("N must be a power of two.");
    if (N > THREADS_PER_GROUP)
        throw std::invalid_argument("N exceeds THREADS_PER_GROUP; increase it.");

    FFTKernel_1D(data, N, dir == FFTDirection::Inverse);
}

void FFT_GPU_Emu_2D(float2* data, uint32_t width, uint32_t height, FFTDirection dir)
{
    if ((width & (width - 1)) != 0 || (height & (height - 1)) != 0)
        throw std::invalid_argument("Width and height must be powers of two.");
    if (width > THREADS_PER_GROUP || height > THREADS_PER_GROUP)
        throw std::invalid_argument("Dimension exceeds THREADS_PER_GROUP; increase it.");

    FFT2D_GPU(data, width, height, dir == FFTDirection::Inverse);
}
