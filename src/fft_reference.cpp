// fft_reference.cpp
// Scalar Cooley-Tukey Butterfly FFT — CPU reference implementation.
// Purpose: correctness baseline and numerical ground truth for validating
//          the GPU emulator (fft_gpu_cpu_emu.cpp) and HLSL shader.
//
// Design rules:
//   - No SIMD, no threading, no external dependencies.
//   - Indexing and loop structure mirror the HLSL compute shader exactly.
//   - All twiddle factors computed the same way as the GPU path.
//   - float2 is used as a plain struct so diffs against the GPU emu are trivial.

#include "fft_reference.h"

#include <cassert>
#include <cmath>
#include <cstring>
#include <stdexcept>
#include <vector>

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

static inline bool IsPowerOfTwo(uint32_t n)
{
    return n > 0 && (n & (n - 1)) == 0;
}

// Compute the number of bits needed to represent indices in an N-element array.
static inline uint32_t Log2(uint32_t n)
{
    uint32_t k = 0;
    while ((1u << k) < n) ++k;
    return k;
}

// Bit-reversal permutation: swap element i with bit-reversed index of i.
// Operates in-place. N must be a power of two.
static void BitReversePermutation(float2* data, uint32_t N)
{
    const uint32_t bits = Log2(N);
    for (uint32_t i = 0; i < N; ++i)
    {
        // Reverse the lower `bits` bits of i.
        uint32_t j = 0;
        uint32_t tmp = i;
        for (uint32_t b = 0; b < bits; ++b)
        {
            j = (j << 1) | (tmp & 1);
            tmp >>= 1;
        }
        if (j > i)
        {
            float2 t = data[i];
            data[i]  = data[j];
            data[j]  = t;
        }
    }
}

// Complex multiply: (a + bi)(c + di) = (ac-bd) + (ad+bc)i
static inline float2 CMul(float2 a, float2 b)
{
    return { a.x * b.x - a.y * b.y,
             a.x * b.y + a.y * b.x };
}

// Complex add
static inline float2 CAdd(float2 a, float2 b)
{
    return { a.x + b.x, a.y + b.y };
}

// Complex subtract
static inline float2 CSub(float2 a, float2 b)
{
    return { a.x - b.x, a.y - b.y };
}

// Twiddle factor W_N^k = e^{-2πi k/N}  (forward transform sign convention).
// Matches twiddle_factor() in fft_gpu_cpu_emu.cpp and fft_gpu.hlsl.
static inline float2 TwiddleFactor(uint32_t k, uint32_t N)
{
    const float angle = -2.0f * FFT_PI * static_cast<float>(k) / static_cast<float>(N);
    return { std::cos(angle), std::sin(angle) };
}

// ---------------------------------------------------------------------------
// Core 1-D FFT (in-place, Cooley-Tukey DIT)
// ---------------------------------------------------------------------------
// Pass structure:
//   pass 0 : butterfly span = 1  (pairs of adjacent elements)
//   pass 1 : butterfly span = 2
//   ...
//   pass p : butterfly span = 2^p
//
// This matches the HLSL compute shader's dispatch loop exactly.

static void FFT1D_InPlace(float2* data, uint32_t N, FFTDirection dir)
{
    assert(IsPowerOfTwo(N));

    BitReversePermutation(data, N);

    const uint32_t numPasses = Log2(N);

    for (uint32_t pass = 0; pass < numPasses; ++pass)
    {
        const uint32_t span       = 1u << pass;        // half-size of butterfly group
        const uint32_t groupSize  = span << 1;         // full butterfly group size

        for (uint32_t groupStart = 0; groupStart < N; groupStart += groupSize)
        {
            for (uint32_t k = 0; k < span; ++k)
            {
                // --- twiddle factor ---
                // Forward: e^{-2πi k/groupSize}
                // Inverse: conjugate of forward
                float2 twiddle = TwiddleFactor(k, groupSize);
                if (dir == FFTDirection::Inverse)
                    twiddle.y = -twiddle.y;

                const uint32_t evenIdx = groupStart + k;
                const uint32_t oddIdx  = evenIdx + span;

                float2 even = data[evenIdx];
                float2 odd  = CMul(twiddle, data[oddIdx]);

                data[evenIdx] = CAdd(even, odd);   // butterfly top
                data[oddIdx]  = CSub(even, odd);   // butterfly bottom
            }
        }
    }

    // Normalize inverse transform.
    if (dir == FFTDirection::Inverse)
    {
        const float invN = 1.0f / static_cast<float>(N);
        for (uint32_t i = 0; i < N; ++i)
        {
            data[i].x *= invN;
            data[i].y *= invN;
        }
    }
}

// ---------------------------------------------------------------------------
// 2-D FFT (row-column decomposition)
// ---------------------------------------------------------------------------

static void FFT2D_InPlace(float2* data, uint32_t width, uint32_t height, FFTDirection dir)
{
    // Transform each row.
    for (uint32_t row = 0; row < height; ++row)
        FFT1D_InPlace(data + row * width, width, dir);

    // Transform each column (requires a temporary column buffer).
    std::vector<float2> col(height);
    for (uint32_t col_idx = 0; col_idx < width; ++col_idx)
    {
        for (uint32_t row = 0; row < height; ++row)
            col[row] = data[row * width + col_idx];

        FFT1D_InPlace(col.data(), height, dir);

        for (uint32_t row = 0; row < height; ++row)
            data[row * width + col_idx] = col[row];
    }
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

void FFT_Reference_1D(float2* data, uint32_t N, FFTDirection dir)
{
    if (!IsPowerOfTwo(N))
        throw std::invalid_argument("FFT size must be a power of two.");
    FFT1D_InPlace(data, N, dir);
}

void FFT_Reference_2D(float2* data, uint32_t width, uint32_t height, FFTDirection dir)
{
    if (!IsPowerOfTwo(width) || !IsPowerOfTwo(height))
        throw std::invalid_argument("FFT dimensions must be powers of two.");
    FFT2D_InPlace(data, width, height, dir);
}
