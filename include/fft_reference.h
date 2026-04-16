// fft_reference.h
// Scalar C++ reference FFT — public interface.
// Include this in test/validation code and in fft_gpu_cpu_emu.cpp for cross-checks.

#pragma once

#include <cstdint>

// ---------------------------------------------------------------------------
// Shared math primitive: 2-component float vector.
// Identical layout to HLSL float2 and the GPU emulator's float2.
// Do NOT replace with std::complex<float> — keeping the plain struct means
// a direct memory-level comparison between CPU and GPU output buffers.
// ---------------------------------------------------------------------------

#ifndef FLOAT2_DEFINED
#define FLOAT2_DEFINED
struct float2
{
    float x, y;   // x = real, y = imaginary
};
#endif

// Pi constant shared with GPU emulator and HLSL.
#ifndef FFT_PI
#define FFT_PI 3.14159265358979323846f
#endif

// ---------------------------------------------------------------------------
// Transform direction
// ---------------------------------------------------------------------------

enum class FFTDirection
{
    Forward,   // DFT  — e^{-2πi}  sign convention (matches HLSL)
    Inverse    // IDFT — e^{+2πi}, normalized by 1/N
};

// ---------------------------------------------------------------------------
// API
// ---------------------------------------------------------------------------

// 1-D in-place FFT/IFFT.
// data : array of N complex samples (float2), modified in-place.
// N    : must be a power of two.
void FFT_Reference_1D(float2* data, uint32_t N, FFTDirection dir);

// 2-D in-place FFT/IFFT via row-column decomposition.
// data : row-major array of width*height complex samples, modified in-place.
// width, height : must each be powers of two.
void FFT_Reference_2D(float2* data, uint32_t width, uint32_t height, FFTDirection dir);
