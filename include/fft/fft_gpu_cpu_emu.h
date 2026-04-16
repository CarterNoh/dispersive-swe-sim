// fft_gpu_cpu_emu.h
// CPU emulator of the HLSL FFT compute shader — public interface.
//
// Call these functions to validate that fft_gpu.hlsl produces the same
// numerical results as fft_reference.cpp before dispatching on the GPU.

#pragma once

#include "fft_reference.h"   // float2, FFTDirection, FFT_PI

// 1-D in-place FFT using GPU-style groupshared butterfly algorithm.
// N must be a power of two and <= THREADS_PER_GROUP (512 by default).
void FFT_GPU_Emu_1D(float2* data, uint32_t N, FFTDirection dir);

// 2-D in-place FFT via two sequential row/column dispatch emulations.
// width and height must each be powers of two and <= THREADS_PER_GROUP.
void FFT_GPU_Emu_2D(float2* data, uint32_t width, uint32_t height, FFTDirection dir);
