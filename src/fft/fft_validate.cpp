// fft_validate.cpp
// Cross-validates fft_reference.cpp against fft_gpu_cpu_emu.cpp.
//
// Build (example):
//   g++ -std=c++17 -O2 fft_validate.cpp fft_reference.cpp fft_gpu_cpu_emu.cpp -o fft_validate
//
// All tests print PASS or FAIL with the maximum absolute error found.

#include "fft_reference.h"
#include "fft_gpu_cpu_emu.h"

#include <cmath>
#include <cstdio>
#include <cstring>
#include <vector>
#include <random>

// ---------------------------------------------------------------------------
// Utilities
// ---------------------------------------------------------------------------

static float MaxAbsError(const float2* a, const float2* b, uint32_t n)
{
    float maxErr = 0.0f;
    for (uint32_t i = 0; i < n; ++i)
    {
        float ex = std::abs(a[i].x - b[i].x);
        float ey = std::abs(a[i].y - b[i].y);
        maxErr = std::max(maxErr, std::max(ex, ey));
    }
    return maxErr;
}

static void FillRandom(float2* data, uint32_t n, uint32_t seed = 42)
{
    std::mt19937 rng(seed);
    std::uniform_real_distribution<float> dist(-1.0f, 1.0f);
    for (uint32_t i = 0; i < n; ++i)
        data[i] = { dist(rng), dist(rng) };
}

static void PrintResult(const char* name, float maxErr, float tol = 1e-4f)
{
    bool pass = maxErr < tol;
    std::printf("  %-40s  max_err=%.2e  %s\n", name, maxErr, pass ? "PASS" : "FAIL");
}

// ---------------------------------------------------------------------------
// Test 1: Forward FFT agreement (1-D)
// ---------------------------------------------------------------------------

static void Test_Forward_1D(uint32_t N)
{
    std::vector<float2> ref(N), gpu(N);
    FillRandom(ref.data(), N);
    std::memcpy(gpu.data(), ref.data(), N * sizeof(float2));

    FFT_Reference_1D(ref.data(), N, FFTDirection::Forward);
    FFT_GPU_Emu_1D (gpu.data(), N, FFTDirection::Forward);

    char label[64];
    std::snprintf(label, sizeof(label), "Forward 1-D N=%u", N);
    PrintResult(label, MaxAbsError(ref.data(), gpu.data(), N));
}

// ---------------------------------------------------------------------------
// Test 2: Inverse FFT agreement (1-D)
// ---------------------------------------------------------------------------

static void Test_Inverse_1D(uint32_t N)
{
    std::vector<float2> ref(N), gpu(N);
    FillRandom(ref.data(), N);
    std::memcpy(gpu.data(), ref.data(), N * sizeof(float2));

    FFT_Reference_1D(ref.data(), N, FFTDirection::Inverse);
    FFT_GPU_Emu_1D (gpu.data(), N, FFTDirection::Inverse);

    char label[64];
    std::snprintf(label, sizeof(label), "Inverse 1-D N=%u", N);
    PrintResult(label, MaxAbsError(ref.data(), gpu.data(), N));
}

// ---------------------------------------------------------------------------
// Test 3: Round-trip (forward then inverse recovers original signal)
// ---------------------------------------------------------------------------

static void Test_RoundTrip_1D(uint32_t N)
{
    std::vector<float2> original(N), data(N);
    FillRandom(original.data(), N);
    std::memcpy(data.data(), original.data(), N * sizeof(float2));

    FFT_GPU_Emu_1D(data.data(), N, FFTDirection::Forward);
    FFT_GPU_Emu_1D(data.data(), N, FFTDirection::Inverse);

    char label[64];
    std::snprintf(label, sizeof(label), "Round-trip 1-D N=%u", N);
    PrintResult(label, MaxAbsError(original.data(), data.data(), N));
}

// ---------------------------------------------------------------------------
// Test 4: Forward FFT agreement (2-D)
// ---------------------------------------------------------------------------

static void Test_Forward_2D(uint32_t W, uint32_t H)
{
    uint32_t n = W * H;
    std::vector<float2> ref(n), gpu(n);
    FillRandom(ref.data(), n);
    std::memcpy(gpu.data(), ref.data(), n * sizeof(float2));

    FFT_Reference_2D(ref.data(), W, H, FFTDirection::Forward);
    FFT_GPU_Emu_2D (gpu.data(), W, H, FFTDirection::Forward);

    char label[64];
    std::snprintf(label, sizeof(label), "Forward 2-D %ux%u", W, H);
    PrintResult(label, MaxAbsError(ref.data(), gpu.data(), n));
}

// ---------------------------------------------------------------------------
// Test 5: Round-trip (2-D)
// ---------------------------------------------------------------------------

static void Test_RoundTrip_2D(uint32_t W, uint32_t H)
{
    uint32_t n = W * H;
    std::vector<float2> original(n), data(n);
    FillRandom(original.data(), n);
    std::memcpy(data.data(), original.data(), n * sizeof(float2));

    FFT_GPU_Emu_2D(data.data(), W, H, FFTDirection::Forward);
    FFT_GPU_Emu_2D(data.data(), W, H, FFTDirection::Inverse);

    char label[64];
    std::snprintf(label, sizeof(label), "Round-trip 2-D %ux%u", W, H);
    PrintResult(label, MaxAbsError(original.data(), data.data(), n));
}

// ---------------------------------------------------------------------------
// Test 6: Known DFT values — impulse at index 0 → flat spectrum
// ---------------------------------------------------------------------------

static void Test_Impulse_Spectrum(uint32_t N)
{
    std::vector<float2> data(N, {0.0f, 0.0f});
    data[0] = {1.0f, 0.0f};   // impulse

    FFT_GPU_Emu_1D(data.data(), N, FFTDirection::Forward);

    // All output bins should equal 1 + 0i.
    float maxErr = 0.0f;
    for (uint32_t i = 0; i < N; ++i)
    {
        maxErr = std::max(maxErr, std::abs(data[i].x - 1.0f));
        maxErr = std::max(maxErr, std::abs(data[i].y - 0.0f));
    }

    char label[64];
    std::snprintf(label, sizeof(label), "Impulse spectrum N=%u", N);
    PrintResult(label, maxErr);
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------

int main()
{
    std::printf("=== FFT Validation Suite ===\n\n");

    std::printf("[ 1-D Forward agreement ]\n");
    Test_Forward_1D(8);
    Test_Forward_1D(64);
    Test_Forward_1D(256);
    Test_Forward_1D(512);

    std::printf("\n[ 1-D Inverse agreement ]\n");
    Test_Inverse_1D(8);
    Test_Inverse_1D(64);
    Test_Inverse_1D(256);
    Test_Inverse_1D(512);

    std::printf("\n[ 1-D Round-trip (GPU emu) ]\n");
    Test_RoundTrip_1D(8);
    Test_RoundTrip_1D(64);
    Test_RoundTrip_1D(256);
    Test_RoundTrip_1D(512);

    std::printf("\n[ 2-D Forward agreement ]\n");
    Test_Forward_2D(16,  16);
    Test_Forward_2D(64,  64);
    Test_Forward_2D(256, 256);
    Test_Forward_2D(512, 512);

    std::printf("\n[ 2-D Round-trip (GPU emu) ]\n");
    Test_RoundTrip_2D(16,  16);
    Test_RoundTrip_2D(64,  64);
    Test_RoundTrip_2D(256, 256);
    Test_RoundTrip_2D(512, 512);

    std::printf("\n[ Known spectral values ]\n");
    Test_Impulse_Spectrum(8);
    Test_Impulse_Spectrum(64);
    Test_Impulse_Spectrum(512);

    std::printf("\nDone.\n");
    return 0;
}
