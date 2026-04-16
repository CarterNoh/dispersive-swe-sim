// Example: How to integrate GPU FFT into Sim2D.cpp
//
// Add this to Sim2D.h:
//     #include "fft/fft_gpu_wrapper.h"
//     std::unique_ptr<GPUFFTCompute> gpuFFT;
//
// Add this to Sim2D constructor:
//     gpuFFT = std::make_unique<GPUFFTCompute>(GRIDSIZE, GRIDSIZE);
//     if (!gpuFFT->Initialize()) {
//         throw std::runtime_error("GPU FFT initialization failed");
//     }
//
// Replace the commented eWaveStep() call in SimStep() with:

void Sim::FFTStep()
{
    // Example: Apply FFT-based wave dispersion correction to surface height

    // 1. Forward FFT: Convert real surface height to frequency domain
    std::vector<float2> htilde_fft = gpuFFT->ComputeFFT(htilde, false);

    // 2. Apply dispersion correction in frequency domain
    // (This is where you'd implement your wave dispersion physics)
    for (size_t i = 0; i < htilde_fft.size(); ++i) {
        // Example: Simple damping based on frequency
        // In practice, you'd compute proper dispersion relations
        float kx = (i % GRIDSIZE) - GRIDSIZE/2;  // wave number in x
        float ky = (i / GRIDSIZE) - GRIDSIZE/2;  // wave number in y
        float k_squared = kx*kx + ky*ky;

        if (k_squared > 0) {
            // Apply some dispersion correction (this is just an example)
            float dispersion_factor = exp(-k_squared * 0.001f);  // damping
            htilde_fft[i].x *= dispersion_factor;
            htilde_fft[i].y *= dispersion_factor;
        }
    }

    // 3. Inverse FFT: Convert back to spatial domain
    htilde = gpuFFT->ComputeFFT(htilde_fft, true);

    // Note: The FFT is unnormalized, so the result will have the correct scale
    // If you need normalized FFT, divide by (GRIDSIZE * GRIDSIZE) after forward FFT
}

// Alternative: If you need to work with both real and imaginary parts
void Sim::FFTStep_Complex()
{
    // Convert htilde to complex format manually
    std::vector<float2> htilde_complex(GRIDSIZE * GRIDSIZE);
    for (int i = 0; i < GRIDSIZE * GRIDSIZE; ++i) {
        htilde_complex[i] = {htilde[i], 0.0f};
    }

    // Forward FFT
    std::vector<float2> htilde_fft = gpuFFT->ComputeFFTComplex(htilde_complex, false);

    // Apply your dispersion correction here...
    // (modify htilde_fft)

    // Inverse FFT
    std::vector<float2> htilde_result = gpuFFT->ComputeFFTComplex(htilde_fft, true);

    // Extract real part back to htilde
    for (int i = 0; i < GRIDSIZE * GRIDSIZE; ++i) {
        htilde[i] = htilde_result[i].x;
    }
}