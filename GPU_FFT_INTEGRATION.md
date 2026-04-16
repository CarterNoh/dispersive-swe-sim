/*
INTEGRATION EXAMPLE: How to use GPU FFT in Sim2D.cpp

In your Sim2D.h, add:
    #include "fft/fft_gpu_wrapper.h"
    std::unique_ptr<GPUFFTCompute> gpuFFT;

In your Sim2D constructor, initialize:
    Sim::Sim(...) {
        // Initialize GPU FFT for 2D grid
        gpuFFT = std::make_unique<GPUFFTCompute>(GRIDSIZE, GRIDSIZE);
        if (!gpuFFT->Initialize()) {
            throw std::runtime_error("Failed to initialize GPU FFT");
        }
    }

In your Sim2D::SimStep() or where you need FFT, call:

    // *** NEW CONVENIENCE API ***
    // Forward FFT: real input -> complex output
    std::vector<float2> htilde_fft = gpuFFT->ComputeFFT(htilde, false);

    // Do stuff to the frequency-domain data (e.g., wave dispersion correction)
    // ... modify htilde_fft ...

    // Inverse FFT: complex input -> real output
    htilde = gpuFFT->ComputeFFT(htilde_fft, true);

    // *** ALTERNATIVE: Complex-to-complex FFT ***
    // For advanced usage where you need both real and imaginary parts
    std::vector<float2> complex_result = gpuFFT->ComputeFFTComplex(complex_input, false);

    // *** LEGACY API (still available) ***
    // Manual step-by-step approach (same as before)
    std::vector<float2> htildeComplex(GRIDSIZE * GRIDSIZE);
    for (int i = 0; i < GRIDSIZE * GRIDSIZE; ++i) {
        htildeComplex[i] = {htilde[i], 0.0f};
    }
    gpuFFT->UploadData(htildeComplex);
    gpuFFT->ComputeFFTRows(false);
    gpuFFT->ComputeFFTColumns(false);
    gpuFFT->DownloadData(htildeComplex);

COMPILATION:
    1. Ensure dxc.exe (DirectX Shader Compiler) is installed
       - Download from: https://github.com/microsoft/DirectXShaderCompiler/releases
       - Or install via: vcpkg install directxshadercompiler
       - Or via Windows SDK

    2. Build the project:
       mkdir build
       cd build
       cmake ..
       cmake --build .

    3. The CMakeLists.txt will:
       - Compile fft_gpu.hlsl to fft_gpu.cso
       - Link against Direct3D 12
       - Include the GPU wrapper in the build

PERFORMANCE NOTES:
    - GPU FFT is efficient for large 2D FFTs (256x256 or larger)
    - For small grids, CPU emulator may be faster (less overhead)
    - Memory upload/download has latency; batch FFT operations if possible
    - Each dispatch waits for GPU completion via fence (synchronous)
    - For async execution, modify WaitForGPU() method

TROUBLESHOOTING:
    - "fft_gpu.cso not found": Ensure dxc.exe compiled the shader successfully
    - "D3D12CreateDevice failed": Check DirectX 12 drivers are installed
    - "shader compilation error": Check HLSL syntax and dxc.exe version (6.4+)
    - GPU memory errors: Reduce GRIDSIZE or increase GPU VRAM
 */
