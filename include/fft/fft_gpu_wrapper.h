#pragma once

#include <windows.h>
#include <d3d12.h>
#include <dxgi1_6.h>
#include <vector>
#include <memory>

struct float2 {
    float x, y;
};

class GPUFFTCompute {
public:
    GPUFFTCompute(uint32_t width, uint32_t height);
    ~GPUFFTCompute();

    // Initialize DirectX 12 device and command queue
    bool Initialize();

    // Compute 1-D FFT row pass
    // N = width, height = number of rows
    bool ComputeFFTRows(bool inverse = false);

    // Compute 1-D FFT column pass
    // N = height, width = number of columns
    bool ComputeFFTColumns(bool inverse = false);

    // Upload data to GPU buffer
    bool UploadData(const std::vector<float2>& data);

    // Download data from GPU buffer
    bool DownloadData(std::vector<float2>& outData);

    // Convenience functions for complete FFT operations
    // Forward FFT: real input -> complex output
    std::vector<float2> ComputeFFT(const std::vector<float>& realData, bool inverse = false);
    
    // Inverse FFT: complex input -> real output  
    std::vector<float> ComputeFFT(const std::vector<float2>& complexData, bool inverse = true);
    
    // Complex-to-complex FFT (for advanced usage)
    std::vector<float2> ComputeFFTComplex(const std::vector<float2>& complexData, bool inverse = false);

    // Get total number of elements
    uint32_t GetElementCount() const { return width * height; }

private:
    // DirectX 12 objects
    Microsoft::WRL::ComPtr<ID3D12Device> device;
    Microsoft::WRL::ComPtr<ID3D12CommandQueue> commandQueue;
    Microsoft::WRL::ComPtr<ID3D12CommandAllocator> commandAllocator;
    Microsoft::WRL::ComPtr<ID3D12GraphicsCommandList> commandList;
    Microsoft::WRL::ComPtr<ID3D12Fence> fence;
    HANDLE fenceEvent;
    uint64_t fenceValue;

    // Shader resources
    Microsoft::WRL::ComPtr<ID3D12PipelineState> pipelineState;
    Microsoft::WRL::ComPtr<ID3D12RootSignature> rootSignature;

    // Buffers
    Microsoft::WRL::ComPtr<ID3D12Resource> dataBuffer;      // UAV for FFT data
    Microsoft::WRL::ComPtr<ID3D12Resource> uploadBuffer;    // Upload heap
    Microsoft::WRL::ComPtr<ID3D12Resource> downloadBuffer;  // Readback heap
    Microsoft::WRL::ComPtr<ID3D12DescriptorHeap> descriptorHeap;

    // Constants
    struct FFTConstants {
        uint32_t N;
        uint32_t Inverse;
        uint32_t RowOffset;
        uint32_t Stride;
    };
    Microsoft::WRL::ComPtr<ID3D12Resource> constantBuffer;

    // Dimensions
    uint32_t width;
    uint32_t height;

    // Helper methods
    bool CreateDevice();
    bool LoadShader();
    bool CreateBuffers();
    bool CreateRootSignature();
    bool WaitForGPU();
    bool UpdateConstants(uint32_t N, uint32_t inverse, uint32_t rowOffset, uint32_t stride);
};
