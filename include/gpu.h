#pragma once
#include <d3d11.h>
#include <d3dcompiler.h>
#include <vector>
#include <string>
#include <iostream>
#include <dxgi.h>
#include <DirectXMath.h> // For the ViewProjection matrix

// The constant buffer must be padded to a multiple of 16 bytes for HLSL
struct SimConstants {
    // Simulation params
    int gridSize;
    float cellSize;
    float timeStep;
    int boundaryType;
    float minWaterHeight;
    // Decomposition params
    int diffusionIterations;
    float deltaT;
    float diffusionPenalty;
    // SWE & Transport Params
    float cflCondition;
    float gammaTransport;
    // Padding for 16-byte alignment
    float buffer[2];
};

struct RenderConstants {
    DirectX::XMMATRIX viewProjection; // 64 bytes
    float gridSize;                   // 4 bytes
    float cellSize;                   // 4 bytes
    float vuffer[2];                  // 8 bytes (Total = 80 bytes, must me multiple of 16)
};

struct FFTConstants {
    uint32_t N;         // Transform size (gridsize) (must be power of two)
    uint32_t Inverse;   // 0 = forward DFT, 1 = inverse DFT
    uint32_t Stride;    // Element stride between samples
    uint32_t Bits;      // log2(N)  
};

struct GPUField {
    ID3D11Texture2D* tex;
    ID3D11ShaderResourceView* srv;   // read
    ID3D11UnorderedAccessView* uav;  // write
};

class GPU {
public:
    GPU();
    ~GPU();

    bool Init(int size); // For headless compute
    bool Init(int size, HWND hwnd);

    // Compute Shaders
    void UpdateConstants(const SimConstants& constants);
    bool CreateGridTexture(GPUField* field, int size, bool isComplex);
    bool UploadToGPU(ID3D11Texture2D* tex, const std::vector<float>& data, int size);
    bool DownloadFromGPU(ID3D11Texture2D* tex, std::vector<float>& data, int size);
    // void ClearUAV(ID3D11UnorderedAccessView* uav, float clearValue);
    bool CompileComputeShader(const std::wstring& file, const std::string& entryPoint, ID3D11ComputeShader** shader);
    void Dispatch(ID3D11ComputeShader* shader, 
                  const std::vector<ID3D11ShaderResourceView*>& srvs, 
                  const std::vector<ID3D11UnorderedAccessView*>& uavs);

    // FFT
    void UpdateFFTConstants(const FFTConstants& constants);
    bool CompileFFTShaders(int size);
    bool UploadToFFT(ID3D11ShaderResourceView* texSRV, ID3D11UnorderedAccessView* fftUAV);
    bool DownloadFromFFT(ID3D11ShaderResourceView* fftSRV, ID3D11UnorderedAccessView* texUAV);
    void ExecuteFFT(ID3D11UnorderedAccessView* fftBufferUAV, int size, bool inverse);

    // Rendering 
    void UpdateRenderConstants(const RenderConstants& constants);
    bool CreateGridVertexBuffer(int gridSize);
    bool CreateGridMesh(int gridSize);
    bool CompileVertexShader(const std::wstring& file, const std::string& entryPoint);
    bool CompilePixelShader(const std::wstring& file, const std::string& entryPoint);
    void Render(ID3D11ShaderResourceView* heightSRV);
    

private:
    // Compute Shaders
    int size;   // simulation grid size
    UINT group; // number of thread groups to dispatch for compute shaders, assuming 16x16 threads per group
    ID3D11Device* device = nullptr;
    ID3D11DeviceContext* context = nullptr;
    ID3D11Buffer* constantBuffer = nullptr;
    ID3D11Texture2D* stagingTex = nullptr; // Used to download data back to CPU
    

    // FFT
    ID3D11Buffer* fftConstantBuffer = nullptr;
    ID3D11ComputeShader* fftShader = nullptr;
    ID3D11ComputeShader* fftUploadShader = nullptr;
    ID3D11ComputeShader* fftDownloadShader = nullptr;

    // Rendering
    IDXGISwapChain* swapChain = nullptr;
    ID3D11RenderTargetView* renderTargetView = nullptr;
    ID3D11VertexShader* vertexShader = nullptr;
    ID3D11PixelShader* pixelShader = nullptr;
    ID3D11Buffer* vertexBuffer = nullptr;
    ID3D11Buffer* renderConstantBuffer = nullptr;
    ID3D11SamplerState* samplerState = nullptr;
    ID3D11InputLayout* inputLayout = nullptr;
    ID3D11RasterizerState* rasterState = nullptr;
    ID3D11Texture2D* depthStencilBuffer = nullptr;
    ID3D11DepthStencilView* depthStencilView = nullptr;
    ID3D11Buffer* indexBuffer = nullptr;
    UINT indexCount = 0; // We need to remember how many indices to draw

    
};