#pragma once
#include <d3d11.h>
#include <d3dcompiler.h>
#include <vector>
#include <complex>
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
    // eWave Params
    int depthNum;
    // Padding for 16-byte alignment
    float buffer;
};

struct FFTConstants {
    int N;         // Transform size (gridsize) (must be power of two)
    int Bits;      // log2(N)  
    int Inverse;   // 0 = forward DFT, 1 = inverse DFT
    int Row;       // Row = 1. Col = 0
};

struct RenderConstants {
    DirectX::XMMATRIX viewProjection; // 64 bytes
    float gridSize;                   // 4 bytes
    float cellSize;                   // 4 bytes
    float buffer[2];                  // 8 bytes (Total = 80 bytes, must be multiple of 16)
};

struct GPUField {
    ID3D11Texture2D* tex;
    ID3D11ShaderResourceView* srv;   // read
    ID3D11UnorderedAccessView* uav;  // write
};

struct GPUBuffer {
    ID3D11Buffer* buf;
    ID3D11ShaderResourceView* srv;
    ID3D11UnorderedAccessView* uav;  // write
};

class GPU {
public:
    GPU();
    ~GPU();

    bool BaseInit(int size);
    bool Init(int size); // For headless compute
    bool Init(int size, HWND hwnd); // For rendering

    // Compute Shaders
    void UpdateConstants(const SimConstants& constants);
    bool CreateGridTexture(GPUField* field, int size, bool isComplex = false, int arraySize = 1);
    bool UploadToGPU(ID3D11Texture2D* tex, const std::vector<float>& data, int gridSize, bool isComplex = false, int arraySize = 1);
    bool DownloadFromGPU(ID3D11Texture2D* tex, std::vector<float>& data, int size);
    bool DownloadFromGPU(ID3D11Texture2D* tex, std::vector<std::complex<float>>& data, int size);
    bool DownloadFromGPU(ID3D11Texture2D* tex, std::vector<std::complex<float>>& data, int size, int arraySize);
    bool CreateBuffer(GPUBuffer* bufferObj, const void* initData, UINT elementCount);
    bool DownloadBuffer(ID3D11Buffer* buf, std::vector<float>& data, int size);
    bool CompileComputeShader(const std::wstring& file, const std::string& entryPoint, ID3D11ComputeShader** shader);
    void BindSRV(int slot, ID3D11ShaderResourceView* srv);
    void Dispatch(ID3D11ComputeShader* shader, 
                  const std::vector<ID3D11ShaderResourceView*>& srvs, 
                  const std::vector<ID3D11UnorderedAccessView*>& uavs,
                  int size = 1);

    // FFT
    void UpdateFFTConstants(const FFTConstants& constants);
    bool CompileFFTShaders(int size);
    void ExecuteFFT(ID3D11UnorderedAccessView* fftBufferUAV, int size, bool inverse, int numLayers = 1);

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
    ID3D11Texture2D* stagingComplexTex = nullptr;
    ID3D11Texture2D* stagingComplexArrayTex = nullptr;
    ID3D11Buffer* stagingBuf = nullptr;
    
    // FFT
    ID3D11Buffer* fftConstantBuffer = nullptr;
    ID3D11ComputeShader* fftShader = nullptr;
    ID3D11ComputeShader* fftArrayShader = nullptr;

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