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
struct alignas(16) SimConstants {
    float time;
    // Sim Params
    int gridSizeX; 
    int gridSizeY; 
    float cellSize;
    float timeStep;
    int spongeThickness;
    float minWaterHeight;
    float surfaceTension;
    float density;
    // Decomposition Params
    int diffusionIterations;
    float deltaT;
    float diffusionPenalty;
    // SWE & Transport Params
    float slopeLimit;
    float cflCondition;
    float gammaTransport;
    // eWave Params
    int depthNum;
    // FFT wave params
    float fetch;
    float windSpeed;
    float windAngle;
    float swell;
    float swellAngle;
    float choppiness;
    float filterSmall;
    float filterBig;
    float filterWidth;
    float filterMin;
    float depthCutoff;
    int paddedGridSizeX;
    int paddedGridSizeY;
    float simConstantPadding[3]; // Align to 16 bytes
};

struct alignas(16) FFTConstants {
    int Nx;        // Transform size X (width)
    int Ny;        // Transform size Y (height)
    int BitsX;     // log2(Nx)
    int BitsY;     // log2(Ny)
    int Inverse;   // 0 = forward DFT, 1 = inverse DFT
    int Row;       // Row = 1. Col = 0
};

struct alignas(16) RenderConstants {
    DirectX::XMMATRIX viewProjection;
    float gridSizeX;
    float gridSizeY;
    float cellSize;
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

    bool BaseInit(int sizeX, int sizeY);
    bool Init(int sizeX, int sizeY); // For headless compute
    bool Init(int sizeX, int sizeY, HWND hwnd); // For rendering

    int GetPaddedSizeX() const { return paddedSizeX; }
    int GetPaddedSizeY() const { return paddedSizeY; }

    // Compute Shaders    
    void UpdateConstants(const SimConstants& constants);
    bool CreateGridTexture(GPUField* field, int width, int height, bool isComplex = false, int arraySize = 1);
    bool UploadToGPU(ID3D11Texture2D* tex, const std::vector<float>& data, int width, int height, bool isComplex = false, int arraySize = 1);
    bool DownloadFromGPU(ID3D11Texture2D* tex, std::vector<float>& data, int width, int height);
    bool DownloadFromGPU(ID3D11Texture2D* tex, std::vector<std::complex<float>>& data, int width, int height);
    bool DownloadFromGPU(ID3D11Texture2D* tex, std::vector<std::complex<float>>& data, int width, int height, int arraySize);
    bool CreateBuffer(GPUBuffer* bufferObj, const void* initData, UINT elementCount);
    bool DownloadBuffer(ID3D11Buffer* buf, std::vector<float>& data, int size);
    bool CompileComputeShader(const std::wstring& file, const std::string& entryPoint, ID3D11ComputeShader** shader);
    void BindSRV(int slot, ID3D11ShaderResourceView* srv);
    void Dispatch(ID3D11ComputeShader* shader, 
                  const std::vector<ID3D11ShaderResourceView*>& srvs, 
                  const std::vector<ID3D11UnorderedAccessView*>& uavs,
                  int layers = 1);
    void DispatchPadded(ID3D11ComputeShader* shader, 
                        const std::vector<ID3D11ShaderResourceView*>& srvs, 
                        const std::vector<ID3D11UnorderedAccessView*>& uavs,
                        int layers = 1);

    // FFT
    void UpdateFFTConstants(const FFTConstants& constants);
    bool CompileFFTShaders(int sizeX, int sizeY);
    void ExecuteFFT(ID3D11UnorderedAccessView* fftBufferUAV, int sizeX, int sizeY, bool inverse, int numLayers = 1);

    // Rendering 
    void UpdateRenderConstants(const RenderConstants& constants);
    bool CreateGridVertexBuffer(int width, int height);
    bool CreateGridMesh(int width, int height);
    bool CompileVertexShader(const std::wstring& file, const std::string& entryPoint);
    bool CompilePixelShader(const std::wstring& file, const std::string& entryPoint);
    // void Render(ID3D11ShaderResourceView* heightSRV);
    void Render(const std::vector<ID3D11ShaderResourceView*>& srvs);
    

private:
    static int NextPowerOf2(int n) {
        if (n <= 0) return 1;
        int p = 1;
        while (p < n) p <<= 1;
        return p;
    }

    // Compute Shaders
    int sizeX;   // simulation grid width
    int sizeY;   // simulation grid height
    int paddedSizeX;
    int paddedSizeY;
    UINT groupX; // number of thread groups in X to dispatch for compute shaders, assuming 16x16 threads per group
    UINT groupY; // number of thread groups in Y
    UINT paddedGroupX;
    UINT paddedGroupY;

    ID3D11Device* device = nullptr;
    ID3D11DeviceContext* context = nullptr;
    ID3D11Buffer* constantBuffer = nullptr;
    ID3D11Texture2D* stagingTex = nullptr; // Used to download data back to CPU
    ID3D11Texture2D* stagingComplexTex = nullptr;
    ID3D11Texture2D* stagingComplexArrayTex = nullptr;
    ID3D11Buffer* stagingBuf = nullptr;
    
    // FFT
    ID3D11Buffer* fftConstantBuffer = nullptr;
    ID3D11ComputeShader* fftShaderX = nullptr;
    ID3D11ComputeShader* fftShaderY = nullptr;
    ID3D11ComputeShader* fftArrayShaderX = nullptr;
    ID3D11ComputeShader* fftArrayShaderY = nullptr;

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