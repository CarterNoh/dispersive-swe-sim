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
    float vuffer[2];                  // 8 bytes (Total = 80 bytes, a multiple of 16!)
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

    bool Init(); // For headless compute
    bool Init(HWND hwnd);
    
    // Memory Management
    bool CreateGridTexture(GPUField* field, int size);
    bool UploadToGPU(ID3D11Texture2D* tex, const std::vector<float>& data, int size);
    bool DownloadFromGPU(ID3D11Texture2D* tex, std::vector<float>& data, int size);
    void ClearUAV(ID3D11UnorderedAccessView* uav, float clearValue);
    bool CreateGridVertexBuffer(int gridSize);

    // Shader Management
    bool CompileComputeShader(const std::wstring& file, const std::string& entryPoint, ID3D11ComputeShader** shader);
    bool CompileVertexShader(const std::wstring& file, const std::string& entryPoint);
    bool CompilePixelShader(const std::wstring& file, const std::string& entryPoint);
    void UpdateConstants(const SimConstants& constants);
    void UpdateRenderConstants(const RenderConstants& constants);
    
    // Execution
    void Dispatch(ID3D11ComputeShader* shader, 
                  const std::vector<ID3D11ShaderResourceView*>& srvs, 
                  const std::vector<ID3D11UnorderedAccessView*>& uavs, 
                  int groupsX, int groupsY);
    bool CreateGridMesh(int gridSize);
    void Render(ID3D11ShaderResourceView* heightSRV); // , int gridVertexCount

private:
    ID3D11Device* device = nullptr;
    ID3D11DeviceContext* context = nullptr;
    ID3D11Buffer* constantBuffer = nullptr;
    ID3D11Texture2D* stagingTex = nullptr; // Used to download data back to CPU

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