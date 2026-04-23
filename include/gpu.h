#pragma once
#include <d3d11.h>
#include <d3dcompiler.h>
#include <vector>
#include <string>
#include <iostream>

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
    float buffer1;
    float buffer2;
    
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

    bool Init();
    
    // Memory Management
    bool CreateGridTexture(GPUField* field, int size);
    bool UploadToGPU(ID3D11Texture2D* tex, const std::vector<float>& data, int size);
    bool DownloadFromGPU(ID3D11Texture2D* tex, std::vector<float>& data, int size);
    void ClearUAV(ID3D11UnorderedAccessView* uav, float clearValue);
    void Flush();

    // Shader Management
    bool CompileShader(const std::wstring& file, const std::string& entryPoint, ID3D11ComputeShader** shader);
    void UpdateConstants(const SimConstants& constants);
    
    // Execution
    void Dispatch(ID3D11ComputeShader* shader, 
                  const std::vector<ID3D11ShaderResourceView*>& srvs, 
                  const std::vector<ID3D11UnorderedAccessView*>& uavs, 
                  int groupsX, int groupsY);

private:
    ID3D11Device* device = nullptr;
    ID3D11DeviceContext* context = nullptr;
    ID3D11Buffer* constantBuffer = nullptr;
    ID3D11Texture2D* stagingTex = nullptr; // Used to download data back to CPU
};