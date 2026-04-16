#pragma comment(lib, "d3d11.lib")
#pragma comment(lib, "d3dcompiler.lib")

#include "gpu.h"

GPU::GPU() {}

GPU::~GPU() {
    if (stagingTex) stagingTex->Release();
    if (constantBuffer) constantBuffer->Release();
    if (context) context->Release();
    if (device) device->Release();
}

bool GPU::Init() {
    // Create a headless D3D11 device
    UINT createDeviceFlags = 0;
#ifdef _DEBUG
    createDeviceFlags |= D3D11_CREATE_DEVICE_DEBUG;
#endif

    D3D_FEATURE_LEVEL featureLevel;
    HRESULT hr = D3D11CreateDevice(nullptr, D3D_DRIVER_TYPE_HARDWARE, nullptr, createDeviceFlags,
                                   nullptr, 0, D3D11_SDK_VERSION, &device, &featureLevel, &context);
    if (FAILED(hr)) return false;

    // Create the Constant Buffer
    if (sizeof(SimConstants) % 16 != 0) {
        std::cerr << "CRITICAL ERROR: SimConstants is " << sizeof(SimConstants) 
                  << " bytes. It MUST be a multiple of 16. Add padding!" << std::endl;
        return false;
    }
    D3D11_BUFFER_DESC cbDesc = {};
    cbDesc.ByteWidth = sizeof(SimConstants);
    cbDesc.Usage = D3D11_USAGE_DYNAMIC;
    cbDesc.BindFlags = D3D11_BIND_CONSTANT_BUFFER;
    cbDesc.CPUAccessFlags = D3D11_CPU_ACCESS_WRITE;
    device->CreateBuffer(&cbDesc, nullptr, &constantBuffer);
    HRESULT hr_cb = device->CreateBuffer(&cbDesc, nullptr, &constantBuffer);
    if (FAILED(hr_cb)) {
        std::cerr << "ERROR: Failed to create Constant Buffer." << std::endl;
        return false;
    }

    return true;
}

bool GPU::CreateGridTexture(ID3D11Texture2D** tex, ID3D11UnorderedAccessView** uav, int size) {
    D3D11_TEXTURE2D_DESC desc = {};
    desc.Width = size;
    desc.Height = size;
    desc.MipLevels = 1;
    desc.ArraySize = 1;
    desc.Format = DXGI_FORMAT_R32_FLOAT; // Standard 32-bit float
    desc.SampleDesc.Count = 1;
    desc.Usage = D3D11_USAGE_DEFAULT;
    desc.BindFlags = D3D11_BIND_UNORDERED_ACCESS | D3D11_BIND_SHADER_RESOURCE;

    if (FAILED(device->CreateTexture2D(&desc, nullptr, tex))) return false;
    
    D3D11_UNORDERED_ACCESS_VIEW_DESC uavDesc = {};
    uavDesc.Format = desc.Format;
    uavDesc.ViewDimension = D3D11_UAV_DIMENSION_TEXTURE2D;
    if (FAILED(device->CreateUnorderedAccessView(*tex, &uavDesc, uav))) return false;

    return true;
}

bool GPU::UploadToGPU(ID3D11Texture2D* tex, const std::vector<float>& data, int size) {
    context->UpdateSubresource(tex, 0, nullptr, data.data(), size * sizeof(float), 0);
    return true;
}

bool GPU::DownloadFromGPU(ID3D11Texture2D* tex, std::vector<float>& data, int size) {
    // Lazy initialize staging texture
    if (!stagingTex) {
        D3D11_TEXTURE2D_DESC stDesc = {};
        tex->GetDesc(&stDesc);
        stDesc.Usage = D3D11_USAGE_STAGING;
        stDesc.BindFlags = 0;
        stDesc.CPUAccessFlags = D3D11_CPU_ACCESS_READ;
        device->CreateTexture2D(&stDesc, nullptr, &stagingTex);
    }

    // Copy from VRAM to System RAM accessible texture
    context->CopyResource(stagingTex, tex);

    D3D11_MAPPED_SUBRESOURCE mapped;
    if (SUCCEEDED(context->Map(stagingTex, 0, D3D11_MAP_READ, 0, &mapped))) {
        float* pData = reinterpret_cast<float*>(mapped.pData);
        for (int y = 0; y < size; ++y) {
            memcpy(&data[y * size], pData + y * (mapped.RowPitch / sizeof(float)), size * sizeof(float));
        }
        context->Unmap(stagingTex, 0);
        return true;
    }
    return false;
}

void GPU::ClearUAV(ID3D11UnorderedAccessView* uav, float clearValue) {
    float clearColor[4] = { clearValue, clearValue, clearValue, clearValue };
    context->ClearUnorderedAccessViewFloat(uav, clearColor);
}

bool GPU::CompileShader(const std::wstring& file, const std::string& entryPoint, ID3D11ComputeShader** shader) {
    ID3DBlob* shaderBlob = nullptr;
    ID3DBlob* errorBlob = nullptr;

    HRESULT hr = D3DCompileFromFile(file.c_str(), nullptr, nullptr, entryPoint.c_str(), "cs_5_0", 0, 0, &shaderBlob, &errorBlob);
    if (FAILED(hr)) {
        if (errorBlob) std::cerr << "Shader Error: " << (char*)errorBlob->GetBufferPointer() << std::endl;
        return false;
    }

    device->CreateComputeShader(shaderBlob->GetBufferPointer(), shaderBlob->GetBufferSize(), nullptr, shader);
    shaderBlob->Release();
    return true;
}

void GPU::UpdateConstants(const SimConstants& constants) {
    if (!context || !constantBuffer) {
        std::cerr << "ERROR: GPU Context or Constant Buffer is null!" << std::endl;
        return;
    }
    D3D11_MAPPED_SUBRESOURCE mapped;
    HRESULT hr = context->Map(constantBuffer, 0, D3D11_MAP_WRITE_DISCARD, 0, &mapped);
    if (SUCCEEDED(hr)) {
        memcpy(mapped.pData, &constants, sizeof(SimConstants));
        context->Unmap(constantBuffer, 0);
        context->CSSetConstantBuffers(0, 1, &constantBuffer);
    } else {
        std::cerr << "ERROR: Failed to map Constant Buffer! HRESULT: " << hr << std::endl;
    }
}

void GPU::Dispatch(ID3D11ComputeShader* shader, const std::vector<ID3D11UnorderedAccessView*>& uavs, int groupsX, int groupsY) {
    context->CSSetShader(shader, nullptr, 0);
    context->CSSetUnorderedAccessViews(0, uavs.size(), uavs.data(), nullptr);
    
    context->Dispatch(groupsX, groupsY, 1);

    // Unbind to prevent read/write hazards
    std::vector<ID3D11UnorderedAccessView*> nullUAVs(uavs.size(), nullptr);
    context->CSSetUnorderedAccessViews(0, nullUAVs.size(), nullUAVs.data(), nullptr);
}

