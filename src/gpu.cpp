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

bool GPU::Init(HWND hwnd) {
    // Create a headless D3D11 device
    UINT createDeviceFlags = 0;
    DXGI_SWAP_CHAIN_DESC sd = {};
    sd.BufferCount = 1;
    sd.BufferDesc.Width = 800;  // Window Width
    sd.BufferDesc.Height = 600; // Window Height
    sd.BufferDesc.Format = DXGI_FORMAT_R8G8B8A8_UNORM;
    sd.BufferUsage = DXGI_USAGE_RENDER_TARGET_OUTPUT;
    sd.OutputWindow = hwnd;
    sd.SampleDesc.Count = 1;
    sd.Windowed = TRUE;

    D3D_FEATURE_LEVEL featureLevel;
    HRESULT hr = D3D11CreateDeviceAndSwapChain(
        nullptr, D3D_DRIVER_TYPE_HARDWARE, nullptr, createDeviceFlags, nullptr, 0, 
        D3D11_SDK_VERSION, &sd, &swapChain, &device, &featureLevel, &context);
    if (FAILED(hr)) return false;
    // Get the back buffer from the swap chain to use as our Render Target
    ID3D11Texture2D* backBuffer = nullptr;
    swapChain->GetBuffer(0, __uuidof(ID3D11Texture2D), (void**)&backBuffer);
    device->CreateRenderTargetView(backBuffer, nullptr, &renderTargetView);
    backBuffer->Release();

    // Create the Depth/Stencil Buffer
    D3D11_TEXTURE2D_DESC descDepth = {};
    descDepth.Width = 800;  // Must match Window Width!
    descDepth.Height = 600; // Must match Window Height!
    descDepth.MipLevels = 1;
    descDepth.ArraySize = 1;
    descDepth.Format = DXGI_FORMAT_D24_UNORM_S8_UINT; // 24 bits for Depth, 8 for Stencil
    descDepth.SampleDesc.Count = 1;
    descDepth.Usage = D3D11_USAGE_DEFAULT;
    descDepth.BindFlags = D3D11_BIND_DEPTH_STENCIL;
    device->CreateTexture2D(&descDepth, nullptr, &depthStencilBuffer);
    device->CreateDepthStencilView(depthStencilBuffer, nullptr, &depthStencilView);
    // Bind BOTH the Render Target (Color) and the Depth View (Z-Buffer)
    // Overwrite your old OMSetRenderTargets line with this one:
    context->OMSetRenderTargets(1, &renderTargetView, depthStencilView);
    // Create the Rasterizer State to disable Backface Culling
    D3D11_RASTERIZER_DESC rsDesc = {};
    rsDesc.FillMode = D3D11_FILL_SOLID; // Change to D3D11_FILL_WIREFRAME to see the grid lines!
    rsDesc.CullMode = D3D11_CULL_NONE;  // This is the magic line that fixes invisibility
    device->CreateRasterizerState(&rsDesc, &rasterState);
    // Bind it to the pipeline
    context->RSSetState(rasterState);

    // Create the Sim Constant Buffer
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
    HRESULT hr_cb = device->CreateBuffer(&cbDesc, nullptr, &constantBuffer);
    if (FAILED(hr_cb)) {
        std::cerr << "ERROR: Failed to create Constant Buffer." << std::endl;
        return false;
    }

    // Create the Render Constant Buffer
    if (sizeof(RenderConstants) % 16 != 0) {
        std::cerr << "CRITICAL ERROR: RenderConstants is " << sizeof(RenderConstants) 
                  << " bytes. It MUST be a multiple of 16. Add padding!" << std::endl;
        return false;
    }
    D3D11_BUFFER_DESC rbDesc = {};
    rbDesc.ByteWidth = sizeof(RenderConstants);
    rbDesc.Usage = D3D11_USAGE_DYNAMIC;
    rbDesc.BindFlags = D3D11_BIND_CONSTANT_BUFFER;
    rbDesc.CPUAccessFlags = D3D11_CPU_ACCESS_WRITE;
    HRESULT hr_rb = device->CreateBuffer(&rbDesc, nullptr, &renderConstantBuffer);
    if (FAILED(hr_rb)) {
        std::cerr << "ERROR: Failed to create Render Constant Buffer." << std::endl;
        return false;
    }

    // Define the Viewport (Must match your Window Size!)
    D3D11_VIEWPORT vp = {};
    vp.Width = 800.0f;
    vp.Height = 600.0f;
    vp.MinDepth = 0.0f;
    vp.MaxDepth = 1.0f;
    vp.TopLeftX = 0;
    vp.TopLeftY = 0;
    context->RSSetViewports(1, &vp);

    // Create the Texture Sampler
    D3D11_SAMPLER_DESC sampDesc = {};
    sampDesc.Filter = D3D11_FILTER_MIN_MAG_MIP_POINT; // Don't blur the data!
    sampDesc.AddressU = D3D11_TEXTURE_ADDRESS_CLAMP;
    sampDesc.AddressV = D3D11_TEXTURE_ADDRESS_CLAMP;
    sampDesc.AddressW = D3D11_TEXTURE_ADDRESS_CLAMP;
    sampDesc.ComparisonFunc = D3D11_COMPARISON_NEVER;
    device->CreateSamplerState(&sampDesc, &samplerState);

    return true;
}

bool GPU::CreateGridTexture(GPUField* field, int size) {
    D3D11_TEXTURE2D_DESC desc = {};
    desc.Width = size;
    desc.Height = size;
    desc.MipLevels = 1;
    desc.ArraySize = 1;
    desc.Format = DXGI_FORMAT_R32_FLOAT; // Standard 32-bit float
    desc.SampleDesc.Count = 1;
    desc.Usage = D3D11_USAGE_DEFAULT;
    desc.BindFlags = D3D11_BIND_UNORDERED_ACCESS | D3D11_BIND_SHADER_RESOURCE;
    if (FAILED(device->CreateTexture2D(&desc, nullptr, &field->tex))) return false;

    D3D11_SHADER_RESOURCE_VIEW_DESC srvDesc = {};
    srvDesc.Format = desc.Format;
    srvDesc.ViewDimension = D3D11_SRV_DIMENSION_TEXTURE2D;
    srvDesc.Texture2D.MipLevels = 1;
    if (FAILED(device->CreateShaderResourceView(field->tex, &srvDesc, &field->srv))) return false;
    
    D3D11_UNORDERED_ACCESS_VIEW_DESC uavDesc = {};
    uavDesc.Format = desc.Format;
    uavDesc.ViewDimension = D3D11_UAV_DIMENSION_TEXTURE2D;
    if (FAILED(device->CreateUnorderedAccessView(field->tex, &uavDesc, &field->uav))) return false;

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

void GPU::Flush() {
    D3D11_QUERY_DESC queryDesc = {};
    queryDesc.Query = D3D11_QUERY_EVENT;
    ID3D11Query* query = nullptr;
    device->CreateQuery(&queryDesc, &query);
    context->End(query);
    BOOL done = FALSE;
    while (context->GetData(query, &done, sizeof(BOOL), 0) == S_FALSE || !done) {}
    query->Release();
}

bool GPU::CreateGridVertexBuffer(int gridSize) {
    std::vector<float> vertices;
    // Generate a flat grid of X, Y coordinates
    for (int y = 0; y < gridSize; ++y) {
        for (int x = 0; x < gridSize; ++x) {
            vertices.push_back((float)x);
            vertices.push_back((float)y);
        }
    }

    // Send this list to the GPU as an immutable Vertex Buffer
    D3D11_BUFFER_DESC bd = {};
    bd.Usage = D3D11_USAGE_IMMUTABLE;
    bd.ByteWidth = sizeof(float) * vertices.size();
    bd.BindFlags = D3D11_BIND_VERTEX_BUFFER;
    
    D3D11_SUBRESOURCE_DATA initData = {};
    initData.pSysMem = vertices.data();
    
    HRESULT hr = device->CreateBuffer(&bd, &initData, &vertexBuffer);
    return SUCCEEDED(hr);
}

bool GPU::CompileComputeShader(const std::wstring& file, const std::string& entryPoint, ID3D11ComputeShader** shader) {
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

bool GPU::CompileVertexShader(const std::wstring& file, const std::string& entryPoint) {
    ID3DBlob* shaderBlob = nullptr;
    ID3DBlob* errorBlob = nullptr;

    HRESULT hr = D3DCompileFromFile(file.c_str(), nullptr, nullptr, entryPoint.c_str(), "vs_5_0", 0, 0, &shaderBlob, &errorBlob);
    if (FAILED(hr)) {
        if (errorBlob) std::cerr << "VS Compile Error: " << (char*)errorBlob->GetBufferPointer() << std::endl;
        return false;
    }

    device->CreateVertexShader(shaderBlob->GetBufferPointer(), shaderBlob->GetBufferSize(), nullptr, &vertexShader);

    D3D11_INPUT_ELEMENT_DESC layoutDesc[] = {
        { "POSITION", 0, DXGI_FORMAT_R32G32_FLOAT, 0, 0, D3D11_INPUT_PER_VERTEX_DATA, 0 }
    };
    device->CreateInputLayout(layoutDesc, 1, shaderBlob->GetBufferPointer(), shaderBlob->GetBufferSize(), &inputLayout);
    
    shaderBlob->Release();
    return true;
}

bool GPU::CompilePixelShader(const std::wstring& file, const std::string& entryPoint) {
    ID3DBlob* shaderBlob = nullptr;
    ID3DBlob* errorBlob = nullptr;

    HRESULT hr = D3DCompileFromFile(file.c_str(), nullptr, nullptr, entryPoint.c_str(), "ps_5_0", 0, 0, &shaderBlob, &errorBlob);
    if (FAILED(hr)) {
        if (errorBlob) std::cerr << "PS Compile Error: " << (char*)errorBlob->GetBufferPointer() << std::endl;
        return false;
    }

    device->CreatePixelShader(shaderBlob->GetBufferPointer(), shaderBlob->GetBufferSize(), nullptr, &pixelShader);
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

void GPU::UpdateRenderConstants(const RenderConstants& constants) {
    if (!context || !renderConstantBuffer) {
        std::cerr << "CRITICAL: Render Constant Buffer is null! Did you create it in Init()?" << std::endl;
        return;
    }
    D3D11_MAPPED_SUBRESOURCE mapped;
    if (SUCCEEDED(context->Map(renderConstantBuffer, 0, D3D11_MAP_WRITE_DISCARD, 0, &mapped))) {
        memcpy(mapped.pData, &constants, sizeof(RenderConstants));
        context->Unmap(renderConstantBuffer, 0);
        
        context->VSSetConstantBuffers(0, 1, &renderConstantBuffer);
        context->PSSetConstantBuffers(0, 1, &renderConstantBuffer);
    }
}

void GPU::Dispatch(ID3D11ComputeShader* shader, 
                   const std::vector<ID3D11ShaderResourceView*>& srvs, 
                   const std::vector<ID3D11UnorderedAccessView*>& uavs, 
                   int groupsX, int groupsY) {
    context->CSSetShader(shader, nullptr, 0);
    context->CSSetShaderResources(0, srvs.size(), srvs.data());
    context->CSSetUnorderedAccessViews(0, uavs.size(), uavs.data(), nullptr);
    context->Dispatch(groupsX, groupsY, 1);

    // Unbind to prevent read/write hazards
    std::vector<ID3D11ShaderResourceView*> nullSRVs(srvs.size(), nullptr);
    std::vector<ID3D11UnorderedAccessView*> nullUAVs(uavs.size(), nullptr);
    context->CSSetShaderResources(0, nullSRVs.size(), nullSRVs.data());
    context->CSSetUnorderedAccessViews(0, nullUAVs.size(), nullUAVs.data(), nullptr);
}

bool GPU::CreateGridMesh(int gridSize) {
    // 1. Generate the Vertices (Same as before)
    std::vector<float> vertices;
    vertices.reserve(gridSize * gridSize * 2);
    for (int y = 0; y < gridSize; ++y) {
        for (int x = 0; x < gridSize; ++x) {
            vertices.push_back((float)x);
            vertices.push_back((float)y);
        }
    }

    // Create the Vertex Buffer
    D3D11_BUFFER_DESC vbd = {};
    vbd.Usage = D3D11_USAGE_IMMUTABLE;
    vbd.ByteWidth = sizeof(float) * vertices.size();
    vbd.BindFlags = D3D11_BIND_VERTEX_BUFFER;
    D3D11_SUBRESOURCE_DATA vInit = {};
    vInit.pSysMem = vertices.data();
    device->CreateBuffer(&vbd, &vInit, &vertexBuffer);

    // 2. Generate the Indices (The "Map" of Triangles)
    std::vector<unsigned int> indices;
    indices.reserve(gridSize * gridSize * 6);
    // We loop to gridSize - 1 because we are building squares BETWEEN the points
    for (int y = 0; y < gridSize - 1; ++y) {
        for (int x = 0; x < gridSize - 1; ++x) {
            // Find the four corners of the current square (quad)
            unsigned int topLeft     = y * gridSize + x;
            unsigned int topRight    = y * gridSize + (x + 1);
            unsigned int bottomLeft  = (y + 1) * gridSize + x;
            unsigned int bottomRight = (y + 1) * gridSize + (x + 1);

            // Triangle 1 (Top-Left, Top-Right, Bottom-Left)
            indices.push_back(topLeft);
            indices.push_back(topRight);
            indices.push_back(bottomLeft);

            // Triangle 2 (Bottom-Left, Top-Right, Bottom-Right)
            indices.push_back(bottomLeft);
            indices.push_back(topRight);
            indices.push_back(bottomRight);
        }
    }

    indexCount = indices.size();

    // Create the Index Buffer
    D3D11_BUFFER_DESC ibd = {};
    ibd.Usage = D3D11_USAGE_IMMUTABLE;
    ibd.ByteWidth = sizeof(unsigned int) * indices.size();
    ibd.BindFlags = D3D11_BIND_INDEX_BUFFER;
    D3D11_SUBRESOURCE_DATA iInit = {};
    iInit.pSysMem = indices.data();
    device->CreateBuffer(&ibd, &iInit, &indexBuffer);

    return true;
}

void GPU::Render(ID3D11ShaderResourceView* heightSRV) { //, int gridVertexCount
    // Clear the screen
    float clearColor[4] = { 1.0f, 1.0f, 1.0f, 1.0f };
    context->ClearRenderTargetView(renderTargetView, clearColor);
    context->ClearDepthStencilView(depthStencilView, D3D11_CLEAR_DEPTH, 1.0f, 0);

    // Bind Graphics Pipeline
    context->OMSetRenderTargets(1, &renderTargetView, depthStencilView);
    context->VSSetShader(vertexShader, nullptr, 0);
    context->PSSetShader(pixelShader, nullptr, 0);

    // Bind the Fluid Data! (As an SRV, NOT a UAV)
    context->VSSetShaderResources(0, 1, &heightSRV);
    context->VSSetSamplers(0, 1, &samplerState);

    // Draw 
    context->IASetInputLayout(inputLayout);
    UINT stride = sizeof(float) * 2; 
    UINT offset = 0;
    context->IASetVertexBuffers(0, 1, &vertexBuffer, &stride, &offset);
    // context->IASetPrimitiveTopology(D3D11_PRIMITIVE_TOPOLOGY_POINTLIST); // Draw points for now!
    // context->Draw(gridVertexCount, 0);
    context->IASetIndexBuffer(indexBuffer, DXGI_FORMAT_R32_UINT, 0);
    context->IASetPrimitiveTopology(D3D11_PRIMITIVE_TOPOLOGY_TRIANGLELIST);
    context->DrawIndexed(indexCount, 0, 0); // 4. Draw using the Indices, not just raw vertices

    // Unbind the SRV
    ID3D11ShaderResourceView* nullSRV = nullptr;
    context->VSSetShaderResources(0, 1, &nullSRV);

    // Present the frame to the Window
    swapChain->Present(0, 0); // (1, 0) = VSync on
}