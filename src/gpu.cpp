#pragma comment(lib, "d3d11.lib")
#pragma comment(lib, "d3dcompiler.lib")

#include "gpu.h"

GPU::GPU() {}

GPU::~GPU() {
    if (stagingTex) stagingTex->Release();
    if (constantBuffer) constantBuffer->Release();
    if (fftConstantBuffer) fftConstantBuffer->Release();
    if (fftShader) fftShader->Release();
    if (fftArrayShader) fftArrayShader->Release();
    if (vertexShader) vertexShader->Release();
    if (pixelShader) pixelShader->Release();
    if (context) context->Release();
    if (device) device->Release();
}

bool GPU::BaseInit(int size) {
    /*** Compute Shaders ***/
    group = (size + 15) / 16; // size for shader dispatch
    // Create Compute Constant Buffer
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
        std::cerr << "ERROR: Failed to create Compute Constant Buffer." << std::endl;
        return false;
    }

    /*** FFT ***/
    // Create FFT Constant Buffer
    if (sizeof(FFTConstants) % 16 != 0) {
        std::cerr << "CRITICAL ERROR: FFTConstants is " << sizeof(FFTConstants) 
                  << " bytes. It MUST be a multiple of 16. Add padding!" << std::endl;
        return false;
    }
    D3D11_BUFFER_DESC fftCbDesc = {};
    fftCbDesc.ByteWidth = sizeof(FFTConstants);
    fftCbDesc.Usage = D3D11_USAGE_DYNAMIC;
    fftCbDesc.BindFlags = D3D11_BIND_CONSTANT_BUFFER;
    fftCbDesc.CPUAccessFlags = D3D11_CPU_ACCESS_WRITE;
    HRESULT hr_fft_cb = device->CreateBuffer(&fftCbDesc, nullptr, &fftConstantBuffer);
    if (FAILED(hr_fft_cb)) {
        std::cerr << "ERROR: Failed to create FFT Constant Buffer." << std::endl;
        return false;
    }
    // Compile FFT Shader
    if (!CompileFFTShaders(size)) {
        std::cerr << "ERROR: Failed to create FFT shaders." << std::endl;
        return false;
    }
}

bool GPU::Init(int size) {
    // Create a headless D3D11 device
    UINT createDeviceFlags = 0;
    D3D_FEATURE_LEVEL featureLevel;
    HRESULT hr = D3D11CreateDevice(nullptr, D3D_DRIVER_TYPE_HARDWARE, nullptr, createDeviceFlags,
                                   nullptr, 0, D3D11_SDK_VERSION, &device, &featureLevel, &context);
    if (FAILED(hr)) return false;

    BaseInit(size);    

    return true;
}

bool GPU::Init(int size, HWND hwnd) {
    /*** Create Device ***/
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

    /*** Rendering ***/
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

    // Define the Viewport (Must match window size in render.hlsl)
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
    sampDesc.Filter = D3D11_FILTER_MIN_MAG_MIP_POINT;
    sampDesc.AddressU = D3D11_TEXTURE_ADDRESS_CLAMP;
    sampDesc.AddressV = D3D11_TEXTURE_ADDRESS_CLAMP;
    sampDesc.AddressW = D3D11_TEXTURE_ADDRESS_CLAMP;
    sampDesc.ComparisonFunc = D3D11_COMPARISON_NEVER;
    device->CreateSamplerState(&sampDesc, &samplerState);

    /*** Everything Else ***/
    BaseInit(size);

    return true;
}


/************ COMPUTE SHADERS **************/
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

bool GPU::CreateGridTexture(GPUField* field, int size, bool isComplex, int arraySize) {
    D3D11_TEXTURE2D_DESC desc = {};
    desc.Width = size;
    desc.Height = size;
    desc.MipLevels = 1;
    desc.ArraySize = arraySize;
    desc.Format = isComplex ? DXGI_FORMAT_R32G32_FLOAT : DXGI_FORMAT_R32_FLOAT;
    desc.ArraySize = arraySize;
    desc.SampleDesc.Count = 1;
    desc.Usage = D3D11_USAGE_DEFAULT;
    desc.BindFlags = D3D11_BIND_UNORDERED_ACCESS | D3D11_BIND_SHADER_RESOURCE;
    if (FAILED(device->CreateTexture2D(&desc, nullptr, &field->tex))) return false;

    D3D11_SHADER_RESOURCE_VIEW_DESC srvDesc = {};
    srvDesc.Format = desc.Format;
    if (arraySize > 1) {
        srvDesc.ViewDimension = D3D11_SRV_DIMENSION_TEXTURE2DARRAY;
        srvDesc.Texture2DArray.MipLevels = 1;
        srvDesc.Texture2DArray.FirstArraySlice = 0;
        srvDesc.Texture2DArray.ArraySize = arraySize;
    } else {
        srvDesc.ViewDimension = D3D11_SRV_DIMENSION_TEXTURE2D;
        srvDesc.Texture2D.MipLevels = 1;
    }
    if (FAILED(device->CreateShaderResourceView(field->tex, &srvDesc, &field->srv))) return false;
    
    D3D11_UNORDERED_ACCESS_VIEW_DESC uavDesc = {};
    uavDesc.Format = desc.Format;
    if (arraySize > 1) {
        uavDesc.ViewDimension = D3D11_UAV_DIMENSION_TEXTURE2DARRAY;
        uavDesc.Texture2DArray.MipSlice = 0;
        uavDesc.Texture2DArray.FirstArraySlice = 0;
        uavDesc.Texture2DArray.ArraySize = arraySize;
    } else {
        uavDesc.ViewDimension = D3D11_UAV_DIMENSION_TEXTURE2D;
        uavDesc.Texture2D.MipSlice = 0;
    }
    if (FAILED(device->CreateUnorderedAccessView(field->tex, &uavDesc, &field->uav))) return false;

    return true;
}

bool GPU::UploadToGPU(ID3D11Texture2D* tex, const std::vector<float>& data, int gridSize, bool isComplex, int arraySize) {
    int floatsPerPixel = isComplex ? 2 : 1; // If complex, there are 2 floats per pixel.
    int rowPitch = gridSize * sizeof(float) * floatsPerPixel; // byte size of a single row
    int layerPitch = gridSize * rowPitch; // byte size of a single 2D layer

    // Upload each layer one by one
    for (int i = 0; i < arraySize; ++i) {
        // Find the starting memory address for this specific layer in the flat vector
        int floatOffset = i * (gridSize * gridSize * floatsPerPixel);
        const float* pLayerData = data.data() + floatOffset;

        // Upload to subresource 'i'
        context->UpdateSubresource(tex, i, nullptr, pLayerData, rowPitch, layerPitch);
    }
    
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

bool GPU::DownloadFromGPU(ID3D11Texture2D* tex, std::vector<std::complex<float>>& data, int size) {
    if (!stagingComplexTex) {
        D3D11_TEXTURE2D_DESC stDesc = {};
        tex->GetDesc(&stDesc);
        stDesc.Usage = D3D11_USAGE_STAGING;
        stDesc.BindFlags = 0;
        stDesc.CPUAccessFlags = D3D11_CPU_ACCESS_READ;
        device->CreateTexture2D(&stDesc, nullptr, &stagingComplexTex);
    }
    context->CopyResource(stagingComplexTex, tex);
    D3D11_MAPPED_SUBRESOURCE mapped;
    if (SUCCEEDED(context->Map(stagingComplexTex, 0, D3D11_MAP_READ, 0, &mapped))) {
        if (data.size() < size * size) data.resize(size * size);
        // Safely cast the raw VRAM bytes to the standard complex type
        std::complex<float>* pData = reinterpret_cast<std::complex<float>*>(mapped.pData);
        for (int y = 0; y < size; ++y) {
            memcpy(&data[y * size], 
                   pData + y * (mapped.RowPitch / sizeof(std::complex<float>)), 
                   size * sizeof(std::complex<float>));
        }
        context->Unmap(stagingComplexTex, 0);
        return true;
    }
    return false;
}

bool GPU::DownloadFromGPU(ID3D11Texture2D* tex, std::vector<std::complex<float>>& data, int size, int arraySize) {
    if (!stagingComplexArrayTex) {
        D3D11_TEXTURE2D_DESC stDesc = {};
        tex->GetDesc(&stDesc);
        stDesc.Usage = D3D11_USAGE_STAGING;
        stDesc.BindFlags = 0;
        stDesc.CPUAccessFlags = D3D11_CPU_ACCESS_READ;
        device->CreateTexture2D(&stDesc, nullptr, &stagingComplexArrayTex);
    }
    context->CopyResource(stagingComplexArrayTex, tex);
    int elementsPerLayer = size * size;
    int totalElements = elementsPerLayer * arraySize;
    if (data.size() < totalElements) {
        data.resize(totalElements);
    }
    for (int slice = 0; slice < arraySize; ++slice) {
        UINT subresource = D3D11CalcSubresource(0, slice, 1);
        D3D11_MAPPED_SUBRESOURCE mapped;
        if (SUCCEEDED(context->Map(stagingComplexArrayTex, subresource, D3D11_MAP_READ, 0, &mapped))) {
            std::complex<float>* pData = reinterpret_cast<std::complex<float>*>(mapped.pData);
            int sliceMemoryOffset = slice * elementsPerLayer;
            for (int y = 0; y < size; ++y) {
                memcpy(&data[sliceMemoryOffset + (y * size)], 
                       pData + y * (mapped.RowPitch / sizeof(std::complex<float>)), 
                       size * sizeof(std::complex<float>));
            }
            context->Unmap(stagingComplexArrayTex, subresource);
        } else {
            return false; // Failed to map a specific slice
        }
    }
    
    return true;
}

bool GPU::CreateBuffer(GPUBuffer* bufferObj, const void* data, UINT elementCount) {
    if (!bufferObj) return false;

    D3D11_BUFFER_DESC desc = {};
    desc.BindFlags = D3D11_BIND_UNORDERED_ACCESS | D3D11_BIND_SHADER_RESOURCE;
    desc.ByteWidth = UINT(elementCount * sizeof(float));
    desc.MiscFlags = D3D11_RESOURCE_MISC_BUFFER_STRUCTURED;
    desc.StructureByteStride = sizeof(float);
    desc.Usage = D3D11_USAGE_DEFAULT;
    D3D11_SUBRESOURCE_DATA initData = {};
    initData.pSysMem = data;
    HRESULT hr = device->CreateBuffer(&desc, &initData, &bufferObj->buf);
    if (FAILED(hr)) return false;

    D3D11_SHADER_RESOURCE_VIEW_DESC srvDesc = {};
    srvDesc.Format = DXGI_FORMAT_UNKNOWN; 
    srvDesc.ViewDimension = D3D11_SRV_DIMENSION_BUFFER;
    srvDesc.Buffer.FirstElement = 0;
    srvDesc.Buffer.NumElements = elementCount;
    hr = device->CreateShaderResourceView(bufferObj->buf, &srvDesc, &bufferObj->srv);
    if (FAILED(hr)) return false;
    
    D3D11_UNORDERED_ACCESS_VIEW_DESC uavDesc = {};
    uavDesc.Format = DXGI_FORMAT_UNKNOWN;
    uavDesc.ViewDimension = D3D11_UAV_DIMENSION_BUFFER;
    uavDesc.Buffer.FirstElement = 0;
    uavDesc.Buffer.NumElements = elementCount;
    hr = device->CreateUnorderedAccessView(bufferObj->buf, &uavDesc, &bufferObj->uav);
    if (FAILED(hr)) return false;

    return true; 
}

bool GPU::DownloadBuffer(ID3D11Buffer* buf, std::vector<float>& data, int size) {
    if (!stagingBuf) {
        D3D11_BUFFER_DESC stDesc = {};
        buf->GetDesc(&stDesc);
        stDesc.Usage = D3D11_USAGE_STAGING;
        stDesc.BindFlags = 0;
        stDesc.MiscFlags = D3D11_RESOURCE_MISC_BUFFER_STRUCTURED;
        stDesc.CPUAccessFlags = D3D11_CPU_ACCESS_READ;
        stDesc.ByteWidth = UINT(size * sizeof(float));
        stDesc.StructureByteStride = sizeof(float);
        device->CreateBuffer(&stDesc, nullptr, &stagingBuf);
    }
    // Copy from VRAM to System RAM accessible buffer
    context->CopyResource(stagingBuf, buf);

    D3D11_MAPPED_SUBRESOURCE mapped;
    if (SUCCEEDED(context->Map(stagingBuf, 0, D3D11_MAP_READ, 0, &mapped))) {
        if (data.size() < size) {
            data.resize(size);
        }
        float* pData = reinterpret_cast<float*>(mapped.pData); // (float*)
        memcpy(data.data(), pData, size * sizeof(float));
        context->Unmap(stagingBuf, 0);
        return true;
    }

    return false;
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

void GPU::BindSRV(int slot, ID3D11ShaderResourceView* srv) {
    context->CSSetShaderResources(slot, 1, &srv);
}

void GPU::Dispatch(ID3D11ComputeShader* shader, 
                   const std::vector<ID3D11ShaderResourceView*>& srvs, 
                   const std::vector<ID3D11UnorderedAccessView*>& uavs, 
                   int size) {
    context->CSSetShader(shader, nullptr, 0);
    context->CSSetShaderResources(0, srvs.size(), srvs.data());
    context->CSSetUnorderedAccessViews(0, uavs.size(), uavs.data(), nullptr);
    context->Dispatch(group, group, size);

    // Unbind to prevent read/write hazards
    std::vector<ID3D11ShaderResourceView*> nullSRVs(srvs.size(), nullptr);
    std::vector<ID3D11UnorderedAccessView*> nullUAVs(uavs.size(), nullptr);
    context->CSSetShaderResources(0, nullSRVs.size(), nullSRVs.data());
    context->CSSetUnorderedAccessViews(0, nullUAVs.size(), nullUAVs.data(), nullptr);
}


/************ FFT **************/
void GPU::UpdateFFTConstants(const FFTConstants& constants) {
    if (!context || !fftConstantBuffer) {
        std::cerr << "ERROR: GPU Context or FFT Constant Buffer is null!" << std::endl;
        return;
    }
    D3D11_MAPPED_SUBRESOURCE mapped;
    HRESULT hr = context->Map(fftConstantBuffer, 0, D3D11_MAP_WRITE_DISCARD, 0, &mapped);
    if (SUCCEEDED(hr)) {
        memcpy(mapped.pData, &constants, sizeof(FFTConstants));
        context->Unmap(fftConstantBuffer, 0);
        context->CSSetConstantBuffers(1, 1, &fftConstantBuffer); // Use slot 1 for FFT constants
    } else {
        std::cerr << "ERROR: Failed to map FFT Constant Buffer! HRESULT: " << hr << std::endl;
    }
}

bool GPU::CompileFFTShaders(int size) {
    // if (!CompileComputeShader(L"shaders/fft.hlsl", "FFTKernel_1D", &fftShader))
    //     return false;

    ID3DBlob* shaderBlob = nullptr;
    ID3DBlob* errorBlob = nullptr;
    D3D_SHADER_MACRO singleMacros[] = {
        { "FFT_SIZE", std::to_string(size).c_str() },
        { NULL, NULL }};
    HRESULT hr = D3DCompileFromFile(L"shaders/fft.hlsl", singleMacros, nullptr, "FFTKernel_1D", "cs_5_0", 0, 0, &shaderBlob, &errorBlob);
    if (FAILED(hr)) {
        if (errorBlob) std::cerr << "Shader Error: " << (char*)errorBlob->GetBufferPointer() << std::endl;
        return false;
    }
    device->CreateComputeShader(shaderBlob->GetBufferPointer(), shaderBlob->GetBufferSize(), nullptr, &fftShader);
    shaderBlob->Release();

    D3D_SHADER_MACRO arrayMacros[] = {
    { "FFT_SIZE", std::to_string(size).c_str() },
    { "IS_ARRAY", "1" },
    { NULL, NULL }};
    HRESULT hr1 = D3DCompileFromFile(L"shaders/fft.hlsl", arrayMacros, nullptr, "FFTKernel_1D", "cs_5_0", 0, 0, &shaderBlob, &errorBlob);
    if (FAILED(hr1)) {
        if (errorBlob) std::cerr << "Shader Error: " << (char*)errorBlob->GetBufferPointer() << std::endl;
        return false;
    }
    device->CreateComputeShader(shaderBlob->GetBufferPointer(), shaderBlob->GetBufferSize(), nullptr, &fftArrayShader);
    shaderBlob->Release();

    return true;
}

void GPU::ExecuteFFT(ID3D11UnorderedAccessView* fftBufferUAV, int size, bool inverse, int numLayers) {
    // Setup
    if (numLayers > 1)
        context->CSSetShader(fftArrayShader, nullptr, 0);
    else
        context->CSSetShader(fftShader, nullptr, 0);
    context->CSSetUnorderedAccessViews(0, 1, &fftBufferUAV, nullptr);

    FFTConstants constants = {};
    constants.N = size;
    constants.Inverse = inverse ? 1 : 0;
    unsigned int bits = 0; // integer Log2
    int temp = size;
    while (temp >>= 1) ++bits;
    constants.Bits = bits;

    // Row pass
    constants.Row = true;
    UpdateFFTConstants(constants);
    context->Dispatch(1, size, numLayers); // One group per row

    // Column pass
    constants.Row = false;
    UpdateFFTConstants(constants);
    context->Dispatch(1, size, numLayers); // One group per column

    // Unbind
    ID3D11UnorderedAccessView* nullUAV = nullptr;
    context->CSSetUnorderedAccessViews(0, 1, &nullUAV, nullptr);
}


/************ RENDERING **************/
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
