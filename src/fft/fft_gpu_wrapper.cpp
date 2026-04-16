#include "fft/fft_gpu_wrapper.h"
#include <stdexcept>
#include <fstream>
#include <wrl.h>
#include <d3dx12.h>

using Microsoft::WRL::ComPtr;

GPUFFTCompute::GPUFFTCompute(uint32_t w, uint32_t h)
    : width(w), height(h), fenceValue(0)
{
    fenceEvent = CreateEvent(nullptr, FALSE, FALSE, nullptr);
    if (fenceEvent == nullptr) {
        throw std::runtime_error("Failed to create fence event");
    }
}

GPUFFTCompute::~GPUFFTCompute()
{
    if (fenceEvent != nullptr) {
        CloseHandle(fenceEvent);
    }
}

bool GPUFFTCompute::Initialize()
{
    try {
        if (!CreateDevice()) return false;
        if (!CreateBuffers()) return false;
        if (!CreateRootSignature()) return false;
        if (!LoadShader()) return false;
        return true;
    }
    catch (const std::exception& e) {
        OutputDebugStringA(("GPUFFT Init Error: " + std::string(e.what()) + "\n").c_str());
        return false;
    }
}

bool GPUFFTCompute::CreateDevice()
{
    // Enable debug layer in debug builds
#ifdef _DEBUG
    ComPtr<ID3D12Debug> debugController;
    if (SUCCEEDED(D3D12GetDebugInterface(IID_PPV_ARGS(&debugController)))) {
        debugController->EnableDebugLayer();
    }
#endif

    // Create DXGI factory
    ComPtr<IDXGIFactory6> factory;
    if (FAILED(CreateDXGIFactory1(IID_PPV_ARGS(&factory)))) {
        OutputDebugStringA("Failed to create DXGI factory\n");
        return false;
    }

    // Enumerate adapters and find discrete GPU
    ComPtr<IDXGIAdapter1> adapter;
    DXGI_ADAPTER_DESC1 desc;
    
    for (UINT i = 0; factory->EnumAdapters1(i, &adapter) != DXGI_ERROR_NOT_FOUND; ++i) {
        adapter->GetDesc1(&desc);
        // Prefer discrete GPU but fall back to any GPU
        if (desc.Flags != DXGI_ADAPTER_FLAG_SOFTWARE) {
            break;
        }
    }

    // Create device
    if (FAILED(D3D12CreateDevice(adapter.Get(), D3D_FEATURE_LEVEL_12_0, IID_PPV_ARGS(&device)))) {
        OutputDebugStringA("Failed to create D3D12 device\n");
        return false;
    }

    // Create command queue
    D3D12_COMMAND_QUEUE_DESC queueDesc = {};
    queueDesc.Type = D3D12_COMMAND_LIST_TYPE_COMPUTE;
    
    if (FAILED(device->CreateCommandQueue(&queueDesc, IID_PPV_ARGS(&commandQueue)))) {
        OutputDebugStringA("Failed to create command queue\n");
        return false;
    }

    // Create command allocator
    if (FAILED(device->CreateCommandAllocator(D3D12_COMMAND_LIST_TYPE_COMPUTE, IID_PPV_ARGS(&commandAllocator)))) {
        OutputDebugStringA("Failed to create command allocator\n");
        return false;
    }

    // Create command list
    if (FAILED(device->CreateCommandList1(0, D3D12_COMMAND_LIST_TYPE_COMPUTE, D3D12_COMMAND_LIST_FLAG_NONE, IID_PPV_ARGS(&commandList)))) {
        OutputDebugStringA("Failed to create command list\n");
        return false;
    }

    // Create fence
    if (FAILED(device->CreateFence(0, D3D12_FENCE_FLAG_NONE, IID_PPV_ARGS(&fence)))) {
        OutputDebugStringA("Failed to create fence\n");
        return false;
    }

    return true;
}

bool GPUFFTCompute::CreateBuffers()
{
    uint32_t elementCount = width * height;
    uint32_t bufferSize = elementCount * sizeof(float2);

    // Create GPU buffer (UAV)
    D3D12_RESOURCE_DESC bufferDesc = {};
    bufferDesc.Dimension = D3D12_RESOURCE_DIMENSION_BUFFER;
    bufferDesc.Alignment = 0;
    bufferDesc.Width = bufferSize;
    bufferDesc.Height = 1;
    bufferDesc.DepthOrArraySize = 1;
    bufferDesc.MipLevels = 1;
    bufferDesc.Format = DXGI_FORMAT_UNKNOWN;
    bufferDesc.SampleDesc.Count = 1;
    bufferDesc.SampleDesc.Quality = 0;
    bufferDesc.Layout = D3D12_TEXTURE_LAYOUT_ROW_MAJOR;
    bufferDesc.Flags = D3D12_RESOURCE_FLAG_ALLOW_UNORDERED_ACCESS;

    D3D12_HEAP_PROPERTIES heapProps = {};
    heapProps.Type = D3D12_HEAP_TYPE_DEFAULT;
    heapProps.CPUPageProperty = D3D12_CPU_PAGE_PROPERTY_UNKNOWN;
    heapProps.MemoryPoolPreference = D3D12_MEMORY_POOL_UNKNOWN;
    heapProps.CreationNodeMask = 1;
    heapProps.VisibleNodeMask = 1;

    if (FAILED(device->CreateCommittedResource(&heapProps, D3D12_HEAP_FLAG_NONE, &bufferDesc,
        D3D12_RESOURCE_STATE_UNORDERED_ACCESS, nullptr, IID_PPV_ARGS(&dataBuffer)))) {
        OutputDebugStringA("Failed to create data buffer\n");
        return false;
    }

    // Create upload buffer
    heapProps.Type = D3D12_HEAP_TYPE_UPLOAD;
    if (FAILED(device->CreateCommittedResource(&heapProps, D3D12_HEAP_FLAG_NONE, &bufferDesc,
        D3D12_RESOURCE_STATE_GENERIC_READ, nullptr, IID_PPV_ARGS(&uploadBuffer)))) {
        OutputDebugStringA("Failed to create upload buffer\n");
        return false;
    }

    // Create readback buffer
    heapProps.Type = D3D12_HEAP_TYPE_READBACK;
    bufferDesc.Flags = D3D12_RESOURCE_FLAG_NONE;
    if (FAILED(device->CreateCommittedResource(&heapProps, D3D12_HEAP_FLAG_NONE, &bufferDesc,
        D3D12_RESOURCE_STATE_COPY_DEST, nullptr, IID_PPV_ARGS(&downloadBuffer)))) {
        OutputDebugStringA("Failed to create readback buffer\n");
        return false;
    }

    // Create constant buffer
    bufferDesc.Width = sizeof(FFTConstants);
    bufferDesc.Flags = D3D12_RESOURCE_FLAG_NONE;
    heapProps.Type = D3D12_HEAP_TYPE_UPLOAD;
    if (FAILED(device->CreateCommittedResource(&heapProps, D3D12_HEAP_FLAG_NONE, &bufferDesc,
        D3D12_RESOURCE_STATE_GENERIC_READ, nullptr, IID_PPV_ARGS(&constantBuffer)))) {
        OutputDebugStringA("Failed to create constant buffer\n");
        return false;
    }

    // Create descriptor heap for UAV
    D3D12_DESCRIPTOR_HEAP_DESC heapDesc = {};
    heapDesc.NumDescriptors = 2;  // One for data, one for constants
    heapDesc.Type = D3D12_DESCRIPTOR_HEAP_TYPE_CBV_SRV_UAV;
    heapDesc.Flags = D3D12_DESCRIPTOR_HEAP_FLAG_SHADER_VISIBLE;

    if (FAILED(device->CreateDescriptorHeap(&heapDesc, IID_PPV_ARGS(&descriptorHeap)))) {
        OutputDebugStringA("Failed to create descriptor heap\n");
        return false;
    }

    // Create UAV descriptor
    D3D12_UNORDERED_ACCESS_VIEW_DESC uavDesc = {};
    uavDesc.Format = DXGI_FORMAT_UNKNOWN;
    uavDesc.ViewDimension = D3D12_UAV_DIMENSION_BUFFER;
    uavDesc.Buffer.FirstElement = 0;
    uavDesc.Buffer.NumElements = elementCount;
    uavDesc.Buffer.StructureByteStride = sizeof(float2);
    uavDesc.Buffer.CounterOffsetInBytes = 0;
    uavDesc.Buffer.Flags = D3D12_BUFFER_UAV_FLAG_NONE;

    CD3DX12_CPU_DESCRIPTOR_HANDLE cpuHandle(descriptorHeap->GetCPUDescriptorHandleForHeapStart());
    device->CreateUnorderedAccessView(dataBuffer.Get(), nullptr, &uavDesc, cpuHandle);

    return true;
}

bool GPUFFTCompute::CreateRootSignature()
{
    // Define root parameters: constant buffer and UAV
    CD3DX12_ROOT_PARAMETER rootParams[2];
    rootParams[0].InitAsConstantBufferView(0);  // CBV at b0
    rootParams[1].InitAsUnorderedAccessView(0); // UAV at u0

    CD3DX12_ROOT_SIGNATURE_DESC rootSigDesc(2, rootParams, 0, nullptr,
        D3D12_ROOT_SIGNATURE_FLAG_NONE);

    ComPtr<ID3DBlob> signature;
    ComPtr<ID3DBlob> error;

    if (FAILED(D3D12SerializeRootSignature(&rootSigDesc, D3D_ROOT_SIGNATURE_VERSION_1, &signature, &error))) {
        if (error) {
            OutputDebugStringA((const char*)error->GetBufferPointer());
        }
        OutputDebugStringA("Failed to serialize root signature\n");
        return false;
    }

    if (FAILED(device->CreateRootSignature(0, signature->GetBufferPointer(), signature->GetBufferSize(),
        IID_PPV_ARGS(&rootSignature)))) {
        OutputDebugStringA("Failed to create root signature\n");
        return false;
    }

    return true;
}

bool GPUFFTCompute::LoadShader()
{
    // Read compiled shader file
    std::ifstream shaderFile("fft_gpu.cso", std::ios::binary | std::ios::ate);
    if (!shaderFile.is_open()) {
        // Try relative to build directory
        shaderFile.open("../../../src/fft/fft_gpu.cso", std::ios::binary | std::ios::ate);
        if (!shaderFile.is_open()) {
            OutputDebugStringA("Failed to find fft_gpu.cso shader file\n");
            return false;
        }
    }

    std::streamsize size = shaderFile.tellg();
    shaderFile.seekg(0, std::ios::beg);
    std::vector<char> buffer(size);
    if (!shaderFile.read(buffer.data(), size)) {
        OutputDebugStringA("Failed to read shader file\n");
        return false;
    }
    shaderFile.close();

    // Create compute pipeline state
    D3D12_COMPUTE_PIPELINE_STATE_DESC psoDesc = {};
    psoDesc.pRootSignature = rootSignature.Get();
    psoDesc.CS = { buffer.data(), (size_t)size };

    if (FAILED(device->CreateComputePipelineState(&psoDesc, IID_PPV_ARGS(&pipelineState)))) {
        OutputDebugStringA("Failed to create compute pipeline state\n");
        return false;
    }

    return true;
}

bool GPUFFTCompute::UploadData(const std::vector<float2>& data)
{
    if (data.size() != width * height) {
        OutputDebugStringA("Data size mismatch\n");
        return false;
    }

    // Copy to upload buffer
    void* uploadPtr;
    if (FAILED(uploadBuffer->Map(0, nullptr, &uploadPtr))) {
        OutputDebugStringA("Failed to map upload buffer\n");
        return false;
    }

    memcpy(uploadPtr, data.data(), data.size() * sizeof(float2));
    uploadBuffer->Unmap(0, nullptr);

    // Record copy command
    commandAllocator->Reset();
    commandList->Reset(commandAllocator.Get(), nullptr);
    commandList->CopyBufferRegion(dataBuffer.Get(), 0, uploadBuffer.Get(), 0, data.size() * sizeof(float2));
    commandList->Close();

    // Execute commands
    ID3D12CommandList* ppCommandLists[] = { commandList.Get() };
    commandQueue->ExecuteCommandLists(1, ppCommandLists);

    return WaitForGPU();
}

bool GPUFFTCompute::DownloadData(std::vector<float2>& outData)
{
    outData.resize(width * height);

    // Record readback command
    commandAllocator->Reset();
    commandList->Reset(commandAllocator.Get(), nullptr);
    commandList->CopyBufferRegion(downloadBuffer.Get(), 0, dataBuffer.Get(), 0, width * height * sizeof(float2));
    commandList->Close();

    // Execute commands
    ID3D12CommandList* ppCommandLists[] = { commandList.Get() };
    commandQueue->ExecuteCommandLists(1, ppCommandLists);

    if (!WaitForGPU()) return false;

    // Map and copy readback buffer
    void* readbackPtr;
    if (FAILED(downloadBuffer->Map(0, nullptr, &readbackPtr))) {
        OutputDebugStringA("Failed to map readback buffer\n");
        return false;
    }

    memcpy(outData.data(), readbackPtr, outData.size() * sizeof(float2));
    downloadBuffer->Unmap(0, nullptr);

    return true;
}

bool GPUFFTCompute::UpdateConstants(uint32_t N, uint32_t inverse, uint32_t rowOffset, uint32_t stride)
{
    FFTConstants constants = { N, inverse, rowOffset, stride };

    void* constPtr;
    if (FAILED(constantBuffer->Map(0, nullptr, &constPtr))) {
        OutputDebugStringA("Failed to map constant buffer\n");
        return false;
    }

    memcpy(constPtr, &constants, sizeof(FFTConstants));
    constantBuffer->Unmap(0, nullptr);

    return true;
}

bool GPUFFTCompute::ComputeFFTRows(bool inverse)
{
    // Upload current data
    std::vector<float2> dummyData(width * height, {0, 0});  // Should be pre-populated via UploadData
    
    commandAllocator->Reset();
    commandList->Reset(commandAllocator.Get(), nullptr);

    commandList->SetPipelineState(pipelineState.Get());
    commandList->SetComputeRootSignature(rootSignature.Get());
    
    // Set constant buffer
    commandList->SetComputeRootConstantBufferView(0, constantBuffer->GetGPUVirtualAddress());
    
    // Set UAV descriptor heap
    ID3D12DescriptorHeap* ppHeaps[] = { descriptorHeap.Get() };
    commandList->SetDescriptorHeaps(1, ppHeaps);
    commandList->SetComputeRootDescriptorTable(1, descriptorHeap->GetGPUDescriptorHandleForHeapStart());

    // Dispatch row passes
    for (uint32_t row = 0; row < height; ++row) {
        UpdateConstants(width, inverse ? 1 : 0, row * width, 1);
        commandList->Dispatch(1, 1, 1);  // One group per row
        
        // UAV barrier between passes
        D3D12_RESOURCE_BARRIER barrier = CD3DX12_RESOURCE_BARRIER::UAV(dataBuffer.Get());
        commandList->ResourceBarrier(1, &barrier);
    }

    commandList->Close();

    ID3D12CommandList* ppCommandLists[] = { commandList.Get() };
    commandQueue->ExecuteCommandLists(1, ppCommandLists);

    return WaitForGPU();
}

bool GPUFFTCompute::ComputeFFTColumns(bool inverse)
{
    commandAllocator->Reset();
    commandList->Reset(commandAllocator.Get(), nullptr);

    commandList->SetPipelineState(pipelineState.Get());
    commandList->SetComputeRootSignature(rootSignature.Get());
    
    // Set constant buffer
    commandList->SetComputeRootConstantBufferView(0, constantBuffer->GetGPUVirtualAddress());
    
    // Set UAV descriptor heap
    ID3D12DescriptorHeap* ppHeaps[] = { descriptorHeap.Get() };
    commandList->SetDescriptorHeaps(1, ppHeaps);
    commandList->SetComputeRootDescriptorTable(1, descriptorHeap->GetGPUDescriptorHandleForHeapStart());

    // Dispatch column passes
    for (uint32_t col = 0; col < width; ++col) {
        UpdateConstants(height, inverse ? 1 : 0, col, width);
        commandList->Dispatch(1, 1, 1);  // One group per column
        
        // UAV barrier between passes
        D3D12_RESOURCE_BARRIER barrier = CD3DX12_RESOURCE_BARRIER::UAV(dataBuffer.Get());
        commandList->ResourceBarrier(1, &barrier);
    }

    commandList->Close();

    ID3D12CommandList* ppCommandLists[] = { commandList.Get() };
    commandQueue->ExecuteCommandLists(1, ppCommandLists);

    return WaitForGPU();
}

bool GPUFFTCompute::WaitForGPU()
{
    const uint64_t signal = ++fenceValue;
    if (FAILED(commandQueue->Signal(fence.Get(), signal))) {
        OutputDebugStringA("Failed to signal fence\n");
        return false;
    }

    if (FAILED(fence->SetEventOnCompletion(signal, fenceEvent))) {
        OutputDebugStringA("Failed to set fence event\n");
        return false;
    }

    WaitForSingleObject(fenceEvent, INFINITE);
    return true;
}

// ---------------------------------------------------------------------------
// Convenience functions for complete FFT operations
// ---------------------------------------------------------------------------

std::vector<float2> GPUFFTCompute::ComputeFFT(const std::vector<float>& realData, bool inverse)
{
    if (realData.size() != width * height) {
        throw std::runtime_error("Input data size mismatch");
    }

    // Convert real data to complex (imaginary part = 0)
    std::vector<float2> complexData(realData.size());
    for (size_t i = 0; i < realData.size(); ++i) {
        complexData[i] = { realData[i], 0.0f };
    }

    // Upload to GPU
    if (!UploadData(complexData)) {
        throw std::runtime_error("Failed to upload data to GPU");
    }

    // Compute FFT (rows then columns)
    if (!ComputeFFTRows(inverse)) {
        throw std::runtime_error("Failed to compute FFT rows");
    }
    if (!ComputeFFTColumns(inverse)) {
        throw std::runtime_error("Failed to compute FFT columns");
    }

    // Download results
    std::vector<float2> result;
    if (!DownloadData(result)) {
        throw std::runtime_error("Failed to download FFT results");
    }

    return result;
}

std::vector<float> GPUFFTCompute::ComputeFFT(const std::vector<float2>& complexData, bool inverse)
{
    if (complexData.size() != width * height) {
        throw std::runtime_error("Input data size mismatch");
    }

    // Upload to GPU
    if (!UploadData(complexData)) {
        throw std::runtime_error("Failed to upload data to GPU");
    }

    // Compute FFT (rows then columns)
    if (!ComputeFFTRows(inverse)) {
        throw std::runtime_error("Failed to compute FFT rows");
    }
    if (!ComputeFFTColumns(inverse)) {
        throw std::runtime_error("Failed to compute FFT columns");
    }

    // Download results
    std::vector<float2> complexResult;
    if (!DownloadData(complexResult)) {
        throw std::runtime_error("Failed to download FFT results");
    }

    // Extract real part
    std::vector<float> realResult(complexResult.size());
    for (size_t i = 0; i < complexResult.size(); ++i) {
        realResult[i] = complexResult[i].x;
    }

    return realResult;
}

std::vector<float2> GPUFFTCompute::ComputeFFTComplex(const std::vector<float2>& complexData, bool inverse)
{
    if (complexData.size() != width * height) {
        throw std::runtime_error("Input data size mismatch");
    }

    // Upload to GPU
    if (!UploadData(complexData)) {
        throw std::runtime_error("Failed to upload data to GPU");
    }

    // Compute FFT (rows then columns)
    if (!ComputeFFTRows(inverse)) {
        throw std::runtime_error("Failed to compute FFT rows");
    }
    if (!ComputeFFTColumns(inverse)) {
        throw std::runtime_error("Failed to compute FFT columns");
    }

    // Download results
    std::vector<float2> result;
    if (!DownloadData(result)) {
        throw std::runtime_error("Failed to download FFT results");
    }

    return result;
}
