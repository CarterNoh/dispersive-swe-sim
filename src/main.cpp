#include <iostream>
#include <vector>
#include <complex>
#include <chrono>
#include <string>
#include <fstream>
#include <filesystem>
#include <windows.h>
#include "Sim2D.h"

extern "C" {
    __declspec(dllexport) DWORD NvOptimusEnablement = 1;      // For NVIDIA GPUs
    __declspec(dllexport) int AmddxGpuControlSelect = 1;      // For AMD GPUs
}

namespace fs = std::filesystem;

// Parameters
int numTicks = 2000; // Number of simulation steps to run
bool render = true; // Whether to render the simulation or just run it headless

void SaveToCSV(const std::vector<float>& h, int size, int tick) {
    if (!fs::exists("data")) {
        fs::create_directory("data");
    }
    std::string filename = "data/frame_" + std::to_string(tick) + ".csv";
    if (tick < 0) {
        filename = "data/terrain.csv";
    }
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << " for writing!" << std::endl;
        return;
    }
    for (int y = 0; y < size; ++y) {
        for (int x = 0; x < size; ++x) {
            file << h[y * size + x] << (x == size - 1 ? "" : ",");
        }
        file << "\n";
    }
    file.close();
}

void SaveToCSV(const std::vector<std::complex<float>>& data, int size, int tick) {
    if (!fs::exists("data")) {
        fs::create_directory("data");
    }
    std::string filename = "data/frame_" + std::to_string(tick) + ".csv";
    if (tick < 0) {
        filename = "data/terrain.csv";
    }
    std::ofstream file(filename);
    if (!file.is_open()) return;

    for (int y = 0; y < size; ++y) {
        for (int x = 0; x < size; ++x) {
            std::complex<float> val = data[y * size + x];
            
            // std::abs() automatically calculates the magnitude (sqrt(r^2 + i^2))
            // std::arg(val) would give you the phase angle!
            float magnitude = std::abs(val);

            file << val.real() << "," << val.imag();// << "," << magnitude;
            
            if (x < size - 1) file << ",";
        }
        file << "\n";
    }
}

void SaveToCSV(const std::vector<std::complex<float>>& data, int size, int arraySize, int tick) {
    if (!fs::exists("data")) {
        fs::create_directory("data");
    }
    std::string filename = "data/frame_" + std::to_string(tick) + ".csv";
    if (tick < 0) {
        filename = "data/terrain.csv";
    }
    std::ofstream file(filename);
    if (!file.is_open()) return;
    // // Header optimized for Python / Pandas DataFrames
    // file << "x,y,layer,real,imag,magnitude\n";
    for (int layer = 0; layer < arraySize; ++layer) {
        int sliceOffset = layer * (size * size);
        for (int y = 0; y < size; ++y) {
            for (int x = 0; x < size; ++x) {
                std::complex<float> val = data[sliceOffset + (y * size) + x];
                float magnitude = std::abs(val);
                file //<< x << "," 
                     //<< y << "," 
                     //<< layer << "," 
                     << val.real() << "," 
                    //  << val.imag() << "," 
                     //<< magnitude 
                     //<< "\n"
                     ;
            }
            file << "\n";
        }
    }
}

LRESULT CALLBACK WindowProc(HWND hwnd, UINT uMsg, WPARAM wParam, LPARAM lParam) {
    // The callback function Windows uses to handle clicks, keys, and closing
    if (uMsg == WM_DESTROY) {
        PostQuitMessage(0);
        return 0;
    }
    return DefWindowProc(hwnd, uMsg, wParam, lParam);
}

RenderConstants InitCamera(int gridSize) {
    // These values must match those in the Vertex Shader
    float maxInitialHeight = 0.2f; 
    float visualScale = 1.f;
    float trueMax = maxInitialHeight * visualScale;

    // Focus the camera on the center point
    float center = gridSize / 2.0f;
    DirectX::XMVECTOR focusPoint = DirectX::XMVectorSet(center, trueMax/2.0f, center, 1.0f);

    // Define camera angles
    float pitch =  25.0f * (DirectX::XM_PI / 180.0f);
    float yaw   = 30.0f * (DirectX::XM_PI / 180.0f);

    // Define the distance to pull the camera back
    float distance = gridSize * 1.5f;

    // Calculate the new Eye Position using a Rotation Matrix
    // We start with the camera pulled straight back along the Z-axis...
    DirectX::XMVECTOR localEye = DirectX::XMVectorSet(0.0f, 0.0f, -distance, 1.0f);
    // ...then we rotate it by Pitch and Yaw...
    DirectX::XMMATRIX rotation = DirectX::XMMatrixRotationRollPitchYaw(pitch, yaw, 0.0f);
    DirectX::XMVECTOR rotatedEye = DirectX::XMVector3Transform(localEye, rotation);
    // ...and finally, move it so it orbits around the center of the grid.
    DirectX::XMVECTOR eyePosition = DirectX::XMVectorAdd(focusPoint, rotatedEye);

    // Build the View Matrix
    DirectX::XMVECTOR upDirection = DirectX::XMVectorSet(0.0f, 1.0f, 0.0f, 0.0f);
    DirectX::XMMATRIX view = DirectX::XMMatrixLookAtLH(eyePosition, focusPoint, upDirection);

    // Use an Orthographic Projection for strict mathematical lines (No perspective distortion)
    float gridDiagonal = gridSize * 1.414f;
    float padding = 1.1f;
    float viewWidth = gridDiagonal * padding;
    float tiltedGridHeight = gridDiagonal * sin(pitch);
    float viewHeight = (tiltedGridHeight + trueMax) * padding;
    DirectX::XMMATRIX projection = DirectX::XMMatrixOrthographicLH(viewWidth, viewHeight, 1.0f, 2000.0f);

    RenderConstants rConsts = {};
    rConsts.viewProjection = DirectX::XMMatrixTranspose(view * projection);
    rConsts.gridSize = (float)gridSize;
    rConsts.cellSize = (float)Sim::CELLSIZE;
    return rConsts;
}

int RunWithRender() {
    HINSTANCE hInstance = GetModuleHandle(NULL);
    // Register the Window Blueprint
    const char CLASS_NAME[] = "SimClass";
    WNDCLASS wc = {};
    wc.lpfnWndProc = WindowProc;
    wc.hInstance = hInstance;
    wc.lpszClassName = CLASS_NAME;
    RegisterClass(&wc);

    // Create the Window
    HWND hwnd = CreateWindowEx(0, CLASS_NAME, "DirectX 11 Simulation",
        WS_OVERLAPPEDWINDOW, CW_USEDEFAULT, CW_USEDEFAULT, 800, 600,
        nullptr, nullptr, hInstance, nullptr);

    if (hwnd == nullptr) return 0;
    ShowWindow(hwnd, SW_SHOW);

    // Initialize Simulation
    Sim sim(hwnd); 
    // Set up a camera looking at the center of the grid
    RenderConstants rConsts = InitCamera(sim.GRIDSIZE);
    sim.gpu->UpdateRenderConstants(rConsts);
    sim.gpu->Render(sim.H.srv);

    // Loop
    MSG msg = {};
    bool running = true;
    int tick = 0;
    
    std::cout << "Starting simulation loop..." << std::endl;
    auto totalStart = std::chrono::high_resolution_clock::now();
    while (running) {
        // Handle Windows OS messages (like moving the window or closing it)
        while (PeekMessage(&msg, nullptr, 0, 0, PM_REMOVE)) {
            if (msg.message == WM_QUIT) running = false;
            TranslateMessage(&msg);
            DispatchMessage(&msg);
        }
        sim.SimStep();
        sim.gpu->Render(sim.H.srv);
        // sim.gpu->Render({sim.H.srv, sim.disp_x.srv, sim.disp_y.srv});
        // sim.gpu->Render({sim.H.srv});
        tick++;
        if (tick >= numTicks) running = false;
    }
    auto totalEnd = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> totalElapsed = totalEnd - totalStart;
    std::cout << "Total Simulation Time for " << numTicks << " ticks: " << totalElapsed.count() << "ms" << std::endl;
    std::cout << "Average Speed: " << ((totalElapsed.count() / tick)) << " ms/tick" << std::endl;
    std::cout << "FPS: " << (1000.0 / (totalElapsed.count() / tick)) << " frames/sec" << std::endl;

    return 0;
}

int RunHeadless() {
    // Clear the data folder before starting the simulation
    if (fs::exists("data")) {  
        for (const auto& entry : fs::directory_iterator("data")) {
            fs::remove(entry.path());
        }
    }

    // Initialize the simulation
    Sim sim;
    int size = sim.GRIDSIZE;
    int arraySize = sim.DEPTH_NUM;
    std::vector<float> terrain(size * size);
    std::vector<float> H(size * size);
    std::vector<std::complex<float>> hHat(size*size);
    std::vector<std::complex<float>> array;
    // sim.gpu->DownloadFromGPU(sim.terrain.tex, terrain, size);
    // sim.gpu->DownloadFromGPU(sim.H.tex, H, size);
    // SaveToCSV(terrain, size, -1); // Save terrain data for reference
    // SaveToCSV(H, size, 0); // Save initial water height data

    std::cout << "Starting simulation loop..." << std::endl;
    auto totalStart = std::chrono::high_resolution_clock::now();
    for (int tick = 1; tick < numTicks; ++tick) {
        auto start = std::chrono::high_resolution_clock::now();
        sim.SimStep();
        sim.gpu->DownloadFromGPU(sim.H.tex, H, size);
        SaveToCSV(H, size, tick);
        // sim.gpu->DownloadFromGPU(sim.HPos.tex, hHat, size);
        // SaveToCSV(hHat, size, tick);
        // sim.gpu->DownloadFromGPU(sim.qHat_x_array.tex, array, size, arraySize);
        // SaveToCSV(array, size, arraySize, tick);

        auto sim_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> elapsed = sim_time - start;
        std::cout << "Tick: " << tick << " | Simulation Time: " << elapsed.count() << "ms" << std::endl;
    }
    auto totalEnd = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> totalElapsed = totalEnd - totalStart;
    std::cout << "Total Simulation Time for " << numTicks << " ticks: " << totalElapsed.count() << "ms" << std::endl;
    std::cout << "Average Speed: " << ((totalElapsed.count() / 1000 / numTicks)) << " sec/tick" << std::endl;
    std::cout << "FPS: " << (1000.0 / (totalElapsed.count() / numTicks)) << " frames/sec" << std::endl;

    return 0;
}

int main() {

    if (render) {
        return RunWithRender();
    }
    else {
        return RunHeadless();
    }
}