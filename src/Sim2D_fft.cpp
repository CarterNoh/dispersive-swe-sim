#define NOMINMAX
#include <windows.h>
#include "Sim2D_fft.h"

// ********************************************************************************************************************
// Init functions
// ********************************************************************************************************************

void Sim::Init(GPU* gpu) {
	// Set up GPU Shaders
	for (int i=0; i < sizeof(shaders) / sizeof(shaders[0]); i++) {
		gpu->CompileComputeShader(L"shaders/fftwaves.hlsl", names[i], shaders[i]);
	}
    gpu->UpdateConstants(constants);

	// Create GPU Textures and Upload Initial Data
	for (int i=0; i < sizeof(fields) / sizeof(fields[0]); i++) {
		gpu->CreateGridTexture(fields[i], GRIDSIZE);
		std::vector<float> temp(GRIDSIZE * GRIDSIZE, 0.0f);
		gpu->UploadToGPU(fields[i]->tex, temp, GRIDSIZE);
	}

    // Complex Textures
	for (int i = 0; i < sizeof(fields_complex) / sizeof(fields_complex[0]); i++) {
		gpu->CreateGridTexture(fields_complex[i], GRIDSIZE, true);
		std::vector<float> temp(GRIDSIZE * GRIDSIZE * 2, 0.0f);
		gpu->UploadToGPU(fields_complex[i]->tex, temp, GRIDSIZE, true);
	}

    gpu->Dispatch(PopulateSpectrum, {}, 
		{omega.uav, HPos.uav, HNeg.uav});
}

Sim::Sim() {
	// Init GPU
	gpu = new GPU();
	if (!gpu->Init(GRIDSIZE)) std::cerr << "GPU INIT FAILED" << std::endl;
	Sim::Init(gpu);
}

Sim::Sim(HWND hwnd = nullptr) {
	gpu = new GPU();
	if (!gpu->Init(GRIDSIZE, hwnd)) std::cerr << "GPU INIT FAILED" << std::endl;
	Sim::Init(gpu);

	// Set up rendering shaders and mesh
	gpu->CompileVertexShader(L"shaders/render.hlsl", "VSMain");
	gpu->CompilePixelShader(L"shaders/render.hlsl", "PSMain");
	gpu->CreateGridMesh(GRIDSIZE);
	gpu->CreateGridVertexBuffer(GRIDSIZE);
}

int Sim::Release(void) {
	return 0;
}


// ********************************************************************************************************************
// Simulation functions
// ********************************************************************************************************************

void Sim::SimStep() {
	// Update timestep
	time += TIMESTEP;
	constants.time = time;
	gpu->UpdateConstants(constants);
	
	// Propagate waves
	gpu->Dispatch(PropagateWaves, 
		{omega.srv, HPos.srv, HNeg.srv}, 
		{HProp.uav, DxProp.uav, DyProp.uav});
	gpu->ExecuteFFT(HProp.uav, GRIDSIZE, true);
	gpu->ExecuteFFT(DxProp.uav, GRIDSIZE, true);
	gpu->ExecuteFFT(DyProp.uav, GRIDSIZE, true);

	// Copy real part of complex textures to H, Dx, Dy
	gpu->Dispatch(ComplexToReal, {HProp.srv}, {H.uav});
	gpu->Dispatch(ComplexToReal, {DxProp.srv}, {Dx.uav});
	gpu->Dispatch(ComplexToReal, {DyProp.srv}, {Dy.uav});
}
