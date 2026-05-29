#define NOMINMAX
#include <windows.h>
#include "Sim2D.h"

// ********************************************************************************************************************
// Init functions
// ********************************************************************************************************************

std::vector<float> Sim::SetTerrain() {
	std::vector<float> terrain(GRIDSIZE * GRIDSIZE, 0.0f);
	float boundary = 2.f * abs(TERRAIN_HEIGHT);
    for (int y = 0; y < GRIDSIZE; y++) {
        for (int x = 0; x < GRIDSIZE; x++) {
            float xf = (float)x / (GRIDSIZE - 1);
            float yf = (float)y / (GRIDSIZE - 1);
            int i = idx(x,y);

            if (TERRAIN_TYPE == 0) { // Flat
                terrain[i] = TERRAIN_HEIGHT; 
            }
            else if (TERRAIN_TYPE == 1) { // Inclined Plane (Ramp)
                terrain[i] = TERRAIN_HEIGHT + (xf * TERRAIN_SCALE);
            }
            else if (TERRAIN_TYPE == 2) { // Bumpy / Natural (Sum of Sines)
                float noise = 0.5f * sin(10.f * xf) * cos(8.f * yf) + 
                              0.2f * sin(25.f * xf + 2.f * yf);
                terrain[i] = TERRAIN_HEIGHT + (noise * TERRAIN_SCALE);
            }
            else if (TERRAIN_TYPE == 3) { // Two Basins with Divider
                // Create two dips with a ridge at xf = 0.5
				float divider = (xf > 0.4f && xf < 0.6f) ? 5.0f : 0.0f;
				terrain[i] = TERRAIN_HEIGHT + divider;
            }
            else if (TERRAIN_TYPE == 4) { // Beach Scene
                // Simple slope with some noise for "sand dunes"
                float dunes = 0.05f * sin(20.f * yf); 
                terrain[i] = TERRAIN_HEIGHT + TERRAIN_SCALE * (xf * (1 + dunes));
            }
			else if (TERRAIN_TYPE == 5) {// 1D Hill in center of scene
				// Bumpy hill near middle of scene
				float hill = -0.9f + 0.1f * xf + 0.03f * sin(20.f * xf) + 0.9f * sin(2.5f * xf);
				terrain[i] = hill * abs(TERRAIN_HEIGHT * 1.2);
			}
			else if (TERRAIN_TYPE == 6) {// 2D Gaussian hill in center
                float dist = sqrt(pow(xf - 0.5f, 2) + pow(yf - 0.5f, 2));
				float hill = (dist < 0.1f) ? TERRAIN_SCALE * cos(dist * PI * 5.0f) : 0.f;
				terrain[i] = TERRAIN_HEIGHT + hill;
			}
			// terrain[idx(x, GRIDSIZE-1)] = boundary;
			// terrain[idx(x, 0)] = boundary;
        }
		// terrain[idx(GRIDSIZE-1, y)] = boundary;
		// terrain[idx(0, y)] = boundary;
    }
	return terrain;
}

std::vector<float> Sim::SetWater(std::vector<float>& terrain) {
    std::vector<float> h(GRIDSIZE * GRIDSIZE, 0.0f);
    for (int y = 0; y < GRIDSIZE; y++) {
        for (int x = 0; x < GRIDSIZE; x++) {
            float xf = (float)x / (GRIDSIZE - 1);
            float yf = (float)y / (GRIDSIZE - 1);
            int i = idx(x,y);
			float waterSurface = WATER_LEVEL;
			if (WATER_TYPE == 0) { // Flat
				// do nothing
			}
			else if (WATER_TYPE == 1) { // Step/Dam Break
                if (xf < 0.3f) 
					waterSurface += WATER_SCALE;
            }
			else if (WATER_TYPE == 2) { // Diagonal slope on 1st half
				waterSurface += 2 * (1 - xf + yf) * WATER_SCALE; // maybe make this only on upper diagonal later or something
			}
            else if (WATER_TYPE == 3) { // Localized splash (Gaussian)
                float dist = sqrt(pow(xf - 0.5f, 2) + pow(yf - 0.5f, 2));
                if (dist < 0.1f) 
					waterSurface += WATER_SCALE * cos(dist * PI * 5.0f);
            }
			else if (WATER_TYPE == 4) { // Surface Ripples
				waterSurface += 0.5f * WATER_SCALE * (cos(2.f * PI * x / 37.f) + cos(2.f * PI * y / 49.f)); // this needs some work
			}
            else if (WATER_TYPE == 5) { // Basin Flood
                // Fill only the left basin (xf < 0.5)
				if (xf < 0.25f)
                    waterSurface = std::max(WATER_LEVEL, terrain[idx(GRIDSIZE/2, y)] + 2.0f);
                else if (xf < 0.5f) 
					waterSurface = std::max(WATER_LEVEL, terrain[idx(GRIDSIZE/2, y)]);
                else 
					waterSurface = terrain[i]; // Start dry
            }
            h[i] = std::max(0.0f, waterSurface - (float)terrain[i]);
        }
    }
	return h;
}

void Sim::Init(GPU* gpu) {
	// Set up GPU Shaders
	for (int i=0; i < sizeof(shaders) / sizeof(shaders[0]); i++) {
		gpu->CompileComputeShader(L"shaders/kernels.hlsl", names[i], shaders[i]);
	}
    gpu->UpdateConstants(constants);

	// Create GPU Textures and Upload Initial Data
	for (int i=0; i < sizeof(fields) / sizeof(fields[0]); i++) {
		gpu->CreateGridTexture(fields[i], GRIDSIZE);
		std::vector<float> temp(GRIDSIZE * GRIDSIZE, 0.0f);
		gpu->UploadToGPU(fields[i]->tex, temp, GRIDSIZE);
	}
	for (int i = 0; i < sizeof(fields_complex) / sizeof(fields_complex[0]); i++) {
		gpu->CreateGridTexture(fields_complex[i], GRIDSIZE, true);
		std::vector<float> temp(GRIDSIZE * GRIDSIZE * 2, 0.0f);
		gpu->UploadToGPU(fields_complex[i]->tex, temp, GRIDSIZE, true);
	}
	for (int i = 0; i < 2; i++) {
		gpu->CreateGridTexture(fields_arrays[i], GRIDSIZE, true, DEPTH_NUM);
		std::vector<float> temp(GRIDSIZE * GRIDSIZE * DEPTH_NUM * 2, 0.0f);
		gpu->UploadToGPU(fields_arrays[i]->tex, temp, GRIDSIZE, true, DEPTH_NUM);
	}
	gpu->CreateBuffer(&depth, depths.data(), DEPTH_NUM);
	gpu->BindSRV(8, depth.srv);
	std::vector<float> terrain_temp = SetTerrain();
	std::vector<float> h_temp = SetWater(terrain_temp);
	std::vector<float> H_temp(GRIDSIZE*GRIDSIZE, 0.0f);
	for (int i = 0; i < GRIDSIZE*GRIDSIZE; i++) {
		H_temp[i] = terrain_temp[i] + h_temp[i];
	}
	gpu->UploadToGPU(terrain.tex, terrain_temp, GRIDSIZE);
	gpu->UploadToGPU(h.tex, h_temp, GRIDSIZE);
	gpu->UploadToGPU(H.tex, H_temp, GRIDSIZE);
	gpu->UploadToGPU(hbar.tex, h_temp, GRIDSIZE);
	gpu->UploadToGPU(hbarOld.tex, h_temp, GRIDSIZE);
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
	DecompositionStep();
	eWaveStep();
	// FFTStep();
	SWEStep();
	TransportStep();
	ComputeValues();
}

void Sim::DecompositionStep() {
	/******* Bulk vs Surface Wave Decomposition ******/
	// Initialize values for diffusion step
    gpu->Dispatch(InitDecomp, 
		{h.srv, q_x.srv, q_y.srv, terrain.srv}, 
		{H.uav, Q_x.uav, Q_y.uav});

	// Calculate diffusion coefficients
    gpu->Dispatch(CalcDiffusionCoeffs, 
		{H.srv, terrain.srv}, 
		{alpha_H.uav, alpha_Q_x.uav, alpha_Q_y.uav});
	// gpu->Dispatch(ApplyBoundaries, {},
	// 	{alpha_H.uav, alpha_Q_x.uav, alpha_Q_y.uav});
	
	// Run diffusion to low-pass filter H and Q
	for (int j = 0; (j < DIFFUSION_ITERATIONS); j++) {
		// Swap H and HPast pointers for ping-ponging
        std::swap(H, HPast); 
        std::swap(Q_x, QPast_x);
		std::swap(Q_y, QPast_y);

		// Diffusion step for H and Q
		gpu->Dispatch(DiffusionStep, 
			{terrain.srv, HPast.srv, QPast_x.srv, QPast_y.srv, alpha_H.srv, alpha_Q_x.srv, alpha_Q_y.srv}, 
			{H.uav, Q_x.uav, Q_y.uav});
		// gpu->Dispatch(ApplyBoundaries, {}, 
		// 	{H.uav, Q_x.uav, Q_y.uav});
	}

	// final conversion to individual solver quantities
	gpu->Dispatch(DecomposeFields, 
		{H.srv, Q_x.srv, Q_y.srv, h.srv, q_x.srv, q_y.srv, terrain.srv}, 
		{hbar.uav, qbar_x.uav, qbar_y.uav, htilde.uav, qtilde_x.uav, qtilde_y.uav});
	// gpu->Dispatch(ApplyBoundaries, {}, 
	// 	{hbar.uav, qbar_x.uav, qbar_y.uav});
	// gpu->Dispatch(ApplyBoundaries, {}, 
	// 	{htilde.uav, qtilde_x.uav, qtilde_y.uav});

	gpu->Dispatch(InitDecomp, 
		{h.srv, nullptr, nullptr, terrain.srv}, 
		{H.uav});
}

void Sim::eWaveStep() {
	// Copy variables to fourier domain & perform FFT
	gpu->Dispatch(TransferToFFT, 
		{htilde.srv, qtilde_x.srv, qtilde_y.srv}, 
		{htildeOld.uav, hHat.uav, qHat_x.uav, qHat_y.uav});
	gpu->ExecuteFFT(hHat.uav, GRIDSIZE, false);
	gpu->ExecuteFFT(qHat_x.uav, GRIDSIZE, false);
	gpu->ExecuteFFT(qHat_y.uav, GRIDSIZE, false);

	// Compute eWave
	gpu->Dispatch(CalcEWave, 
		{hHat.srv, qHat_x.srv, qHat_y.srv},
		{qHat_x_array.uav, qHat_y_array.uav}, DEPTH_NUM);

	// Inverse FFT fourier variables
	gpu->ExecuteFFT(qHat_x_array.uav, GRIDSIZE, true, DEPTH_NUM);
	gpu->ExecuteFFT(qHat_y_array.uav, GRIDSIZE, true, DEPTH_NUM);

	// Interpolate between depths to get qtilde
	gpu->Dispatch(InterpQ,
		{hbar.srv, qHat_x_array.srv, qHat_y_array.srv},
		{qtilde_x.uav, qtilde_y.uav}, DEPTH_NUM);
	// gpu->Dispatch(ApplyBoundaries, {}, 
	// 	{qtilde_x.uav, qtilde_y.uav});
}

void Sim::SWEStep() {
	// SWE bulk simulation using [Stelling03]

	// qbar to ubar using hbar from last timestep	
	gpu->Dispatch(CalcUbar, 
		{qbar_x.srv, qbar_y.srv, hbarOld.srv,}, 
		{ubar_x.uav, ubar_y.uav});
		
	// Compute time derivative of u_bar and integrate to get new u_bar, then 
	// transfer back to flow rate using upwinding on most recent hbar
	gpu->Dispatch(CalcSWE,
		{ubar_x.srv, ubar_y.srv, hbar.srv, H.srv}, 
		{ubarNew_x.uav, ubarNew_y.uav, qbar_x.uav, qbar_y.uav});
	// gpu->Dispatch(ApplyBoundaries, {}, 
	// 	{ubarNew_x.uav, ubarNew_y.uav, qbar_x.uav, qbar_y.uav});

	// store current hbar for next timestep
	std::swap(hbar, hbarOld);
}

void Sim::TransportStep() {
	/****** Advect high-frequency wave height and flow rate through bulk velocity ******/ 

	// Adjust qtilde to account for advection by ubar, using cubic sampling to get better accuracy.
	std::swap(qtilde_x, qtildePast_x);
	std::swap(qtilde_y, qtildePast_y);
	gpu->Dispatch(UpdateTilde, 
		{ubarNew_x.srv, ubar_x.srv, ubarNew_y.srv, ubar_y.srv, 
			qtildePast_x.srv, qtildePast_y.srv, h.srv, htilde.srv},
		{htilde.uav, qtilde_x.uav, qtilde_y.uav});	
	// gpu->Dispatch(ApplyBoundaries, {}, 
	// 	{htilde.uav, qtilde_x.uav, qtilde_y.uav});
	
	// Advection of h through ubar:
	// Construct q_advect = ubar * htilde sampled at cell edges using cubic sampling
	gpu->Dispatch(CalcQAdvect, 
		{ubarNew_x.srv, ubarNew_y.srv, htilde.srv},
		{qAdvect_x.uav, qAdvect_y.uav});

	std::swap(h, hPast);
}
	
void Sim::ComputeValues() {
	gpu->Dispatch(IntegrateH, 
		{qbar_x.srv, qtilde_x.srv, qAdvect_x.srv, qbar_y.srv, qtilde_y.srv, qAdvect_y.srv, 
			hPast.srv, terrain.srv}, 
		{h.uav, q_x.uav, q_y.uav});
	// gpu->Dispatch(ApplyBoundaries, {}, 
	// 	{h.uav, q_x.uav, q_y.uav});

	gpu->Dispatch(InitDecomp, 
		{h.srv, nullptr, nullptr, terrain.srv}, 
		{H.uav});
}