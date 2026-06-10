#define NOMINMAX
#include <windows.h>
#include "Sim2D_fft.h"

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
		gpu->CompileComputeShader(L"shaders/fftwaves.hlsl", names[i], shaders[i]);
	}
    gpu->UpdateConstants(constants);

	// Create GPU Textures and Upload Initial Data
	// Real Textures
	for (int i=0; i < sizeof(fields_r) / sizeof(fields_r[0]); i++) {
		gpu->CreateGridTexture(fields_r[i], GRIDSIZE);
		std::vector<float> temp(GRIDSIZE * GRIDSIZE, 0.0f);
		gpu->UploadToGPU(fields_r[i]->tex, temp, GRIDSIZE);
	}
    // // Complex Textures
	// for (int i = 0; i < sizeof(fields_c) / sizeof(fields_c[0]); i++) {
	// 	gpu->CreateGridTexture(fields_c[i], GRIDSIZE, true);
	// 	std::vector<float> temp(GRIDSIZE * GRIDSIZE * 2, 0.0f);
	// 	gpu->UploadToGPU(fields_c[i]->tex, temp, GRIDSIZE, true);
	// }
	// // Real Texture Arrays
	// for (int i = 0; i < sizeof(fields_arrays_r) / sizeof(fields_arrays_r[0]); i++) {
	// 	gpu->CreateGridTexture(fields_arrays_r[i], GRIDSIZE, false, DEPTH_NUM);
	// 	std::vector<float> temp(GRIDSIZE * GRIDSIZE * DEPTH_NUM, 0.0f);
	// 	gpu->UploadToGPU(fields_arrays_r[i]->tex, temp, GRIDSIZE, false, DEPTH_NUM);
	// }
	// Complex Texture Arrays
	for (int i = 0; i < sizeof(fields_arrays_c) / sizeof(fields_arrays_c[0]); i++) {
		gpu->CreateGridTexture(fields_arrays_c[i], GRIDSIZE, true, DEPTH_NUM);
		std::vector<float> temp(GRIDSIZE * GRIDSIZE * DEPTH_NUM * 2, 0.0f);
		gpu->UploadToGPU(fields_arrays_c[i]->tex, temp, GRIDSIZE, true, DEPTH_NUM);
	}

	// terrain
	std::vector<float> terrain_temp = SetTerrain();
	gpu->UploadToGPU(terrain.tex, terrain_temp, GRIDSIZE);
	// nothing for H, assume FFT waves fluctuate around 0. Will need to adjust for tide possibly at some point
	gpu->CreateBuffer(&depth, depths.data(), DEPTH_NUM);
	gpu->BindSRV(8, depth.srv);

	// Initialize FFT Wave Spectrum
    gpu->Dispatch(PopulateSpectrum, {}, 
		{HPos.uav, HNeg.uav}, DEPTH_NUM);
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
		{HPos.srv, HNeg.srv},
		{HProp.uav, DxProp.uav, DyProp.uav, UxProp.uav, UyProp.uav, DelSx.uav, DelSy.uav}, DEPTH_NUM);
	gpu->ExecuteFFT(HProp.uav,  GRIDSIZE, true, DEPTH_NUM);
	gpu->ExecuteFFT(DxProp.uav, GRIDSIZE, true, DEPTH_NUM);
	gpu->ExecuteFFT(DyProp.uav, GRIDSIZE, true, DEPTH_NUM);
	gpu->ExecuteFFT(UxProp.uav, GRIDSIZE, true, DEPTH_NUM);
	gpu->ExecuteFFT(UyProp.uav, GRIDSIZE, true, DEPTH_NUM);
	gpu->ExecuteFFT(DelSx.uav,  GRIDSIZE, true, DEPTH_NUM);
	gpu->ExecuteFFT(DelSy.uav,  GRIDSIZE, true, DEPTH_NUM);

	// Interpolate outputs between depths
	gpu->Dispatch(Interp, 
		{HProp.srv, DxProp.srv, DyProp.srv, UxProp.srv, UyProp.srv, DelSx.srv, DelSy.srv, terrain.srv}, 
		{H.uav, Dx.uav, Dy.uav, Ux.uav, Uy.uav, Sx.uav, Sy.uav});
}
