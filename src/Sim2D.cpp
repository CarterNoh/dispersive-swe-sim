#define NOMINMAX
#include <windows.h>
#include "Sim2D.h"

/*
Note: Every loop through the grid only covers interior cells and a separate loop covers boundaries. This was done 
to increase simulation speed, but added a lot of complexity. I'm not positive that all the boundary checks are 
right, and one thing to do eventually would be compare the speed of the sim with and without parts to see where 
it matters. It would cut down on complexity and total length of code to do it all in the loop oinstead of separate 
boundary handling. 
*/

// ********************************************************************************************************************
// Helper functions
// ********************************************************************************************************************

inline float Sim::LimitFlowRate(float flow_rate_in, float waterDepth_left, float waterDepth_right)
{
	if (flow_rate_in >= 0.f)
		return std::min(flow_rate_in, CFL_CONDITION * waterDepth_left * CELLSIZE / TIMESTEP);  // 0.25 since other neighbor might take from this source cell as well
	else
		return std::max(flow_rate_in, -CFL_CONDITION * waterDepth_right * CELLSIZE / TIMESTEP);
}

inline float Sim::LimitVelocity(float velocity_in)
{
	if (velocity_in >= 0.f)
		return std::min(velocity_in, CFL_CONDITION * CELLSIZE / TIMESTEP);   // 0.25 since other neighbors might take from this source cell as well
	else
		return std::max(velocity_in, -CFL_CONDITION * CELLSIZE / TIMESTEP);
}

bool Sim::StopFlowOnTerrainBoundary(int x, int y, std::vector<float>& h, std::vector<float>& terrain, bool isYDirection = false)
{
	// test if the terrain boundary stops any flow across x+0.5
	// Key: 1 = stop in x, 2 = stop in y, 3 = stop in both, 0 = no stop
	bool result_x = 0;
	bool result_y = 0;

	// Test x boundary
	if (!isYDirection)
	{
		if ((h[idx(x,y)] <= MIN_WATER_HEIGHT) && (terrain[idx(x,y)] >= terrain[idx(x+1,y)] + h[idx(x+1,y)])) // positive q_x
			return true;
		if ((h[idx(x+1,y)] <= MIN_WATER_HEIGHT) && (terrain[idx(x+1,y)] > terrain[idx(x,y)] + h[idx(x,y)])) // negative q_x
			return true;
		return false;
	}
	else
	{
		if ((h[idx(x,y)] <= MIN_WATER_HEIGHT) && (terrain[idx(x,y)] >= terrain[idx(x,y+1)] + h[idx(x,y+1)])) // positive q_y
			return true;
		if ((h[idx(x,y+1)] <= MIN_WATER_HEIGHT) && (terrain[idx(x,y+1)] > terrain[idx(x,y)] + h[idx(x,y)])) // negative q_y
			return true;
		return false;
	}
}



// Boundary Condition Functions
void Sim::HandleWallBoundary(std::vector<float>& field) 
{
	// Handle Top and Bottom Edges (Horizontal)
	for (int i = 0; i < GRIDSIZE; ++i) {
		field[idx(i, 0)] = field[idx(i, 1)]; // Bottom Edge
		field[idx(i, GRIDSIZE-1)] = field[idx(i, GRIDSIZE-2)]; // Top Edge
		field[idx(0, i)] = field[idx(1, i)]; // Left Edge
		field[idx(GRIDSIZE-1, i)] = field[idx(GRIDSIZE-2, i)]; // Right Edge
	}
	// Corners
	field[idx(0, 0)] = field[idx(1, 1)]; // Bottom Left Corner
	field[idx(0, GRIDSIZE-1)] = field[idx(1, GRIDSIZE-2)]; // Top Left Corner
	field[idx(GRIDSIZE-1, 0)] = field[idx(GRIDSIZE-2, 1)]; // Bottom Right Corner
	field[idx(GRIDSIZE-1, GRIDSIZE-1)] = field[idx(GRIDSIZE-2, GRIDSIZE-2)]; // Top Right Corner
}

void Sim::HandleFreeBoundary(std::vector<float>& field)
{
	// Free boundary: extrapolate values from interior to boundaries using linear extrapolation
	for (int i = 0; i < GRIDSIZE; ++i) {
		field[idx(i, 0)] = 2.f * field[idx(i, 1)] - field[idx(i, 2)]; // Top Edge
		field[idx(i, GRIDSIZE-1)] = 2.f * field[idx(i, GRIDSIZE-2)] - field[idx(i, GRIDSIZE-3)]; // Bottom Edge
		field[idx(0, i)] = 2.f * field[idx(1, i)] - field[idx(2, i)]; // Left Edge
		field[idx(GRIDSIZE-1, i)] = 2.f * field[idx(GRIDSIZE-2, i)] - field[idx(GRIDSIZE-3, i)]; // Right Edge
	}
	// Corners
	field[idx(0, 0)] = 2.f * field[idx(1, 1)] - field[idx(2, 2)]; // Bottom Left Corner
	field[idx(0, GRIDSIZE-1)] = 2.f * field[idx(1, GRIDSIZE-2)] - field[idx(2, GRIDSIZE-3)]; // Top Left Corner
	field[idx(GRIDSIZE-1, 0)] = 2.f * field[idx(GRIDSIZE-2, 1)] - field[idx(GRIDSIZE-3, 2)]; // Bottom Right Corner
	field[idx(GRIDSIZE-1, GRIDSIZE-1)] = 2.f * field[idx(GRIDSIZE-2, GRIDSIZE-2)] - field[idx(GRIDSIZE-3, GRIDSIZE-3)]; // Top Right Corner
}

void Sim::HandleZeroBoundary(std::vector<float>& field)
{
	for (int i = 0; i < GRIDSIZE; ++i) {
		field[idx(i, 0)] = WATER_LEVEL; // Top Edge
		field[idx(i, GRIDSIZE - 1)] = WATER_LEVEL; // Bottom Edge
		field[idx(0, i)] = WATER_LEVEL; // Left Edge
		field[idx(GRIDSIZE - 1, i)] = WATER_LEVEL; // Right Edge
	}
}

void Sim::ApplyBoundaries(std::vector<float>& field, int type)
{
	if (type == 0)
		Sim::HandleWallBoundary(field);
	else if (type == 1)
		Sim::HandleFreeBoundary(field);
	else if (type == 2)
		Sim::HandleZeroBoundary(field);
}


// ********************************************************************************************************************
// Init functions
// ********************************************************************************************************************

std::vector<float> Sim::SetTerrain(int type) {
	std::vector<float> terrain(GRIDSIZE * GRIDSIZE, 0.0f);
    for (int y = 0; y < GRIDSIZE; y++) {
        for (int x = 0; x < GRIDSIZE; x++) {
            float xf = (float)x / (GRIDSIZE - 1);
            float yf = (float)y / (GRIDSIZE - 1);
            int i = idx(x, y);

            if (type == 0) { // Flat
                terrain[i] = TERRAIN_HEIGHT; 
            }
            else if (type == 1) { // Inclined Plane (Ramp)
                terrain[i] = TERRAIN_HEIGHT + (xf * TERRAIN_SCALE);
            }
            else if (type == 2) { // Bumpy / Natural (Sum of Sines)
                float noise = 0.5f * sin(10.f * xf) * cos(8.f * yf) + 
                              0.2f * sin(25.f * xf + 2.f * yf);
                terrain[i] = TERRAIN_HEIGHT + (noise * TERRAIN_SCALE);
            }
            else if (type == 3) { // Two Basins with Divider
                // Create two dips with a ridge at xf = 0.5
                float divider = (xf > 0.48f && xf < 0.52f) ? 5.0f : 0.0f;
                terrain[i] = TERRAIN_HEIGHT + divider;
            }
            else if (type == 4) { // Beach Scene
                // Simple slope with some noise for "sand dunes"
                float slope = xf * 15.0f;
                float dunes = 0.5f * sin(30.f * yf) * xf; 
                terrain[i] = TERRAIN_HEIGHT + slope + dunes;
            }
        }
    }
	return terrain;
}

std::vector<float> Sim::SetWater(int type, float level, std::vector<float>& terrain) {
    std::vector<float> h(GRIDSIZE * GRIDSIZE, 0.0f);
    for (int y = 0; y < GRIDSIZE; y++) {
        for (int x = 0; x < GRIDSIZE; x++) {
            float xf = (float)x / (GRIDSIZE - 1);
            float yf = (float)y / (GRIDSIZE - 1);
            int i = idx(x, y);

            float waterSurface = level;

            if (type == 0) { // Localized splash (Gaussian pill)
                float dist = sqrt(pow(xf - 0.5f, 2) + pow(yf - 0.5f, 2));
                if (dist < 0.1f) 
					waterSurface += 10.0f * cos(dist * PI * 5.0f);
            }
            else if (type == 1) { // Step/Dam Break
                if (xf < 0.3f) 
					waterSurface += 10.0f;
            }
            else if (type == 2) { // Basin Flood
                // Fill only the left basin (xf < 0.5)
				if (xf < 0.25f)
                    waterSurface = std::max(level, terrain[idx(GRIDSIZE/2, y)] + 2.0f);
                else if (xf < 0.5f) 
					waterSurface = std::max(level, terrain[idx(GRIDSIZE/2, y)]);
                else 
					waterSurface = terrain[i]; // Start dry
            }

            // Standardize height as (Surface - Terrain)
            h[i] = std::max(0.0f, waterSurface - (float)terrain[i]);
        }
    }
	return h;
}

Sim::Sim()
{
	// Init GPU
	gpu = new GPU();
	if (!gpu->Init()) std::cerr << "GPU INIT FAILED" << std::endl;

	// Set up GPU Shaders
	for (int i=0; i < sizeof(shaders) / sizeof(shaders[0]); i++) {
		gpu->CompileComputeShader(L"shaders/kernels.hlsl", names[i], shaders[i]);
	}
    gpu->UpdateConstants(constants);

	// Create GPU Textures and Upload Initial Data
	for (int i=0; i < sizeof(fields) / sizeof(fields[0]); i++) {
		gpu->CreateGridTexture(fields[i], GRIDSIZE);
		std::vector<float> temp(GRIDSIZE*GRIDSIZE, 0.0f);
		gpu->UploadToGPU(fields[i]->tex, temp, GRIDSIZE);
	}
	std::vector<float> terrain_temp = SetTerrain(TERRAIN_TYPE);
	std::vector<float> h_temp = SetWater(WATER_TYPE, WATER_LEVEL, terrain_temp);
	gpu->UploadToGPU(terrain.tex, terrain_temp, GRIDSIZE);
	gpu->UploadToGPU(h.tex, h_temp, GRIDSIZE);
	gpu->UploadToGPU(hbar.tex, h_temp, GRIDSIZE);
	gpu->UploadToGPU(hbarOld.tex, h_temp, GRIDSIZE);
}

Sim::Sim(HWND hwnd = nullptr)
{
	gpu = new GPU();
	if (!gpu->Init(hwnd)) std::cerr << "GPU INIT FAILED" << std::endl;

	// Set up GPU Shaders
	for (int i=0; i < sizeof(shaders) / sizeof(shaders[0]); i++) {
		gpu->CompileComputeShader(L"shaders/kernels.hlsl", names[i], shaders[i]);
	}
    gpu->UpdateConstants(constants);

	// Create GPU Textures and Upload Initial Data
	for (int i=0; i < sizeof(fields) / sizeof(fields[0]); i++) {
		gpu->CreateGridTexture(fields[i], GRIDSIZE);
		std::vector<float> temp(GRIDSIZE*GRIDSIZE, 0.0f);
		gpu->UploadToGPU(fields[i]->tex, temp, GRIDSIZE);
	}
	std::vector<float> terrain_temp = SetTerrain(TERRAIN_TYPE);
	std::vector<float> h_temp = SetWater(WATER_TYPE, WATER_LEVEL, terrain_temp);
	std::vector<float> H_temp(GRIDSIZE*GRIDSIZE, 0.0f);
	for (int i = 0; i < GRIDSIZE*GRIDSIZE; i++) {
		H_temp[i] = terrain_temp[i] + h_temp[i];
	}
	gpu->UploadToGPU(terrain.tex, terrain_temp, GRIDSIZE);
	gpu->UploadToGPU(h.tex, h_temp, GRIDSIZE);
	gpu->UploadToGPU(H.tex, H_temp, GRIDSIZE);
	gpu->UploadToGPU(hbar.tex, h_temp, GRIDSIZE);
	gpu->UploadToGPU(hbarOld.tex, h_temp, GRIDSIZE);

	// Set up rendering shaders and mesh
	gpu->CompileVertexShader(L"shaders/render.hlsl", "VSMain");
	gpu->CompilePixelShader(L"shaders/render.hlsl", "PSMain");
	gpu->CreateGridMesh(GRIDSIZE);
	gpu->CreateGridVertexBuffer(GRIDSIZE);
    gpu->UpdateConstants(constants);
}

int Sim::Release(void)
{
	return 0;
}


// ********************************************************************************************************************
// Simulation functions
// ********************************************************************************************************************

void Sim::SimStep()
{
	DecompositionStep();
	// eWaveStep();
	// 	FFTStep();
	SWEStep();
	TransportStep();
	ComputeValues();
	// time += TIMESTEP;
}

void Sim::DecompositionStep()
{
	/******* Bulk vs Surface Wave Decomposition ******/
	// Initialize values for diffusion step
    gpu->Dispatch(shader_InitDecomp, {h.srv, q_x.srv, q_y.srv, terrain.srv}, 
		{H.uav, Q_x.uav, Q_y.uav}, group, group);

	// Calculate diffusion coefficients
    gpu->Dispatch(shader_CalcDiff, {H.srv}, 
		{alpha_H.uav, alpha_Q_x.uav, alpha_Q_y.uav}, group, group);
	gpu->Dispatch(shader_Boundaries, {}, 
		{alpha_H.uav, alpha_Q_x.uav, alpha_Q_y.uav}, group, group);
	
	// Run diffusion to low-pass filter H and Q
	for (int j = 0; (j < DIFFUSION_ITERATIONS); j++)
	{
		// Swap H and HPast pointers for ping-ponging
        std::swap(H, HPast); 
        std::swap(Q_x, QPast_x);
		std::swap(Q_y, QPast_y);

		// Diffusion step for H and Q
		gpu->Dispatch(shader_Diffusion, 
			{terrain.srv, HPast.srv, QPast_x.srv, QPast_y.srv, alpha_H.srv, alpha_Q_x.srv, alpha_Q_y.srv}, 
			{H.uav, Q_x.uav, Q_y.uav}, group, group);
		gpu->Dispatch(shader_Boundaries, {}, 
			{H.uav, Q_x.uav, Q_y.uav}, group, group);
	}
	// gpu->Dispatch(shader_Boundaries,	{}, 
	// 	{H.uav, Q_x.uav, Q_y.uav}, group, group);

	// final conversion to individual solver quantities
	gpu->Dispatch(shader_Decompose, 
		{H.srv, Q_x.srv, Q_y.srv, h.srv, q_x.srv, q_y.srv, terrain.srv}, 
		{hbar.uav, qbar_x.uav, qbar_y.uav, htilde.uav, qtilde_x.uav, qtilde_y.uav}, 
		group, group);
	gpu->Dispatch(shader_Boundaries, {}, 
		{hbar.uav, htilde.uav}, group, group);
	gpu->Dispatch(shader_Boundaries, {}, 
		{qbar_x.uav, qbar_y.uav, qtilde_x.uav, qtilde_y.uav}, group, group);
}

// void Sim::eWaveStep()
// {
// 	// surface velocity update using eWave
// 	for (int x = 0; x < GRIDSIZE; x++)
// 	{
// 		htildehat[x].x = 0.5f * (htilde[x] + htildeOld[x]);
// 		htildeOld[x] = htilde[x];
// 		htildehat[x].y = 0.;
// 		qtildehat[x].x = qtilde[x];
// 		qtildehat[x].y = 0.;
// 	}
// 	fftc1d(htildehat);   //https://www.alglib.net/download.php#cpp
// 	fftc1d(qtildehat);
	
// 	for (int x = 0; x < GRIDSIZE; x++)
// 	{
// 		// physical k from grid position
// 		double kx = GRIDSIZE / 2. - abs(GRIDSIZE / 2. - x);  // this gives [0,..,m_gridSizeX / 2.f-1, m_gridSizeX / 2.f, .. 1]
// 		double k = 2. * PI * fabs(kx) / GRIDSIZE / CELLSIZE;
// 		double kNonZero = std::max(0.01, k);
// 		double kS = k;  // signed k
// 		if (x > (double)(GRIDSIZE) / 2.f)
// 			kS = -k;
// 		// Fourier gradient: multiply by -i k
// 		double real = htildehat[x].x;
// 		double imag = htildehat[x].y;
// 		htildehat[x].x = -kS * imag;
// 		htildehat[x].y = kS * real;
// 		// phase shift to translate function to cell boundaries
// 		real = htildehat[x].x;
// 		imag = htildehat[x].y;
// 		double beta = 0.5 * CELLSIZE * kS;
// 		htildehat[x].x = cos(beta) * real - sin(beta) * imag;
// 		htildehat[x].y = sin(beta) * real + cos(beta) * imag;
// 		for (int depth = 0; depth < DEPTH_NUM; depth++)
// 		{
// 			double k2 = std::max(0.0001, 2. * kx / GRIDSIZE);  //k2 = 0..1
// 			double omega = sqrtf(GRAVITY * k * tanhf(k * Depth[depth]));
// 			omega *= 1.f / sqrt(2.0 / (k2 * PI) * sin(k2 * PI / 2.0));  // grid dispersion correction
// 			qtildehat_depth[depth][x].x = qtildehat[x].x * cos(omega * TIMESTEP) - omega / (kNonZero * kNonZero) * htildehat[x].x * sin(omega * TIMESTEP);
// 			qtildehat_depth[depth][x].y = qtildehat[x].y * cos(omega * TIMESTEP) - omega / (kNonZero * kNonZero) * htildehat[x].y * sin(omega * TIMESTEP);
// 		}
// 	}
// 	for (int depth = 0; depth < DEPTH_NUM; depth++)
// 		fftc1dinv(qtildehat_depth[depth]); // Back transform
// 	// interpolate surface velocity from the two closest water depth solutions
// 	for (int x = 0; x < GRIDSIZE; x++)
// 	{
// 		float waterDepth = std::max(hbar[x], hbar[x_plus]);
// 		int depth1 = 0;
// 		for (int depth = 0; depth < DEPTH_NUM; depth++)
// 			if (waterDepth >= Depth[depth])
// 				depth1 = depth;
// 		int depth2 = std::min(DEPTH_NUM - 1, depth1 + 1);
// 		float s = 0.f;
// 		if (depth1 != depth2)
// 			s = (Depth[depth2] - waterDepth) / (Depth[depth2] - Depth[depth1]);
// 		qtilde[x] = s * qtildehat_depth[depth1][x].x + (1.f - s) * qtildehat_depth[depth2][x].x;
// 	}
// }

void Sim::SWEStep()
{
	// SWE bulk simulation using [Stelling03]

	// qbar to ubar using hbar from LAST timestep	
	gpu->Dispatch(shader_Ubar, 
		{qbar_x.srv, qbar_y.srv, hbarOld.srv,}, 
		{ubar_x.uav, ubar_y.uav}, 
		group, group);
	gpu->Dispatch(shader_Boundaries, {}, 
		{ubar_x.uav, ubar_y.uav}, 
		group, group);
		
	// Compute time derivative of u_bar and integrate to get new u_bar, then 
	// transfer back to flow rate using upwinding on *most recent* hbar
	gpu->Dispatch(shader_SWE,
		{ubar_x.srv, ubar_y.srv, hbar.srv, H.srv}, 
		{ubarNew_x.uav, ubarNew_y.uav, qbar_x.uav, qbar_y.uav}, 
		group, group);
	gpu->Dispatch(shader_Boundaries, {}, 
		{ubarNew_x.uav, ubarNew_y.uav, qbar_x.uav, qbar_y.uav}, 
		group, group);

	// store current hbar for next timestep
	std::swap(hbar, hbarOld);
}

void Sim::TransportStep()
{
	/****** Advect high-frequency wave height and flow rate through bulk velocity ******/ 

	// Adjust qtilde to account for advection by ubar, using cubic sampling to get better accuracy.
	std::swap(qtilde_x, qtildePast_x);  // store current qtilde for sampling, we will update qtilde in place
	std::swap(qtilde_y, qtildePast_y);
	gpu->Dispatch(shader_UpdateTilde, 
		{ubarNew_x.srv, ubar_x.srv, ubarNew_y.srv, ubar_y.srv, 
			qtildePast_x.srv, qtildePast_y.srv, h.srv, htilde.srv},
		{qtilde_x.uav, qtilde_y.uav, htilde.uav}, group, group);	
	gpu->Dispatch(shader_Boundaries, {}, 
		{qtilde_x.uav, qtilde_y.uav, htilde.uav}, group, group);
	
	// Advection of h through ubar:
	// Construct q_advect = ubar * htilde sampled at cell edges using cubic sampling
	gpu->Dispatch(shader_QAdvect, 
		{ubarNew_x.srv, ubarNew_y.srv, htilde.srv},
		{qAdvect_x.uav, qAdvect_y.uav}, group, group);
	// Use q_advect to update h using finite volume update: h_new = h_old + dt * (Del . (q + q_advect))
	std::swap(h, hPast);
	gpu->Dispatch(shader_HAdvect, 
		{qAdvect_x.srv, qAdvect_y.srv, hPast.srv, terrain.srv}, 
		{h.uav}, group, group);
	gpu->Dispatch(shader_Boundaries, {}, {h.uav}, group, group);
}
	
void Sim::ComputeValues()
{
	gpu->Dispatch(shader_IntegrateH, 
		{qbar_x.srv, qtilde_x.srv, qAdvect_x.srv, qbar_y.srv, qtilde_y.srv, qAdvect_y.srv, 
			hPast.srv, terrain.srv}, 
		{h.uav, q_x.uav, q_y.uav}, group, group);
	gpu->Dispatch(shader_Boundaries, {}, 
		{h.uav, q_x.uav, q_y.uav}, group, group);
}