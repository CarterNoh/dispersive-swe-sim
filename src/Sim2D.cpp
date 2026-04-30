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

void Sim::SetTerrain(int type) {
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
}

void Sim::SetWater(int type, float level) {
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
            
            // Sync all buffers to prevent NaN/jitter on first frame
            hbar[i] = h[i];
            hbarOld[i] = h[i];
            q_x[i] = q_y[i] = qbar_x[i] = qbar_y[i] = 0.0f;
            htilde[i] = 0.0f;
        }
    }
}

Sim::Sim()
{
	// Define the total number of cells
    size_t totalCells = GRIDSIZE * GRIDSIZE;
	time = 0.f;
	// Allocate memory for all member vectors
	for (int i=0; i < sizeof(fields) / sizeof(fields[0]); i++) {
		// std::vector<float>& field = *fields[i];
		fields[i]->resize(totalCells, 0.0);
	}
    // Initialize terrain and water
	SetTerrain(TERRAIN_TYPE);
	SetWater(WATER_TYPE, WATER_LEVEL);

	// Set up GPU Shaders
	gpu = new GPU();
	if (!gpu->Init()) std::cerr << "GPU INIT FAILED" << std::endl;
	for (int i=0; i < sizeof(gpu_fields) / sizeof(gpu_fields[0]); i++) {
		gpu->CreateGridTexture(gpu_fields[i], GRIDSIZE);
		gpu->UploadToGPU(gpu_fields[i]->tex, *fields[i], GRIDSIZE);
	}
	for (int i=0; i < sizeof(compute_shaders) / sizeof(compute_shaders[0]); i++) {
		gpu->CompileComputeShader(L"shaders/kernels.hlsl", compute_names[i], compute_shaders[i]);
	}
    gpu->UpdateConstants(compute_constants);
	gpu->UploadToGPU(gpu_terrain.tex, terrain, GRIDSIZE);
}

Sim::Sim(HWND hwnd = nullptr)
{
	// Define the total number of cells
    size_t totalCells = GRIDSIZE * GRIDSIZE;
	time = 0.f;
	// Allocate memory for all member vectors
	for (int i=0; i < sizeof(fields) / sizeof(fields[0]); i++) {
		fields[i]->resize(totalCells, 0.0);
	}
    // Initialize terrain and water
	SetTerrain(TERRAIN_TYPE);
	SetWater(WATER_TYPE, WATER_LEVEL);

	// Set up GPU Shaders
	gpu = new GPU();
	if (!gpu->Init(hwnd)) std::cerr << "GPU INIT FAILED" << std::endl;
	for (int i=0; i < sizeof(gpu_fields) / sizeof(gpu_fields[0]); i++) {
		gpu->CreateGridTexture(gpu_fields[i], GRIDSIZE);
		gpu->UploadToGPU(gpu_fields[i]->tex, *fields[i], GRIDSIZE);
	}
	for (int i=0; i < sizeof(compute_shaders) / sizeof(compute_shaders[0]); i++) {
		gpu->CompileComputeShader(L"shaders/kernels.hlsl", compute_names[i], compute_shaders[i]);
	}
	gpu->CompileVertexShader(L"shaders/render.hlsl", "VSMain");
	gpu->CompilePixelShader(L"shaders/render.hlsl", "PSMain");
	gpu->CreateGridMesh(GRIDSIZE);
	gpu->CreateGridVertexBuffer(GRIDSIZE);
    gpu->UpdateConstants(compute_constants);
	gpu->UploadToGPU(gpu_terrain.tex, terrain, GRIDSIZE);
}

int Sim::Release(void)
{
	return 0;
}


// ********************************************************************************************************************
// Simulation functions
// ********************************************************************************************************************

void Sim::SimStep(bool SWEonly)
{
	DecompositionStep(SWEonly);
	// if (!SWEonly) {
	// eWaveStep();
	// 	FFTStep();
	// }
	SWEStep();
	TransportStep();
	ComputeValues();
	time += TIMESTEP;
}

void Sim::DecompositionStep(bool SWEonly)
{
	/******* Bulk vs Surface Wave Decomposition ******/

	if (SWEonly)
	{
		// If we're only doing SWE, then we skip the decomposition and just set the bulk values to be the total values (and surface values to 0)
		for (int i = 0; i < GRIDSIZE*GRIDSIZE; i++)
		{
			hbar[i] = h[i];
			qbar_x[i] = q_x[i];
			qbar_y[i] = q_y[i];
			htilde[i] = 0.f;
			qtilde_x[i] = 0.f;
			qtilde_y[i] = 0.f;
		}
		return;
	}

	// Initialize values for diffusion step
    gpu->Dispatch(shader_InitDecomp, {gpu_h.srv, gpu_q_x.srv, gpu_q_y.srv, gpu_terrain.srv}, 
		{gpu_H.uav, gpu_Q_x.uav, gpu_Q_y.uav}, group, group);

	// Calculate diffusion coefficients
    gpu->Dispatch(shader_CalcDiff, {gpu_H.srv}, 
		{gpu_alpha_H.uav, gpu_alpha_Q_x.uav, gpu_alpha_Q_y.uav}, group, group);
	gpu->Dispatch(shader_Boundaries, {}, 
		{gpu_alpha_H.uav, gpu_alpha_Q_x.uav, gpu_alpha_Q_y.uav}, group, group);
	
	// Run diffusion to low-pass filter H and Q
	for (int j = 0; (j < DIFFUSION_ITERATIONS); j++)
	{
		// Swap H and HPast pointers for ping-ponging
        std::swap(gpu_H, gpu_HPast); 
        std::swap(gpu_Q_x, gpu_QPast_x);
		std::swap(gpu_Q_y, gpu_QPast_y);

		// Diffusion step for H and Q
		gpu->Dispatch(shader_Diffusion, 
			{gpu_terrain.srv, gpu_HPast.srv, gpu_QPast_x.srv, gpu_QPast_y.srv, gpu_alpha_H.srv, gpu_alpha_Q_x.srv, gpu_alpha_Q_y.srv}, 
			{gpu_H.uav, gpu_Q_x.uav, gpu_Q_y.uav}, group, group);
		gpu->Dispatch(shader_Boundaries, {}, 
			{gpu_H.uav, gpu_Q_x.uav, gpu_Q_y.uav}, group, group);
	}
	// gpu->Dispatch(shader_Boundaries,	{}, 
	// 	{gpu_H.uav, gpu_Q_x.uav, gpu_Q_y.uav}, group, group);

	// final conversion to individual solver quantities
	gpu->Dispatch(shader_Decompose, 
		{gpu_H.srv, gpu_Q_x.srv, gpu_Q_y.srv, gpu_h.srv, gpu_q_x.srv, gpu_q_y.srv, gpu_terrain.srv}, 
		{gpu_hbar.uav, gpu_qbar_x.uav, gpu_qbar_y.uav, gpu_htilde.uav, gpu_qtilde_x.uav, gpu_qtilde_y.uav}, 
		group, group);
	gpu->Dispatch(shader_Boundaries, {}, 
		{gpu_hbar.uav, gpu_htilde.uav}, group, group);
	gpu->Dispatch(shader_Boundaries, {}, 
		{gpu_qbar_x.uav, gpu_qbar_y.uav, gpu_qtilde_x.uav, gpu_qtilde_y.uav}, 
		group, group);
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
		{gpu_qbar_x.srv, gpu_qbar_y.srv, gpu_hbarOld.srv,}, 
		{gpu_ubar_x.uav, gpu_ubar_y.uav}, 
		group, group);
	gpu->Dispatch(shader_Boundaries, {}, 
		{gpu_ubar_x.uav, gpu_ubar_y.uav}, 
		group, group);
		
	// Compute time derivative of u_bar and integrate to get new u_bar, then 
	// transfer back to flow rate using upwinding on *most recent* hbar
	gpu->Dispatch(shader_SWE,
		{gpu_ubar_x.srv, gpu_ubar_y.srv, gpu_hbar.srv, gpu_H.srv}, 
		{gpu_ubarNew_x.uav, gpu_ubarNew_y.uav, gpu_qbar_x.uav, gpu_qbar_y.uav}, 
		group, group);
	gpu->Dispatch(shader_Boundaries, {}, 
		{gpu_ubarNew_x.uav, gpu_ubarNew_y.uav, gpu_qbar_x.uav, gpu_qbar_y.uav}, 
		group, group);

	// store current hbar for next timestep
	std::swap(gpu_hbar, gpu_hbarOld);
}

void Sim::TransportStep()
{
	// Advect high-frequency wave height and flow rate through bulk velocity

	// Adjust qtilde to account for advection by ubar, using cubic sampling to get better accuracy.
	std::swap(gpu_qtilde_x, gpu_qtildePast_x);  // store current qtilde for sampling, we will update qtilde in place
	std::swap(gpu_qtilde_y, gpu_qtildePast_y);
	gpu->Dispatch(shader_UpdateTilde, 
		{gpu_ubarNew_x.srv, gpu_ubar_x.srv, gpu_ubarNew_y.srv, gpu_ubar_y.srv, 
			gpu_qtildePast_x.srv, gpu_qtildePast_y.srv, gpu_h.srv, gpu_htilde.srv},
		{gpu_qtilde_x.uav, gpu_qtilde_y.uav, gpu_htilde.uav}, 
		group, group);	
	gpu->Dispatch(shader_Boundaries, {}, 
		{gpu_qtilde_x.uav, gpu_qtilde_y.uav, gpu_htilde.uav}, 
		group, group);
	
	// Advection of h through ubar:
	// Construct q_advect = ubar * htilde sampled at cell edges using cubic sampling
	gpu->Dispatch(shader_QAdvect, 
		{gpu_ubarNew_x.srv, gpu_ubarNew_y.srv, gpu_htilde.srv},
		{gpu_qAdvect_x.uav, gpu_qAdvect_y.uav}, 
		group, group);
	// Use q_advect to update h using finite volume update: h_new = h_old + dt * (Del . (q + q_advect))
	std::swap(gpu_h, gpu_hPast);
	gpu->Dispatch(shader_HAdvect, 
		{gpu_qAdvect_x.srv, gpu_qAdvect_y.srv, gpu_hPast.srv, gpu_terrain.srv},
		{gpu_h.uav}, 
		group, group);
	gpu->Dispatch(shader_Boundaries, {}, {gpu_h.uav}, group, group);
}
	
void Sim::ComputeValues()
{
	gpu->Dispatch(shader_IntegrateH, 
		{gpu_qbar_x.srv, gpu_qtilde_x.srv, gpu_qAdvect_x.srv, gpu_qbar_y.srv, gpu_qtilde_y.srv, gpu_qAdvect_y.srv, 
			gpu_hPast.srv, gpu_terrain.srv}, 
		{gpu_h.uav, gpu_q_x.uav, gpu_q_y.uav}, 
		group, group);
	gpu->Dispatch(shader_Boundaries, {}, 
		{gpu_h.uav, gpu_q_x.uav, gpu_q_y.uav}, 
		group, group);
	// gpu->DownloadFromGPU(gpu_h.tex, h, GRIDSIZE);



	// gpu->DownloadFromGPU(gpu_h.tex, h, GRIDSIZE);
	// gpu->DownloadFromGPU(gpu_qbar_x.tex, qbar_x, GRIDSIZE);
	// gpu->DownloadFromGPU(gpu_qbar_y.tex, qbar_y, GRIDSIZE);
	// gpu->DownloadFromGPU(gpu_qtilde_x.tex, qtilde_x, GRIDSIZE);
	// gpu->DownloadFromGPU(gpu_qtilde_y.tex, qtilde_y, GRIDSIZE);
	// gpu->DownloadFromGPU(gpu_qAdvect_x.tex, qAdvect_x, GRIDSIZE);
	// gpu->DownloadFromGPU(gpu_qAdvect_y.tex, qAdvect_y, GRIDSIZE);

	// // Recombine bulk and surface flow
	// for (int i = 0; i < GRIDSIZE*GRIDSIZE; i++)
	// {
	// 	q_x[i] = qbar_x[i] + qtilde_x[i];
	// 	q_y[i] = qbar_y[i] + qtilde_y[i];
	// }
	// for (int y = 0; y < GRIDSIZE-1; y++)
	// {
	// 	for (int x = 0; x < GRIDSIZE-1; x++)
	// 	{
	// 		q_x[idx(x,y)] = LimitFlowRate(q_x[idx(x,y)], h[idx(x,y)], h[idx(std::min(x+1,GRIDSIZE-1),y)]);
	// 		q_y[idx(x,y)] = LimitFlowRate(q_y[idx(x,y)], h[idx(x,y)], h[idx(x,std::min(y+1,GRIDSIZE-1))]);
	// 		if ( (StopFlowOnTerrainBoundary(x, y, h, terrain, false)) || (x == 0) || (x >= GRIDSIZE - 2) )
	// 			q_x[idx(x,y)] = 0.f;
	// 		if ( (StopFlowOnTerrainBoundary(x, y, h, terrain, true)) || (y == 0) || (y >= GRIDSIZE - 2) )
	// 			q_y[idx(x,y)] = 0.f;
	// 		// NOTE: making edges a reflective boundary, revisit later
	// 	}
	// }
	// ApplyBoundaries(q_x, BOUNDARY_TYPE);
	// ApplyBoundaries(q_y, BOUNDARY_TYPE);
	
	// // Height integration 
	// for (int y = 1; y < GRIDSIZE; y++)
	// {
	// 	for (int x = 1; x < GRIDSIZE; x++)
	// 	{
	// 		float div_q = (q_x[idx(x,y)] - q_x[idx(x-1,y)] + q_y[idx(x,y)] - q_y[idx(x,y-1)]) / CELLSIZE;
	// 		h[idx(x,y)] = std::max(0.f, h[idx(x,y)] - TIMESTEP * div_q);
	// 	}
	// }
	// ApplyBoundaries(h, BOUNDARY_TYPE);

	// for (int y = 1; y < GRIDSIZE-1; y++)
	// {
	// 	for (int x = 1; x < GRIDSIZE-1; x++)
	// 	{

	// 		float fq_x = qbar_x[idx(x,y)] + qtilde_x[idx(x,y)] + qAdvect_x[idx(x,y)];
	// 		float fq_xm = qbar_x[idx(x-1,y)] + qtilde_x[idx(x-1,y)] + qAdvect_x[idx(x-1,y)];
	// 		float fq_y = qbar_y[idx(x,y)] + qtilde_y[idx(x,y)] + qAdvect_y[idx(x,y)];
	// 		float fq_ym = qbar_y[idx(x,y-1)] + qtilde_y[idx(x,y-1)] + qAdvect_y[idx(x,y-1)];

	// 		fq_x  = LimitFlowRate(fq_x, h[idx(x,y)], h[idx(x+1,y)]);
	// 		fq_xm = LimitFlowRate(fq_xm, h[idx(x-1,y)], h[idx(x,y)]);
	// 		fq_y  = LimitFlowRate(fq_y, h[idx(x,y)], h[idx(x,y+1)]);
	// 		fq_ym = LimitFlowRate(fq_ym, h[idx(x,y-1)], h[idx(x,y)]);

	// 		// if ((StopFlowOnTerrainBoundary(in4, in5, id.xy, false)) || (id.x == 0) || (id.x == gridSize - 2))
	// 		//     out0[id.xy] = 0.f;
	// 		// else
	// 		//     out0[id.xy] = LimitFlowRate(q_x, in4[id.xy], in4[id.xy + uint2(1, 0)]);
			
	// 		// if ((StopFlowOnTerrainBoundary(in4, in5, id.xy, true)) || (id.y == 0) || (id.y == gridSize - 2))
	// 		//     out1[id.xy] = 0.f;
	// 		// else        
	// 		//     out1[id.xy] = LimitFlowRate(q_y, in4[id.xy], in4[id.xy + uint2(0, 1)]);

	// 		float div_q = (fq_x - fq_xm + fq_y - fq_ym) / CELLSIZE;
	// 		h[idx(x,y)] = std::max(0.f, h[idx(x,y)] - TIMESTEP * div_q);
	// 		q_x[idx(x,y)] = fq_x;// - qAdvect_x[idx(x,y)]; // qbar + qtilde, removing qAdvect
	// 		q_y[idx(x,y)] = fq_y;// - qAdvect_y[idx(x,y)];

	// 		// q_x[idx(x,y)] = LimitFlowRate(q_x[idx(x,y)], h[idx(x,y)], h[idx(std::min(x+1,GRIDSIZE-1),y)]);
	// 		// q_y[idx(x,y)] = LimitFlowRate(q_y[idx(x,y)], h[idx(x,y)], h[idx(x,std::min(y+1,GRIDSIZE-1))]);
	// 		// if ( (StopFlowOnTerrainBoundary(x, y, h, terrain, false)) || (x == 0) || (x >= GRIDSIZE - 2) )
	// 		// 	q_x[idx(x,y)] = 0.f;
	// 		// if ( (StopFlowOnTerrainBoundary(x, y, h, terrain, true)) || (y == 0) || (y >= GRIDSIZE - 2) )
	// 		// 	q_y[idx(x,y)] = 0.f;
	// 		// NOTE: making edges a reflective boundary, revisit later
	// 	}
	// }
	// ApplyBoundaries(h, BOUNDARY_TYPE);	
	// ApplyBoundaries(q_x, BOUNDARY_TYPE);
	// ApplyBoundaries(q_y, BOUNDARY_TYPE);


	// // stability measure to not drag too much water from a cell in a single timestep (important for extreme initial conditions)
	// for (int y = 0; y < GRIDSIZE - 1; y++)
	// {
	// 	for (int x = 0; x < GRIDSIZE - 1; x++)
	// 	{
	// 		q_x[idx(x,y)] = LimitFlowRate(q_x[idx(x,y)], h[idx(x,y)], h[idx(x+1,y)]);
	// 		q_y[idx(x,y)] = LimitFlowRate(q_y[idx(x,y)], h[idx(x,y)], h[idx(x,y+1)]);
	// 	}
	// }
	// ApplyBoundaries(q_x, BOUNDARY_TYPE);
	// ApplyBoundaries(q_y, BOUNDARY_TYPE);

	// gpu->UploadToGPU(gpu_h.tex, h, GRIDSIZE);
	// gpu->UploadToGPU(gpu_q_x.tex, q_x, GRIDSIZE);
	// gpu->UploadToGPU(gpu_q_y.tex, q_y, GRIDSIZE);

}