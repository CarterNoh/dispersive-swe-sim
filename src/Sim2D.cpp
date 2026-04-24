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

void Sim::ExtractLocal(std::vector<float>& target, std::vector<float>& dataField, int index, bool isYDirection)
{
	if (isYDirection)
	{
		int x = index % GRIDSIZE;
		for (int y = 0; y < GRIDSIZE; y++)
			target[y] = dataField[idx(x,y)];
	}
	else
	{
		int y = floor(index / GRIDSIZE);
		for (int x = 0; x < GRIDSIZE; x++)
			target[x] = dataField[idx(x,y)];
	}
}

float Sim::SampleCubicClamped(float samplePos, std::vector<float>& dataField)
{
	// cubic interpolation with Catmull-Rom Spline https://en.wikipedia.org/wiki/Cubic_Hermite_spline#Interpolating_a_data_set
	int id_start = floor(samplePos) - 1;
	int id0 = std::max(0, std::min(id_start + 0, GRIDSIZE - 1));
	int id1 = std::max(0, std::min(id_start + 1, GRIDSIZE - 1));
	int id2 = std::max(0, std::min(id_start + 2, GRIDSIZE - 1));
	int id3 = std::max(0, std::min(id_start + 3, GRIDSIZE - 1));
	float fx = std::max(0.f, std::min(1.f, static_cast<float>(samplePos - floor(samplePos))));
	float x2 = fx * fx;
	float x3 = x2 * fx;
	float xcubicX = -x3 + 2.f * x2 - fx;
	float xcubicY =  3.f * x3 - 5.f * x2 + 2.f;
	float xcubicZ = -3.f * x3 + 4.f * x2 + fx;
	float xcubicW =  x3 - x2;
	float out = 0.5f * (xcubicX * dataField[id0] + xcubicY * dataField[id1] + xcubicZ * dataField[id2] + xcubicW * dataField[id3]);
	out = std::min(std::max(dataField[id1], dataField[id2]), out);  // value-limiting
	out = std::max(std::min(dataField[id1], dataField[id2]), out);  // value-limiting
	return out;
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
					waterSurface += 5.0f * cos(dist * PI * 5.0f);
            }
            else if (type == 1) { // Step/Dam Break
                if (xf < 0.3f) 
					waterSurface += 5.0f;
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
		fields[i]->resize(totalCells, 0.0);
	}
}

Sim::Sim(int terrainType = 0, int waterType = 0, float waterLevel = 5.0f)
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
	SetTerrain(terrainType);
	SetWater(waterType, WATER_LEVEL);
	// Set up GPU Shaders
	gpu = new GPU();
	if (!gpu->Init()) std::cerr << "GPU INIT FAILED" << std::endl;
	for (int i=0; i < sizeof(gpu_fields) / sizeof(gpu_fields[0]); i++) {
		gpu->CreateGridTexture(gpu_fields[i], GRIDSIZE);
		gpu->UploadToGPU(gpu_fields[i]->tex, *fields[i], GRIDSIZE);
	}
	for (int i=0; i < sizeof(shaders) / sizeof(shaders[0]); i++) {
		gpu->CompileShader(L"shaders/kernels.hlsl", names[i], shaders[i]);
	}
    gpu->UpdateConstants(constants);
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
	gpu->UploadToGPU(gpu_h.tex, h, GRIDSIZE);
	gpu->UploadToGPU(gpu_q_x.tex, q_x, GRIDSIZE);
	gpu->UploadToGPU(gpu_q_y.tex, q_y, GRIDSIZE);
    gpu->Dispatch(shader_InitDecomp, {gpu_h.srv, gpu_q_x.srv, gpu_q_y.srv, gpu_terrain.srv}, 
		{gpu_H.uav, gpu_Q_x.uav, gpu_Q_y.uav}, group, group);

	// Calculate diffusion coefficients
    gpu->Dispatch(shader_CalcDiff, {gpu_H.srv}, 
		{gpu_alpha_H.uav, gpu_alpha_Q_x.uav, gpu_alpha_Q_y.uav}, group, group);
	gpu->Dispatch(shader_Boundaries, {gpu_alpha_H.srv, gpu_alpha_Q_x.srv, gpu_alpha_Q_y.srv}, 
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
			{gpu_H.uav, gpu_Q_x.uav, gpu_Q_y.uav}, 
			group, group);
		gpu->Dispatch(shader_Boundaries, 
			{gpu_H.srv, gpu_Q_x.srv, gpu_Q_y.srv}, 
			{gpu_H.uav, gpu_Q_x.uav, gpu_Q_y.uav}, 
			group, group);
	}
	
	
	// final conversion to individual solver quantities
	gpu->DownloadFromGPU(gpu_H.tex, H, GRIDSIZE);
	gpu->DownloadFromGPU(gpu_Q_x.tex, Q_x, GRIDSIZE);
	gpu->DownloadFromGPU(gpu_Q_y.tex, Q_y, GRIDSIZE);
	for (int i = 0; i < GRIDSIZE*GRIDSIZE; i++)
	{
		hbar[i] = std::max(0.f, H[i] - terrain[i]);
		qbar_x[i] = Q_x[i];
		qbar_y[i] = Q_y[i];
		htilde[i] = h[i] - hbar[i];
		qtilde_x[i] = q_x[i] - qbar_x[i];
		qtilde_y[i] = q_y[i] - qbar_y[i];
	}
	// Enforce no-flow conditions at terrain boundaries
	for (int y = 0; y < GRIDSIZE-1; y++)
	{
		for (int x = 0; x < GRIDSIZE-1; x++)
		{
			if (StopFlowOnTerrainBoundary(x, y, h, terrain, false))  // stop flow in x direction
			{
				qbar_x[idx(x,y)] = 0.f;
				qtilde_x[idx(x,y)] = 0.f;
			}
			if (StopFlowOnTerrainBoundary(x, y, h, terrain, true))  // stop flow in y direction
			{
				qbar_y[idx(x,y)] = 0.f;
				qtilde_y[idx(x,y)] = 0.f;
			}
		}
	}
	ApplyBoundaries(qbar_x, BOUNDARY_TYPE);
	ApplyBoundaries(qbar_y, BOUNDARY_TYPE);
	ApplyBoundaries(qtilde_x, BOUNDARY_TYPE);
	ApplyBoundaries(qtilde_y, BOUNDARY_TYPE);

	// // final conversion to individual solver quantities
	// gpu->Dispatch(shader_Decompose, 
	// 	{gpu_H.srv, gpu_Q_x.srv, gpu_Q_y.srv, gpu_h.srv, gpu_q_x.srv, gpu_q_y.srv, gpu_terrain.srv}, 
	// 	{gpu_hbar.uav, gpu_qbar_x.uav, gpu_qbar_y.uav, gpu_htilde.uav, gpu_qtilde_x.uav, gpu_qtilde_y.uav}, 
	// 	group, group);
	// gpu->Dispatch(shader_StopFlow, 
	// 	{gpu_h.srv, gpu_terrain.srv}, 
	// 	{gpu_qbar_x.uav, gpu_qbar_y.uav, gpu_qtilde_x.uav, gpu_qtilde_y.uav}, 
	// 	group, group);
	// gpu->Dispatch(shader_Boundaries, 
	// 	{gpu_qbar_x.srv, gpu_qbar_y.srv, gpu_qtilde_x.srv, gpu_qtilde_y.srv}, 
	// 	{gpu_qbar_x.uav, gpu_qbar_y.uav, gpu_qtilde_x.uav, gpu_qtilde_y.uav}, 
	// 	group, group);	
	// gpu->DownloadFromGPU(gpu_H.tex, H, GRIDSIZE);
	// gpu->DownloadFromGPU(gpu_hbar.tex, hbar, GRIDSIZE);
	// gpu->DownloadFromGPU(gpu_qbar_x.tex, qbar_x, GRIDSIZE);
	// gpu->DownloadFromGPU(gpu_qbar_y.tex, qbar_y, GRIDSIZE);
	// gpu->DownloadFromGPU(gpu_htilde.tex, htilde, GRIDSIZE);
	// gpu->DownloadFromGPU(gpu_qtilde_x.tex, qtilde_x, GRIDSIZE);
	// gpu->DownloadFromGPU(gpu_qtilde_y.tex, qtilde_y, GRIDSIZE);
	
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

	// // qbar to ubar using hbar from LAST timestep
	// for (int y = 0; y < GRIDSIZE; y++)
	// {
	// 	for (int x = 0; x < GRIDSIZE; x++)
	// 	{
	// 		ubar_x[idx(x,y)] = qbar_x[idx(x,y)];
	// 		ubar_y[idx(x,y)] = qbar_y[idx(x,y)];

	// 		// First-Order Up-Winding
	// 		// SOMEDAY: Try interpolating h or using higher-order upwinding for better accuracy?
	// 		// Technically u = q / H, not h?? Different derivations differ here
	// 		if (ubar_x[idx(x,y)] >= 0.f || x == GRIDSIZE-1)
	// 			ubar_x[idx(x,y)] /= std::max(MIN_WATER_HEIGHT, hbarOld[idx(x,y)]);
	// 		else
	// 			ubar_x[idx(x,y)] /= std::max(MIN_WATER_HEIGHT, hbarOld[idx(x+1,y)]);
	// 		if (ubar_y[idx(x,y)] >= 0.f || y == GRIDSIZE-1)
	// 			ubar_y[idx(x,y)] /= std::max(MIN_WATER_HEIGHT, hbarOld[idx(x,y)]);
	// 		else
	// 			ubar_y[idx(x,y)] /= std::max(MIN_WATER_HEIGHT, hbarOld[idx(x,y+1)]);

	// 		// Enforcing CFL condition for later surface waves advection
	// 		ubar_x[idx(x,y)] = LimitVelocity(ubar_x[idx(x,y)]);  
	// 		ubar_y[idx(x,y)] = LimitVelocity(ubar_y[idx(x,y)]);  
	// 	}
	// }

	// qbar to ubar using hbar from LAST timestep
	// gpu->UploadToGPU(gpu_H.tex, H, GRIDSIZE);
	gpu->UploadToGPU(gpu_hbar.tex, hbar, GRIDSIZE);
	gpu->UploadToGPU(gpu_qbar_x.tex, qbar_x, GRIDSIZE);
	gpu->UploadToGPU(gpu_qbar_y.tex, qbar_y, GRIDSIZE);
	gpu->UploadToGPU(gpu_hbarOld.tex, hbarOld, GRIDSIZE);
	gpu->Dispatch(shader_Ubar, 
		{gpu_qbar_x.srv, gpu_qbar_y.srv, gpu_hbarOld.srv,}, 
		{gpu_ubar_x.uav, gpu_ubar_y.uav}, 
		group, group);	
	gpu->DownloadFromGPU(gpu_ubar_x.tex, ubar_x, GRIDSIZE);
	gpu->DownloadFromGPU(gpu_ubar_y.tex, ubar_y, GRIDSIZE);

	gpu->Dispatch(shader_SWE,
		{gpu_ubar_x.srv, gpu_ubar_y.srv, gpu_hbar.srv, gpu_H.srv}, 
		{gpu_ubarNew_x.uav, gpu_ubarNew_y.uav}, 
		group, group);
	gpu->Dispatch(shader_Boundaries, 
		{gpu_ubarNew_x.srv, gpu_ubarNew_y.srv}, 
		{gpu_ubarNew_x.uav, gpu_ubarNew_y.uav}, 
		group, group);
	gpu->DownloadFromGPU(gpu_ubarNew_x.tex, ubarNew_x, GRIDSIZE);
	gpu->DownloadFromGPU(gpu_ubarNew_y.tex, ubarNew_y, GRIDSIZE);

	// // Compute time derivative of u_bar and integrate to get new u_bar
	// for (int y = 1; y < GRIDSIZE-1; y++)
	// {
	// 	for (int x = 1; x < GRIDSIZE-1; x++)
	// 	{
	// 		// Compute intermediate values needed for du/dt calculations
	// 		// Need: q_x/y_ij, q_x_i-0.5_j, q_y_i_j-0.5, 
	// 		// 	     h_i+0.5_j, h_i_j+0.5, h_ij=hbar[idx(x,y)], h_i+1_j=hbar[idx(x+1,y)], h_i_j+1=hbar[idx(x,y+1)],
	// 		//       u_x_i+0.5_j=ubar_x[idx(x,y)], u_x_i-0.5_j=ubar_x[idx(x-1,y)], 
	// 		//       u_y_i_j+0.5=ubar_y[idx(x,y)], u_y_i_j-0.5=ubar_y[idx(x,y-1)], 
	// 		//       u_x_i+0.5_j-1=ubar_x[idx(x,y-1)], u_y_i-1_j+0.5=ubar_y[idx(x-1,y)]
	// 		// Note: The 1D implementation uses far more intermeidate values that the paper, and I don't know why. 
	// 		// The commented sections below are the conversion of their 1D code with all its complexity, and the 
	// 		// uncommented sections are my simplified version that tries to stay more true to the paper. We'll see 
	// 		// if it causes stability/accuracy issues.

	// 		////// X DIRECTION //////
	// 		// Use upwinding to evaluate q_(i-0.5,j), q_(i+0.5,j), q_(i+1.5,j) 
	// 		// for x direction to get q_x_(i,j), q_x_(i+1,j)
	// 		float q_x_m05 = ubar_x[idx(x-1,y)];
	// 		if (q_x_m05 >= 0.f)
	// 			q_x_m05 *= hbar[idx(x-1,y)];
	// 		else
	// 			q_x_m05 *= hbar[idx(x,y)];
	// 		float q_x_p05 = ubar_x[idx(x,y)];  
	// 		if (q_x_p05 >= 0.f)
	// 			q_x_p05 *= hbar[idx(x,y)];
	// 		else
	// 			q_x_p05 *= hbar[idx(x+1,y)];
	// 		float q_x_0 = 0.5f * (q_x_m05 + q_x_p05);
	// 		// float q_x_p15 = ubar_x[idx(x+1,y)];  //q_(i+1.5,j) = hfr at position x
	// 		// if (q_x_p15 >= 0.f)
	// 		// 	q_x_p15 *= hbar[idx(x+1,y)];
	// 		// else
	// 		// 	q_x_p15 *= hbar[std::min(idx(x+1,y) + 1, GRIDSIZE*GRIDSIZE-1)];
	// 		// float q_x_p1 = 0.5f * (q_x_p05 + q_x_p15);
	// 		// Calculate corresponding vaules for u_x_(i,j) using upwinding
	// 		// (why do we use upwinding here instead of averaging like q?)
	// 		// float u_star_x_0 = 0.f;
	// 		// if (q_x_0 >= 0.f)
	// 		// 	u_star_x_0 = ubar_x[idx(x-1,y)];
	// 		// else
	// 		// 	u_star_x_0 = ubar_x[idx(x,y)];
	// 		// float u_star_x_p1 = 0.f;
	// 		// if (q_x_p1 > 0.f)
	// 		// 	u_star_x_p1 = ubar_x[idx(x,y)];
	// 		// else
	// 		// 	u_star_x_p1 = ubar_x[idx(x+1,y)];

	// 		// Calculate h_(i+0.5,j) and h_(i-0.5,j) using upwinding
	// 		// float h_avg_x_p05 = (hbar[idx(x,y)] + hbar[idx(x+1,y)]) / 2.f; // averaging, for some reason the old paper used this
	// 		float h_x_p05 = 0.f;
	// 		if (ubar_x[idx(x,y)] >= 0.f)
	// 			h_x_p05 = hbar[idx(x,y)];
	// 		else
	// 			h_x_p05 = hbar[idx(x+1,y)];


	// 		/////// Y DIRECTION //////
	// 		// Use upwinding to evaluate q_(i,j-0.5), q_(i,j+0.5), q_(i,j+1.5) 
	// 		// for y direction to get q_y_(i,j), q_y_(i,j+1)
	// 		float q_y_m05 = ubar_y[idx(x,y-1)];
	// 		if (q_y_m05 >= 0.f)
	// 			q_y_m05 *= hbar[idx(x,y-1)];
	// 		else
	// 			q_y_m05 *= hbar[idx(x,y)];
	// 		float q_y_p05 = ubar_y[idx(x,y)];  
	// 		if (q_y_p05 >= 0.f)
	// 			q_y_p05 *= hbar[idx(x,y)];
	// 		else
	// 			q_y_p05 *= hbar[idx(x,y+1)];
	// 		float q_y_0 = 0.5f * (q_y_m05 + q_y_p05);

	// 		// Calculate h_(i,j+0.5) and h_(i,j-0.5) using upwinding
	// 		float h_y_p05 = 0.f;
	// 		if (ubar_y[idx(x,y)] >= 0.f)
	// 			h_y_p05 = hbar[idx(x,y)];
	// 		else
	// 			h_y_p05 = hbar[idx(x,y+1)];


	// 		// Compute dux_dt and duy_dt
	// 		// X DIRECTION
	// 		float dux_dt = - (1/CELLSIZE) * ((q_x_0/h_x_p05) * (ubar_x[idx(x,y)] - ubar_x[idx(x-1,y)]) + (q_y_m05/h_x_p05) * (ubar_x[idx(x,y)] - ubar_x[idx(x,y-1)]) + GRAVITY * (H[idx(x+1,y)] - H[idx(x,y)]));
	// 		ubarNew_x[idx(x,y)] = LimitVelocity(ubar_x[idx(x,y)] + TIMESTEP * dux_dt);  // Enforcing CFL condition
	// 		// Y DIRECTION
	// 		float duy_dt = - (1/CELLSIZE) * ((q_y_0/h_y_p05) * (ubar_y[idx(x,y)] - ubar_y[idx(x,y-1)]) + (q_x_m05/h_y_p05) * (ubar_y[idx(x,y)] - ubar_y[idx(x-1,y)]) + GRAVITY * (H[idx(x,y+1)] - H[idx(x,y)]));
	// 		ubarNew_y[idx(x,y)] = LimitVelocity(ubar_y[idx(x,y)] + TIMESTEP * duy_dt);  // Enforcing CFL condition
	// 	}
	// }
	// ApplyBoundaries(ubarNew_x, BOUNDARY_TYPE);
	// ApplyBoundaries(ubarNew_y, BOUNDARY_TYPE);

	// transfer back to flow rate using upwinding on *most recent* hbar
	for (int y = 0; y < GRIDSIZE; y++)
	{
		for (int x = 0; x < GRIDSIZE; x++)
		{
			if (ubarNew_x[idx(x,y)] >= 0.f || x == GRIDSIZE-1)
				qbar_x[idx(x,y)] = ubarNew_x[idx(x,y)] * hbar[idx(x,y)];
			else
				qbar_x[idx(x,y)] = ubarNew_x[idx(x,y)] * hbar[idx(x+1,y)];
			if (ubarNew_y[idx(x,y)] >= 0.f || y == GRIDSIZE-1)
				qbar_y[idx(x,y)] = ubarNew_y[idx(x,y)] * hbar[idx(x,y)];
			else
				qbar_y[idx(x,y)] = ubarNew_y[idx(x,y)] * hbar[idx(x,y+1)];
		}
	}

	hbarOld = hbar;   // store current hbar for next timestep
}

void Sim::TransportStep()
{
	// Advect high-frequency wave height and flow rate through bulk velocity

	// Adjust qtilde to account for advection by ubar, using cubic sampling to get better accuracy.
	// qtildePast_x = qtilde_x;  // store current qtilde for sampling, we will update qtilde in place
	// qtildePast_y = qtilde_y; 
	std::swap(qtilde_x, qtildePast_x);  // store current qtilde for sampling, we will update qtilde in place
	std::swap(qtilde_y, qtildePast_y);
	for (int y = 1; y < GRIDSIZE-1; y++)
	{
		for (int x = 1; x < GRIDSIZE-1; x++)
		{
			// Extract row and column for cubic sampling
			std::vector<float> local_x(GRIDSIZE);
			std::vector<float> local_y(GRIDSIZE);
			ExtractLocal(local_x, qtildePast_x, idx(x,y), false);
			ExtractLocal(local_y, qtildePast_y, idx(x,y), true);
			// Adjust qtilde to account for bulk flow
			float bulkVelocity_x = 0.5f * (ubarNew_x[idx(x,y)] + ubar_x[idx(x,y)]); // ubar is on same timestep as h, need to get back to timestep of q 
			float bulkVelocity_y = 0.5f * (ubarNew_y[idx(x,y)] + ubar_y[idx(x,y)]);
			float step_x = - bulkVelocity_x * TIMESTEP / CELLSIZE; // unitless (cells)
			float step_y = - bulkVelocity_y * TIMESTEP / CELLSIZE; // unitless (cells)
			qtilde_x[idx(x,y)] = SampleCubicClamped(x + step_x, local_x);
			qtilde_y[idx(x,y)] = SampleCubicClamped(y + step_y, local_y);
			if (((bulkVelocity_x >= 0.f) && (h[idx(x,y)] < MIN_WATER_HEIGHT)) ||
				((bulkVelocity_x < 0.f) && (h[idx(x+1,y)] < MIN_WATER_HEIGHT)))
				qtilde_x[idx(x,y)] = 0.f;
			if (((bulkVelocity_y >= 0.f) && (h[idx(x,y)] < MIN_WATER_HEIGHT)) ||
				((bulkVelocity_y < 0.f) && (h[idx(x,y+1)] < MIN_WATER_HEIGHT)))
				qtilde_y[idx(x,y)] = 0.f;
			
			// Update qtilde from ubar divergence: dq/dt = -q * div(ubar)
			float ubar_x_m1 = 0.5f * (ubarNew_x[idx(x-1,y)] + ubar_x[idx(x-1,y)]);
			float ubar_x_p1 = 0.5f * (ubarNew_x[idx(x+1,y)] + ubar_x[idx(x+1,y)]);
			float ubar_y_m1 = 0.5f * (ubarNew_y[idx(x,y-1)] + ubar_y[idx(x,y-1)]);
			float ubar_y_p1 = 0.5f * (ubarNew_y[idx(x,y+1)] + ubar_y[idx(x,y+1)]);
			// central difference to get divergence of ubar (using central because q and u are on same grid)
			float div_ubar_x = (ubar_x_p1 - ubar_x_m1) / (2.f * CELLSIZE);
			float div_ubar_y = (ubar_y_p1 - ubar_y_m1) / (2.f * CELLSIZE);
			float div_ubar = div_ubar_x + div_ubar_y;
			// dampen if converging to avoid breaking waves
			if (div_ubar < 0.f)
				div_ubar *= GAMMA_TRANSPORT;	
					
			qtilde_x[idx(x,y)] *= exp(-div_ubar * TIMESTEP);
			qtilde_y[idx(x,y)] *= exp(-div_ubar * TIMESTEP);

			// Update htilde from ubar divergence: dh/dt = -h * div(ubar)
			// backward difference to get divergence of ubar for h (using backward because u is on staggered grid from h)
			div_ubar_x = (ubarNew_x[idx(x,y)] - ubarNew_x[idx(x-1,y)]) / CELLSIZE;
			div_ubar_y = (ubarNew_y[idx(x,y)] - ubarNew_y[idx(x,y-1)]) / CELLSIZE;
			div_ubar = div_ubar_x + div_ubar_y;
			// dampen if converging to avoid breaking waves
			if (div_ubar < 0.f)
				div_ubar *= GAMMA_TRANSPORT;
			htilde[idx(x,y)] *= exp(-div_ubar * TIMESTEP);
		}
	}
	ApplyBoundaries(qtilde_x, BOUNDARY_TYPE);
	ApplyBoundaries(qtilde_y, BOUNDARY_TYPE);
	ApplyBoundaries(htilde, BOUNDARY_TYPE);	
	

	// Advection of h through ubar
	// First, construct q_advect = ubar * htilde sampled at cell edges using cubic sampling
	for (int y = 0; y < GRIDSIZE; y++)
	{
		for (int x = 0; x < GRIDSIZE; x++)
		{
			// extract row and column for cubic sampling
			std::vector<float> local_x(GRIDSIZE);
			std::vector<float> local_y(GRIDSIZE);
			ExtractLocal(local_x, htilde, idx(x,y), false);
			ExtractLocal(local_y, htilde, idx(x,y), true);

			// cubic reconstruction: 1/2 cell accounts for staggered grid, 0.5 * dt accounts for h and u at different 1/2 times
			float step_x = 0.5f - ubarNew_x[idx(x,y)] * (0.5f * TIMESTEP) / CELLSIZE;
			float step_y = 0.5f - ubarNew_y[idx(x,y)] * (0.5f * TIMESTEP) / CELLSIZE; 
			qAdvect_x[idx(x,y)] = ubarNew_x[idx(x,y)] * SampleCubicClamped(x + step_x, local_x);  
			qAdvect_y[idx(x,y)] = ubarNew_y[idx(x,y)] * SampleCubicClamped(y + step_y, local_y);  
		}
	}
	// Next, use q_advect to update h using finite volume update: h_new = h_old + dt * (Del . (q + q_advect))
	hPast = h;
	for (int y = 1; y < GRIDSIZE-1; y++)
	{
		for (int x = 1; x < GRIDSIZE-1; x++)
		{
			float q_x_l = LimitFlowRate(qAdvect_x[idx(x-1,y)], h[idx(x-1,y)], h[idx(x,y)]);
			float q_x_r = LimitFlowRate(qAdvect_x[idx(x,y)], h[idx(x,y)], h[idx(x+1,y)]);
			float q_y_l = LimitFlowRate(qAdvect_y[idx(x,y-1)], h[idx(x,y-1)], h[idx(x,y)]);
			float q_y_r = LimitFlowRate(qAdvect_y[idx(x,y)], h[idx(x,y)], h[idx(x,y+1)]);
			if ( ((h[idx(x-1,y)] == 0.f) && (h[idx(x,y)] == 0.f)) || (StopFlowOnTerrainBoundary(x-1, y, h, terrain, false)) )
				q_x_l = 0.f;
			if ( ((h[idx(x,y)] == 0.f) && (h[idx(x+1,y)] == 0.f)) || (StopFlowOnTerrainBoundary(x, y, h, terrain, false)) )
				q_x_r = 0.f;
			if ( ((h[idx(x,y-1)] == 0.f) && (h[idx(x,y)] == 0.f)) || (StopFlowOnTerrainBoundary(x, y-1, h, terrain, true)) )
				q_y_l = 0.f;
			if ( ((h[idx(x,y)] == 0.f) && (h[idx(x,y+1)] == 0.f)) || (StopFlowOnTerrainBoundary(x, y, h, terrain, true)) )
				q_y_r = 0.f;

			// update h using htilde advected through ubar
			float div_q = (q_x_l - q_x_r + q_y_l - q_y_r) / CELLSIZE;
			h[idx(x,y)] = std::max(0.f, hPast[idx(x,y)] + TIMESTEP * div_q);
		}
	}
	Sim::ApplyBoundaries(h, BOUNDARY_TYPE);
}
	
void Sim::ComputeValues()
{
	// Recombine bulk and surface flow
	for (int i = 0; i < GRIDSIZE*GRIDSIZE; i++)
	{
		q_x[i] = qbar_x[i] + qtilde_x[i];
		q_y[i] = qbar_y[i] + qtilde_y[i];
	}
	for (int y = 0; y < GRIDSIZE-1; y++)
	{
		for (int x = 0; x < GRIDSIZE-1; x++)
		{
			q_x[idx(x,y)] = LimitFlowRate(q_x[idx(x,y)], h[idx(x,y)], h[idx(std::min(x+1,GRIDSIZE-1),y)]);
			q_y[idx(x,y)] = LimitFlowRate(q_y[idx(x,y)], h[idx(x,y)], h[idx(x,std::min(y+1,GRIDSIZE-1))]);
			if ( (StopFlowOnTerrainBoundary(x, y, h, terrain, false)) || (x == 0) || (x >= GRIDSIZE - 2) )
				q_x[idx(x,y)] = 0.f;
			if ( (StopFlowOnTerrainBoundary(x, y, h, terrain, true)) || (y == 0) || (y >= GRIDSIZE - 2) )
				q_y[idx(x,y)] = 0.f;
			// NOTE: making edges a reflective boundary, revisit later
		}
	}
	ApplyBoundaries(q_x, BOUNDARY_TYPE);
	ApplyBoundaries(q_y, BOUNDARY_TYPE);
	
	// Height integration 
	for (int y = 1; y < GRIDSIZE; y++)
	{
		for (int x = 1; x < GRIDSIZE; x++)
		{
			float div_q = (q_x[idx(x-1,y)] - q_x[idx(x,y)] + q_y[idx(x,y-1)] - q_y[idx(x,y)]) / CELLSIZE;
			h[idx(x,y)] = std::max(0.f, h[idx(x,y)] + TIMESTEP * div_q);
		}
	}
	ApplyBoundaries(h, BOUNDARY_TYPE);

	// stability measure to not drag too much water from a cell in a single timestep (important for extreme initial conditions)
	for (int y = 0; y < GRIDSIZE - 1; y++)
	{
		for (int x = 0; x < GRIDSIZE - 1; x++)
		{
			q_x[idx(x,y)] = LimitFlowRate(q_x[idx(x,y)], h[idx(x,y)], h[idx(x+1,y)]);
			q_y[idx(x,y)] = LimitFlowRate(q_y[idx(x,y)], h[idx(x,y)], h[idx(x,y+1)]);
		}
	}
	ApplyBoundaries(q_x, BOUNDARY_TYPE);
	ApplyBoundaries(q_y, BOUNDARY_TYPE);
}