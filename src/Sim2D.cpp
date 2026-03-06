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

inline int idx(int x, int y)
{
	return y * GRIDSIZE + x;
}

inline float LimitFlowRate(float flow_rate_in, float waterDepth_left, float waterDepth_right)
{
	if (flow_rate_in >= 0.f)
		return min(flow_rate_in, CFL_CONDITION * waterDepth_left * CELLSIZE / TIMESTEP);  // 0.25 since other neighbor might take from this source cell as well
	else
		return max(flow_rate_in, -CFL_CONDITION * waterDepth_right * CELLSIZE / TIMESTEP);
}

inline float LimitVelocity(float velocity_in)
{
	if (velocity_in >= 0.f)
		return min(velocity_in, CFL_CONDITION * CELLSIZE / TIMESTEP);   // 0.25 since other neighbors might take from this source cell as well
	else
		return max(velocity_in, -CFL_CONDITION * CELLSIZE / TIMESTEP);
}

bool StopFlowOnTerrainBoundary(int x, int y, std::vector<float>& h, std::vector<float>& terrain, bool isYDirection = false)
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

float ExtractLocal(std::vector<float>& target, std::vector<float>& dataField, int index, bool isYDirection)
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

float SampleCubicClamped(float samplePos, std::vector<float>& dataField)
{
	// cubic interpolation with Catmull-Rom Spline https://en.wikipedia.org/wiki/Cubic_Hermite_spline#Interpolating_a_data_set
	int id_start = floor(samplePos) - 1;
	int id0 = max(0, min(id_start + 0, GRIDSIZE - 1));
	int id1 = max(0, min(id_start + 1, GRIDSIZE - 1));
	int id2 = max(0, min(id_start + 2, GRIDSIZE - 1));
	int id3 = max(0, min(id_start + 3, GRIDSIZE - 1));
	float fx = max(0.f, min(1.f, samplePos - floor(samplePos)));
	float x2 = fx * fx;
	float x3 = x2 * fx;
	// pretty sure theirs is wrong, see wiki page formula
	// const float s = 0.5f;
	// float xcubicX = -s * x3 + 2.f * s * x2 - s * fx;
	// float xcubicY = (2.f - s) * x3 + (s - 3.f) * x2 + 1.f;
	// float xcubicZ = (s - 2.f) * x3 + (3.f - 2.f * s) * x2 + s * fx;
	// float xcubicW = s * x3 - s * x2;
	// float out = xcubicX * dataField[id0] + xcubicY * dataField[id1] + xcubicZ * dataField[id2] + xcubicW * dataField[id3];
	float xcubicX = -x3 + 2.f * x2 - fx;
	float xcubicY =  3.f * x3 - 5.f * x2 + 2.f;
	float xcubicZ = -3.f * x3 + 4.f * x2 + fx;
	float xcubicW =  x3 - x2;
	float out = 0.5f * (xcubicX * dataField[id0] + xcubicY * dataField[id1] + xcubicZ * dataField[id2] + xcubicW * dataField[id3]);
	out = min(max(dataField[id1], dataField[id2]), out);  // value-limiting
	out = max(min(dataField[id1], dataField[id2]), out);  // value-limiting
	return out;
}

// Boundary Condition Functions
void HandleWallBoundary(std::vector<float>& field) 
{
	// Handle Top and Bottom Edges (Horizontal)
	for (int i = 0; i < GRIDSIZE; ++i) {
		field[idx(i, 0)] = field[idx(i, 1)]; // Top Edge
		field[idx(i, GRIDSIZE-1)] = field[idx(i, GRIDSIZE-2)]; // Bottom Edge
		field[idx(0, i)] = field[idx(1, i)]; // Left Edge
		field[idx(GRIDSIZE-1, i)] = field[idx(GRIDSIZE-2, i)]; // Right Edge
	}
}

void HandleFreeBoundary(std::vector<float>& field)
{
	// Free boundary: extrapolate values from interior to boundaries using linear extrapolation
	for (int i = 0; i < GRIDSIZE; ++i) {
		field[idx(i, 0)] = 2.f * field[idx(i, 1)] - field[idx(i, 2)]; // Top Edge
		field[idx(i, GRIDSIZE-1)] = 2.f * field[idx(i, GRIDSIZE-2)] - field[idx(i, GRIDSIZE-3)]; // Bottom Edge
		field[idx(0, i)] = 2.f * field[idx(1, i)] - field[idx(2, i)]; // Left Edge
		field[idx(GRIDSIZE-1, i)] = 2.f * field[idx(GRIDSIZE-2, i)] - field[idx(GRIDSIZE-3, i)]; // Right Edge
	}
}

void HandleZeroBoundary(std::vector<float>& field)
{
	for (int i = 0; i < GRIDSIZE; ++i) {
		field[idx(i, 0)] = 0.f; // Top Edge
		field[idx(i, GRIDSIZE - 1)] = 0.f; // Bottom Edge
		field[idx(0, i)] = 0.f; // Left Edge
		field[idx(GRIDSIZE - 1, i)] = 0.f; // Right Edge
	}
}

void ApplyBoundaries(std::vector<float>& field, int type)
{
	if (type == 0)
		HandleWallBoundary(field);
	else if (type == 1)
		HandleFreeBoundary(field);
	else if (type == 2)
		HandleZeroBoundary(field);
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
					waterSurface += 10.0f;
            }
            else if (type == 2) { // Basin Flood
                // Fill only the left basin (xf < 0.5)
				if (xf < 0.25f)
                    waterSurface = max(level, terrain[idx(GRIDSIZE/2, y)] + 2.0f);
                else if (xf < 0.5f) 
					waterSurface = max(level, terrain[idx(GRIDSIZE/2, y)]);
                else 
					waterSurface = terrain[i]; // Start dry
            }

            // Standardize height as (Surface - Terrain)
            h[i] = max(0.0f, waterSurface - (float)terrain[i]);
            
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
	time = 0.0f;

    // Allocate memory for all member vectors
    terrain.resize(totalCells, 0.0);
    h.resize(totalCells, 0.0);
    q_x.resize(totalCells, 0.0);
    q_y.resize(totalCells, 0.0);
    hbarOld.resize(totalCells, 0.0);
    htildeOld.resize(totalCells, 0.0);
    hbar.resize(totalCells, 0.0);
    qbar_x.resize(totalCells, 0.0);
    qbar_y.resize(totalCells, 0.0);
	htilde.resize(totalCells, 0.0);
    qtilde_x.resize(totalCells, 0.0);
    qtilde_y.resize(totalCells, 0.0);
    ubar_x.resize(totalCells, 0.0);
    ubar_y.resize(totalCells, 0.0);
    ubarNew_x.resize(totalCells, 0.0);
    ubarNew_y.resize(totalCells, 0.0);
	alpha_H.resize(totalCells, 0.0);
	alpha_Q_x.resize(totalCells, 0.0);
	alpha_Q_y.resize(totalCells, 0.0);
	H.resize(totalCells, 0.0);
	Q_x.resize(totalCells, 0.0);
	Q_y.resize(totalCells, 0.0);
	HPast.resize(totalCells, 0.0);
	QPast_x.resize(totalCells, 0.0);
	QPast_y.resize(totalCells, 0.0);
	qtildePast_x.resize(totalCells, 0.0);
	qtildePast_y.resize(totalCells, 0.0);
	qAdvect_x.resize(totalCells, 0.0);
	qAdvect_y.resize(totalCells, 0.0);
	hPast.resize(totalCells, 0.0);
}

Sim::Sim(int terrainType = 0, int waterType = 0, float waterLevel = 5.0f)
{
	// Define the total number of cells
    size_t totalCells = GRIDSIZE * GRIDSIZE;
	time = 0.0f;

    // Allocate memory for all member vectors
    terrain.resize(totalCells, 0.0);
    h.resize(totalCells, 0.0);
    q_x.resize(totalCells, 0.0);
    q_y.resize(totalCells, 0.0);
    hbarOld.resize(totalCells, 0.0);
    htildeOld.resize(totalCells, 0.0);
    hbar.resize(totalCells, 0.0);
    qbar_x.resize(totalCells, 0.0);
    qbar_y.resize(totalCells, 0.0);
	htilde.resize(totalCells, 0.0);
    qtilde_x.resize(totalCells, 0.0);
    qtilde_y.resize(totalCells, 0.0);
    ubar_x.resize(totalCells, 0.0);
    ubar_y.resize(totalCells, 0.0);
    ubarNew_x.resize(totalCells, 0.0);
    ubarNew_y.resize(totalCells, 0.0);
	alpha_H.resize(totalCells, 0.0);
	alpha_Q_x.resize(totalCells, 0.0);
	alpha_Q_y.resize(totalCells, 0.0);
	H.resize(totalCells, 0.0);
	Q_x.resize(totalCells, 0.0);
	Q_y.resize(totalCells, 0.0);
	HPast.resize(totalCells, 0.0);
	QPast_x.resize(totalCells, 0.0);
	QPast_y.resize(totalCells, 0.0);
	qtildePast_x.resize(totalCells, 0.0);
	qtildePast_y.resize(totalCells, 0.0);
	qAdvect_x.resize(totalCells, 0.0);
	qAdvect_y.resize(totalCells, 0.0);
	hPast.resize(totalCells, 0.0);

    
    // Initialize terrain and water to default states (flat terrain, no water)
	SetTerrain(terrainType);
	SetWater(waterType, waterLevel);
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
	// if (!SWEonly)
		// eWaveStep();
		// FFTStep();
	SWEStep();
	TransportStep();
	ComputeValues();
	time += TIMESTEP;
}

void Sim::DecompositionStep(bool SWEonly)
{
	/******* Bulk vs Surface Wave Decomposition ******/

	// if (SWEonly)
	// {
	// 	// If we're only doing SWE, then we skip the decomposition and just set the bulk values to be the total values (and surface values to 0)
	// 	for (int i = 0; i < GRIDSIZE*GRIDSIZE; i++)
	// 	{
	// 		hbar[i] = h[i];
	// 		qbar_x[i] = q_x[i];
	// 		qbar_y[i] = q_y[i];
	// 		htilde[i] = 0.f;
	// 		qtilde_x[i] = 0.f;
	// 		qtilde_y[i] = 0.f;
	// 	}
	// 	return;
	// }

	// Calculate diffusion coefficient (alpha) at every location
	for (int i = 0; i < GRIDSIZE*GRIDSIZE; i++)
	{
		alpha_H[i] = 0.f;
		alpha_Q_x[i] = 0.f;
		alpha_Q_y[i] = 0.f;
		H[i] = terrain[i] + h[i];
		Q_x[i] = q_x[i];
		Q_y[i] = q_y[i];
	}

	/*	NOTE: Their implementation uses an average height across neighbor cells (I assume to improve stability, 
	but reduces accuracy?). I am choosing to use only local cell values to align with eqn from paper, we'll 
	see how it does. 

	Their code for initial calculation:
	// Identify the correct height (sigma) to use for diffusivity calculation
	alpha_H[idx(x,y)] = 0.f;
	float maxGround = max(terrain[idx(x,y)], terrain[idx(x+1,y)], terrain[idx(x,y+1)]); // Why do this?
	float minWaterlevel = (H[idx(x,y)] + H[idx(x+1,y)] + H[idx(x,y+1)]) / 3.f; // Why average here?
	if ((h[idx(x,y)] > 0.f) && (h[idx(x+1,y)] > 0.f) && (h[idx(x,y+1)] > 0.f))
	{
		static const float sigma_max = 8.f;
		// they limit diffusion coefficient to between 0 and 1, maybe for stability?
		float sigma = min(sigma_max, max(0.f, minWaterlevel - maxGround));
		alpha_H[idx(x,y)] = sigma * sigma / (2*DELTA_T*DIFFUSION_ITERATIONS);
	} 
	*/

	// Loop through main grid to calculate diffusion coefficients
	for (int y = 0; y < GRIDSIZE-1; y++)
	{
		for (int x = 0; x < GRIDSIZE-1; x++)
		{
			// Alpha_H
			float denom = 2*DELTA_T*DIFFUSION_ITERATIONS;
			alpha_H[idx(x,y)] = h[idx(x,y)] * h[idx(x,y)] / denom;
			// Penalize steep gradients
			// NOTE: they used H, I switched it to h to stay strict with the paper. We'll see if this causes bugs.
			float gradient_x = (h[idx(x+1,y)] - h[idx(x,y)]) / CELLSIZE; // could use central difference here
			float gradient_y = (h[idx(x,y+1)] - h[idx(x,y)]) / CELLSIZE;
			float penalty = - DIFFUSION_PENALTY * (gradient_x * gradient_x + gradient_y * gradient_y);
			alpha_H[idx(x,y)] *= exp(penalty);

			// Alpha_Q
			float avg_h_x = 0.5f * (h[idx(x,y)] + h[idx(x+1,y)]);
			float avg_h_y = 0.5f * (h[idx(x,y)] + h[idx(x,y+1)]);
			alpha_Q_x[idx(x,y)] = avg_h_x * avg_h_x / denom;
			alpha_Q_y[idx(x,y)] = avg_h_y * avg_h_y / denom;
			alpha_Q_x[idx(x,y)] *= exp(penalty);
			alpha_Q_y[idx(x,y)] *= exp(penalty);
		}
	}
	// Handle boundaries
	for (int i = 0; i < GRIDSIZE; i++)
	{
		float denom = 2*DELTA_T*DIFFUSION_ITERATIONS;
		// Right Edge
		alpha_H[idx(GRIDSIZE-1,i)] = h[idx(GRIDSIZE-1,i)] * h[idx(GRIDSIZE-1,i)] / denom; // Right Edge
		float grad_x = (h[idx(GRIDSIZE-1,i)] - h[idx(GRIDSIZE-2,i)]) / CELLSIZE;
		float grad_y = (h[idx(GRIDSIZE-1, min(i+1, GRIDSIZE-1))] - h[idx(GRIDSIZE-1,i)]) / CELLSIZE;
		float penalty = - DIFFUSION_PENALTY * (grad_x * grad_x + grad_y * grad_y);
		alpha_H[idx(GRIDSIZE-1,i)] *= exp(penalty);
		float h_next_y = h[idx(GRIDSIZE-1, min(i+1, GRIDSIZE-1))];
		float avg_h_y = 0.5f * (h[idx(GRIDSIZE-1,i)] + h_next_y);
		alpha_Q_x[idx(GRIDSIZE-1,i)] = alpha_H[idx(GRIDSIZE-1,i)];
		alpha_Q_y[idx(GRIDSIZE-1,i)] = avg_h_y * avg_h_y / denom * exp(penalty);
		// Top Edge
		alpha_H[idx(i,GRIDSIZE-1)] = h[idx(i,GRIDSIZE-1)] * h[idx(i,GRIDSIZE-1)] / denom;
		grad_x = (h[idx(min(i+1, GRIDSIZE-1),GRIDSIZE-1)] - h[idx(i,GRIDSIZE-1)]) / CELLSIZE;
		grad_y = (h[idx(i,GRIDSIZE-1)] - h[idx(i,GRIDSIZE-2)]) / CELLSIZE;
		penalty = - DIFFUSION_PENALTY * (grad_x * grad_x + grad_y * grad_y);
		alpha_H[idx(i, GRIDSIZE-1)] *= exp(penalty);
		float h_next_x = h[idx(min(i+1, GRIDSIZE-1),GRIDSIZE-1)];
		float avg_h_x = 0.5f * (h[idx(i,GRIDSIZE-1)] + h_next_x);
		alpha_Q_x[idx(i,GRIDSIZE-1)] = avg_h_x * avg_h_x / denom * exp(penalty);
		alpha_Q_y[idx(i,GRIDSIZE-1)] = alpha_H[idx(i,GRIDSIZE-1)];
	}
	
	// Run diffusion to low-pass filter H and Q
	// SOMEDAY: Improve this implementation of diffusion by replacing Euler integration with something better
	for (int j = 0; (j < DIFFUSION_ITERATIONS); j++)
	{
		HPast = H;
		QPast_x = Q_x;
		QPast_y = Q_y;
		for (int y = 1; y < GRIDSIZE-1; y++) // one diffusion iteration
		{
			for (int x = 1; x < GRIDSIZE - 1; x++)
			{
				// Diffusion step for H: dH/dt = Del * ( alpha_H * Del H )
				float dH_x = alpha_H[idx(x,y)] * (HPast[idx(x+1,y)] - HPast[idx(x,y)]) - alpha_H[idx(x-1,y)] * (HPast[idx(x,y)] - HPast[idx(x-1,y)]);
				float dH_y = alpha_H[idx(x,y)] * (HPast[idx(x,y+1)] - HPast[idx(x,y)]) - alpha_H[idx(x,y-1)] * (HPast[idx(x,y)] - HPast[idx(x,y-1)]);
				float dHdT = (dH_x + dH_y) / (CELLSIZE*CELLSIZE);
				H[idx(x,y)] = HPast[idx(x,y)] + DELTA_T * dHdT;
				H[idx(x,y)] = max(terrain[idx(x,y)], H[idx(x,y)]); // ensure water surface is above terrain

				// Diffusion step for Q: dQ/dt = Del * ( alpha_Q * Del Q )
				// Q has two components, so we do them separately
				// This could be redundant, maybe only need x component for Q_x and y component for Q_y? Leave for now unless looks weird. 
				float dQ_x_x = alpha_Q_x[idx(x,y)] * (QPast_x[idx(x+1,y)] - QPast_x[idx(x,y)]) - alpha_Q_x[idx(x-1,y)] * (QPast_x[idx(x,y)] - QPast_x[idx(x-1,y)]);
				float dQ_x_y = alpha_Q_y[idx(x,y)] * (QPast_x[idx(x,y+1)] - QPast_x[idx(x,y)]) - alpha_Q_y[idx(x,y-1)] * (QPast_x[idx(x,y)] - QPast_x[idx(x,y-1)]); // = 0.f ? 
				float dQdT_x = (dQ_x_x + dQ_x_y) / (CELLSIZE*CELLSIZE);
				Q_x[idx(x,y)] = QPast_x[idx(x,y)] + DELTA_T * dQdT_x;
				float dQ_y_x = alpha_Q_y[idx(x,y)] * (QPast_y[idx(x,y+1)] - QPast_y[idx(x,y)]) - alpha_Q_y[idx(x,y-1)] * (QPast_y[idx(x,y)] - QPast_y[idx(x,y-1)]); // = 0.f ?
				float dQ_y_y = alpha_Q_x[idx(x,y)] * (QPast_y[idx(x+1,y)] - QPast_y[idx(x,y)]) - alpha_Q_x[idx(x-1,y)] * (QPast_y[idx(x,y)] - QPast_y[idx(x-1,y)]);
				float dQdT_y = (dQ_y_x + dQ_y_y) / (CELLSIZE*CELLSIZE);
				Q_y[idx(x,y)] = QPast_y[idx(x,y)] + DELTA_T * dQdT_y;
			}
		}
		// Handle boundaries for H and Q after each diffusion iteration
		ApplyBoundaries(H, BOUNDARY_TYPE);
		ApplyBoundaries(Q_x, BOUNDARY_TYPE);
		ApplyBoundaries(Q_y, BOUNDARY_TYPE);
	}

	// final conversion to individual solver quantities
	for (int i = 0; i < GRIDSIZE*GRIDSIZE; i++)
	{
		hbar[i] = max(0.f, H[i] - terrain[i]);
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
	// Apply to boundaries
	for (int i = 0; i < GRIDSIZE; i++)
	{
		// Right Edge
		if (qbar_x[idx(GRIDSIZE-1,i)] == 0.f)  // if flow is going out of domain
		{
			qbar_x[idx(GRIDSIZE-1,i)] = 0.f;
			qtilde_x[idx(GRIDSIZE-1,i)] = 0.f;
		}
		// Top Edge
		if (qbar_y[idx(i,GRIDSIZE-1)] == 0.f)  // if flow is going out of domain
		{
			qbar_y[idx(i,GRIDSIZE-1)] = 0.f;
			qtilde_y[idx(i,GRIDSIZE-1)] = 0.f;
		}
	}
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
// 		double kNonZero = max(0.01, k);
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
// 			double k2 = max(0.0001, 2. * kx / GRIDSIZE);  //k2 = 0..1
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
// 		float waterDepth = max(hbar[x], hbar[x_plus]);
// 		int depth1 = 0;
// 		for (int depth = 0; depth < DEPTH_NUM; depth++)
// 			if (waterDepth >= Depth[depth])
// 				depth1 = depth;
// 		int depth2 = min(DEPTH_NUM - 1, depth1 + 1);
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
	for (int y = 0; y < GRIDSIZE; y++)
	{
		for (int x = 0; x < GRIDSIZE; x++)
		{
			ubar_x[idx(x,y)] = qbar_x[idx(x,y)];
			ubar_y[idx(x,y)] = qbar_y[idx(x,y)];

			// First-Order Up-Winding
			// SOMEDAY: Try interpolating h or using higher-order upwinding for better accuracy?
			// Technically u = q / H, not h?? Different derivations differ here
			if (ubar_x[idx(x,y)] >= 0.f || x == GRIDSIZE-1)
				ubar_x[idx(x,y)] /= max(MIN_WATER_HEIGHT, hbarOld[idx(x,y)]);
			else
				ubar_x[idx(x,y)] /= max(MIN_WATER_HEIGHT, hbarOld[idx(x+1,y)]);
			if (ubar_y[idx(x,y)] >= 0.f || y == GRIDSIZE-1)
				ubar_y[idx(x,y)] /= max(MIN_WATER_HEIGHT, hbarOld[idx(x,y)]);
			else
				ubar_y[idx(x,y)] /= max(MIN_WATER_HEIGHT, hbarOld[idx(x,y+1)]);

			// Enforcing CFL condition for later surface waves advection
			ubar_x[idx(x,y)] = LimitVelocity(ubar_x[idx(x,y)]);  
			ubar_y[idx(x,y)] = LimitVelocity(ubar_y[idx(x,y)]);  
		}
	}
	hbarOld = hbar;   // store current hbar for next timestep

	// Compute time derivative of u_bar and integrate to get new u_bar
	for (int y = 1; y < GRIDSIZE-1; y++)
	{
		for (int x = 1; x < GRIDSIZE-1; x++)
		{
			// Compute intermediate values needed for du/dt calculations
			// Need: q_x/y_ij, q_x_i-0.5_j, q_y_i_j-0.5, 
			// 	     h_i+0.5_j, h_i_j+0.5, h_ij=hbar[idx(x,y)], h_i+1_j=hbar[idx(x+1,y)], h_i_j+1=hbar[idx(x,y+1)],
			//       u_x_i+0.5_j=ubar_x[idx(x,y)], u_x_i-0.5_j=ubar_x[idx(x-1,y)], 
			//       u_y_i_j+0.5=ubar_y[idx(x,y)], u_y_i_j-0.5=ubar_y[idx(x,y-1)], 
			//       u_x_i+0.5_j-1=ubar_x[idx(x,y-1)], u_y_i-1_j+0.5=ubar_y[idx(x-1,y)]
			// Note: The 1D implementation uses far more intermeidate values that the paper, and I don't know why. 
			// The commented sections below are the conversion of their 1D code with all its complexity, and the 
			// uncommented sections are my simplified version that tries to stay more true to the paper. We'll see 
			// if it causes stability/accuracy issues.

			////// X DIRECTION //////
			// Use upwinding to evaluate q_(i-0.5,j), q_(i+0.5,j), q_(i+1.5,j) 
			// for x direction to get q_x_(i,j), q_x_(i+1,j)
			float q_x_m05 = ubar_x[idx(x-1,y)];
			if (q_x_m05 >= 0.f)
				q_x_m05 *= hbar[idx(x-1,y)];
			else
				q_x_m05 *= hbar[idx(x,y)];
			float q_x_p05 = ubar_x[idx(x,y)];  
			if (q_x_p05 >= 0.f)
				q_x_p05 *= hbar[idx(x,y)];
			else
				q_x_p05 *= hbar[idx(x+1,y)];
			float q_x_0 = 0.5f * (q_x_m05 + q_x_p05);
			// float q_x_p15 = ubar_x[idx(x+1,y)];  //q_(i+1.5,j) = hfr at position x
			// if (q_x_p15 >= 0.f)
			// 	q_x_p15 *= hbar[idx(x+1,y)];
			// else
			// 	q_x_p15 *= hbar[min(idx(x+1,y) + 1, GRIDSIZE*GRIDSIZE-1)];
			// float q_x_p1 = 0.5f * (q_x_p05 + q_x_p15);
			// Calculate corresponding vaules for u_x_(i,j) using upwinding
			// (why do we use upwinding here instead of averaging like q?)
			// float u_star_x_0 = 0.f;
			// if (q_x_0 >= 0.f)
			// 	u_star_x_0 = ubar_x[idx(x-1,y)];
			// else
			// 	u_star_x_0 = ubar_x[idx(x,y)];
			// float u_star_x_p1 = 0.f;
			// if (q_x_p1 > 0.f)
			// 	u_star_x_p1 = ubar_x[idx(x,y)];
			// else
			// 	u_star_x_p1 = ubar_x[idx(x+1,y)];

			// Calculate h_(i+0.5,j) and h_(i-0.5,j) using upwinding
			// float h_avg_x_p05 = (hbar[idx(x,y)] + hbar[idx(x+1,y)]) / 2.f; // averaging, for some reason the old paper used this
			float h_x_p05 = 0.f;
			if (ubar_x[idx(x,y)] >= 0.f)
				h_x_p05 = hbar[idx(x,y)];
			else
				h_x_p05 = hbar[idx(x+1,y)];


			/////// Y DIRECTION //////
			// Use upwinding to evaluate q_(i,j-0.5), q_(i,j+0.5), q_(i,j+1.5) 
			// for y direction to get q_y_(i,j), q_y_(i,j+1)
			float q_y_m05 = ubar_y[idx(x,y-1)];
			if (q_y_m05 >= 0.f)
				q_y_m05 *= hbar[idx(x,y-1)];
			else
				q_y_m05 *= hbar[idx(x,y)];
			float q_y_p05 = ubar_y[idx(x,y)];  
			if (q_y_p05 >= 0.f)
				q_y_p05 *= hbar[idx(x,y)];
			else
				q_y_p05 *= hbar[idx(x,y+1)];
			float q_y_0 = 0.5f * (q_y_m05 + q_y_p05);
			// float q_y_p15 = ubar_y[idx(x,y+1)];  //q_(i,j+1.5) = hfr at position y
			// if (q_y_p15 >= 0.f)
			// 	q_y_p15 *= hbar[idx(x,y+1)];
			// else
			// 	q_y_p15 *= hbar[min(idx(x,y+2), GRIDSIZE*GRIDSIZE-1)];
			// float q_y_p1 = 0.5f * (q_y_p05 + q_y_p15);
			// Calculate corresponding vaules for u_y_(i,j) using upwinding
			// float u_star_y_0 = 0.f;
			// if (q_y_0 >= 0.f)
			// 	u_star_y_0 = ubar_y[idx(x,y-1)];
			// else
			// 	u_star_y_0 = ubar_y[idx(x,y+1)];
			// float u_star_y_p1 = 0.f;
			// if (q_y_p1 > 0.f)
			// 	u_star_y_p1 = ubar_y[idx(x,y)];
			// else
			// 	u_star_y_p1 = ubar_y[idx(x,y+1)];

			// Calculate h_(i,j+0.5) and h_(i,j-0.5) using upwinding
			float h_y_p05 = 0.f;
			if (ubar_y[idx(x,y)] >= 0.f)
				h_y_p05 = hbar[idx(x,y)];
			else
				h_y_p05 = hbar[idx(x,y+1)];


			// Compute dux_dt and duy_dt
			// X DIRECTION
			float dux_dt = - (1/CELLSIZE) * ((q_x_0/h_x_p05) * (ubar_x[idx(x,y)] - ubar_x[idx(x-1,y)]) + (q_y_m05/h_x_p05) * (ubar_x[idx(x,y)] - ubar_x[idx(x,y-1)]) + GRAVITY * (H[idx(x+1,y)] - H[idx(x,y)]));
			ubarNew_x[idx(x,y)] = LimitVelocity(ubar_x[idx(x,y)] + TIMESTEP * dux_dt);  // Enforcing CFL condition
			// Y DIRECTION
			float duy_dt = - (1/CELLSIZE) * ((q_y_0/h_y_p05) * (ubar_y[idx(x,y)] - ubar_y[idx(x,y-1)]) + (q_x_m05/h_y_p05) * (ubar_y[idx(x,y)] - ubar_y[idx(x-1,y)]) + GRAVITY * (H[idx(x,y+1)] - H[idx(x,y)]));
			ubarNew_y[idx(x,y)] = LimitVelocity(ubar_y[idx(x,y)] + TIMESTEP * duy_dt);  // Enforcing CFL condition
		}
	}
	ApplyBoundaries(ubarNew_x, BOUNDARY_TYPE);
	ApplyBoundaries(ubarNew_y, BOUNDARY_TYPE);

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
}

void Sim::TransportStep()
{
	// Advect high-frequency wave height and flow rate through bulk velocity

	// Adjust qtilde to account for advection by ubar, using cubic sampling to get better accuracy.
	qtildePast_x = qtilde_x;  // store current qtilde for sampling, we will update qtilde in place
	qtildePast_y = qtilde_y; 
	for (int y = 0; y < GRIDSIZE-1; y++)
	{
		for (int x = 0; x < GRIDSIZE-1; x++)
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
		}
	}
	// handle boundaries
	for (int i = 0; i < GRIDSIZE; i++)
	{
		std::vector<float> local_x(GRIDSIZE);
		std::vector<float> local_y(GRIDSIZE);

		// Right Edge
		int idx_r = idx(GRIDSIZE-1,i);
		int idx_rplus = idx(GRIDSIZE-1, min(i+1, GRIDSIZE-1));
		ExtractLocal(local_x, qtildePast_x, idx_r, false);
		ExtractLocal(local_y, qtildePast_y, idx_r, true);
		float bulkVelocity_x = 0.5f * (ubarNew_x[idx_r] + ubar_x[idx_r]);
		float bulkVelocity_y = 0.5f * (ubarNew_y[idx_r] + ubar_y[idx_r]);
		float step_x = - bulkVelocity_x * TIMESTEP / CELLSIZE; // unitless (cells)
		float step_y = - bulkVelocity_y * TIMESTEP / CELLSIZE; // unitless (cells)
		qtilde_x[idx_r] = SampleCubicClamped(GRIDSIZE-1 + step_x, local_x);
		qtilde_y[idx_r] = SampleCubicClamped(i + step_y, local_y);
		if (((bulkVelocity_x >= 0.f) && (h[idx_r] < MIN_WATER_HEIGHT)) ||
			((bulkVelocity_x < 0.f) && (h[idx_r] < MIN_WATER_HEIGHT)))
			qtilde_x[idx_r] = 0.f;
		if (((bulkVelocity_y >= 0.f) && (h[idx_r] < MIN_WATER_HEIGHT)) ||
			((bulkVelocity_y < 0.f) && (h[idx_rplus] < MIN_WATER_HEIGHT)))
			qtilde_y[idx_r] = 0.f;

		// Top Edge
		int idx_t = idx(i,GRIDSIZE-1);
		int idx_tplus = idx(min(i+1, GRIDSIZE-1), GRIDSIZE-1);
		ExtractLocal(local_x, qtildePast_x, idx_t, false);
		ExtractLocal(local_y, qtildePast_y, idx_t, true);
		float bulkVelocity_x = 0.5f * (ubarNew_x[idx_t] + ubar_x[idx_t]);
		float bulkVelocity_y = 0.5f * (ubarNew_y[idx_t] + ubar_y[idx_t]);
		float step_x = - bulkVelocity_x * TIMESTEP / CELLSIZE; // unitless (cells)
		float step_y = - bulkVelocity_y * TIMESTEP / CELLSIZE; // unitless (cells)
		qtilde_x[idx_t] = SampleCubicClamped(i + step_x, local_x);
		qtilde_y[idx_t] = SampleCubicClamped(GRIDSIZE-1 + step_y, local_y);
		if (((bulkVelocity_x >= 0.f) && (h[idx_t] < MIN_WATER_HEIGHT)) ||
			((bulkVelocity_x < 0.f) && (h[idx_t] < MIN_WATER_HEIGHT)))
			qtilde_x[idx_t] = 0.f;
		if (((bulkVelocity_y >= 0.f) && (h[idx_t] < MIN_WATER_HEIGHT)) ||
			((bulkVelocity_y < 0.f) && (h[idx_tplus] < MIN_WATER_HEIGHT)))
			qtilde_y[idx_t] = 0.f;
	}


	// Update qtilde from ubar divergence: dq/dt = -q * div(ubar)
	for (int y = 1; y < GRIDSIZE-1; y++)
	{
		for (int x = 1; x < GRIDSIZE-1; x++)
		{
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
				div_ubar *= GAMMA;	
					
			qtilde_x[idx(x,y)] *= exp(-div_ubar * TIMESTEP);
			qtilde_y[idx(x,y)] *= exp(-div_ubar * TIMESTEP);
		}
	}
	// handle boundaries
	for (int i = 0; i < GRIDSIZE; i++)
	{
		// Left Edge
		int idx_lplusy = idx(0, min(i+1, GRIDSIZE-1));
		int idx_lminy = idx(0, max(i-1, 0));
		float ubar_x_m1 = 0.5f * (ubarNew_x[idx(0,i)] + ubar_x[idx(0,i)]);
		float ubar_x_p1 = 0.5f * (ubarNew_x[idx(1,i)] + ubar_x[idx(1,i)]);
		float ubar_y_m1 = 0.5f * (ubarNew_y[idx_lminy] + ubar_y[idx_lminy]);
		float ubar_y_p1 = 0.5f * (ubarNew_y[idx_lplusy] + ubar_y[idx_lplusy]);
		float div_ubar_x = (ubar_x_p1 - ubar_x_m1) / (2.f * CELLSIZE);
		float div_ubar_y = (ubar_y_p1 - ubar_y_m1) / (2.f * CELLSIZE);
		float div_ubar = div_ubar_x + div_ubar_y;
		if (div_ubar < 0.f)
			div_ubar *= GAMMA;
		qtilde_x[idx(0,i)] *= exp(-div_ubar * TIMESTEP);
		qtilde_y[idx(0,i)] *= exp(-div_ubar * TIMESTEP);

		// Right Edge
		int idx_rplusy = idx(GRIDSIZE-1, min(i+1, GRIDSIZE-1));
		int idx_rminy = idx(GRIDSIZE-1, max(i-1, 0));
		ubar_x_m1 = 0.5f * (ubarNew_x[idx(GRIDSIZE-2,i)] + ubar_x[idx(GRIDSIZE-2,i)]);
		ubar_x_p1 = 0.5f * (ubarNew_x[idx(GRIDSIZE-1,i)] + ubar_x[idx(GRIDSIZE-1,i)]);
		ubar_y_m1 = 0.5f * (ubarNew_y[idx_rminy] + ubar_y[idx_rminy]);
		ubar_y_p1 = 0.5f * (ubarNew_y[idx_rplusy] + ubar_y[idx_rplusy]);
		div_ubar_x = (ubar_x_p1 - ubar_x_m1) / (2.f * CELLSIZE);
		div_ubar_y = (ubar_y_p1 - ubar_y_m1) / (2.f * CELLSIZE);
		div_ubar = div_ubar_x + div_ubar_y;
		if (div_ubar < 0.f)
			div_ubar *= GAMMA;
		qtilde_x[idx(GRIDSIZE-1,i)] *= exp(-div_ubar * TIMESTEP);
		qtilde_y[idx(GRIDSIZE-1,i)] *= exp(-div_ubar * TIMESTEP);		

		// Top Edge
		int idx_tplusx = idx(min(i+1,GRIDSIZE-1),GRIDSIZE-1);
		int idx_tminx = idx(max(i-1,0),GRIDSIZE-1);
		ubar_x_m1 = 0.5f * (ubarNew_x[idx_tminx] + ubar_x[idx_tminx]);
		ubar_x_p1 = 0.5f * (ubarNew_x[idx_tplusx] + ubar_x[idx_tplusx]);
		ubar_y_m1 = 0.5f * (ubarNew_y[idx(i,GRIDSIZE-2)] + ubar_y[idx(i,GRIDSIZE-2)]);
		ubar_y_p1 = 0.5f * (ubarNew_y[idx(i,GRIDSIZE-1)] + ubar_y[idx(i,GRIDSIZE-1)]);
		div_ubar_x = (ubar_x_p1 - ubar_x_m1) / (2.f * CELLSIZE);
		div_ubar_y = (ubar_y_p1 - ubar_y_m1) / (2.f * CELLSIZE);
		div_ubar = div_ubar_x + div_ubar_y;
		if (div_ubar < 0.f)
			div_ubar *= GAMMA;
		qtilde_x[idx(i,GRIDSIZE-1)] *= exp(-div_ubar * TIMESTEP);
		qtilde_y[idx(i,GRIDSIZE-1)] *= exp(-div_ubar * TIMESTEP);

		// Bottom Edge
		int idx_bplusx = idx(min(i+1,GRIDSIZE-1),0);
		int idx_bminx = idx(max(i-1,0),0);
		ubar_x_m1 = 0.5f * (ubarNew_x[idx_bminx] + ubar_x[idx_bminx]);
		ubar_x_p1 = 0.5f * (ubarNew_x[idx_bplusx] + ubar_x[idx_bplusx]);
		ubar_y_m1 = 0.5f * (ubarNew_y[idx(i,0)] + ubar_y[idx(i,0)]);
		ubar_y_p1 = 0.5f * (ubarNew_y[idx(i,1)] + ubar_y[idx(i,1)]);
		div_ubar_x = (ubar_x_p1 - ubar_x_m1) / (2.f * CELLSIZE);
		div_ubar_y = (ubar_y_p1 - ubar_y_m1) / (2.f * CELLSIZE);
		div_ubar = div_ubar_x + div_ubar_y;
		if (div_ubar < 0.f)
			div_ubar *= GAMMA;
		qtilde_x[idx(i,0)] *= exp(-div_ubar * TIMESTEP);
		qtilde_y[idx(i,0)] *= exp(-div_ubar * TIMESTEP);
	}
	

	// Update htilde from ubar divergence: dh/dt = -h * div(ubar)
	for (int y = 1; y < GRIDSIZE; y++)
	{
		for (int x = 1; x < GRIDSIZE; x++)
		{
			// backward difference to get divergence of ubar (using backward because u is on staggered grid from h)
			float div_ubar_x = (ubarNew_x[idx(x,y)] - ubarNew_x[idx(x-1,y)]) / CELLSIZE;
			float div_ubar_y = (ubarNew_y[idx(x,y)] - ubarNew_y[idx(x,y-1)]) / CELLSIZE;
			float div_ubar = div_ubar_x + div_ubar_y;
			// dampen if converging to avoid breaking waves
			if (div_ubar < 0.f)
				div_ubar *= GAMMA;
			htilde[idx(x,y)] *= exp(-div_ubar * TIMESTEP);
		}
	}
	// handle boundaries: left and bottom
	for (int i = 0; i < GRIDSIZE; i++)
	{
		// Assume div_u_bar is same as adjacent cell, but h_tilde is different, so need to recalculate fully
		// Left Edge
		int idx_lminy = idx(0, max(i-1,0));
		float div_ubar_x = (ubarNew_x[idx(1,i)] - ubarNew_x[idx(0,i)]) / CELLSIZE; 
		float div_ubar_y = (ubarNew_y[idx_lminy + GRIDSIZE] - ubarNew_y[idx_lminy]) / CELLSIZE;
		float div_ubar = div_ubar_x + div_ubar_y;
		if (div_ubar < 0.f)
			div_ubar *= GAMMA;
		htilde[idx(0,i)] *= exp(-div_ubar * TIMESTEP);

		// Bottom Edge
		int idx_bminx = idx(max(i-1,0),0);
		float div_ubar_x = (ubarNew_x[idx_bminx + 1] - ubarNew_x[idx_bminx]) / CELLSIZE;
		float div_ubar_y = (ubarNew_y[idx(i,1)] - ubarNew_y[idx(i,0)]) / CELLSIZE;
		float div_ubar = div_ubar_x + div_ubar_y;
		if (div_ubar < 0.f)
			div_ubar *= GAMMA;
		htilde[idx(i,0)] *= exp(-div_ubar * TIMESTEP);
	}


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
			h[idx(x,y)] = max(0.f, hPast[idx(x,y)] + TIMESTEP * div_q);
		}
	}
	// handle boundaries
	for (int i = 0; i < GRIDSIZE; i++)
	{
		// Left Edge
		int idx_lplusy = idx(0, min(i+1, GRIDSIZE-1));
		int idx_lminy = idx(0, max(i-1, 0));
		float q_x_l = LimitFlowRate(qAdvect_x[idx(0,i)], h[idx(0,i)], h[idx(0,i)]);
		float q_x_r = LimitFlowRate(qAdvect_x[idx(0,i)], h[idx(0,i)], h[idx(1,i)]);
		float q_y_l = LimitFlowRate(qAdvect_y[idx_lminy], h[idx_lminy], h[idx_lminy + GRIDSIZE]);
		float q_y_r = LimitFlowRate(qAdvect_y[idx_lplusy - GRIDSIZE], h[idx_lplusy - GRIDSIZE], h[idx_lplusy]);
		if ( ((h[idx(0,i)] == 0.f)) || (StopFlowOnTerrainBoundary(0, i, h, terrain, false)) )
			q_x_l = 0.f;
		if ( ((h[idx(0,i)] == 0.f) && (h[idx(1,i)] == 0.f)) || (StopFlowOnTerrainBoundary(0, i, h, terrain, false)) )
			q_x_r = 0.f;
		if ( ((h[idx_lminy] == 0.f) && (h[idx_lminy+GRIDSIZE] == 0.f)) || (StopFlowOnTerrainBoundary(0, max(i-1,0), h, terrain, true)) )
			q_y_l = 0.f;
		if ( ((h[idx_lplusy] == 0.f) && (h[idx(0,i)] == 0.f)) || (StopFlowOnTerrainBoundary(0, min(i+1,GRIDSIZE-1), h, terrain, true)) )
			q_y_r = 0.f;
		float div_q = (q_x_l - q_x_r + q_y_l - q_y_r) / CELLSIZE;
		h[idx(0,i)] = max(0.f, hPast[idx(0,i)] + TIMESTEP * div_q);

		// Right Edge
		int idx_rplusy = idx(GRIDSIZE-1, min(i+1,GRIDSIZE-1));
		int idx_rminy = idx(GRIDSIZE-1, max(i-1,0));
		float q_x_l = LimitFlowRate(qAdvect_x[idx(GRIDSIZE-2,i)], h[idx(GRIDSIZE-2,i)], h[idx(GRIDSIZE-1,i)]);
		float q_x_r = LimitFlowRate(qAdvect_x[idx(GRIDSIZE-1,i)], h[idx(GRIDSIZE-1,i)], h[idx(GRIDSIZE-1,i)]);
		float q_y_l = LimitFlowRate(qAdvect_y[idx_rminy], h[idx_rminy], h[idx_rminy + GRIDSIZE]);
		float q_y_r = LimitFlowRate(qAdvect_y[idx_rplusy - GRIDSIZE], h[idx_rplusy - GRIDSIZE], h[idx_rplusy]);
		if ( ((h[idx(GRIDSIZE-2,i)] == 0.f) && (h[idx(GRIDSIZE-1,i)] == 0.f)) || (StopFlowOnTerrainBoundary(GRIDSIZE-2, i, h, terrain, false)) )
			q_x_l = 0.f;
		if ( ((h[idx(GRIDSIZE-1,i)] == 0.f)) || (StopFlowOnTerrainBoundary(GRIDSIZE-1, i, h, terrain, false)) )
			q_x_r = 0.f;
		if ( ((h[idx_rminy] == 0.f) && (h[idx_rminy+GRIDSIZE] == 0.f)) || (StopFlowOnTerrainBoundary(GRIDSIZE-1, max(i-1,0), h, terrain, true)) )
			q_y_l = 0.f;
		if ( ((h[idx_rplusy] == 0.f) && (h[idx(GRIDSIZE-1,i)] == 0.f)) || (StopFlowOnTerrainBoundary(GRIDSIZE-1, min(i+1,GRIDSIZE-1), h, terrain, true)) )
			q_y_r = 0.f;
		float div_q = (q_x_l - q_x_r + q_y_l - q_y_r) / CELLSIZE;
		h[idx(GRIDSIZE-1,i)] = max(0.f, hPast[idx(GRIDSIZE-1,i)] + TIMESTEP * div_q);

		// Top Edge
		int idx_tplusx = idx(min(i+1,GRIDSIZE-1),GRIDSIZE-1);
		int idx_tminx = idx(max(i-1,0),GRIDSIZE-1);
		float q_x_l = LimitFlowRate(qAdvect_x[idx_tminx], h[idx_tminx], h[idx_tminx + 1]);
		float q_x_r = LimitFlowRate(qAdvect_x[idx_tplusx - 1], h[idx_tplusx - 1], h[idx_tplusx]);
		float q_y_l = LimitFlowRate(qAdvect_y[idx(i,GRIDSIZE-2)], h[idx(i,GRIDSIZE-2)], h[idx(i,GRIDSIZE-1)]);
		float q_y_r = LimitFlowRate(qAdvect_y[idx(i,GRIDSIZE-1)], h[idx(i,GRIDSIZE-1)], h[idx(i,GRIDSIZE-1)]);	
		if ( ((h[idx_tminx] == 0.f) && (h[idx_tminx+1] == 0.f)) || (StopFlowOnTerrainBoundary(max(i-1,0), GRIDSIZE-1, h, terrain, false)) )
			q_x_l = 0.f;
		if ( ((h[idx_tplusx-1] == 0.f) && (h[idx_tplusx] == 0.f)) || (StopFlowOnTerrainBoundary(min(i,GRIDSIZE-2), GRIDSIZE-1, h, terrain, false)) )
			q_x_r = 0.f;
		if ( ((h[idx(i,GRIDSIZE-2)] == 0.f) && (h[idx(i,GRIDSIZE-1)] == 0.f)) || (StopFlowOnTerrainBoundary(i, GRIDSIZE-2, h, terrain, true)) )
			q_y_l = 0.f;
		if ( ((h[idx(i,GRIDSIZE-1)] == 0.f)) || (StopFlowOnTerrainBoundary(i, GRIDSIZE-1, h, terrain, true)) )
			q_y_r = 0.f;
		float div_q = (q_x_l - q_x_r + q_y_l - q_y_r) / CELLSIZE;
		h[idx(i,GRIDSIZE-1)] = max(0.f, hPast[idx(i,GRIDSIZE-1)] + TIMESTEP * div_q);

		// Bottom Edge
		int idx_bplusx = idx(min(i+1,GRIDSIZE-1),0);
		int idx_bminx = idx(max(i-1,0),0);
		float q_x_l = LimitFlowRate(qAdvect_x[idx_bminx], h[idx_bminx], h[idx_bminx + 1]);
		float q_x_r = LimitFlowRate(qAdvect_x[idx_bplusx - 1], h[idx_bplusx - 1], h[idx_bplusx]);
		float q_y_l = LimitFlowRate(qAdvect_y[idx(i,0)], h[idx(i,0)], h[idx(i,0)]);
		float q_y_r = LimitFlowRate(qAdvect_y[idx(i,0)], h[idx(i,0)], h[idx(i,1)]);
		if ( ((h[idx_bminx] == 0.f) && (h[idx_bminx+1] == 0.f)) || (StopFlowOnTerrainBoundary(max(i-1,0), 0, h, terrain, false)) )
			q_x_l = 0.f;
		if ( ((h[idx_bplusx-1] == 0.f) && (h[idx_bplusx] == 0.f)) || (StopFlowOnTerrainBoundary(min(i,GRIDSIZE-2), 0, h, terrain, false)) )
			q_x_r = 0.f;
		if ( ((h[idx(i,0)] == 0.f)) || (StopFlowOnTerrainBoundary(i, 0, h, terrain, true)) )
			q_y_l = 0.f;
		if ( ((h[idx(i,0)] == 0.f) && (h[idx(i,1)] == 0.f)) || (StopFlowOnTerrainBoundary(i, 0, h, terrain, true)) )
			q_y_r = 0.f;
		float div_q = (q_x_l - q_x_r + q_y_l - q_y_r) / CELLSIZE;
		h[idx(i,0)] = max(0.f, hPast[idx(i,0)] + TIMESTEP * div_q);
	}
}
	
void Sim::ComputeValues()
{
	// Recombine bulk and surface flow
	for (int i = 0; i < GRIDSIZE*GRIDSIZE; i++)
	{
		q_x[i] = qbar_x[i] + qtilde_x[i];
		q_y[i] = qbar_y[i] + qtilde_y[i];
	}
	for (int y = 0; y < GRIDSIZE; y++)
	{
		for (int x = 0; x < GRIDSIZE; x++)
		{
			q_x[idx(x,y)] = LimitFlowRate(q_x[idx(x,y)], h[idx(x,y)], h[idx(min(x+1,GRIDSIZE-1),y)]);
			q_y[idx(x,y)] = LimitFlowRate(q_y[idx(x,y)], h[idx(x,y)], h[idx(x,min(y+1,GRIDSIZE-1))]);
			if ( (StopFlowOnTerrainBoundary(x, y, h, terrain, false)) || (x == 0) || (x >= GRIDSIZE - 2) )
				q_x[idx(x,y)] = 0.f;
			if ( (StopFlowOnTerrainBoundary(x, y, h, terrain, true)) || (y == 0) || (y >= GRIDSIZE - 2) )
				q_y[idx(x,y)] = 0.f;
			// NOTE: making edges a reflective boundary, revisit later
		}
	}
	
	// Height integration 
	for (int y = 1; y < GRIDSIZE; y++)
	{
		for (int x = 1; x < GRIDSIZE; x++)
		{
			h[idx(x,y)] = max(0.f, h[idx(x,y)] + TIMESTEP * -(q_x[idx(x,y)] - q_x[idx(x-1,y)] + q_y[idx(x,y)] - q_y[idx(x,y-1)]) / CELLSIZE);
		}
	}
	// handle boundaries (left and bottom only)
	for (int i = 0; i < GRIDSIZE; i++)
	{
		h[idx(0,i)] = max(0.f, h[idx(0,i)] + TIMESTEP * -(q_x[idx(0,i)] - 0.f + q_y[idx(0,i)] - q_y[idx(0,min(i-1,0))]) / CELLSIZE); // Left
		h[idx(i,0)] = max(0.f, h[idx(i,0)] + TIMESTEP * -(q_x[idx(i,0)] - q_x[idx(min(i-1,0),0)] + q_y[idx(i,0)] - 0.f) / CELLSIZE); // Bottom
	}
	
	// stability measure to not drag too much water from a cell in a single timestep (important for extreme initial conditions)
	for (int y = 0; y < GRIDSIZE - 1; y++)
	{
		for (int x = 0; x < GRIDSIZE - 1; x++)
		{
			q_x[idx(x,y)] = LimitFlowRate(q_x[idx(x,y)], h[idx(x,y)], h[idx(x+1,y)]);
			q_y[idx(x,y)] = LimitFlowRate(q_y[idx(x,y)], h[idx(x,y)], h[idx(x,y+1)]);
		}
	}
	// handle boundaries (top and right only)
	for (int i = 0; i < GRIDSIZE; i++)
	{
		q_x[idx(GRIDSIZE-1,i)] = LimitFlowRate(q_x[idx(GRIDSIZE-1,i)], h[idx(GRIDSIZE-1,i)], h[idx(GRIDSIZE-1,i)]); // Right
		q_y[idx(i,GRIDSIZE-1)] = LimitFlowRate(q_y[idx(i,GRIDSIZE-1)], h[idx(i,GRIDSIZE-1)], h[idx(i,GRIDSIZE-1)]); // Top
	}
}