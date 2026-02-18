#include <windows.h>
#include "Sim2D.h"

 

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

// test if the terrain boundary stops any flow across x+0.5
int StopFlowOnTerrainBoundary(int x, int y, float* h, float* terrain)
{
	// Key: 1 = stop in x, 2 = stop in y, 3 = stop in both, 0 = no stop
	bool result_x = 0;
	bool result_y = 0;
	int result = 0;

	// Test x boundary
	if ((h[idx(x,y)] <= MIN_WATER_HEIGHT) && (terrain[idx(x,y)] >= terrain[idx(x+1,y)] + h[idx(x+1,y)])) // positive q_x
		result_x = 1;
	if ((h[idx(x+1,y)] <= MIN_WATER_HEIGHT) && (terrain[idx(x+1,y)] > terrain[idx(x,y)] + h[idx(x,y)])) // negative q_x
		result_x = 1;

	// Test y boundary
	if ((h[idx(x,y)] <= MIN_WATER_HEIGHT) && (terrain[idx(x,y)] >= terrain[idx(x,y+1)] + h[idx(x,y+1)])) // positive q_y
		result_y = 1;
	if ((h[idx(x,y+1)] <= MIN_WATER_HEIGHT) && (terrain[idx(x,y+1)] > terrain[idx(x,y)] + h[idx(x,y)])) // negative q_y
		result_y = 1;

	// Combine results
	if (result_x && result_y)
		result = 3;
	else if (result_x)
		result = 1;
	else if (result_y)
		result = 2;
	return result;
}


// // xCoord and size in (0..1), factor determines how much to add/subtract
// void Sim::EditWaterLocal(float xCoord, float yCoord, float size, float factor)
// {
// 	for (int y = 0; y < GRIDSIZE; y++)
// 	{
// 		for (int x = 0; x < GRIDSIZE; x++)
// 			if (fabs((float)(x) / GRIDSIZE - xCoord) < size)
// 					h[idx(x,y)] = max(0.f, h[idx(x,y)] + factor * 1.f);
// 		h[idx(0,y)] = 0.0f;
// 		h[idx(GRIDSIZE - 1,y)] = 0.0f;
// 	}
// }


float ExtractLocal(float* target, float* dataField, float index, bool isYDirection)
{
	if (isYDirection)
	{
		int x = index % GRIDSIZE;
		for (int y = 0; y < GRIDSIZE; y++)
			target[y] = dataField[y * GRIDSIZE + x];
	}
	else
	{
		int y = floor(index / GRIDSIZE);
		for (int x = 0; x < GRIDSIZE; x++)
			target[x] = dataField[y * GRIDSIZE + x];
	}
}

// cubic interpolation with Catmull-Rom Spline https://en.wikipedia.org/wiki/Cubic_Hermite_spline#Interpolating_a_data_set
float SampleCubicClamped(float samplePos, float* dataField)
{
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
	out = min(max(dataField[id1], dataField[id2]), out);  // value-limiting (same as in BFECC)
	out = max(min(dataField[id1], dataField[id2]), out);  // value-limiting (same as in BFECC)
	return out;
}


// Boundary Condition Functions
void HandleWallBoundary(float* field) 
{
	// Handle Top and Bottom Edges (Horizontal)
	for (int x = 0; x < GRIDSIZE; ++x) {
		// Top Edge
		field[idx(x, 0)] = field[idx(x, 1)];
		// Bottom Edge
		field[idx(x, GRIDSIZE - 1)] = field[idx(x, GRIDSIZE - 2)];
	}
	// Handle Left and Right Edges (Vertical)
	for (int y = 0; y < GRIDSIZE; ++y) {
		// Left Edge
		field[idx(0, y)] = field[idx(1, y)];
		// Right Edge
		field[idx(GRIDSIZE - 1, y)] = field[idx(GRIDSIZE - 2, y)];
	}
}

void HandleFreeBoundary(float* field)
{
	// Free boundary: extrapolate values from interior to boundaries using linear extrapolation

	// Handle Top and Bottom Edges (Horizontal)
	for (int x = 0; x < GRIDSIZE; ++x) {
		// Top Edge
		field[idx(x, 0)] = 2.f * field[idx(x, 1)] - field[idx(x, 2)];
		// Bottom Edge
		field[idx(x, GRIDSIZE - 1)] = 2.f * field[idx(x, GRIDSIZE - 2)] - field[idx(x, GRIDSIZE - 3)];
	}
	// Handle Left and Right Edges (Vertical)
	for (int y = 0; y < GRIDSIZE; ++y) {
		// Left Edge
		field[idx(0, y)] = 2.f * field[idx(1, y)] - field[idx(2, y)];
		// Right Edge
		field[idx(GRIDSIZE - 1, y)] = 2.f * field[idx(GRIDSIZE - 2, y)] - field[idx(GRIDSIZE - 3, y)];
	}
}

void HandleZeroBoundary(float* field)
{
	// Handle Top and Bottom Edges (Horizontal)
	for (int x = 0; x < GRIDSIZE; ++x) {
		// Top Edge
		field[idx(x, 0)] = 0.f;
		// Bottom Edge
		field[idx(x, GRIDSIZE - 1)] = 0.f;
	}
	// Handle Left and Right Edges (Vertical)
	for (int y = 0; y < GRIDSIZE; ++y) {
		// Left Edge
		field[idx(0, y)] = 0.f;
		// Right Edge
		field[idx(GRIDSIZE - 1, y)] = 0.f;
	}
}

void HandlePeriodicBoundary(float* field)
{
	// Handle Top and Bottom Edges (Horizontal)
	for (int x = 0; x < GRIDSIZE; ++x) {
		// Top Edge
		field[idx(x, 0)] = field[idx(x, GRIDSIZE - 2)];
		// Bottom Edge
		field[idx(x, GRIDSIZE - 1)] = field[idx(x, 1)];
	}
	// Handle Left and Right Edges (Vertical)
	for (int y = 0; y < GRIDSIZE; ++y) {
		// Left Edge
		field[idx(0, y)] = field[idx(GRIDSIZE - 2, y)];
		// Right Edge
		field[idx(GRIDSIZE - 1, y)] = field[idx(1, y)];
	}
}

void ApplyBoundaries(float* field, int type)
{
	if (type == 0)
		HandleWallBoundary(field);
	else if (type == 1)
		HandleFreeBoundary(field);
	else if (type == 2)
		HandleZeroBoundary(field);
	else if (type == 3)
		HandlePeriodicBoundary(field);
}

// ********************************************************************************************************************
// Init functions (TODO: TRANSITION TO 2D)
// ********************************************************************************************************************

// type: 0=flat, 1=hill
void Sim::ResetTerrain(int type)
{
	for (int x = 0; x < GRIDSIZE; x++)
		if (type == 0)  // flat
			terrain[x] = -abs(TERRAIN_HEIGHT_SHIFT_INIT);
		else if (type == 1) // hill
			terrain[x] = (-1.f + 0.1f + 0.1f * x / GRIDSIZE + 0.03f * sin(20.f * x / GRIDSIZE) + 0.9f * sin(2.5f * x / GRIDSIZE)) * abs(TERRAIN_HEIGHT_SHIFT_INIT);
			//terrain[x] = (-1.f + 0.1f + 0.5f * (0.1f * x / GRIDSIZE + 0.03f * sin(20.f * x / GRIDSIZE) + 0.9f * sin(2.5f * x / GRIDSIZE))) * abs(TERRAIN_HEIGHT_SHIFT_INIT);
	terrain[0] = 1.8f * abs(TERRAIN_HEIGHT_SHIFT_INIT);
	terrain[GRIDSIZE - 1] = 1.8f * abs(TERRAIN_HEIGHT_SHIFT_INIT);
}


// type: 0=constant level, 1=dam break, 2=sloped  level = y coordinate of water in domain
void Sim::ResetWater(int type, float level)
{
	for (int x = 0; x < GRIDSIZE; x++)
	{
		if (type == 0) //constant level
			h[x] = max(0.f, (level - terrain[x]));
		if (type == 1)  // dam break
			if (x <= GRIDSIZE / 2)
				h[x] = 0.f;
			else
				h[x] = max(0.f, (level - terrain[x]));
		if (type == 2)  // sloped water
			h[x] = max(0.f, level + (2.f * (-0.5f + (float)(x) / GRIDSIZE) * fabs(-0.5f + (float)(x) / GRIDSIZE)) * abs(0.5f * TERRAIN_HEIGHT_SHIFT_INIT) - terrain[x]);
		if (type == 3)  // flat with initial surface waves
		{
			float lambda = 10.f;
			h[x] = max(0.f, (level + 0.5f * cos(2.f * PI * (x / lambda)) - terrain[x]));
		}
		hbar[x] = h[x];
		hbar_past[x] = h[x];
		htilde[x] = 0.f;
		htilde_past[x] = 0.f;
		qbar_x[x] = 0.f;
		qbar_y[x] = 0.f;
		q_x[x] = 0.f;
		q_y[x] = 0.f;
	}
	// clear left and right boundaries
	
	h[0] = 0.f;
	h[GRIDSIZE - 1] = 0.f;
	time = 0.f;
}


Sim::Sim()
{
	htildehat.setlength(GRIDSIZE*GRIDSIZE);
	qtildehat_x.setlength(GRIDSIZE*GRIDSIZE);
	qtildehat_y.setlength(GRIDSIZE*GRIDSIZE);
	for (int i=0; i < DEPTH_NUM; i++)
		qtildehat_depth_x[i].setlength(GRIDSIZE*GRIDSIZE);
		qtildehat_depth_y[i].setlength(GRIDSIZE*GRIDSIZE);
	ResetTerrain(1);
	ResetWater(2, 0.f);
}


int Sim::Release(void)
{
	return 0;
}


// ********************************************************************************************************************
// simulation functions
// ********************************************************************************************************************


// simulate a timestep
void Sim::SimStep(bool SWEonly)
{
	DecompositionStep(SWEonly);
	if (!SWEonly)
		eWaveStep();
		// FFTStep();
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

	// Calculate diffusion coefficient (alpha) at every location
	static float alpha_H[GRIDSIZE*GRIDSIZE]; // = zeros
	static float alpha_Q[GRIDSIZE*GRIDSIZE];
	static float H[GRIDSIZE*GRIDSIZE];
	static float Q_x[GRIDSIZE*GRIDSIZE];
	static float Q_y[GRIDSIZE*GRIDSIZE];
	H = terrain + h;  // start off with the current water surface
	Q_x = q_x;
	Q_y = q_y;
	// Loop through main grid, avoid boundaries
	for (int y = 1; y < GRIDSIZE-1; y++)
	{
		for (int x = 1; x < GRIDSIZE-1; x++)
		{
			// Identify the correct height (sigma) to use for diffusivity calculation
			alpha_H[idx(x,y)] = 0.f;

			// // Their implementation: using averages (I assume to improve stability, but reduces accuracy?)
			// float maxGround = max(terrain[idx(x,y)], terrain[idx(x+1,y)], terrain[idx(x,y+1)]); // Why do this?
			// float minWaterlevel = (H[idx(x,y)] + H[idx(x+1,y)] + H[idx(x,y+1)]) / 3.f; // Why average here?
			// if ((h[idx(x,y)] > 0.f) && (h[idx(x+1,y)] > 0.f) && (h[idx(x,y+1)] > 0.f))
			// {
			// 	static const float sigma_max = 8.f;
			// // they limit diffusion coefficient to between 0 and 1, but that isn't requred by the math = maybe for stability?
			// 	float sigma = min(sigma_max, max(0.f, minWaterlevel - maxGround));
			// 	alpha_H[idx(x,y)] = sigma * sigma / (2*DELTA_T*DIFFUSION_ITERATIONS);
			// } 

			// My implementation: using local cell values only to align with eqn from paper
			if (h[idx(x,y)] > 0.f)
			{
				float denom = 2*DELTA_T*DIFFUSION_ITERATIONS;
				alpha_H[idx(x,y)] = h[idx(x,y)] * h[idx(x,y)] / denom;
				alpha_H[idx(x,y)] = min(std::sqrt(denom), alpha_H[idx(x,y)]); // clamp to improve stability; limits max depth to 
			}
			
			// Extra gradient filter
			// NOTE: they used H, I switched it to h to stay strict with the paper. We'll see if this causes bugs. 
			float gradient_x = (h[idx(x+1,y)] - h[idx(x,y)]) / CELLSIZE; // could use central difference here
			float gradient_y = (h[idx(x,y+1)] - h[idx(x,y)]) / CELLSIZE;
			alpha_H[idx(x,y)] *= exp(- DIFFUSION_PENALTY * (gradient_x * gradient_x + gradient_y * gradient_y));
			alpha_Q[idx(x,y)] = 0.5f * (alpha_H[idx(x-1,y)] + alpha_H[idx(x,y)]); // Where does this come from??
		}
	}
	ApplyBoundaries(alpha_H, BOUNDARY_TYPE);
	ApplyBoundaries(alpha_Q, BOUNDARY_TYPE);
	
	// Run diffusion to low-pass filter H and Q
	// SOMEDAY: Improve this implementation of diffusion by replacing Euler integration with FFT or something
	static float H_past[GRIDSIZE*GRIDSIZE];
	static float Q_past_x[GRIDSIZE*GRIDSIZE];
	static float Q_past_y[GRIDSIZE*GRIDSIZE];
	for (int j = 0; (j < DIFFUSION_ITERATIONS); j++)  // 64 diffusion iterations
	{
		memcpy(H_past, H, GRIDSIZE * GRIDSIZE * sizeof(float));
		memcpy(Q_past_x, Q_x, GRIDSIZE * GRIDSIZE * sizeof(float));
		memcpy(Q_past_y, Q_y, GRIDSIZE * GRIDSIZE * sizeof(float));
		for (int y = 1; y < GRIDSIZE-1; y++) // one diffusion iteration
		{
			for (int x = 1; x < GRIDSIZE - 1; x++)
			{
				// Diffusion step for H: dH/dt = Del * ( alpha_H * Del H )
				float dH_x = (alpha_H[idx(x,y)] * (H_past[idx(x+1,y)] - H_past[idx(x,y)]) - alpha_H[idx(x-1,y)] * (H_past[idx(x,y)] - H_past[idx(x-1,y)]));
				float dH_y = (alpha_H[idx(x,y)] * (H_past[idx(x,y+1)] - H_past[idx(x,y)]) - alpha_H[idx(x,y-1)] * (H_past[idx(x,y)] - H_past[idx(x,y-1)]));
				float dHdT = (dH_x + dH_y) / (CELLSIZE*CELLSIZE);
				H[idx(x,y)] = H_past[idx(x,y)] + DELTA_T * dHdT;
				H[idx(x,y)] = max(terrain[idx(x,y)], H[idx(x,y)]); // ensure water surface is above terrain

				// Diffusion step for Q: dQ/dt = Del * ( alpha_Q * Del Q )
				// Q has two components, so we do them separately
				float dQ_x_x = (alpha_Q[idx(x,y)] * (Q_past_x[idx(x+1,y)] - Q_past_x[idx(x,y)]) - alpha_Q[idx(x-1,y)] * (Q_past_x[idx(x,y)] - Q_past_x[idx(x-1,y)]));
				float dQ_x_y = (alpha_Q[idx(x,y)] * (Q_past_x[idx(x,y+1)] - Q_past_x[idx(x,y)]) - alpha_Q[idx(x,y-1)] * (Q_past_x[idx(x,y)] - Q_past_x[idx(x,y-1)]));
				float dQdT_x = (dQ_x_x + dQ_x_y) / (CELLSIZE*CELLSIZE);
				Q_x[idx(x,y)] = Q_past_x[idx(x,y)] + DELTA_T * dQdT_x;
				float dQ_y_x = (alpha_Q[idx(x,y)] * (Q_past_y[idx(x,y+1)] - Q_past_y[idx(x,y)]) - alpha_Q[idx(x,y-1)] * (Q_past_y[idx(x,y)] - Q_past_y[idx(x,y-1)]));
				float dQ_y_y = (alpha_Q[idx(x,y)] * (Q_past_y[idx(x+1,y)] - Q_past_y[idx(x,y)]) - alpha_Q[idx(x-1,y)] * (Q_past_y[idx(x,y)] - Q_past_y[idx(x-1,y)]));
				float dQdT_y = (dQ_y_x + dQ_y_y) / (CELLSIZE*CELLSIZE);
				Q_y[idx(x,y)] = Q_past_y[idx(x,y)] + DELTA_T * dQdT_y;
			}
		}
	}
	ApplyBoundaries(H, BOUNDARY_TYPE);
	ApplyBoundaries(Q_x, BOUNDARY_TYPE);
	ApplyBoundaries(Q_y, BOUNDARY_TYPE);

	// final conversion to individual solver quantities
	hbar = max(0.f, H - terrain);
	qbar_x = Q_x;
	qbar_y = Q_y;
	htilde = h - hbar;
	qtilde_x = q_x - qbar_x;
	qtilde_y = q_y - qbar_y;

	// Enforce no-flow conditions at terrain boundaries
	for (int y = 1; y < GRIDSIZE-1; y++)
	{
		for (int x = 1; x < GRIDSIZE-1; x++)
		{
			int stop_flow = StopFlowOnTerrainBoundary(x, y, h, terrain);
			if (stop_flow & 1)  // stop flow in x direction
			{
				qbar_x[idx(x,y)] = 0.f;
				qtilde_x[idx(x,y)] = 0.f;
			}
			if (stop_flow & 2)  // stop flow in y direction
			{
				qbar_y[idx(x,y)] = 0.f;
				qtilde_y[idx(x,y)] = 0.f;
			}
		}
	}
}

void Sim::eWaveStep()
{
	// surface velocity update using eWave
	for (int x = 0; x < GRIDSIZE; x++)
	{
		htildehat[x].x = 0.5f * (htilde[x] + htildeOld[x]);
		htildeOld[x] = htilde[x];
		htildehat[x].y = 0.;
		qtildehat[x].x = qtilde[x];
		qtildehat[x].y = 0.;
	}
	fftc1d(htildehat);   //https://www.alglib.net/download.php#cpp
	fftc1d(qtildehat);
	
	for (int x = 0; x < GRIDSIZE; x++)
	{
		// physical k from grid position
		double kx = GRIDSIZE / 2. - abs(GRIDSIZE / 2. - x);  // this gives [0,..,m_gridSizeX / 2.f-1, m_gridSizeX / 2.f, .. 1]
		double k = 2. * PI * fabs(kx) / GRIDSIZE / CELLSIZE;
		double kNonZero = max(0.01, k);
		double kS = k;  // signed k
		if (x > (double)(GRIDSIZE) / 2.f)
			kS = -k;
		// Fourier gradient: multiply by -i k
		double real = htildehat[x].x;
		double imag = htildehat[x].y;
		htildehat[x].x = -kS * imag;
		htildehat[x].y = kS * real;
		// phase shift to translate function to cell boundaries
		real = htildehat[x].x;
		imag = htildehat[x].y;
		double beta = 0.5 * CELLSIZE * kS;
		htildehat[x].x = cos(beta) * real - sin(beta) * imag;
		htildehat[x].y = sin(beta) * real + cos(beta) * imag;
		for (int depth = 0; depth < DEPTH_NUM; depth++)
		{
			double k2 = max(0.0001, 2. * kx / GRIDSIZE);  //k2 = 0..1
			double omega = sqrtf(GRAVITY * k * tanhf(k * Depth[depth]));
			omega *= 1.f / sqrt(2.0 / (k2 * PI) * sin(k2 * PI / 2.0));  // grid dispersion correction
			qtildehat_depth[depth][x].x = qtildehat[x].x * cos(omega * TIMESTEP) - omega / (kNonZero * kNonZero) * htildehat[x].x * sin(omega * TIMESTEP);
			qtildehat_depth[depth][x].y = qtildehat[x].y * cos(omega * TIMESTEP) - omega / (kNonZero * kNonZero) * htildehat[x].y * sin(omega * TIMESTEP);
		}
	}
	for (int depth = 0; depth < DEPTH_NUM; depth++)
		fftc1dinv(qtildehat_depth[depth]); // Back transform
	// interpolate surface velocity from the two closest water depth solutions
	for (int x = 0; x < GRIDSIZE; x++)
	{
		float waterDepth = max(hbar[x], hbar[x_plus]);
		int depth1 = 0;
		for (int depth = 0; depth < DEPTH_NUM; depth++)
			if (waterDepth >= Depth[depth])
				depth1 = depth;
		int depth2 = min(DEPTH_NUM - 1, depth1 + 1);
		float s = 0.f;
		if (depth1 != depth2)
			s = (Depth[depth2] - waterDepth) / (Depth[depth2] - Depth[depth1]);
		qtilde[x] = s * qtildehat_depth[depth1][x].x + (1.f - s) * qtildehat_depth[depth2][x].x;
	}
}

void Sim::SWEStep()
{
	// SWE bulk simulation using [Stelling03]

	// qbar to ubar using hbar from last timestep
	static float ubar_x[GRIDSIZE*GRIDSIZE];
	static float ubar_y[GRIDSIZE*GRIDSIZE];
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
	memcpy(hbarOld, hbar, GRIDSIZE * GRIDSIZE * sizeof(float));   // store current hbar for next timestep

	// Compute time derivative of u_bar
	static float ubarNew_x[GRIDSIZE*GRIDSIZE];
	static float ubarNew_y[GRIDSIZE*GRIDSIZE];
	float H[GRIDSIZE*GRIDSIZE] = terrain + hbar;
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
			// Note: The 1D implementaiton uses far more intermeidate values that the paper, and I don't know why. 
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
			// 	q_x_p15 *= hbar[min(idx(x+1,y) + 1, GRIDSIZE - 1)];
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
			// 	q_y_p15 *= hbar[min(idx(x,y+2), GRIDSIZE * GRIDSIZE - 1)];
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
			ubarNew_x[idx(x,y)] = ubar_x[idx(x,y)] + TIMESTEP * dux_dt;
			ubarNew_x[idx(x,y)] = LimitVelocity(ubarNew_x[idx(x,y)]);  // Enforcing CFL condition
			// Y DIRECTION
			float duy_dt = - (1/CELLSIZE) * ((q_y_0/h_y_p05) * (ubar_y[idx(x,y)] - ubar_y[idx(x,y-1)]) + (q_x_m05/h_y_p05) * (ubar_y[idx(x,y)] - ubar_y[idx(x-1,y)]) + GRAVITY * (H[idx(x,y+1)] - H[idx(x,y)]));
			ubarNew_y[idx(x,y)] = ubar_y[idx(x,y)] + TIMESTEP * duy_dt;
			ubarNew_y[idx(x,y)] = LimitVelocity(ubarNew_y[idx(x,y)]);  // Enforcing CFL condition
		}
	}
	ApplyBoundaries(ubarNew_x, BOUNDARY_TYPE);
	ApplyBoundaries(ubarNew_y, BOUNDARY_TYPE);

	// transfer back to flow rate using *most recent* hbar
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
	static float qtilde_x_dummy[GRIDSIZE*GRIDSIZE];
	static float qtilde_y_dummy[GRIDSIZE*GRIDSIZE];
	memcpy(qtilde_x_dummy, qtilde_x, GRIDSIZE * GRIDSIZE * sizeof(float));
	memcpy(qtilde_y_dummy, qtilde_y, GRIDSIZE * GRIDSIZE * sizeof(float));
	for (int y = 1; y < GRIDSIZE-1; y++)
	{
		for (int x = 1; x < GRIDSIZE-1; x++)
		{
			// extract row and column for cubic sampling
			float local_x[GRIDSIZE];
			float local_y[GRIDSIZE];
			ExtractLocal(local_x, qtilde_x_dummy, idx(x,y), false);
			ExtractLocal(local_y, qtilde_y_dummy, idx(x,y), true);

			// 
			float bulkVelocity_x = 0.5f * (ubarNew_x[idx(x,y)] + ubar_x[idx(x,y)]);
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
	// handle boundaries!!

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
	// handle boundaries!!
	

	// Update htilde from ubar divergence: dh/dt = -h * div(ubar)
	for (int y = 1; y < GRIDSIZE-1; y++)
	{
		for (int x = 1; x < GRIDSIZE-1; x++)
		{
			// backward difference to get divergence of ubar (using backward because h is on staggered grid)
			float div_ubar_x = (ubarNew[idx(x,y)] - ubarNew[idx(x-1,y)]) / CELLSIZE;
			float div_ubar_y = (ubarNew[idx(x,y)] - ubarNew[idx(x,y-1)]) / CELLSIZE;
			float div_ubar = div_ubar_x + div_ubar_y;
			// dampen if converging to avoid breaking waves
			if (div_ubar < 0.f)
				div_ubar *= GAMMA;
			htilde[idx(x,y)] *= exp(-div_ubar * TIMESTEP);
		}
	}

	// advection of h through ubar
	// first, construct q_advect = ubar * htilde sampled at cell edges using cubic sampling
	static float q_x_advect[GRIDSIZE*GRIDSIZE];
	static float q_y_advect[GRIDSIZE*GRIDSIZE];
	for (int y = 0; y < GRIDSIZE; y++)
	{
		for (int x = 0; x < GRIDSIZE; x++)
		{
			// extract row and column for cubic sampling
			float local_x[GRIDSIZE];
			float local_y[GRIDSIZE];
			ExtractLocal(local_x, htilde, idx(x,y), false);
			ExtractLocal(local_y, htilde, idx(x,y), true);

			// cubic reconstruction: 1/2 cell accounts for staggered grid, 0.5 * dt accounts for h and u at different 1/2 times
			float step_x = 0.5f*CELLSIZE - ubarNew_x[idx(x,y)] * 0.5f * TIMESTEP / CELLSIZE; // unitless (cells)
			float step_y = 0.5f*CELLSIZE - ubarNew_y[idx(x,y)] * 0.5f * TIMESTEP / CELLSIZE; 
			q_x_advect[idx(x,y)] = ubarNew_x[idx(x,y)] * SampleCubicClamped(x + step_x, local_x);  
			q_y_advect[idx(x,y)] = ubarNew_y[idx(x,y)] * SampleCubicClamped(y + step_y, local_y);  
		}
	}
	// use q_advect to update h using finite volume update: h_new = h_old + dt * (Del . (q + q_advect))
	float h_dummy[GRIDSIZE*GRIDSIZE];
	memcpy(h_dummy, h, GRIDSIZE * GRIDSIZE * sizeof(float));
	for (int y = 1; y < GRIDSIZE-1; y++)
	{
		for (int x = 1; x < GRIDSIZE-1; x++)
		{
			float q_x_l = LimitFlowRate(q_x_advect[idx(x-1,y)], h[idx(x-1,y)], h[idx(x,y)]);
			float q_x_r = LimitFlowRate(q_x_advect[idx(x,y)], h[idx(x,y)], h[idx(x+1,y)]);
			float q_y_l = LimitFlowRate(q_y_advect[idx(x,y-1)], h[idx(x,y-1)], h[idx(x,y)]);
			float q_y_r = LimitFlowRate(q_y_advect[idx(x,y)], h[idx(x,y)], h[idx(x,y+1)]);
			if ( ((h[idx(x-1,y)] == 0.f) && (h[idx(x,y)] == 0.f)) || (StopFlowOnTerrainBoundary(x-1, y, h, terrain)) )
				q_x_l = 0.f;
			if ( ((h[idx(x,y)] == 0.f) && (h[idx(x+1,y)] == 0.f)) || (StopFlowOnTerrainBoundary(x, y, h, terrain)) )
				q_x_r = 0.f;
			if ( ((h[idx(x,y-1)] == 0.f) && (h[idx(x,y)] == 0.f)) || (StopFlowOnTerrainBoundary(x, y-1, h, terrain)) )
				q_y_l = 0.f;
			if ( ((h[idx(x,y)] == 0.f) && (h[idx(x,y+1)] == 0.f)) || (StopFlowOnTerrainBoundary(x, y, h, terrain)) )
				q_y_r = 0.f;

			// update h using htilde advected through ubar
			float div_q = (q_x_l - q_x_r + q_y_l - q_y_r) / CELLSIZE;
			h[idx(x,y)] = max(0.f, h_dummy[idx(x,y)] + TIMESTEP * div_q);
		}
	}
	// handle boundaries!!
}
	
void Sim::ComputeValues()
{
	// Recombine bulk and surface flow
	for (int y = 1; y < GRIDSIZE - 1; y++)
	{
		for (int x = 1; x < GRIDSIZE - 1; x++)
		{
			q_x[idx(x,y)] = LimitFlowRate(qbar_x[idx(x,y)] + qtilde_x[idx(x,y)], h[idx(x,y)], h[idx(x+1,y)]);
			q_y[idx(x,y)] = LimitFlowRate(qbar_y[idx(x,y)] + qtilde_y[idx(x,y)], h[idx(x,y)], h[idx(x,y+1)]);
			if ( (StopFlowOnTerrainBoundary(x, y, h, terrain)) || (x == 0) || (x >= GRIDSIZE - 2) )
				q_x[idx(x,y)] = 0.f;
			if ( (StopFlowOnTerrainBoundary(x, y, h, terrain)) || (y == 0) || (y >= GRIDSIZE - 2) )
				q_y[idx(x,y)] = 0.f;
		}
	}
	// handle boundaries!!
	
	// height integration 
	for (int y = 1; y < GRIDSIZE - 1; y++)
	{
		for (int x = 1; x < GRIDSIZE - 1; x++)
		{
			h[idx(x,y)] = max(0.f, h[idx(x,y)] + TIMESTEP * -(q_x[idx(x,y)] - q_x[idx(x-1,y)]) / CELLSIZE);
		}
	}
	// handle boundaries!!
	
	// stability measure to not drag too much water from a cell in a single timestep (important for extreme initial conditions)
	for (int y = 1; y < GRIDSIZE - 1; y++)
	{
		for (int x = 1; x < GRIDSIZE - 1; x++)
		{
			q_x[idx(x,y)] = LimitFlowRate(q_x[idx(x,y)], h[idx(x,y)], h[idx(x+1,y)]);
			q_y[idx(x,y)] = LimitFlowRate(q_y[idx(x,y)], h[idx(x,y)], h[idx(x,y+1)]);
		}
	}
	// handle boundaries!!
}