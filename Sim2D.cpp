#include <windows.h>
#include "Sim2D.h"

// This is the main simulation code for the 2D solver. For now it's just the 1D code that I am using as reference. 

// ********************************************************************************************************************
// Helper functions (TODO: TRANSITION TO 2D)
// ********************************************************************************************************************

// 1D: 0.5 of the source cell content per timestep guarantees volume conservation
inline float LimitFlowRate(float flow_rate_in, float waterDepth_left, float waterDepth_right)
{
	if (flow_rate_in >= 0.f)
		return min(flow_rate_in, 0.25f * waterDepth_left * CELLSIZE / TIMESTEP);  // 0.25 since other neighbor might take from this source cell as well..
	else
		return max(flow_rate_in, -0.25f * waterDepth_right * CELLSIZE / TIMESTEP);
}


// 1D: 0.5 guarantees CFL condition
inline float LimitVelocity(float velocity_in)
{
	if (velocity_in >= 0.f)
		return min(velocity_in, 0.25f * CELLSIZE / TIMESTEP);   // 0.25 since other neighbors might take from this source cell as well..
	else
		return max(velocity_in, -0.25f * CELLSIZE / TIMESTEP);
}


// // xCoord and size in (0..1), factor determines how much to add/subtract
// void Sim::EditWaterLocal(float xCoord, float yCoord, float size, float factor)
// {
// 	for (int y = 0; y < GRIDSIZE; y++)
// 	{
// 		for (int x = 0; x < GRIDSIZE; x++)
// 			if (fabs((float)(x) / GRIDSIZE - xCoord) < size)
// 					h[idx] = max(0.f, h[idx] + factor * 1.f);
// 		h[0] = 0.0f;
// 		h[GRIDSIZE - 1] = 0.0f;
// 	}
// }


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
	const float s = 0.5f;
	float xcubicX = -s * x3 + 2.f * s * x2 - s * fx;
	float xcubicY = (2.f - s) * x3 + (s - 3.f) * x2 + 1.f;
	float xcubicZ = (s - 2.f) * x3 + (3.f - 2.f * s) * x2 + s * fx;
	float xcubicW = s * x3 - s * x2;
	float out = xcubicX * dataField[id0] + xcubicY * dataField[id1] + xcubicZ * dataField[id2] + xcubicW * dataField[id3];
	out = min(max(dataField[id1], dataField[id2]), out);  // value-limiting (same as in BFECC)
	out = max(min(dataField[id1], dataField[id2]), out);  // value-limiting (same as in BFECC)
	return out;
}


// test if the terrain boundary stops any flow across x+0.5
bool StopFlowOnTerrainBoundary(int x, int y, float* h, float* terrain)
{
	float epsilon = 0.01f;
	if ((h[idx] <= epsilon) && ((terrain[idx] >= terrain[idx_xplus] + h[idx_xplus]) || (terrain[idx] >= terrain[idx_yplus] + h[idx_yplus])))
		return true;
	if ((h[idx_xplus] <= epsilon) && ((terrain[idx_xplus] > terrain[idx] + h[idx]) || (terrain[idx_yplus] > terrain[idx] + h[idx])))
		return true;
	return false;
}


// Boundary Condition Functions
void HandleWallBoundary(float* field1, float* field2) 
{
	// Handle Top and Bottom Edges (Horizontal)
	for (int x = 0; x < GRIDSIZE; ++x) {
		// Top Edge
		field1[idx(x, 0)] = field1[idx(x, 1)];
		field2[idx(x, 0)] = field2[idx(x, 1)];
		// Bottom Edge
		field1[idx(x, GRIDSIZE - 1)] = field1[idx(x, GRIDSIZE - 2)];
		field2[idx(x, GRIDSIZE - 1)] = field2[idx(x, GRIDSIZE - 2)];
	}
	// Handle Left and Right Edges (Vertical)
	for (int y = 0; y < GRIDSIZE; ++y) {
		// Left Edge
		field1[idx(0, y)] = field1[idx(1, y)];
		field2[idx(0, y)] = field2[idx(1, y)];
		// Right Edge
		field1[idx(GRIDSIZE - 1, y)] = field1[idx(GRIDSIZE - 2, y)];
		field2[idx(GRIDSIZE - 1, y)] = field2[idx(GRIDSIZE - 2, y)];
	}
}

void HandleZeroBoundary(float* field1, float* field2)
{
	// Handle Top and Bottom Edges (Horizontal)
	for (int x = 0; x < GRIDSIZE; ++x) {
		// Top Edge
		field1[idx(x, 0)] = 0.f;
		field2[idx(x, 0)] = 0.f;
		// Bottom Edge
		field1[idx(x, GRIDSIZE - 1)] = 0.f;
		field2[idx(x, GRIDSIZE - 1)] = 0.f;
	}
	// Handle Left and Right Edges (Vertical)
	for (int y = 0; y < GRIDSIZE; ++y) {
		// Left Edge
		field1[idx(0, y)] = 0.f;
		field2[idx(0, y)] = 0.f;
		// Right Edge
		field1[idx(GRIDSIZE - 1, y)] = 0.f;
		field2[idx(GRIDSIZE - 1, y)] = 0.f;
	}
}

void HandleFreeBoundary(float* field1, float* field2)
{
	// Handle Top and Bottom Edges (Horizontal)
	for (int x = 0; x < GRIDSIZE; ++x) {
		// Top Edge
		field1[idx(x, 0)] = 2.f * field1[idx(x, 1)] - field1[idx(x, 2)];
		field2[idx(x, 0)] = 2.f * field2[idx(x, 1)] - field2[idx(x, 2)];
		// Bottom Edge
		field1[idx(x, GRIDSIZE - 1)] = 2.f * field1[idx(x, GRIDSIZE - 2)] - field1[idx(x, GRIDSIZE - 3)];
		field2[idx(x, GRIDSIZE - 1)] = 2.f * field2[idx(x, GRIDSIZE - 2)] - field2[idx(x, GRIDSIZE - 3)];
	}
	// Handle Left and Right Edges (Vertical)
	for (int y = 0; y < GRIDSIZE; ++y) {
		// Left Edge
		field1[idx(0, y)] = 2.f * field1[idx(1, y)] - field1[idx(2, y)];
		field2[idx(0, y)] = 2.f * field2[idx(1, y)] - field2[idx(2, y)];
		// Right Edge
		field1[idx(GRIDSIZE - 1, y)] = 2.f * field1[idx(GRIDSIZE - 2, y)] - field1[idx(GRIDSIZE - 3, y)];
		field2[idx(GRIDSIZE - 1, y)] = 2.f * field2[idx(GRIDSIZE - 2, y)] - field2[idx(GRIDSIZE - 3, y)];
	}
}

void ApplyBoundaries(float* field1, float* field2, bool isWall)
{
	if (isWall)
		HandleWallBoundary(field1, field2);
	else
		HandleFreeBoundary(field1, field2);
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
	SWEStep();
	TransportStep();
	ComputeValues();
	time += TIMESTEP;
}

void Sim::DecompositionStep(bool SWEonly)
{
	/******* Bulk vs Surface Wave Decomposition ******/

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
			alpha_H[idx] = 0.f;

			// // Their implementation: using averages (I assume to improve stability, but reduces accuracy?)
			// float maxGround = max(terrain[idx], terrain[idx_xplus], terrain[idx_yplus]); // Why do this?
			// float minWaterlevel = (H[idx] + H[idx_xplus] + H[idx_yplus]) / 3.f; // Why average here?
			// if ((h[idx] > 0.f) && (h[idx_xplus] > 0.f) && (h[idx_yplus] > 0.f))
			// {
			// 	static const float sigma_max = 8.f;
			// // they limit diffusion coefficient to between 0 and 1, but that isn't requred by the math = maybe for stability
			// 	float sigma = min(sigma_max, max(0.f, minWaterlevel - maxGround));
			// 	alpha_H[idx] = sigma * sigma / (2*DELTA_T*DIFFUSION_ITERATIONS);
			// } 

			// My implementation: using local cell values only to align with eqn from paper
			if (h[idx] > 0.f)
			{
				float denom = 2*DELTA_T*DIFFUSION_ITERATIONS;
				alpha_H[idx] = h[idx] * h[idx] / denom;
				alpha_H[idx] = min(std::sqrt(denom), alpha_H[idx]); // clamp to improve stability; limits max depth to 
			}
			
			// Extra gradient filter
			// NOTE: they used H, I switched it to h to stay strict with the paper. We'll see if this causes bugs. 
			float gradient_x = (h[idx_xplus] - h[idx]) / CELLSIZE; // could use central difference here
			float gradient_y = (h[idx_yplus] - h[idx]) / CELLSIZE;
			alpha_H[idx] *= exp(- DIFFUSION_PENALTY * (gradient_x * gradient_x + gradient_y * gradient_y));
			alpha_Q[idx] = 0.5f * (alpha_H[idx_xminus] + alpha_H[idx]); // Where does this come from??
		}
	}
	ApplyBoundaries(alpha_H, alpha_Q, true);
	
	// Run diffusion to low-pass filter H and Q
	// SOMEDAY: Improve this implementation of diffusion by replacing Euler integration with FFT or something
	static float H_past[GRIDSIZE*GRIDSIZE];
	static float Q_past_x[GRIDSIZE*GRIDSIZE];
	static float Q_past_y[GRIDSIZE*GRIDSIZE];
	for (int j = 0; (j < DIFFUSION_ITERATIONS) && (!SWEonly); j++)  // 64 diffusion iterations
	{
		memcpy(H_past, H, GRIDSIZE * GRIDSIZE * sizeof(float));
		memcpy(Q_past_x, Q_x, GRIDSIZE * GRIDSIZE * sizeof(float));
		memcpy(Q_past_y, Q_y, GRIDSIZE * GRIDSIZE * sizeof(float));
		for (int y = 1; y < GRIDSIZE-1; y++) // one diffusion iteration
		{
			for (int x = 1; x < GRIDSIZE - 1; x++)
			{
				// Diffusion step for H: dH/dt = Del * ( alpha_H * Del H )
				float dH_x = (alpha_H[idx] * (H_past[idx_xplus] - H_past[idx]) - alpha_H[idx_xminus] * (H_past[idx] - H_past[idx_xminus]));
				float dH_y = (alpha_H[idx] * (H_past[idx_yplus] - H_past[idx]) - alpha_H[idx_yminus] * (H_past[idx] - H_past[idx_yminus]));
				float dHdT = (dH_x + dH_y) / (CELLSIZE*CELLSIZE);
				H[idx] = H_past[idx] + DELTA_T * dHdT;
				H[idx] = max(terrain[idx], H[idx]); // ensure water surface is above terrain

				// Diffusion step for Q: dQ/dt = Del * ( alpha_Q * Del Q )
				// Q has two components, so we do them separately
				float dQ_x_x = (alpha_Q[idx] * (Q_past_x[idx_xplus] - Q_past_x[idx]) - alpha_Q[idx_xminus] * (Q_past_x[idx] - Q_past_x[idx_xminus]));
				float dQ_x_y = (alpha_Q[idx] * (Q_past_x[idx_yplus] - Q_past_x[idx]) - alpha_Q[idx_yminus] * (Q_past_x[idx] - Q_past_x[idx_yminus]));
				float dQdT_x = (dQ_x_x + dQ_x_y) / (CELLSIZE*CELLSIZE);
				Q_x[idx] = Q_past_x[idx] + DELTA_T * dQdT_x;
				float dQ_y_x = (alpha_Q[idx] * (Q_past_y[idx_yplus] - Q_past_y[idx]) - alpha_Q[idx_yminus] * (Q_past_y[idx] - Q_past_y[idx_yminus]));
				float dQ_y_y = (alpha_Q[idx] * (Q_past_y[idx_xplus] - Q_past_y[idx]) - alpha_Q[idx_xminus] * (Q_past_y[idx] - Q_past_y[idx_xminus]));
				float dQdT_y = (dQ_y_x + dQ_y_y) / (CELLSIZE*CELLSIZE);
				Q_y[idx] = Q_past_y[idx] + DELTA_T * dQdT_y;
			}
		}
	}
	ApplyBoundaries(H, Q, true);

	// final conversion to individual solver quantities
	hbar = max(0.f, H - terrain);
	qbar_x = Q_x;
	qbar_y = Q_y;
	htilde = h - hbar;
	qtilde_x = q_x - qbar_x;
	qtilde_y = q_y - qbar_y;

	// Enforce no-flow conditions at terrain boundaries
	for (int y = 0; y < GRIDSIZE; y++)
	{
		for (int x = 0; x < GRIDSIZE; x++)
		{
			if (StopFlowOnTerrainBoundary(x, y, h, terrain))
			{
				qbar_x[idx] = 0.f;
				qbar_y[idx] = 0.f;
				qtilde_x[idx] = 0.f;
				qtilde_y[idx] = 0.f;
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
	fftc1d(htildehat);   https://www.alglib.net/download.php#cpp
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
	static float ubar[GRIDSIZE];
	static float ubarNew[GRIDSIZE];
	for (int x = 0; x < GRIDSIZE; x++)
	{
		ubar[x] = qbar[x];
		if (ubar[x] >= 0.f)
			ubar[x] /= max(0.01f, hbarOld[x]);
		else
			ubar[x] /= max(0.01f, hbarOld[x_plus]);
		ubar[x] = LimitVelocity(ubar[x]);  // CFL<0.25 will be important later for surface waves advection
	}
	memcpy(hbarOld, hbar, GRIDSIZE * sizeof(float));   // store current hbar for next timestep
	for (int x = 0; x < GRIDSIZE; x++)
	{
		float q_m05 = ubar[x_minus];  
		if (q_m05 >= 0.f)
			q_m05 *= hbar[x_minus];
		else
			q_m05 *= hbar[x];
		float q_p05 = ubar[x];  
		if (q_p05 >= 0.f)
			q_p05 *= hbar[x];
		else
			q_p05 *= hbar[x_plus];
		float q_p15 = ubar[x_plus];  //q_(i+0.5) = hfr at position x
		if (q_p15 >= 0.f)
			q_p15 *= hbar[x_plus];
		else
			q_p15 *= hbar[min(x + 2, GRIDSIZE - 1)];
		float q_bar_0 = 0.5f * (q_m05 + q_p05);
		float q_bar_p1 = 0.5f * (q_p05 + q_p15);
		float u_star_0 = 0.f;
		if (q_bar_0 >= 0.f)
			u_star_0 = ubar[max(x - 1, 0)];
		else
			u_star_0 = ubar[x];
		float u_star_p1 = 0.f;
		if (q_bar_p1 > 0.f)
			u_star_p1 = ubar[x];
		else
			u_star_p1 = ubar[x_plus];
		float uu_x = 2.f / max(0.01f, hbar[x] + hbar[x_plus]) * ((q_bar_p1 * u_star_p1 - q_bar_0 * u_star_0) / CELLSIZE - ubar[x] * (q_bar_p1 - q_bar_0) / CELLSIZE);
		ubarNew[x] = ubar[x] - TIMESTEP * uu_x;  // self-advection
		ubarNew[x] += -GRAVITY * TIMESTEP * (terrain[x_plus] + hbar[x_plus] - terrain[x] - hbar[x]) / CELLSIZE;  // GRAVITY force
		ubarNew[x] = LimitVelocity(ubarNew[x]);
	}
	// transfer back to flow rate using *most recent* hbar
	for (int x = 0; x < GRIDSIZE; x++)
		if (ubarNew[x] >= 0.f)
			qbar[x] = ubarNew[x] * hbar[x];
		else
			qbar[x] = ubarNew[x] * hbar[x_plus];
}

void Sim::TransportStep()
{
	// advect surface HFR through bulk velocity
	static float qtilde_dummy[GRIDSIZE];
	memcpy(qtilde_dummy, qtilde, GRIDSIZE * sizeof(float));
	for (int x = 0; x < GRIDSIZE; x++)
	{
		float bulkVelocity = 0.5f * (ubarNew[x] + ubar[x]);
		qtilde[x] = SampleCubicClamped(x - TIMESTEP * bulkVelocity, qtilde_dummy);
		if (((bulkVelocity >= 0.f) && (h[x] < 0.01f)) ||
			((bulkVelocity < 0.f) && (h[min(x + 1, GRIDSIZE - 1)] < 0.01f)))
			qtilde[x] = 0.f;
	}

	// div(ubar) qtilde
	for (int x = 0; x < GRIDSIZE; x++)
	{
		float ubar_m1 = 0.5f * (ubarNew[x_minus] + ubar[x_minus]);
		float ubar_p1 = 0.5f * (ubarNew[x_plus] + ubar[x_plus]);
		float div_ubar = (ubar_p1 - ubar_m1) / (2.f * CELLSIZE);
		if (div_ubar < 0.f)
			div_ubar *= 0.25f;
		qtilde[x] *= exp(-div_ubar * TIMESTEP);
	}

	// div(ubar) htilde
	for (int x = 0; x < GRIDSIZE; x++)
	{
		float div_ubar = (ubarNew[x] - ubarNew[x_minus]) / CELLSIZE;
		if (div_ubar < 0.f)
			div_ubar *= 0.25f;
		htilde[x] *= exp(-div_ubar * TIMESTEP);
	}
}

void Sim::ComputeValues()
{
	// bulk advecting surface displacements
	static float advectHFR[GRIDSIZE];
	for (int x = 0; x < GRIDSIZE; x++)
		advectHFR[x] = ubarNew[x] * SampleCubicClamped(x + 0.5f - 0.5f * TIMESTEP * ubarNew[x], htilde);  // cubic reconstruction: 0.5 * dt accounts for h and u at different 1/2 times
	float h_dummy[GRIDSIZE];
	memcpy(h_dummy, h, GRIDSIZE * sizeof(float));
	for (int x = 0; x < GRIDSIZE; x++)
	{
		float q_l = Limit_flow_rate(advectHFR[x_minus], h_dummy[x_minus], h_dummy[x]);
		float q_r = Limit_flow_rate(advectHFR[x], h_dummy[x], h_dummy[x_plus]);
		if ( ((h_dummy[x_minus] == 0.f) && (h_dummy[x] == 0.f)) || (StopFlowOnTerrainBoundary(x_minus, h_dummy, terrain)) )
			q_l = 0.f;
		if ( ((h_dummy[x] == 0.f) && (h_dummy[x_plus] == 0.f)) || (StopFlowOnTerrainBoundary(x, h_dummy, terrain)) )
			q_r = 0.f;
		h[x] = max(0.f, h_dummy[x] - TIMESTEP / CELLSIZE * (q_r - q_l));
	}

	// Recombine bulk and surface HFR
	for (int x = 0; x < GRIDSIZE; x++)
	{
		q[x] = Limit_flow_rate(qbar[x] + qtilde[x], h[x], h[x_plus]);
		if ( (StopFlowOnTerrainBoundary(x, h, terrain)) || (x == 0) || (x >= GRIDSIZE - 2) )
			q[x] = 0.f;
	}

	// height integration 
	for (int x = 0; x < GRIDSIZE; x++)
		h[x] = max(0.f, h[x] + TIMESTEP * -(q[x] - q[x_minus]) / CELLSIZE);

	// stability measure to not drag too much water from a cell in a single timestep (important for extreme initial conditions)
	for (int x = 0; x < GRIDSIZE; x++)
		q[x] = Limit_flow_rate(q[x], h[x], h[x_plus]);
}