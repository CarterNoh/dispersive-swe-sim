
#include <windows.h>
#include "Sim.h"


// ********************************************************************************************************************
// Helper functions
// ********************************************************************************************************************

// 1D: 0.5 of the source cell content per timestep guarantees volume conservation
inline float Limit_flow_rate(float flow_rate_in, float waterDepth_left, float waterDepth_right)
{
	if (flow_rate_in >= 0.f)
		return min(flow_rate_in, 0.5f * waterDepth_left * GRIDCELLSIZE / TIMESTEP);  // 0.5 since other neighbor might take from this source cell as well..
	else
		return max(flow_rate_in, -0.5f * waterDepth_right * GRIDCELLSIZE / TIMESTEP);
}


// 1D: 0.5 guarantees CFL condition
inline float LimitVelocity(float velocity_in)
{
	if (velocity_in >= 0.f)
		return min(velocity_in, 0.5f * GRIDCELLSIZE / TIMESTEP);   // 0.5 since other neighbor might take from this source cell as well..
	else
		return max(velocity_in, -0.5f * GRIDCELLSIZE / TIMESTEP);
}


// xCoord and size in (0..1), factor determines how much to add/subtract
void Sim::EditWaterLocal(float xCoord, float size, float factor)
{
	for (int x = 0; x < GRIDRESOLUTION; x++)
		if (fabs((float)(x) / GRIDRESOLUTION - xCoord) < size)
				h[x] = max(0.f, h[x] + factor * 1.f);
	h[0] = 0.0f;
	h[GRIDRESOLUTION - 1] = 0.0f;
}


// cubic interpolation with Catmull-Rom Spline https://en.wikipedia.org/wiki/Cubic_Hermite_spline#Interpolating_a_data_set
float SampleCubicClamped(float samplePos, float* dataField)
{
	int id_start = floor(samplePos) - 1;
	int id0 = max(0, min(id_start + 0, GRIDRESOLUTION - 1));
	int id1 = max(0, min(id_start + 1, GRIDRESOLUTION - 1));
	int id2 = max(0, min(id_start + 2, GRIDRESOLUTION - 1));
	int id3 = max(0, min(id_start + 3, GRIDRESOLUTION - 1));
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
bool StopFlowOnTerrainBoundary(int x, float* h, float* terrain)
{
	float epsilon = 0.01f;
	if ((h[x] <= epsilon) && (terrain[x] >= terrain[min(GRIDRESOLUTION - 1, x + 1)] + h[min(GRIDRESOLUTION - 1, x + 1)]))
		return true;
	if ((h[min(GRIDRESOLUTION - 1, x + 1)] <= epsilon) && (terrain[min(GRIDRESOLUTION - 1, x + 1)] > terrain[x] + h[x]))
		return true;
	return false;
}

// ********************************************************************************************************************
// Init functions 
// ********************************************************************************************************************

// type: 0=flat, 1=hill
void Sim::ResetTerrain(int type)
{
	for (int x = 0; x < GRIDRESOLUTION; x++)
		if (type == 0)  // flat
			terrain[x] = -abs(TERRAIN_HEIGHT_SHIFT_INIT);
		else if (type == 1) // hill
			terrain[x] = (-1.f + 0.1f + 0.1f * x / GRIDRESOLUTION + 0.03f * sin(20.f * x / GRIDRESOLUTION) + 0.9f * sin(2.5f * x / GRIDRESOLUTION)) * abs(TERRAIN_HEIGHT_SHIFT_INIT);
			//terrain[x] = (-1.f + 0.1f + 0.5f * (0.1f * x / GRIDRESOLUTION + 0.03f * sin(20.f * x / GRIDRESOLUTION) + 0.9f * sin(2.5f * x / GRIDRESOLUTION))) * abs(TERRAIN_HEIGHT_SHIFT_INIT);
	terrain[0] = 1.8f * abs(TERRAIN_HEIGHT_SHIFT_INIT);
	terrain[GRIDRESOLUTION - 1] = 1.8f * abs(TERRAIN_HEIGHT_SHIFT_INIT);
}


// type: 0=constant level, 1=dam break, 2=sloped  level = y coordinate of water in domain
void Sim::ResetWater(int type, float level)
{
	for (int x = 0; x < GRIDRESOLUTION; x++)
	{
		if (type == 0) //constant level
			h[x] = max(0.f, (level - terrain[x]));
		if (type == 1)  // dam break
			if (x <= GRIDRESOLUTION / 2)
				h[x] = 0.f;
			else
				h[x] = max(0.f, (level - terrain[x]));
		if (type == 2)  // sloped water
			h[x] = max(0.f, level + (2.f * (-0.5f + (float)(x) / GRIDRESOLUTION) * fabs(-0.5f + (float)(x) / GRIDRESOLUTION)) * abs(0.5f * TERRAIN_HEIGHT_SHIFT_INIT) - terrain[x]);
		if (type == 3)  // flat with initial surface waves
		{
			float lambda = 10.f;
			h[x] = max(0.f, (level + 0.5f * cos(2.f * PI * (x / lambda)) - terrain[x]));
		}
		hbar[x] = h[x];
		hbarOld[x] = h[x];
		htilde[x] = 0.f;
		htildeOld[x] = 0.f;
		qbar[x] = 0.f;
		q[x] = 0.f;
	}
	// clear left and right boundaries
	h[0] = 0.f;
	h[GRIDRESOLUTION - 1] = 0.f;
	time = 0.f;
}


Sim::Sim()
{
	htildehat.setlength(GRIDRESOLUTION);
	qtildehat.setlength(GRIDRESOLUTION);
	for (int i=0; i < DEPTH_NUM; i++)
		qtildehat_depth[i].setlength(GRIDRESOLUTION);
	ResetTerrain(1);
	ResetWater(2, 0.f);
}


int Sim::Release(void)
{
	return 0;
}


// ********************************************************************************************************************
// simulation function
// ********************************************************************************************************************


// simulate a timestep
void Sim::SimStep(bool SWEonly)
{
	//************************
	// bulk vs surface decomposition
	static float alpha_H[GRIDRESOLUTION]; // = zeros
	static float alpha_Q[GRIDRESOLUTION];
	static float H[GRIDRESOLUTION];
	static float Q[GRIDRESOLUTION];
	// h,q diffusivity intialization
	H = terrain + h;  // start off with the current water surface
	Q = q;
	for (int x = 0; x < GRIDRESOLUTION; x++)
	{
		alpha_H[x] = 0.f;
		float maxGround = max(terrain[x], terrain[x_plus])
		float minWaterlevel = 0.5f * (H[x] + H[x_plus]);
		if ((h[x] > 0.f) && (h[x_plus] > 0.f))
		{
			static const float sigma_max = 8.f;
			float sigma = min(sigma_max, max(0.f, minWaterlevel - maxGround));
			alpha_H[x] = sigma * sigma / (2*DELTA_T*DIFFUSION_ITERATIONS);
		} 
		// gradient filter
		float gradient = abs(H[x] - H[x_plus]);
		alpha_H[x] *= exp(- 0.01f * gradient * gradient);
		alpha_Q[x] = 0.5f * (alpha_H[max(0, x - 1)] + alpha_H[x]);
	}
	// run diffusion
	static float H_dummy[GRIDRESOLUTION];
	static float Q_dummy[GRIDRESOLUTION];
	for (int j = 0; (j < DIFFUSION_ITERATIONS) && (!SWEonly); j++)  // 64 diffusion iterations
	{
		memcpy(H_dummy, H, GRIDRESOLUTION * sizeof(float));
		memcpy(Q_dummy, Q, GRIDRESOLUTION * sizeof(float));
		for (int x = 0; x < GRIDRESOLUTION - 1; x++) // one diffusion iteration
		{
			H[x] = max(terrain[x], H_dummy[x] + DELTA_T * (alpha_H[x] * (H_dummy[x_plus] - H_dummy[x]) - alpha_H[x_minus] * (H_dummy[x] - H_dummy[x_minus])));
			Q[x] = Q_dummy[x] + DELTA_T * (alpha_Q[x_plus] * (Q_dummy[x_plus] - Q_dummy[x]) - alpha_Q[x] * (Q_dummy[x] - Q_dummy[x_minus]));
		}
	}
	// final conversion to individual solver quantities
	hbar = max(0.f, H - terrain);
	qbar = Q;
	htilde = h - hbar;
	qtilde = q - qbar;
	for (int x = 0; x < GRIDRESOLUTION; x++)
	{
		if (StopFlowOnTerrainBoundary(x, h, terrain))
		{
			qbar[x] = 0.f;
			qtilde[x] = 0.f;
		}
	}

	//************************
	// surface velocity update using eWave
	for (int x = 0; x < GRIDRESOLUTION; x++)
	{
		htildehat[x].x = 0.5f * (htilde[x] + htildeOld[x]);
		htildeOld[x] = htilde[x];
		htildehat[x].y = 0.;
		qtildehat[x].x = qtilde[x];
		qtildehat[x].y = 0.;
	}
	fftc1d(htildehat);   https://www.alglib.net/download.php#cpp
	fftc1d(qtildehat);
	for (int x = 0; x < GRIDRESOLUTION; x++)
	{
		// physical k from grid position
		double kx = GRIDRESOLUTION / 2. - abs(GRIDRESOLUTION / 2. - x);  // this gives [0,..,m_gridSizeX / 2.f-1, m_gridSizeX / 2.f, .. 1]
		double k = 2. * PI * fabs(kx) / GRIDRESOLUTION / GRIDCELLSIZE;
		double kNonZero = max(0.01, k);
		double kS = k;  // signed k
		if (x > (double)(GRIDRESOLUTION) / 2.f)
			kS = -k;
		// Fourier gradient: multiply by -i k
		double real = htildehat[x].x;
		double imag = htildehat[x].y;
		htildehat[x].x = -kS * imag;
		htildehat[x].y = kS * real;
		// phase shift to translate function to cell boundaries
		real = htildehat[x].x;
		imag = htildehat[x].y;
		double beta = 0.5 * GRIDCELLSIZE * kS;
		htildehat[x].x = cos(beta) * real - sin(beta) * imag;
		htildehat[x].y = sin(beta) * real + cos(beta) * imag;
		for (int depth = 0; depth < DEPTH_NUM; depth++)
		{
			double k2 = max(0.0001, 2. * kx / GRIDRESOLUTION);  //k2 = 0..1
			double omega = sqrtf(GRAVITY * k * tanhf(k * Depth[depth]));
			omega *= 1.f / sqrt(2.0 / (k2 * PI) * sin(k2 * PI / 2.0));  // grid dispersion correction
			qtildehat_depth[depth][x].x = qtildehat[x].x * cos(omega * TIMESTEP) - omega / (kNonZero * kNonZero) * htildehat[x].x * sin(omega * TIMESTEP);
			qtildehat_depth[depth][x].y = qtildehat[x].y * cos(omega * TIMESTEP) - omega / (kNonZero * kNonZero) * htildehat[x].y * sin(omega * TIMESTEP);
		}
	}
	for (int depth = 0; depth < DEPTH_NUM; depth++)
		fftc1dinv(qtildehat_depth[depth]); // Back transform
	// interpolate surface velocity from the two closest water depth solutions
	for (int x = 0; x < GRIDRESOLUTION; x++)
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

	//************************
	// SWE bulk simulation using [Stelling03]
	// qbar to ubar using hbar from last timestep
	static float ubar[GRIDRESOLUTION];
	static float ubarNew[GRIDRESOLUTION];
	for (int x = 0; x < GRIDRESOLUTION; x++)
	{
		ubar[x] = qbar[x];
		if (ubar[x] >= 0.f)
			ubar[x] /= max(0.01f, hbarOld[x]);
		else
			ubar[x] /= max(0.01f, hbarOld[x_plus]);
		ubar[x] = LimitVelocity(ubar[x]);  // CFL<0.5 will be important later for surface waves advection
	}
	memcpy(hbarOld, hbar, GRIDRESOLUTION * sizeof(float));   // store current hbar for next timestep
	for (int x = 0; x < GRIDRESOLUTION; x++)
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
			q_p15 *= hbar[min(x + 2, GRIDRESOLUTION - 1)];
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
		float uu_x = 2.f / max(0.01f, hbar[x] + hbar[x_plus]) * ((q_bar_p1 * u_star_p1 - q_bar_0 * u_star_0) / GRIDCELLSIZE - ubar[x] * (q_bar_p1 - q_bar_0) / GRIDCELLSIZE);
		ubarNew[x] = ubar[x] - TIMESTEP * uu_x;  // self-advection
		ubarNew[x] += -GRAVITY * TIMESTEP * (terrain[x_plus] + hbar[x_plus] - terrain[x] - hbar[x]) / GRIDCELLSIZE;  // GRAVITY force
		ubarNew[x] = LimitVelocity(ubarNew[x]);
	}
	// transfer back to flow rate using *most recent* hbar
	for (int x = 0; x < GRIDRESOLUTION; x++)
		if (ubarNew[x] >= 0.f)
			qbar[x] = ubarNew[x] * hbar[x];
		else
			qbar[x] = ubarNew[x] * hbar[x_plus];

	//************************
	// advect surface HFR through bulk velocity
	static float qtilde_dummy[GRIDRESOLUTION];
	memcpy(qtilde_dummy, qtilde, GRIDRESOLUTION * sizeof(float));
	for (int x = 0; x < GRIDRESOLUTION; x++)
	{
		float bulkVelocity = 0.5f * (ubarNew[x] + ubar[x]);
		qtilde[x] = SampleCubicClamped(x - TIMESTEP * bulkVelocity, qtilde_dummy);
		if (((bulkVelocity >= 0.f) && (h[x] < 0.01f)) ||
			((bulkVelocity < 0.f) && (h[min(x + 1, GRIDRESOLUTION - 1)] < 0.01f)))
			qtilde[x] = 0.f;
	}

	//************************
	// div(ubar) qtilde
	for (int x = 0; x < GRIDRESOLUTION; x++)
	{
		float ubar_m1 = 0.5f * (ubarNew[x_minus] + ubar[x_minus]);
		float ubar_p1 = 0.5f * (ubarNew[x_plus] + ubar[x_plus]);
		float div_ubar = (ubar_p1 - ubar_m1) / (2.f * GRIDCELLSIZE);
		if (div_ubar < 0.f)
			div_ubar *= 0.25f;
		qtilde[x] *= exp(-div_ubar * TIMESTEP);
	}

	// **********************************
	// div(ubar) htilde
	for (int x = 0; x < GRIDRESOLUTION; x++)
	{
		float div_ubar = (ubarNew[x] - ubarNew[x_minus]) / GRIDCELLSIZE;
		if (div_ubar < 0.f)
			div_ubar *= 0.25f;
		htilde[x] *= exp(-div_ubar * TIMESTEP);
	}

	// **********************************
	// bulk advecting surface displacements
	static float advectHFR[GRIDRESOLUTION];
	for (int x = 0; x < GRIDRESOLUTION; x++)
		advectHFR[x] = ubarNew[x] * SampleCubicClamped(x + 0.5f - 0.5f * TIMESTEP * ubarNew[x], htilde);  // cubic reconstruction: 0.5 * dt accounts for h and u at different 1/2 times
	float h_dummy[GRIDRESOLUTION];
	memcpy(h_dummy, h, GRIDRESOLUTION * sizeof(float));
	for (int x = 0; x < GRIDRESOLUTION; x++)
	{
		float q_l = Limit_flow_rate(advectHFR[x_minus], h_dummy[x_minus], h_dummy[x]);
		float q_r = Limit_flow_rate(advectHFR[x], h_dummy[x], h_dummy[x_plus]);
		if ( ((h_dummy[x_minus] == 0.f) && (h_dummy[x] == 0.f)) || (StopFlowOnTerrainBoundary(x_minus, h_dummy, terrain)) )
			q_l = 0.f;
		if ( ((h_dummy[x] == 0.f) && (h_dummy[x_plus] == 0.f)) || (StopFlowOnTerrainBoundary(x, h_dummy, terrain)) )
			q_r = 0.f;
		h[x] = max(0.f, h_dummy[x] - TIMESTEP / GRIDCELLSIZE * (q_r - q_l));
	}

	// **********************************
	// Recombine bulk and surface HFR
	for (int x = 0; x < GRIDRESOLUTION; x++)
	{
		q[x] = Limit_flow_rate(qbar[x] + qtilde[x], h[x], h[x_plus]);
		if ( (StopFlowOnTerrainBoundary(x, h, terrain)) || (x == 0) || (x >= GRIDRESOLUTION - 2) )
			q[x] = 0.f;
	}

	// **********************************
	// height integration 
	for (int x = 0; x < GRIDRESOLUTION; x++)
		h[x] = max(0.f, h[x] + TIMESTEP * -(q[x] - q[x_minus]) / GRIDCELLSIZE);

	// **********************************
	// stability measure to not drag too much water from a cell in a single timestep (important for extreme initial conditions)
	for (int x = 0; x < GRIDRESOLUTION; x++)
		q[x] = Limit_flow_rate(q[x], h[x], h[x_plus]);

	time += TIMESTEP;
}


