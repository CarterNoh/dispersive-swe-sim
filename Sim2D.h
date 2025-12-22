#pragma once

#include "alglib\fasttransforms.h"  // https://www.alglib.net/download.php#cpp

#define GRAVITY 9.80665
#define PI 3.14159265359  // used in FFT stencil code

// sim parameters
#define GRIDRESOLUTION 256
#define GRIDCELLSIZE 1		// this should not be changed in this implementation!
#define TIMESTEP (1.f/60.f)
#define DEPTH_NUM 4
const float Depth[DEPTH_NUM] = { 1.f, 4.f, 16.f, 64.f };
#define TERRAIN_HEIGHT_SHIFT_INIT -10.f // -20.f //-10.f
#define TERRAIN_HEIGHT_SCALE_INIT 20.f //40.f  //20.f
// diffusion parameters
#define DIFFUSION_ITERATIONS 64
#define DELTA_T 0.5f

// helpful shortcuts
#define x_plus min(GRIDRESOLUTION - 1, x + 1)
#define x_minus max(0, x - 1)


class Sim
{
public:
	// variables carried from one timestep to the next
	float	terrain[GRIDRESOLUTION];						// terrain
	float   h[GRIDRESOLUTION];								// overall water height
	float	q[GRIDRESOLUTION];								// overall flow rate
	float	hbarOld[GRIDRESOLUTION];						// last timestep hbar, used for resampling in time 
	float	htildeOld[GRIDRESOLUTION];						// last timestep htilde, used for resampling in time 

	// variables that could be allocated locally but for potential visualizations we store them globally
	float	hbar[GRIDRESOLUTION]; 							// bulk height
	float	qbar[GRIDRESOLUTION];							// bulk flow rate
	float	htilde[GRIDRESOLUTION];							// surface displacement
	float	qtilde[GRIDRESOLUTION];							// surface flow rate
	alglib::complex_1d_array htildehat, qtildehat;			// eWave inputs
	alglib::complex_1d_array qtildehat_depth[DEPTH_NUM];	// eWave outputs

	// time is exclusively used for video recording
	float time;

	// functions
	Sim();
	int Sim::Release(void);
	void ResetTerrain(int type);
	void ResetWater(int type, float level);
	void SimStep(bool SWEonly);													// advects the simulation by one timestep
	void EditWaterLocal(float xCoord, float size, float factor);	// add or subtract water locally.
};
