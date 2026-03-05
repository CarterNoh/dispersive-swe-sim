#pragma once

#include <vector>
// #include "alglib\fasttransforms.h"  // https://www.alglib.net/download.php#cpp


// sim parameters
#define GRIDSIZE 256 	// grid size in one dimension (# cells)
#define CELLSIZE 1		// cell size in one dimension (meters?/cell)
#define TIMESTEP (1.f/60.f)
#define DEPTH_NUM 4
const float Depth[DEPTH_NUM] = { 1.f, 4.f, 16.f, 64.f };
#define TERRAIN_HEIGHT_SHIFT_INIT -10.f // -10 or -20 
#define TERRAIN_HEIGHT_SCALE_INIT 20.f // 20 or 40
#define CFL_CONDITION 0.25f  // max allowed CFL condition for stability
#define MIN_WATER_HEIGHT 0.01f  // minimum water height for stability
#define BOUNDARY_TYPE 1 // 0 = wall, 1 = free, 2 = zero
#define TERRAIN_TYPE 1 // 0 = flat, 1 = hill

// diffusion parameters
#define DIFFUSION_ITERATIONS 128
#define DELTA_T 0.25f
#define DIFFUSION_PENALTY 0.01f

// transport parameters
#define GAMMA 0.25f

// helpful shortcuts
#define GRAVITY 9.80665
#define PI 3.14159265359  // used in FFT stencil code

class Sim
{
public:
	// Functions
	Sim();
	int Release(void);
	void SimStep(bool SWEonly);	// ticks the simulation by one timestep using the following substeps:
	// void ResetTerrain(int type);
	// void ResetWater(int type, float level);
	// void EditWaterLocal(float xCoord, float yCoord, float size, float factor);	// add or subtract water locally.

	// Variables carried from one timestep to the next
	std::vector<float> h;		// overall water height
	std::vector<float> terrain;	// terrain
	std::vector<float> q_x;		// overall flow rate
	std::vector<float> q_y;		// overall flow rate
	float time;
	
private:
	// Functions
	void DecompositionStep(bool SWEonly); 	// bulk vs surface decomposition
	// void eWaveStep(bool SWEonly);			// surface wave simulation step
	void SWEStep();							// SWE bulk simulation step
	void TransportStep();					// transport of bulk and surface quantities
	void ComputeValues();					// compute final h and q values

	// variables carried from one timestep to the next
	std::vector<float> hbarOld;		// last timestep hbar, used for resampling in time
	std::vector<float> htildeOld;	// last timestep htilde, used for resampling in time

	// variables that could be allocated locally but we make gloabl so all functions can modify
	std::vector<float> hbar;		// bulk height
	std::vector<float> qbar_x;		// bulk flow rate
	std::vector<float> qbar_y;		// bulk flow rate
	std::vector<float> htilde;		// surface height
	std::vector<float> qtilde_x;	// surface flow rate												
	std::vector<float> qtilde_y;	// surface flow rate
	std::vector<float> ubar_x;		// bulk velocity
	std::vector<float> ubar_y;		// bulk velocity
	std::vector<float> ubarNew_x;	// bulk velocity after SWE step but before transport step
	std::vector<float> ubarNew_y;	// bulk velocity after SWE step but before transport step

	// Decomposition Variables
	std::vector<float> alpha_H; 	// diffusion coefficient for height
	std::vector<float> alpha_Q_x; 	// diffusion coefficient for flow rate in x direction
	std::vector<float> alpha_Q_y; 	// diffusion coefficient for flow rate in y direction
	std::vector<float> H; 			// Combined water + terrain height used for diffusion
	std::vector<float> Q_x; 		// current flow rate in x direction used for diffusion
	std::vector<float> Q_y; 		// current flow rate in y direction used for diffusion
	std::vector<float> HPast; 		// Combined water + terrain height from last diffusion iteration
	std::vector<float> QPast_x; 	// flow rate in x direction from last diffusion iteration
	std::vector<float> QPast_y; 	// flow rate in y direction from last diffusion iteration

	// eWave Variables
	// alglib::complex_1d_array htildehat, qtildehat_x, qtildehat_y;	// eWave inputs
	// alglib::complex_1d_array qtildehat_depth_x[DEPTH_NUM];			// eWave outputs
	// alglib::complex_1d_array qtildehat_depth_y[DEPTH_NUM];			// eWave outputs

	// SWE Variables

	// Transport Variables
	std::vector<float> qtildePast_x;// flow rate change in x direction due to advection
	std::vector<float> qtildePast_y;// flow rate change in y direction due to advection
	std::vector<float> qAdvect_x; 	// flow rate change in x direction due to advection
	std::vector<float> qAdvect_y; 	// flow rate change in y direction due to advection
	std::vector<float> hPast;		// bulk height from last timestep




};
