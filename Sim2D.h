#pragma once

#include <vector>
#include "alglib\fasttransforms.h"  // https://www.alglib.net/download.php#cpp


// sim parameters
#define GRIDSIZE 258 	// grid size in one dimension (meters?)
#define CELLSIZE 1		// cell size in one dimension (meters?/cell)
#define TIMESTEP (1.f/60.f)
#define DEPTH_NUM 4
const float Depth[DEPTH_NUM] = { 1.f, 4.f, 16.f, 64.f };
#define TERRAIN_HEIGHT_SHIFT_INIT -10.f // -20.f //-10.f
#define TERRAIN_HEIGHT_SCALE_INIT 20.f //40.f  //20.f
// diffusion parameters
#define DIFFUSION_ITERATIONS 128
#define DELTA_T 0.25f
#define DIFFUSION_PENALTY 0.01f


// helpful shortcuts
#define GRAVITY 9.80665
#define PI 3.14159265359  // used in FFT stencil code
#define idx y * GRIDSIZE + x
#define idx_xplus y * GRIDSIZE + x + 1
#define idx_xminus y * GRIDSIZE + x - 1
#define idx_yplus (y + 1) * GRIDSIZE + x
#define idx_yminus (y - 1) * GRIDSIZE + x

class Sim
{
public:
	// variables carried from one timestep to the next
	std::vector<double> terrain;	// terrain
	std::vector<double> h;			// overall water height
	std::vector<double> q;			// overall flow rate
	std::vector<double> hbarOld;	// last timestep hbar, used for resampling in time
	std::vector<double> htildeOld;	// last timestep htilde, used for resampling in time

	// variables that could be allocated locally but for potential visualizations we store them globally
	std::vector<double> hbar;		// bulk height
	std::vector<double> qbar;		// bulk flow rate
	std::vector<double> htilde;		// surface height
	std::vector<double> qtilde;		// surface flow rate												
	alglib::complex_1d_array htildehat, qtildehat;	// eWave inputs
	alglib::complex_1d_array qtildehat_depth[DEPTH_NUM];			// eWave outputs

	// time is exclusively used for video recording
	float time;

	// functions
	Sim(): {};
	// Sim();
	int Sim::Release(void);

	void ResetTerrain(int type);
	void ResetWater(int type, float level);
	void EditWaterLocal(float xCoord, float yCoord, float size, float factor);	// add or subtract water locally.

	void SimStep(bool SWEonly);				// advects the simulation by one timestep
	void DecompositionStep(bool SWEonly); 	// bulk vs surface decomposition
	void eWaveStep(bool SWEonly);			// surface wave simulation step
	void SWEStep();							// SWE bulk simulation step
	void TransportStep();					// transport of bulk and surface quantities
	void ComputeValues();					// compute final h and q values
	
															
};
