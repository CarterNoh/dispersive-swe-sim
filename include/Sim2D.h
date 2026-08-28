#pragma once
#include <iostream>
#include <vector>
#include <algorithm>
#include "gpu.h"

// constants
#define GRAVITY 9.80665
#define PI 3.14159265359  // used in FFT stencil code

class Sim {
public:
	// Simulation parameters
	static constexpr int GRIDSIZE_X = 512; // grid size in X dimension (# cells)
	static constexpr int GRIDSIZE_Y = 512; // grid size in Y dimension (# cells)
	static constexpr float DOMAIN_SIZE_X = 252.f; // domain size in meters
	static constexpr float CELLSIZE = DOMAIN_SIZE_X / GRIDSIZE_X;	// cell size in one dimension (meters/cell)
	static constexpr float TIMESTEP = 1.f / 60.f;
	static constexpr float CFL_CONDITION = 0.25f;  // max allowed CFL condition for stability
	static constexpr float MIN_WATER_HEIGHT = 0.001f; // minimum water height for stability
    static constexpr float SURFACE_TENSION = 0.001f;
    static constexpr float DENSITY = 999.f;

	// Terrain & Water Parameters
	static constexpr int TERRAIN_TYPE = 4; 		   // 0 = flat, 1 = ramp, 2 = bumps, 3 = basins, 4 = beach, 5 = 1D hill, 6 = 2D hill
	static constexpr int WATER_TYPE   = 1; 		   // 0 = flat, 1 = step/dam break, 2 = diagonal slope, 3 = splash, 4 = ripples, 5 = basin flood
	static constexpr float TERRAIN_HEIGHT = -13.f; // base height of terrain features (meters)
	static constexpr float TERRAIN_SCALE = 20.f;    // scale of terrain features (meters)
	static constexpr float WATER_LEVEL = 0.f; 	   // level of water free surface at start (H)
	static constexpr float WATER_SCALE = 0.f;     // scale of water height features
	
	// Decomposition Parameters
	static constexpr int DIFFUSION_ITERATIONS = 256;   // number of iterations for diffusion step, more iterations means more stable but also more expensive
	static constexpr int MAX_DIFFUSION_CELLS = 8; 	  // maximum height of diffusion stencil in cells, higher means more diffusion in deep water
	static constexpr float DIFFUSION_PENALTY = 0.001f; // penalty factor for diffusion, higher means more diffusion and more stability but also more damping of waves
	
	// eWave Parameters
	std::vector<float> depths = { 1.0f, 2.0f, 4.0f, 16.0f, 64.0f };
	int DEPTH_NUM = depths.size(); // number of discrete water depth solutions to compute for eWave dispersion correction

	// SWE & Transport Parameters
	static constexpr float SLOPE_LIMIT = 1.f; 	    // max slope of buld flow surface before energy dissipation occurs
	static constexpr float GAMMA_TRANSPORT = 0.25f; // blending factor for transport step, lower means more damping of waves
	static constexpr int   SPONGE_THICKNESS = 8; // thickness in cells of the sponge layer used to absorb waves at the boundaries
	static constexpr float LAPLACIAN_DAMPING = 0.01f; // damping factor for Laplacian smoothing to reduce spikes and unstable grid-scale ripples

    // FFT Wave Parameters
    float time = 0.f;
    static constexpr float FETCH        = 200.f;    // kilometers
    static constexpr float WIND_SPEED   = 14.f;     // m/s, at 10 meters above surface
    static constexpr float WIND_ANGLE   = 45.f;    // degrees from x-axis
    static constexpr float SWELL        = 0.3f;     // [0, 1] // this has issues at 0 but it shouldn't, I can't find the bug, hunt down later
    static constexpr float SWELL_ANGLE  = 45.f;    // degrees from x-axis
    static constexpr float CHOPPINESS   = 1.f;      // Amount of horizontal displacement in waves
    static constexpr float FILTER_SMALL = 0.f;      // Set to really wide, not really using this right now 
    static constexpr float FILTER_BIG   = 10000.f;  // 
    static constexpr float FILTER_WIDTH = 1.f;      // 
    static constexpr float FILTER_MIN   = 0.01f;    // 
	static constexpr float DEPTH_CUTOFF = 4.f; 		// depth to start attenuating FFT waves

	// Simulation variables
	GPUField terrain, H, Q_x, Q_y, h, q_x, q_y,
			 HOrig, QOrig_x, QOrig_y, HPast, QPast_x, QPast_y, 
			 alpha_H, alpha_Q_x, alpha_Q_y, 
			 hbar, qbar_x, qbar_y, htilde, qtilde_x, qtilde_y,
			 ubar_x, ubar_y, ubarNew_x, ubarNew_y,
			 htildePast, qtildePast_x, qtildePast_y, qAdvect_x, qAdvect_y, 
			 hPast, hbarOld, htildeOld, htildeOldNext,
			 hHat, qHat_x, qHat_y, qHat_x_array, qHat_y_array;
	
    GPUField HPos, HNeg, HProp, DelH_x, DelH_y, Disp_x, Disp_y, FlowX, FlowY, // Outputs of PopulateSpectrum, PropagateWaves (complex arrays)
             hFFT, delH_x, delH_y, disp_x, disp_y; // iFFT'd variables after interpolation
	GPUField* fields[40] = {
		&terrain, &H, &Q_x, &Q_y, &h, &q_x, &q_y,
		&HOrig, &QOrig_x, &QOrig_y, &HPast, &QPast_x, &QPast_y, 
		&alpha_H, &alpha_Q_x, &alpha_Q_y,
		&hbar, &qbar_x, &qbar_y, &htilde, &qtilde_x, &qtilde_y,
		&ubar_x, &ubar_y, &ubarNew_x, &ubarNew_y,
		&htildePast, &qtildePast_x, &qtildePast_y, &qAdvect_x, &qAdvect_y,
		&hPast, &hbarOld, &htildeOld, &htildeOldNext,
        &hFFT, &delH_x, &delH_y, &disp_x, &disp_y,
	};
	GPUField* fields_complex[3] = {&hHat, &qHat_x, &qHat_y};
	GPUField* fields_arrays[11] = {
		&qHat_x_array, &qHat_y_array, &HPos, &HNeg, &HProp, 
        &DelH_x, &DelH_y, &Disp_x, &Disp_y, &FlowX, &FlowY};
	GPUBuffer depth;
	GPU* gpu;

	// Functions
	Sim();	// default constructor
	Sim(HWND hwnd);	// constructor with handle for rendering
	int Release(void);
	void SimStep();	// ticks the simulation by one timestep
	void DecompositionStep();	// bulk vs surface decomposition
    void FFTStep();				// propagate FFT wave simulation
	void eWaveStep();			// surface wave simulation step
	void SWEStep();				// SWE bulk simulation step
	void TransportStep();		// transport of bulk and surface quantities
	void ComputeValues();		// compute final h and q values

	// Constants
	SimConstants constants = {time, GRIDSIZE_X, GRIDSIZE_Y, CELLSIZE, TIMESTEP, MIN_WATER_HEIGHT, SURFACE_TENSION, DENSITY, // Sim Params
							  DIFFUSION_ITERATIONS, MAX_DIFFUSION_CELLS, DIFFUSION_PENALTY, // Diffusion Params
							  SLOPE_LIMIT, CFL_CONDITION, GAMMA_TRANSPORT, SPONGE_THICKNESS, LAPLACIAN_DAMPING, DEPTH_NUM, // SWE & eWave Params
							  FETCH, WIND_SPEED, WIND_ANGLE, SWELL, SWELL_ANGLE, CHOPPINESS, // FFT Params
							  FILTER_SMALL, FILTER_BIG, FILTER_WIDTH, FILTER_MIN, DEPTH_CUTOFF,
							  0, 0, 0.0f, {0.0f}};    
	RenderConstants render_constants = {DirectX::XMMatrixIdentity(), (float)GRIDSIZE_X, (float)GRIDSIZE_Y, CELLSIZE};

private:
	int paddedSizeX;
	int paddedSizeY;

	// Functions
	void Init(GPU* gpu);
    std::vector<float> SetTerrain();
	std::vector<float> SetWater(std::vector<float>& terrain);	

	// Helper functions
	inline int idx(int x, int y) const {return y * GRIDSIZE_X + x;}

	// Compute Shaders
	ID3D11ComputeShader *InitDecomp, *CalcDiffusionCoeffs, *DiffusionStep, *DecomposeFields, 
						*CalcUbar, *CalcSWE, *UpdateTilde, *CalcQAdvect, *IntegrateH, 
						*TransferToFFT, *CalcEWave, *InterpQ;
	ID3D11ComputeShader** shaders[12] = {
						&InitDecomp, &CalcDiffusionCoeffs, &DiffusionStep, &DecomposeFields, 
						&CalcUbar, &CalcSWE, &UpdateTilde, &CalcQAdvect, &IntegrateH, 
						&TransferToFFT, &CalcEWave, &InterpQ};
	char* names[12] = {"InitDecomp", "CalcDiffusionCoeffs", "DiffusionStep", "DecomposeFields", 
					   "CalcUbar", "CalcSWE", "UpdateTilde", "CalcQAdvect", "IntegrateH", 
					   "TransferToFFT", "CalcEWave", "InterpQ"};
    // FFT Wave Compute Shaders
	ID3D11ComputeShader *PopulateSpectrum, *PropagateWaves, *Interp;
	ID3D11ComputeShader** waveShaders[3] = {&PopulateSpectrum, &PropagateWaves, &Interp};
	char* waveNames[3] = {"PopulateSpectrum", "PropagateWaves", "Interp"};	
};
