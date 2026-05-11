#pragma once
#include <iostream>
#include <vector>
#include <algorithm>
#include "gpu.h"

// constants
#define GRAVITY 9.80665
#define PI 3.14159265359  // used in FFT stencil code

class Sim
{
public:
	// Simulation parameters
	static constexpr int GRIDSIZE = 512;	// grid size in one dimension (# cells)
	static constexpr float CELLSIZE = 1.f;	// cell size in one dimension (meters/cell)
	static constexpr float TIMESTEP = 1.f / 60.f;
	static constexpr float TERRAIN_HEIGHT = -10.f; // base height of terrain features (meters)
	static constexpr float TERRAIN_SCALE = 15.f; // scale of terrain features (meters)
	static constexpr int TERRAIN_TYPE = 0; // 0 = flat, 1 = ramp, 2 = bumps, 3 = basins, 4 = beach
	static constexpr int WATER_TYPE = 0; // 0 = localized splash, 1 = step/dam break, 2 = basin flood
	static constexpr float WATER_LEVEL = 10.f; 
	static constexpr float WATER_SCALE = 10.f; 
	static constexpr int BOUNDARY_TYPE = 0; // 0 = wall, 1 = free
	
	static constexpr float MIN_WATER_HEIGHT = 0.01f;  // minimum water height for stability
	
	// Decomposition Parameters
	static constexpr int DIFFUSION_ITERATIONS = 128; // number of iterations for diffusion step, more iterations means more stable but also more expensive
	static constexpr float DELTA_T = 0.25f;
	static constexpr float DIFFUSION_PENALTY = 0.01f; // penalty factor for diffusion, higher means more diffusion and more stability but also more damping of waves
	
	// eWave Parameters
	std::vector<float> depths = { 1.0f, 2.0f, 4.0f, 16.0f, 64.0f };
	int DEPTH_NUM = depths.size(); // number of discrete water depth solutions to compute for eWave dispersion correction

	// Transport Parameters
	static constexpr float CFL_CONDITION = 0.25f;  // max allowed CFL condition for stability of SWE step, can be higher than overall CFL condition since diffusion and transport steps handle stability as well
	static constexpr float GAMMA_TRANSPORT = 0.25f; // blending factor for transport step, higher means more stable but also more damping of waves

	// Simulation variables
	GPUField terrain, H, Q_x, Q_y, h, q_x, q_y, 
			 HPast, QPast_x, QPast_y, alpha_H, alpha_Q_x, alpha_Q_y, 
			 hbar, qbar_x, qbar_y, htilde, qtilde_x, qtilde_y,
			 ubar_x, ubar_y, ubarNew_x, ubarNew_y,
			 qtildePast_x, qtildePast_y, qAdvect_x, qAdvect_y, 
			 hPast, hbarOld, htildeOld, 
			 hHat, qHat_x, qHat_y, wavenum, dhHat_dx, dhHat_dy,
			 qHat_x_array, qHat_y_array, qtilde_x_array, qtilde_y_array;
	GPUField* fields[30] = {
		&terrain, &H, &Q_x, &Q_y, &h, &q_x, &q_y,
		&HPast, &QPast_x, &QPast_y, &alpha_H, &alpha_Q_x, &alpha_Q_y,
		&hbar, &qbar_x, &qbar_y, &htilde, &qtilde_x, &qtilde_y,
		&ubar_x, &ubar_y, &ubarNew_x, &ubarNew_y,
		&qtildePast_x, &qtildePast_y, &qAdvect_x, &qAdvect_y,
		&hPast, &hbarOld, &htildeOld};
	GPUField* fields_complex[6] = {&hHat, &qHat_x, &qHat_y, &dhHat_dx, &dhHat_dy, &wavenum};
	GPUField* q_arrays[2] = {&qHat_x_array, &qHat_y_array};
	GPUBuffer depth;

	GPU* gpu;

	// Functions
	Sim();	// default constructor
	Sim(HWND hwnd);	// constructor with handle for rendering
	int Release(void);
	void SimStep();	// ticks the simulation by one timestep using the following substeps:
	std::vector<float> SetTerrain();
	std::vector<float> SetWater(std::vector<float>& terrain);

private:
	void Init(GPU* gpu);
	// Functions
	void DecompositionStep();	// bulk vs surface decomposition
	void eWaveStep();		// surface wave simulation step
	// void FFTStep();
	void SWEStep();				// SWE bulk simulation step
	void TransportStep();		// transport of bulk and surface quantities
	void ComputeValues();		// compute final h and q values

	// Helper functions
	inline int idx(int x, int y) const {return y * GRIDSIZE + x;}

	// GPU Compute Shaders
	ID3D11ComputeShader *ApplyBoundaries, *InitDecomp, *CalcDiffusionCoeffs, *DiffusionStep, *DecomposeFields, 
						*CalcUbar, *CalcSWE, *UpdateTilde, *CalcQAdvect, *IntegrateH, 
						*TransferToFFT, *CalcEWave, *InterpQ;
	ID3D11ComputeShader** shaders[13] = {
						&ApplyBoundaries, &InitDecomp, &CalcDiffusionCoeffs, &DiffusionStep, &DecomposeFields, 
						&CalcUbar, &CalcSWE, &UpdateTilde, &CalcQAdvect, &IntegrateH, 
						&TransferToFFT, &CalcEWave, &InterpQ};
	char* names[13] = {"ApplyBoundaries", "InitDecomp", "CalcDiffusionCoeffs", "DiffusionStep", "DecomposeFields", 
					   "CalcUbar", "CalcSWE", "UpdateTilde", "CalcQAdvect", "IntegrateH", 
					   "TransferToFFT", "CalcEWave", "InterpQ"};
	SimConstants constants = {GRIDSIZE, CELLSIZE, TIMESTEP, BOUNDARY_TYPE, MIN_WATER_HEIGHT,// Sim Params
							  DIFFUSION_ITERATIONS, DELTA_T, DIFFUSION_PENALTY, 			// Diffusion Params
							  CFL_CONDITION, GAMMA_TRANSPORT, DEPTH_NUM, 0.0f}; 				// SWE & eWave Params
	RenderConstants render_constants = {DirectX::XMMatrixIdentity(), (float)GRIDSIZE, CELLSIZE, {0.f, 0.f}};

};
