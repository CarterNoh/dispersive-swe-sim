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
	static constexpr int GRIDSIZE = 256;	// grid size in one dimension (# cells)
	static constexpr float CELLSIZE = 0.2f;	// cell size in one dimension (meters/cell)
	static constexpr float TIMESTEP = 1.f / 30.f;
	static constexpr int BOUNDARY_TYPE = 0; // not in use any more, need to remove probably
	static constexpr float MIN_WATER_HEIGHT = 0.001f; // minimum water height for stability
    static constexpr float SURFACE_TENSION = 0.001f;
    static constexpr float DENSITY = 999.f;

	// Terrain & Water Parameters
	static constexpr int TERRAIN_TYPE = 4; 		   // 0 = flat, 1 = ramp, 2 = bumps, 3 = basins, 4 = beach, 5 = 1D hill, 6 = 2D hill
	static constexpr int WATER_TYPE   = 2; 		   // 0 = flat, 1 = step/dam break, 2 = diagonal slope, 3 = splash, 4 = ripples, 5 = basin flood
	static constexpr float TERRAIN_HEIGHT = -0.6f; // base height of terrain features (meters)
	static constexpr float TERRAIN_SCALE = 1.f;    // scale of terrain features (meters)
	static constexpr float WATER_LEVEL = 0.f; 	   // level of water free surface at start (H)
	static constexpr float WATER_SCALE = 0.1f;     // scale of water height features
	
	// Decomposition Parameters
	static constexpr int DIFFUSION_ITERATIONS = 128;  // number of iterations for diffusion step, more iterations means more stable but also more expensive
	static constexpr float DELTA_T = 0.25f; 		  // difffusion timestep, "seconds" per subtick of the diffusion loop. Smaller is less diffusion, larger is more
	static constexpr float DIFFUSION_PENALTY = 0.01f; // penalty factor for diffusion, higher means more diffusion and more stability but also more damping of waves
	
	// eWave Parameters
	std::vector<float> depths = { 1.0f, 2.0f, 4.0f, 16.0f, 64.0f };
	int DEPTH_NUM = depths.size(); // number of discrete water depth solutions to compute for eWave dispersion correction

	// Transport Parameters
	static constexpr float CFL_CONDITION = 0.25f;  // max allowed CFL condition for stability of SWE step, can be higher than overall CFL condition since diffusion and transport steps handle stability as well
	static constexpr float GAMMA_TRANSPORT = 0.25f; // blending factor for transport step, higher means more stable but also more damping of waves

    // FFT Parameters
    static float time;
    static constexpr float DEPTH        = 2.5f;     // meters
    static constexpr float FETCH        = 200.f;     // kilometers
    static constexpr float WIND_SPEED   = 14.f;     // m/s
    static constexpr float WIND_ANGLE   = 180.f;     // degrees from x-axis
    static constexpr float SWELL        = 0.9f;    // [0, 1] // this has issues at 0 but it shouldn't, I can't find the bug, hunt down later
    static constexpr float SWELL_ANGLE  = 180.f;    // degrees from x-axis
    static constexpr float CHOPPINESS   = 1.f;      // 
    static constexpr float FILTER_SMALL = 0.f;      // Set to really wide, not really using this right now 
    static constexpr float FILTER_BIG   = 10000.f;  // 
    static constexpr float FILTER_WIDTH = 1.f;      // 
    static constexpr float FILTER_MIN   = 0.01f;    // 

	// Simulation variables
	GPUField terrain, H, Q_x, Q_y, h, q_x, q_y, 
			 HPast, QPast_x, QPast_y, alpha_H, alpha_Q_x, alpha_Q_y, 
			 hbar, qbar_x, qbar_y, htilde, qtilde_x, qtilde_y,
			 ubar_x, ubar_y, ubarNew_x, ubarNew_y,
			 qtildePast_x, qtildePast_y, qAdvect_x, qAdvect_y, 
			 hPast, hbarOld, htildeOld, 
			 hHat, qHat_x, qHat_y, qHat_x_array, qHat_y_array,
             HPos, HNeg, HHigh, UHigh_x, UHigh_y, DelS_x, DelS_y, // Outputs of PopulateSpectrum, PropagateWaves (complex arrays)
             htildeFFT, utildeFFT_x, utildeFFT_y, delS_x, delS_y, // iFFT'd variables after interpolation
             // two options: can interpolate then feed back in, or not interpolate before feeding in. For now just gonna not interpolate, see how it works. 
             // if interpolate: output of Propagate is complex array, iFFT in place, interpolate to real, mult real to make complex, FFT back
             // if NOT interp : output of Propagate is complex array, iFFT in place, mult (use real parts) to make complex array, FFT back
             qHatFFT_x, qHatFFT_y, qHatFFT_x_array, qHatFFT_y_array, disp_x, disp_y; // clean this up later to remove the one I don't use
	GPUField* fields[37] = {
		&terrain, &H, &Q_x, &Q_y, &h, &q_x, &q_y,
		&HPast, &QPast_x, &QPast_y, &alpha_H, &alpha_Q_x, &alpha_Q_y,
		&hbar, &qbar_x, &qbar_y, &htilde, &qtilde_x, &qtilde_y,
		&ubar_x, &ubar_y, &ubarNew_x, &ubarNew_y,
		&qtildePast_x, &qtildePast_y, &qAdvect_x, &qAdvect_y,
		&hPast, &hbarOld, &htildeOld, 
        &htildeFFT, &utildeFFT_x, &utildeFFT_y, &delS_x, &delS_y, &disp_x, &disp_y
        };
	GPUField* fields_complex[5] = {&hHat, &qHat_x, &qHat_y, &qHatFFT_x, &qHatFFT_y};
	GPUField* fields_arrays[13] = {&qHat_x_array, &qHat_y_array, 
        &HPos, &HNeg, &HHigh, &UHigh_x, &UHigh_y, &DelS_x, &DelS_y, &qHatFFT_x_array, &qHatFFT_y_array};
	GPUBuffer depth;

	GPU* gpu;

	// Functions
	Sim();	// default constructor
	Sim(HWND hwnd);	// constructor with handle for rendering
	int Release(void);
	void SimStep();	// ticks the simulation by one timestep using the following substeps:

private:
	void Init(GPU* gpu);
    std::vector<float> SetTerrain();
	std::vector<float> SetWater(std::vector<float>& terrain);
	// Functions
	void DecompositionStep();	// bulk vs surface decomposition
    void FFTStep();
	void eWaveStep();		// surface wave simulation step
	void SWEStep();				// SWE bulk simulation step
	void TransportStep();		// transport of bulk and surface quantities
	void ComputeValues();		// compute final h and q values

    

	// Helper functions
	inline int idx(int x, int y) const {return y * GRIDSIZE + x;}

	// SWE/Airy Compute Shaders
	ID3D11ComputeShader *ApplyBoundaries, *InitDecomp, *CalcDiffusionCoeffs, *DiffusionStep, *DecomposeFields, 
						*CalcQFFT, *CalcUbar, *CalcSWE, *UpdateTilde, *CalcQAdvect, *IntegrateH, 
						*TransferToFFT, *CalcEWave, *InterpQ;
	ID3D11ComputeShader** shaders[14] = {
						&ApplyBoundaries, &InitDecomp, &CalcDiffusionCoeffs, &DiffusionStep, &DecomposeFields, 
						&CalcQFFT, &CalcUbar, &CalcSWE, &UpdateTilde, &CalcQAdvect, &IntegrateH, 
						&TransferToFFT, &CalcEWave, &InterpQ};
	char* names[14] = {"ApplyBoundaries", "InitDecomp", "CalcDiffusionCoeffs", "DiffusionStep", "DecomposeFields", 
					   "CalcQFFT", "CalcUbar", "CalcSWE", "UpdateTilde", "CalcQAdvect", "IntegrateH", 
					   "TransferToFFT", "CalcEWave", "InterpQ"};
	SimConstants constants = {GRIDSIZE, CELLSIZE, TIMESTEP, BOUNDARY_TYPE, MIN_WATER_HEIGHT, SURFACE_TENSION, DENSITY, // Sim Params
							  DIFFUSION_ITERATIONS, DELTA_T, DIFFUSION_PENALTY, 		    // Diffusion Params
							  CFL_CONDITION, GAMMA_TRANSPORT, DEPTH_NUM, 0.f, 0.f, 0.f};    // SWE & eWave Params
    // FFT Wave Compute Shaders
	ID3D11ComputeShader *PopulateSpectrum, *PropagateWaves, *Interp;
	ID3D11ComputeShader** shaders_fft[3] = {&PopulateSpectrum, &PropagateWaves, &Interp};
	char* names_fft[3] = {"PopulateSpectrum", "PropagateWaves", "Interp"};
    FFTWaveConstants fftConstants = {0.f, FETCH, WIND_SPEED, WIND_ANGLE, SWELL, SWELL_ANGLE, CHOPPINESS, 
                                  FILTER_SMALL, FILTER_BIG, FILTER_WIDTH, FILTER_MIN, 0.f,};	
    // Render Shaders
	RenderConstants render_constants = {DirectX::XMMatrixIdentity(), (float)GRIDSIZE, CELLSIZE, {0.f, 0.f}};
};
