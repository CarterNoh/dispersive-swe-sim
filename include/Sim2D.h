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
	// Variables carried from one timestep to the next
	std::vector<float> terrain;	// terrain
	std::vector<float> h;		// overall water height
	std::vector<float> q_x;		// overall flow rate
	std::vector<float> q_y;		// overall flow rate
	int time;					// current simulation time in seconds

	// Simulation parameters
	static constexpr int GRIDSIZE = 512;	// grid size in one dimension (# cells)
	static constexpr float CELLSIZE = 1.f;	// cell size in one dimension (meters/cell)
	static constexpr float TIMESTEP = 1.f / 30.f;
	static constexpr float TERRAIN_HEIGHT = -10.f; // base height of terrain features (meters)
	static constexpr float TERRAIN_SCALE = 15.f; // scale of terrain features (meters)
	static constexpr int TERRAIN_TYPE = 1; // 0 = flat, 1 = ramp, 2 = bumps, 3 = basins, 4 = beach
	static constexpr int WATER_TYPE = 1; // 0 = localized splash, 1 = step/dam break, 2 = basin flood
	static constexpr int BOUNDARY_TYPE = 0; // 0 = wall, 1 = free, 2 = zero
	static constexpr float WATER_LEVEL = 10.f; 
	static constexpr float MIN_WATER_HEIGHT = 0.01f;  // minimum water height for stability
	
	// Decomposition Parameters
	static constexpr int DIFFUSION_ITERATIONS = 128; // number of iterations for diffusion step, more iterations means more stable but also more expensive
	static constexpr float DELTA_T = 0.25f;
	static constexpr float DIFFUSION_PENALTY = 0.01f; // penalty factor for diffusion, higher means more diffusion and more stability but also more damping of waves
	
	// eWave, SWE, & Transport Parameters
	static constexpr int DEPTH_NUM = 4;		// number of discrete water depth solutions to compute for eWave dispersion correction
	static constexpr float Depth[4] = { 1.f, 4.f, 16.f, 64.f };

	// Transport Parameters
	static constexpr float CFL_CONDITION = 0.5f;  // max allowed CFL condition for stability of SWE step, can be higher than overall CFL condition since diffusion and transport steps handle stability as well
	static constexpr float GAMMA_TRANSPORT = 0.25f; // blending factor for transport step, higher means more stable but also more damping of waves

	
	// variables carried from one timestep to the next
	std::vector<float> hbarOld;		// last timestep hbar, used for resampling in time
	std::vector<float> htildeOld;	// last timestep htilde, used for resampling in time

	// variables that could be allocated locally but we make global so all functions can modify
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

	std::vector<float>* fields[30] = {
		&terrain, &H, &Q_x, &Q_y, &h, &q_x, &q_y, 
		&HPast, &QPast_x, &QPast_y, &alpha_H, &alpha_Q_x, &alpha_Q_y,
		&hbar, &qbar_x, &qbar_y, &htilde, &qtilde_x, &qtilde_y, 
		&ubar_x, &ubar_y, &ubarNew_x, &ubarNew_y,
		&qtildePast_x, &qtildePast_y, &qAdvect_x, &qAdvect_y,
		&hPast, &hbarOld, &htildeOld
	};


	// GPU Stuff
	GPU* gpu;
	GPUField gpu_terrain, gpu_H, gpu_Q_x, gpu_Q_y, gpu_h, gpu_q_x, gpu_q_y, 
			 gpu_HPast, gpu_QPast_x, gpu_QPast_y, gpu_alpha_H, gpu_alpha_Q_x, gpu_alpha_Q_y, 
			 gpu_hbar, gpu_qbar_x, gpu_qbar_y, gpu_htilde, gpu_qtilde_x, gpu_qtilde_y,
			 gpu_ubar_x, gpu_ubar_y, gpu_ubarNew_x, gpu_ubarNew_y,
			 gpu_qtildePast_x, gpu_qtildePast_y, gpu_qAdvect_x, gpu_qAdvect_y, 
			 gpu_hPast, gpu_hbarOld, gpu_htildeOld;
	GPUField* gpu_fields[30] = {
		&gpu_terrain, &gpu_H, &gpu_Q_x, &gpu_Q_y, &gpu_h, &gpu_q_x, &gpu_q_y,
		&gpu_HPast, &gpu_QPast_x, &gpu_QPast_y, &gpu_alpha_H, &gpu_alpha_Q_x, &gpu_alpha_Q_y,
		&gpu_hbar, &gpu_qbar_x, &gpu_qbar_y, &gpu_htilde, &gpu_qtilde_x, &gpu_qtilde_y,
		&gpu_ubar_x, &gpu_ubar_y, &gpu_ubarNew_x, &gpu_ubarNew_y,
		&gpu_qtildePast_x, &gpu_qtildePast_y, &gpu_qAdvect_x, &gpu_qAdvect_y,
		&gpu_hPast, &gpu_hbarOld, &gpu_htildeOld
	};
    ID3D11ComputeShader *shader_Boundaries, *shader_InitDecomp, *shader_CalcDiff, *shader_Diffusion, 
						*shader_Decompose, *shader_Ubar, *shader_SWE, *shader_UpdateTilde, 
						*shader_QAdvect, *shader_HAdvect, 
						// *shader_CombineQ, 
						*shader_IntegrateH;
	ID3D11ComputeShader** compute_shaders[11] = {
		&shader_Boundaries, &shader_InitDecomp, &shader_CalcDiff, &shader_Diffusion, &shader_Decompose, 
		&shader_Ubar, &shader_SWE, &shader_UpdateTilde, &shader_QAdvect, &shader_HAdvect, 
		// &shader_CombineQ, 
		&shader_IntegrateH
	};
	char* compute_names[11] = {"ApplyBoundaries", "InitDecomp", "CalcDiffusionCoeffs", "DiffusionStep", "DecomposeFields", 
					   "CalcUbar", "CalcSWE", "UpdateTilde", "CalcQAdvect", "CalcHAdvect", //"CombineQ", 
					   "IntegrateH"};
	SimConstants compute_constants = {GRIDSIZE, CELLSIZE, TIMESTEP, BOUNDARY_TYPE, MIN_WATER_HEIGHT,
							  DIFFUSION_ITERATIONS, DELTA_T, DIFFUSION_PENALTY, CFL_CONDITION, GAMMA_TRANSPORT, {0.f, 0.f}};
	int group = GRIDSIZE / 16; // number of thread groups to dispatch for compute shaders, assuming 16x16 threads per group
	RenderConstants render_constants = {DirectX::XMMatrixIdentity(), (float)GRIDSIZE, CELLSIZE, {0.f, 0.f}};

	// Functions
	Sim();	// default constructor
	Sim(HWND hwnd);	// constructor with handle for rendering
	int Release(void);
	void SimStep(bool SWEonly);	// ticks the simulation by one timestep using the following substeps:
	void SetTerrain(int type);
	void SetWater(int type, float level);

private:
	// Functions
	void DecompositionStep(bool SWEonly); 	// bulk vs surface decomposition
	// void eWaveStep(bool SWEonly);			// surface wave simulation step
	void SWEStep();							// SWE bulk simulation step
	void TransportStep();					// transport of bulk and surface quantities
	void ComputeValues();					// compute final h and q values

	// Helper functions
	inline int idx(int x, int y) const {return y * GRIDSIZE + x;}
	inline float LimitFlowRate(float flow_rate_in, float waterDepth_left, float waterDepth_right);
	inline float LimitVelocity(float velocity_in);
	bool StopFlowOnTerrainBoundary(int x, int y, std::vector<float>& h, std::vector<float>& terrain, bool isYDirection);
	void ExtractLocal(std::vector<float>& target, std::vector<float>& dataField, int index, bool isYDirection);
	float SampleCubicClamped(float samplePos, std::vector<float>& dataField);
	void HandleWallBoundary(std::vector<float>& field);
	void HandleFreeBoundary(std::vector<float>& field);
	void HandleZeroBoundary(std::vector<float>& field);
	void ApplyBoundaries(std::vector<float>& field, int type);




};
