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
	static constexpr float TIMESTEP = 1.f / 30.f;
	static constexpr float TERRAIN_HEIGHT = -10.f; // base height of terrain features (meters)
	static constexpr float TERRAIN_SCALE = 15.f; // scale of terrain features (meters)
	static constexpr int TERRAIN_TYPE = 0; // 0 = flat, 1 = ramp, 2 = bumps, 3 = basins, 4 = beach
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

	// Simulation variables
	GPUField terrain, H, Q_x, Q_y, h, q_x, q_y, 
			 HPast, QPast_x, QPast_y, alpha_H, alpha_Q_x, alpha_Q_y, 
			 hbar, qbar_x, qbar_y, htilde, qtilde_x, qtilde_y,
			 ubar_x, ubar_y, ubarNew_x, ubarNew_y,
			 qtildePast_x, qtildePast_y, qAdvect_x, qAdvect_y, 
			 hPast, hbarOld, htildeOld;
	GPUField* fields[30] = {
		&terrain, &H, &Q_x, &Q_y, &h, &q_x, &q_y,
		&HPast, &QPast_x, &QPast_y, &alpha_H, &alpha_Q_x, &alpha_Q_y,
		&hbar, &qbar_x, &qbar_y, &htilde, &qtilde_x, &qtilde_y,
		&ubar_x, &ubar_y, &ubarNew_x, &ubarNew_y,
		&qtildePast_x, &qtildePast_y, &qAdvect_x, &qAdvect_y,
		&hPast, &hbarOld, &htildeOld
	};

	// GPU Stuff
	GPU* gpu;
    ID3D11ComputeShader *shader_Boundaries, *shader_InitDecomp, *shader_CalcDiff, *shader_Diffusion, 
						*shader_Decompose, *shader_Ubar, *shader_SWE, *shader_UpdateTilde, 
						*shader_QAdvect, *shader_HAdvect, *shader_IntegrateH;
	ID3D11ComputeShader** shaders[11] = {
		&shader_Boundaries, &shader_InitDecomp, &shader_CalcDiff, &shader_Diffusion, &shader_Decompose, 
		&shader_Ubar, &shader_SWE, &shader_UpdateTilde, &shader_QAdvect, &shader_HAdvect, 
		// &shader_CombineQ, 
		&shader_IntegrateH
	};
	char* names[11] = {"ApplyBoundaries", "InitDecomp", "CalcDiffusionCoeffs", "DiffusionStep", "DecomposeFields", 
					   "CalcUbar", "CalcSWE", "UpdateTilde", "CalcQAdvect", "CalcHAdvect", "IntegrateH"};
	SimConstants constants = {GRIDSIZE, CELLSIZE, TIMESTEP, BOUNDARY_TYPE, MIN_WATER_HEIGHT,
							  DIFFUSION_ITERATIONS, DELTA_T, DIFFUSION_PENALTY, CFL_CONDITION, GAMMA_TRANSPORT, {0.f, 0.f}};
	int group = GRIDSIZE / 16; // number of thread groups to dispatch for compute shaders, assuming 16x16 threads per group
	RenderConstants render_constants = {DirectX::XMMatrixIdentity(), (float)GRIDSIZE, CELLSIZE, {0.f, 0.f}};

	// Functions
	Sim();	// default constructor
	Sim(HWND hwnd);	// constructor with handle for rendering
	int Release(void);
	void SimStep();	// ticks the simulation by one timestep using the following substeps:
	std::vector<float> SetTerrain(int type);
	std::vector<float> SetWater(int type, float level, std::vector<float>& terrain);

private:
	// Functions
	void DecompositionStep();	// bulk vs surface decomposition
	// void eWaveStep();		// surface wave simulation step
	// void FFTStep();
	void SWEStep();				// SWE bulk simulation step
	void TransportStep();		// transport of bulk and surface quantities
	void ComputeValues();		// compute final h and q values

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
