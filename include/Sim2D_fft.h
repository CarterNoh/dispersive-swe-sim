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
    static constexpr float CELLSIZE = 1.f;	// physical size of each cell (meters)
	static constexpr float TIMESTEP = 1.f / 30.f;
    static constexpr float SURFACE_TENSION = 0.001f; // minimum water height for stability
    static constexpr float DENSITY = 999.f; // minimum water height for stability

    static constexpr float DEPTH        = 100.f;     // meters
    static constexpr float FETCH        = 200.f;     // kilometers
    static constexpr float WIND_SPEED   = 15.f;     // m/s
    static constexpr float WIND_ANGLE   = 0.f;     // degrees from x-axis
    static constexpr float SWELL        = 0.5f;     // [0, 1]
    static constexpr float SWELL_ANGLE  = 0.f;     // degrees from x-axis
    static constexpr float CHOPPINESS   = 0.f;     // 
    static constexpr float FILTER_SMALL = 0.f;    // 
    static constexpr float FILTER_BIG   = 10000.f;  // 
    static constexpr float FILTER_WIDTH = 1.f;      // 
    static constexpr float FILTER_MIN   = 0.01f;    // 
	
	// eWave Parameters
	std::vector<float> depths = { 1.0f, 2.0f, 4.0f, 16.0f, 64.0f };
	int DEPTH_NUM = depths.size(); // number of discrete water depth solutions to compute for eWave dispersion correction

	// Simulation variables
	GPUField omega, HPos, HNeg, HProp, DxProp, DyProp, H, Dx, Dy;
	GPUField* fields[4] = {&omega, &H, &Dx, &Dy};
	GPUField* fields_complex[5] = {&HPos, &HNeg, &HProp, &DxProp, &DyProp};

	GPU* gpu;

	// Functions
	Sim();	// default constructor
	Sim(HWND hwnd);	// constructor with handle for rendering
	int Release(void);
	void SimStep();	// ticks the simulation by one timestep using the following substeps:

    float time = 0.f; // simulation time in seconds

private:
	void Init(GPU* gpu);

	// Helper functions
	inline int idx(int x, int y) const {return y * GRIDSIZE + x;}

	// GPU Compute Shaders
	ID3D11ComputeShader *PopulateSpectrum, *PropagateWaves, *ComplexToReal;
	ID3D11ComputeShader** shaders[3] = {&PopulateSpectrum, &PropagateWaves, &ComplexToReal};
	char* names[3] = {"PopulateSpectrum", "PropagateWaves", "ComplexToReal"};
	SimConstants constants = {time, GRIDSIZE, CELLSIZE, TIMESTEP, DEPTH_NUM, SURFACE_TENSION, DENSITY, 0.f,
    						  DEPTH, FETCH, WIND_SPEED, WIND_ANGLE, SWELL, SWELL_ANGLE, CHOPPINESS, 
							  FILTER_SMALL, FILTER_BIG, FILTER_WIDTH, FILTER_MIN, 0.f};			
	RenderConstants render_constants = {DirectX::XMMatrixIdentity(), (float)GRIDSIZE, CELLSIZE, {0.f, 0.f}};

};
