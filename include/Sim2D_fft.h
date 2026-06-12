// #pragma once
// #include <iostream>
// #include <vector>
// #include <algorithm>
// #include "gpu.h"

// // constants
// #define GRAVITY 9.80665
// #define PI 3.14159265359  // used in FFT stencil code

// class Sim
// {
// public:
// 	// Simulation parameters
// 	static constexpr int GRIDSIZE = 256;	// grid size in one dimension (# cells)
//     static constexpr float CELLSIZE = 1.f;	// physical size of each cell (meters)
// 	static constexpr float TIMESTEP = 1.f / 30.f;
//     static constexpr float SURFACE_TENSION = 0.001f;
//     static constexpr float DENSITY = 999.f;

// 	// FFT Parameters
//     static constexpr float DEPTH        = 2.5f;     // meters
//     static constexpr float FETCH        = 200.f;     // kilometers
//     static constexpr float WIND_SPEED   = 14.f;     // m/s
//     static constexpr float WIND_ANGLE   = 180.f;     // degrees from x-axis
//     static constexpr float SWELL        = 0.9f;    // [0, 1] // this has issues at 0 but it shouldn't, I can't find the bug, hunt down later
//     static constexpr float SWELL_ANGLE  = 180.f;    // degrees from x-axis
//     static constexpr float CHOPPINESS   = 1.f;      // 
//     static constexpr float FILTER_SMALL = 0.f;      // Set to really wide, not really using this right now 
//     static constexpr float FILTER_BIG   = 10000.f;  // 
//     static constexpr float FILTER_WIDTH = 1.f;      // 
//     static constexpr float FILTER_MIN   = 0.01f;    // 

// 	// Terrain & Water Parameters
// 	static constexpr int TERRAIN_TYPE = 4; 		   // 0 = flat, 1 = ramp, 2 = bumps, 3 = basins, 4 = beach, 5 = 1D hill, 6 = 2D hill
// 	static constexpr float TERRAIN_HEIGHT = -15.f; // base height of terrain features (meters)
// 	static constexpr float TERRAIN_SCALE = 20.f;    // scale of terrain features (meters)
// 	static constexpr int WATER_TYPE   = 2; 		   // 0 = flat, 1 = step/dam break, 2 = diagonal slope, 3 = splash, 4 = ripples, 5 = basin flood
// 	static constexpr float WATER_LEVEL = 0.f; 	   // level of water free surface at start (H)
// 	static constexpr float WATER_SCALE = 0.1f;     // scale of water height features
	
// 	// eWave Parameters
// 	std::vector<float> depths = { 1.0f, 2.0f, 4.0f, 16.0f, 64.0f };
// 	int DEPTH_NUM = depths.size(); // number of discrete water depth solutions to compute for eWave dispersion correction

// 	// Simulation variables
// 	GPUField terrain, HPos, HNeg, HProp, DxProp, DyProp, UxProp, UyProp, DelSx, DelSy, 
// 										H, Dx, Dy, Ux, Uy, Sx, Sy;
// 	GPUField* fields_r[8] = {&terrain, &H, &Dx, &Dy, &Ux, &Uy, &Sx, &Sy};
// 	GPUField* fields_c[5] = {};
// 	GPUField* fields_arrays_r[1] = {};
// 	GPUField* fields_arrays_c[9] = {&HPos, &HNeg, &HProp, &DxProp, &DyProp, &UxProp, &UyProp, &DelSx, &DelSy};
// 	GPUBuffer depth;

// 	GPU* gpu;

// 	// Functions
// 	Sim();	// default constructor
// 	Sim(HWND hwnd);	// constructor with handle for rendering
// 	int Release(void);
// 	void SimStep();	// ticks the simulation by one timestep using the following substeps:

//     float time = 0.f; // simulation time in seconds

// private:
// 	void Init(GPU* gpu);

// 	std::vector<float> SetTerrain();
// 	std::vector<float> SetWater(std::vector<float>& terrain);

// 	// Helper functions
// 	inline int idx(int x, int y) const {return y * GRIDSIZE + x;}

// 	// GPU Compute Shaders
// 	ID3D11ComputeShader *PopulateSpectrum, *PropagateWaves, *Interp;
// 	ID3D11ComputeShader** shaders[3] = {&PopulateSpectrum, &PropagateWaves, &Interp};
// 	char* names[3] = {"PopulateSpectrum", "PropagateWaves", "Interp"};
// 	SimConstants constants = {time, GRIDSIZE, CELLSIZE, TIMESTEP, DEPTH_NUM, SURFACE_TENSION, DENSITY, 0.f,
//     						  DEPTH, FETCH, WIND_SPEED, WIND_ANGLE, SWELL, SWELL_ANGLE, CHOPPINESS, 
// 							  FILTER_SMALL, FILTER_BIG, FILTER_WIDTH, FILTER_MIN, 0.f};			
// 	RenderConstants render_constants = {DirectX::XMMatrixIdentity(), (float)GRIDSIZE, CELLSIZE, {0.f, 0.f}};

// };
