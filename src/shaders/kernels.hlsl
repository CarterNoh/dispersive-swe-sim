// Constant buffer matching the C++ struct
cbuffer Constants : register(b0) {
    // Sim params
    int gridSize; 
    float cellSize;
    float timeStep;
    int boundaryType;
    float minWaterHeight;

    // Decomposition params
    int diffusionIterations;
    float deltaT;
    float diffusionPenalty;

    // SWE & Transport Params
    float cflCondition;
    float gammaTransport;

    // eWave params
    int depthNum;

    // Padding for 16-byte alignment 
    float buffer;
};

#define GRAVITY 9.80665f
#define PI 3.14159265358979323846f

// Directional helpers
#define CURR  id.xy
#define RIGHT uint2(min(id.x + 1, gridSize - 1), id.y)
#define LEFT  uint2(max(id.x - 1, 0), id.y)
#define UP    uint2(id.x, min(id.y + 1, gridSize - 1))
#define DOWN  uint2(id.x, max(id.y - 1, 0))
#define RIGHT2 uint2(min(id.x + 2, gridSize - 1), id.y)
#define UP2    uint2(id.x, min(id.y + 2, gridSize - 1))

// Texture Registers
Texture2D<float>   in0 : register(t0);
Texture2D<float>   in1 : register(t1);
Texture2D<float>   in2 : register(t2);
Texture2D<float>   in3 : register(t3);
Texture2D<float>   in4 : register(t4);
Texture2D<float>   in5 : register(t5);
Texture2D<float>   in6 : register(t6);
Texture2D<float>   in7 : register(t7);
StructuredBuffer<float> in8 : register(t8);
RWTexture2D<float> out0: register(u0);
RWTexture2D<float> out1: register(u1);
RWTexture2D<float> out2: register(u2);
RWTexture2D<float> out3: register(u3);
RWTexture2D<float> out4: register(u4);
RWTexture2D<float> out5: register(u5);


//////////////////// HELPER FUNCTIONS /////////////////////////
float CalcDiffusion(Texture2D<float> f, Texture2D<float> a, uint2 curr) {
    uint2 right = uint2(min(curr.x + 1, gridSize - 1), curr.y);
    uint2 left  = uint2(max(curr.x - 1, 0), curr.y);
    uint2 up    = uint2(curr.x, min(curr.y + 1, gridSize - 1));
    uint2 down  = uint2(curr.x, max(curr.y - 1, 0));
    float dF_x = a[curr] * (f[right] - f[curr]) - a[left] * (f[curr] - f[left]);
    float dF_y = a[curr] * (f[up] - f[curr]) - a[down] * (f[curr] - f[down]);
    float dFdX = (dF_x + dF_y) / (cellSize * cellSize);
    return dFdX;
}

bool StopFlowOnTerrainBoundary(Texture2D<float> h, Texture2D<float> terrain, uint2 curr, bool isYDirection=false) {
    uint2 right = uint2(min(curr.x + 1, gridSize - 1), curr.y);
    uint2 up    = uint2(curr.x, min(curr.y + 1, gridSize - 1));
    // test if the terrain boundary stops any flow across x+0.5
	if (!isYDirection) {
		if ((h[curr] <= minWaterHeight) && (terrain[curr] >= terrain[right] + h[right])) // positive q_x
			return true;
		if ((h[right] <= minWaterHeight) && (terrain[right] > terrain[curr] + h[curr])) // negative q_x
			return true;
		return false;
	}
	else {
		if ((h[curr] <= minWaterHeight) && (terrain[curr] >= terrain[up] + h[up])) // positive q_y
			return true;
		if ((h[up] <= minWaterHeight) && (terrain[up] > terrain[curr] + h[curr])) // negative q_y
			return true;
		return false;
	}
}

float LimitVelocity(float velocity_in) {
	if (velocity_in >= 0.f)
		return min(velocity_in, cflCondition * cellSize / timeStep);
	else
		return max(velocity_in, -cflCondition * cellSize / timeStep);
}

float LimitFlowRate(float flow_rate_in, float waterDepth_left, float waterDepth_right) {
	if (flow_rate_in >= 0.f)
		return min(flow_rate_in, cflCondition * waterDepth_left * cellSize / timeStep);
	else
		return max(flow_rate_in, -cflCondition * waterDepth_right * cellSize / timeStep);
}

float SampleCubicClamped2D(Texture2D<float> dataField, float2 samplePos) {
    // 1. Calculate integer base and fractional coordinates
    int2 id_start = (int2)floor(samplePos) - 1;
    float2 frac = samplePos - floor(samplePos);
    
    // 2. Safely clamp all 4x4 grid indices to the domain boundaries
    int maxGrid = (int)gridSize - 1;
    int2 id0 = clamp(id_start + 0, 0, maxGrid);
    int2 id1 = clamp(id_start + 1, 0, maxGrid);
    int2 id2 = clamp(id_start + 2, 0, maxGrid);
    int2 id3 = clamp(id_start + 3, 0, maxGrid);

    // 3. Precompute cubic weights for X and Y
    float2 frac2 = frac * frac;
    float2 frac3 = frac2 * frac;
    float4 wX = float4(
        -frac3.x + 2.f * frac2.x - frac.x,
         3.f * frac3.x - 5.f * frac2.x + 2.f,
        -3.f * frac3.x + 4.f * frac2.x + frac.x,
         frac3.x - frac2.x
    ) * 0.5f;
    float4 wY = float4(
        -frac3.y + 2.f * frac2.y - frac.y,
         3.f * frac3.y - 5.f * frac2.y + 2.f,
        -3.f * frac3.y + 4.f * frac2.y + frac.y,
         frac3.y - frac2.y
    ) * 0.5f;

    // 4. Sample the 4x4 grid and interpolate along X for each row
    int4 xIndices = int4(id0.x, id1.x, id2.x, id3.x);
    float row0 = wX.x * dataField[int2(xIndices.x, id0.y)] + wX.y * dataField[int2(xIndices.y, id0.y)] + 
                 wX.z * dataField[int2(xIndices.z, id0.y)] + wX.w * dataField[int2(xIndices.w, id0.y)];
    float row1 = wX.x * dataField[int2(xIndices.x, id1.y)] + wX.y * dataField[int2(xIndices.y, id1.y)] + 
                 wX.z * dataField[int2(xIndices.z, id1.y)] + wX.w * dataField[int2(xIndices.w, id1.y)];
    float row2 = wX.x * dataField[int2(xIndices.x, id2.y)] + wX.y * dataField[int2(xIndices.y, id2.y)] + 
                 wX.z * dataField[int2(xIndices.z, id2.y)] + wX.w * dataField[int2(xIndices.w, id2.y)];
    float row3 = wX.x * dataField[int2(xIndices.x, id3.y)] + wX.y * dataField[int2(xIndices.y, id3.y)] + 
                 wX.z * dataField[int2(xIndices.z, id3.y)] + wX.w * dataField[int2(xIndices.w, id3.y)];

    // 5. Interpolate the 4 row results along Y
    float result = wY.x * row0 + wY.y * row1 + wY.z * row2 + wY.w * row3;

    // 6. Monotonicity clamp: limit result to the bounding box of the 4 immediate neighbors (id1, id2)
    float c00 = dataField[int2(id1.x, id1.y)];
    float c10 = dataField[int2(id2.x, id1.y)];
    float c01 = dataField[int2(id1.x, id2.y)];
    float c11 = dataField[int2(id2.x, id2.y)];
    float minVal = min(min(c00, c10), min(c01, c11));
    float maxVal = max(max(c00, c10), max(c01, c11));

    return clamp(result, minVal, maxVal);
}

float2 ComplexMul(float2 a, float2 b) {
    return float2(a.x * b.x - a.y * b.y,
                  a.x * b.y + a.y * b.x);
}

uint4 CheckSource(uint2 id, uint type) {
    if (id.x >= (uint)(gridSize) || id.y >= (uint)(gridSize)) return uint4(0, 0, 0, 0);

    // type: 0 = cell centers, 1 = x faces, 2 = y faces

    bool left   = (id.x == 0);
    bool right  = (id.x == (uint)(gridSize) - 1);
    bool bottom = (id.y == 0);
    bool top    = (id.y == (uint)(gridSize) - 1);
    // return uint4(left, right, bottom, top);
    if (!(left || right || top || bottom)) return uint4(0, 0, 0, 0); // interior cells, skip
    // if (type == 1 && !left && !right) return uint4(0, 0, 0, 0); // X faces: only update if adjacent to left or right boundary
    // if (type == 2 && !top && !bottom) return uint4(0, 0, 0, 0); // Y faces: only update if adjacent to top or bottom boundary

  
    uint2 src; // Pick a source interior cell depending on boundary side
    // uint2 dir = uint2(0, 0); // for interpolation, pick another source cell one step inward
    // if (type == 0) {
        src = uint2(left ? 1 : right ? gridSize-2 : id.x,
                    bottom ? 1 : top ? gridSize-2 : id.y);
    // }
    // if (type == 1) { // x faces
    //     // src = uint2(left ? 1 : right ? gridSize-2 : id.x, id.y);
    //     dir.x = (left ? 1 : right ? -1 : 0);
    // }
    // if (type == 2) {// y faces
    //     // src = uint2(id.x, bottom ? 1 : top ? gridSize-2 : id.y);
    //     dir.y = (bottom ? 1 : top ? -1 : 0);
    // }
    uint2 dir = uint2((left ? 1 : right ? -1 : 0),
                      (bottom ? 1 : top ? -1 : 0));
    return uint4(src, dir);
}

//////////////////// COMPUTE SHADERS /////////////////////////

// Apply Boundary Conditions
[numthreads(16, 16, 1)]
void ApplyBoundariesCenter(uint3 id : SV_DispatchThreadID) {
    uint4 src = CheckSource(id.xy, 0);
    if (src.x == 0 && src.y == 0) return; // not a boundary cell, skip
    // out0[CURR] = out0[src.xy];
    // out1[CURR] = out1[src.xy];
    out0[CURR] = 2.0f * out0[src.xy] - out0[src.xy + src.zw];
    out1[CURR] = 2.0f * out1[src.xy] - out1[src.xy + src.zw];
}

[numthreads(16, 16, 1)]
void ApplyBoundariesXFaces(uint3 id : SV_DispatchThreadID) {
    uint4 src = CheckSource(id.xy, 1);
    if (src.x == 0 && src.y == 0) return; // not a boundary cell, skip

    bool bottom = (id.y == 0);
    bool top    = (id.y == (uint)(gridSize) - 1);

    // for faces adjacent to top or bottom boundaries (including corners for stability), use simple copy since edge is perpendicular to boundary
    // if (top || bottom) {
    //     out0[CURR] = out0[src.xy];
    //     out1[CURR] = out1[src.xy];
    //     return;
    // }
    // for faces adjacent to left or right boundaries, use linear extrapolation to better transmit flow across boundary
    out0[CURR] = 2.0f * out0[src.xy] - out0[src.xy + src.zw];
    out1[CURR] = 2.0f * out1[src.xy] - out1[src.xy + src.zw];
    // out0[CURR] = out0[src.xy] - timeStep * sqrt(GRAVITY * in0[id.xy]) * (out0[src.xy] - out0[src.xy + src.zw]) / cellSize;
    // out1[CURR] = out1[src.xy] - timeStep * sqrt(GRAVITY * in0[id.xy]) * (out1[src.xy] - out1[src.xy + src.zw]) / cellSize;
}

[numthreads(16, 16, 1)]
void ApplyBoundariesYFaces(uint3 id : SV_DispatchThreadID) {
    uint4 src = CheckSource(id.xy, 2);
    if (src.x == 0 && src.y == 0) return; // not a boundary cell, skip

    bool left   = (id.x == 0);
    bool right  = (id.x == (uint)(gridSize) - 1);

    
    // // for faces adjacent to left or right boundaries, use simple copy since edge is perpendicular to boundary
    // if (left || right) {
    //     out0[CURR] = out0[src.xy];
    //     out1[CURR] = out1[src.xy];
    //     return;
    // }
    // for faces adjacent to top or bottom boundaries, use linear extrapolation to better transmit flow across boundary
    out0[CURR] = 2.0f * out0[src.xy] - out0[src.xy + src.zw];
    out1[CURR] = 2.0f * out1[src.xy] - out1[src.xy + src.zw];
    // out0[CURR] = out0[src.xy] - timeStep * sqrt(GRAVITY * in0[id.xy]) * (out0[src.xy] - out0[src.xy + src.zw]) / cellSize;
    // out1[CURR] = out1[src.xy] - timeStep * sqrt(GRAVITY * in0[id.xy]) * (out1[src.xy] - out1[src.xy + src.zw]) / cellSize;
}


/////////// Decomposition ///////////

[numthreads(16, 16, 1)]
void InitDecomp(uint3 id : SV_DispatchThreadID) {
    // Inputs: in0 = h, in1 = q_x, in2 = q_y, in3 = terrain
    // Outputs: out0 = H, out1 = Q_x, out2 = Q_y
    if (id.x >= (uint)(gridSize) || id.y >= (uint)(gridSize)) return;
    out0[CURR] = in0[CURR] + in3[CURR];
    out1[CURR] = in1[CURR];
    out2[CURR] = in2[CURR];
}

[numthreads(16, 16, 1)]
void CalcDiffusionCoeffs(uint3 id : SV_DispatchThreadID) {
    // Inputs: in0 = H, in1 = terrain
    // Outputs: out0 = alpha_H, out1 = alpha_Q_x, out2 = alpha_Q_y
    if (id.x >= (uint)(gridSize) || id.y >= (uint)(gridSize)) return;

    // float max_ground = max(in1[CURR], in1[RIGHT], in1[UP]);
    // float min_water = min(in0[CURR], in0[RIGHT], in0[UP]); // or average of these 3
    float h = in0[CURR] - in1[CURR]; // min_water - max_ground;
    float hr = in0[RIGHT] - in1[RIGHT];
    float hu = in0[UP] - in1[UP];
    hr = (h + hr) / 2;
    hu = (h + hu) / 2;
    float denom = 2.0f * deltaT * diffusionIterations;
    float grad_x = (in0[RIGHT] - in0[CURR]) / cellSize;
    float grad_y = (in0[UP] - in0[CURR]) / cellSize;
    float penalty = exp(- diffusionPenalty * (grad_x * grad_x + grad_y * grad_y));
    float max_alpha = (cellSize * cellSize) / (4.0f * deltaT); // Von Neumann Stability Condition
    out0[CURR] = min(max_alpha, h * h   / denom) * penalty;
    out1[CURR] = min(max_alpha, hr * hr / denom) * penalty;
    out2[CURR] = min(max_alpha, hu * hu / denom) * penalty;
}

[numthreads(16, 16, 1)]
void DiffusionStep(uint3 id : SV_DispatchThreadID) {
    // Inputs: in0 = terrain, in1 = HPast, in2 = QPast_x, in3 = QPast_y, 
    //         in4 = alpha_H, in5 = alpha_Q_x, in6 = alpha_Q_y
    // Outputs: out0 = H, out1 = Q_x, out3 = Q_y
    if (id.x >= (uint)(gridSize) || id.y >= (uint)(gridSize)) return;
    // if (id.x <= 0 || id.x >= (uint)(gridSize-1) || id.y <= 0 || id.y >= (uint)(gridSize-1)) return;

    float newH = in1[CURR] + deltaT * CalcDiffusion(in1, in4, CURR);
    out0[CURR] = max(in0[CURR], newH);
    out1[CURR] = in2[CURR] + deltaT * CalcDiffusion(in2, in5, CURR);
    out2[CURR] = in3[CURR] + deltaT * CalcDiffusion(in3, in6, CURR);
}

[numthreads(16, 16, 1)]
void DecomposeFields(uint3 id : SV_DispatchThreadID) {
    // Inputs: in0 = H, in1 = Q_x, in2 = Q_y
    //         in3 = h, in4 = q_x, in5 = q_y, in6 = terrain
    // Outputs: out0 = hbar, out1 = qbar_x, out2 = qbar_y, 
    //          out3 = htilde, out4 = qtilde_x, out5 = qtilde_y
    if (id.x >= (uint)(gridSize) || id.y >= (uint)(gridSize)) return;
    
    float hbar = max(0.f, in0[CURR] - in6[CURR]);
    out0[CURR] = hbar;
    out1[CURR] = in1[CURR];
    out2[CURR] = in2[CURR];
    out3[CURR] = in3[CURR] - hbar;
    out4[CURR] = in4[CURR] - in1[CURR];
    out5[CURR] = in5[CURR] - in2[CURR];

    // // Airy Only
    // out0[CURR] = max(0.f, in0[CURR] - in6[CURR]);
    // out1[CURR] = 0.f;
    // out2[CURR] = 0.f;
    // out3[CURR] = in3[CURR] - out0[CURR];
    // out4[CURR] = in4[CURR];
    // out5[CURR] = in5[CURR];

    // // SWE Only
    // out0[CURR] = in3[CURR];
    // out1[CURR] = in4[CURR];
    // out2[CURR] = in5[CURR];
    // out3[CURR] = 0;
    // out4[CURR] = 0;
    // out5[CURR] = 0;

    if (StopFlowOnTerrainBoundary(in3, in6, CURR, false)) { // stop flow in x direction
        out1[CURR] = 0.f;
        out4[CURR] = 0.f;
    }
    if (StopFlowOnTerrainBoundary(in3, in6, CURR, true)) { // stop flow in y direction
        out2[CURR] = 0.f;
        out5[CURR] = 0.f;
    }
}


/////////// SWE ///////////

[numthreads(16, 16, 1)]
void CalcUbar(uint3 id : SV_DispatchThreadID) {
    // Inputs: in0 = qbar_x, in1 = qbar_y, in2 = hbarOld
    // Outputs: out0 = ubar_x, out1 = ubar_y
    if (id.x >= (uint)(gridSize) || id.y >= (uint)(gridSize)) return;
    // if (id.x >= (uint)(gridSize-1) || id.y >= (uint)(gridSize-1)) return;

    // First-Order Up-Winding
    float hx = (in0[CURR] >= 0.f || id.x == gridSize-1) ? in2[CURR] : in2[RIGHT];
    float hy = (in1[CURR] >= 0.f || id.y == gridSize-1) ? in2[CURR] : in2[UP];
    float ubar_x = in0[CURR] / max(minWaterHeight, hx);
    float ubar_y = in1[CURR] / max(minWaterHeight, hy);

    // Enforcing CFL condition for later surface waves advection
    out0[CURR] = LimitVelocity(ubar_x);  
    out1[CURR] = LimitVelocity(ubar_y); 
}

[numthreads(16, 16, 1)]
void CalcSWE(uint3 id : SV_DispatchThreadID) {
    // Inputs: in0 = ubar_x, in1 = ubar_y, in2 = hbar, in3 = H
    // Outputs: out0 = ubarNew_x, out1 = ubarNew_y, out2 = qbar_x, out3 = qbar_y
    if (id.x >= (uint)(gridSize) || id.y >= (uint)(gridSize)) return;
    // if (id.x <= 0 || id.x >= (uint)(gridSize-1) || id.y <= 0 || id.y >= (uint)(gridSize-1)) return;

    // Compute intermediate values needed for du/dt calculations
    // Need: q_x/y_ij, q_x_i-0.5_j, q_y_i_j-0.5, 
    // 	     h_i+0.5_j, h_i_j+0.5, h_ij=hbar[x,y], h_i+1_j=hbar[x+1,y], h_i_j+1=hbar[x,y+1],
    //       u_x_i+0.5_j=ubar_x[x,y], u_x_i-0.5_j=ubar_x[x-1,y], 
    //       u_y_i_j+0.5=ubar_y[x,y], u_y_i_j-0.5=ubar_y[x,y-1], 
    //       u_x_i+0.5_j-1=ubar_x[x,y-1], u_y_i-1_j+0.5=ubar_y[x-1,y]

    // Use upwinding to evaluate q_(i-0.5,j), q_(i+0.5,j), q_(i+1.5,j) 
    // for x direction to get q_x_(i,j), q_x_(i+1,j)
    float q_x_m05 = in0[LEFT] * ((in0[LEFT] >= 0.f) ? in2[LEFT] : in2[CURR]);
    float q_x_p05 = in0[CURR] * ((in0[CURR] >= 0.f) ? in2[CURR] : in2[RIGHT]);
    float q_x_p15 = in0[RIGHT]* ((in0[RIGHT]>= 0.f) ? in2[RIGHT]: in2[RIGHT2]);
    float q_x_0 = 0.5f * (q_x_m05 + q_x_p05);
    float q_x_1 = 0.5f * (q_x_p05 + q_x_p15);
    float q_y_m05 = in1[DOWN] * ((in1[DOWN] >= 0.f) ? in2[DOWN] : in2[CURR]);
    float q_y_p05 = in1[CURR] * ((in1[CURR] >= 0.f) ? in2[CURR] : in2[UP]);
    float q_y_p15 = in1[UP]   * ((in1[UP]   >= 0.f) ? in2[UP]   : in2[UP2]);
    float q_y_0 = 0.5f * (q_y_m05 + q_y_p05);
    float q_y_1 = 0.5f * (q_y_p05 + q_y_p15);
    
    // Calculate h_(i+0.5,j)
    // float h_x_p05 = (in0[CURR] >= 0.f) ? in2[CURR] : in2[RIGHT]; // upwinding
    // float h_y_p05 = (in1[CURR] >= 0.f) ? in2[CURR] : in2[UP];
    float h_x_p05 = (in2[CURR] + in2[RIGHT]) / 2.f;  
    float h_y_p05 = (in2[CURR] + in2[UP]) / 2.f;

    // Calculate corresponding values for u_x_(i,j) using upwinding
    // Why not average here like we average the q's? 
    float u_x_0 = (q_x_0 >= 0.f) ? in0[LEFT] : in0[CURR];
    float u_x_1 = (q_x_1 >= 0.f) ? in0[CURR] : in0[RIGHT]; 
    float u_y_0 = (q_y_0 >= 0.f) ? in1[DOWN] : in1[CURR];
    float u_y_1 = (q_y_1 >= 0.f) ? in1[CURR] : in1[UP];       

    // Compute dux_dt and duy_dt 
    // float dux_dt = - (1/cellSize) * ((q_x_0/h_x_p05) * (in0[CURR] - in0[LEFT]) + (q_y_m05/h_x_p05) * (in0[CURR] - in0[DOWN]));
    // float duy_dt = - (1/cellSize) * ((q_y_0/h_y_p05) * (in1[CURR] - in1[DOWN]) + (q_x_m05/h_y_p05) * (in1[CURR] - in1[LEFT]));
    float dux_dt = - (1/(cellSize * h_x_p05)) * ((q_x_1 * u_x_1 - q_x_0 * u_x_0) - in0[CURR] * (q_x_1 - q_x_0));
    float duy_dt = - (1/(cellSize * h_y_p05)) * ((q_y_1 * u_y_1 - q_y_0 * u_y_0) - in1[CURR] * (q_y_1 - q_y_0));
    dux_dt -= (1/cellSize) * GRAVITY * (in3[RIGHT] - in3[CURR]);
    duy_dt -= (1/cellSize) * GRAVITY * (in3[UP]    - in3[CURR]);
    float ubarNew_x = LimitVelocity(in0[CURR] + timeStep * dux_dt);
    float ubarNew_y = LimitVelocity(in1[CURR] + timeStep * duy_dt);
    
    out0[CURR] = ubarNew_x;
    out1[CURR] = ubarNew_y;
    out2[CURR] = ubarNew_x * ((ubarNew_x >= 0.f) ? in2[CURR] : in2[RIGHT]);
    out3[CURR] = ubarNew_y * ((ubarNew_y >= 0.f) ? in2[CURR] : in2[UP]);
}


/////////// Transport & Combine ///////////

[numthreads(16, 16, 1)]
void UpdateTilde(uint3 id : SV_DispatchThreadID) {
    // Inputs: in0 = ubarNew_x, in1 = ubar_x, in2 = ubarNew_y, in3 = ubar_y, 
    //         in4 = qtildePast_x, in5 = qtildePast_y, in6 = h, in7 = htilde
    // Outputs: out0 = htilde, out1 = qtilde_x, out2 = qtilde_y
    if (id.x >= (uint)(gridSize) || id.y >= (uint)(gridSize)) return;
    // if (id.x <= 0 || id.x >= (uint)(gridSize-1) || id.y <= 0 || id.y >= (uint)(gridSize-1)) return;

    uint2 RDOWN = uint2(min(id.x + 1, gridSize - 1), max(id.y - 1, 0));
    uint2 ULEFT = uint2(max(id.x - 1, 0), min(id.y + 1, gridSize - 1));

    // ubar is on same timestep as h, need to get back to timestep of q 
    float ubar_x_avg = 0.5f * (in0[CURR] + in1[CURR]); 
    float ubar_x_m1  = 0.5f * (in0[LEFT] + in1[LEFT]);
    float ubar_x_p1  = 0.5f * (in0[RIGHT] + in1[RIGHT]);
    float ubar_y_avg = 0.5f * (in2[CURR] + in3[CURR]);
    float ubar_y_m1  = 0.5f * (in2[DOWN] + in3[DOWN]);
    float ubar_y_p1  = 0.5f * (in2[UP] + in3[UP]);

    // ubar divergence
    float div_ubar  = (in0[CURR]  - in0[LEFT] + in2[CURR]  - in2[DOWN])  / cellSize;  // at cell center
    float div_right = (in0[RIGHT] - in0[CURR] + in2[RIGHT] - in2[RDOWN]) / cellSize; // at center of right cell
    float div_up    = (in0[UP]    - in0[ULEFT] + in2[UP]   - in2[CURR])  / cellSize;  // at center of up cell
    float div_ux = 0.5f * (div_ubar + div_right); // at right boundary
    float div_uy = 0.5f * (div_ubar + div_up); // at up boundary
    div_ubar *= (div_ubar < 0.f) ? gammaTransport : 1; // Dampen if converging to avoid breaking waves
    div_ux   *= (div_ux   < 0.f) ? gammaTransport : 1;
    div_uy   *= (div_uy   < 0.f) ? gammaTransport : 1;

    // Update qtilde using bulk flow and ubar divergence: dq/dt = -q * div(ubar)
    float step_x = - ubar_x_avg * timeStep / cellSize; // unitless (cells)
    float step_y = - ubar_y_avg * timeStep / cellSize; // unitless (cells)
    float2 samplePos = float2(id.x + step_x, id.y + step_y);
    out1[CURR] = SampleCubicClamped2D(in4, samplePos) * exp(-div_ux * timeStep);
    out2[CURR] = SampleCubicClamped2D(in5, samplePos) * exp(-div_uy * timeStep);

    // // ubar divergence using central difference
    // float div_ubar = (ubar_x_p1 - ubar_x_m1 + ubar_y_p1  - ubar_y_m1) / (2.f * cellSize);
    // div_ubar *= (div_ubar < 0.f) ? gammaTransport : 1; // dampen if converging to avoid breaking waves
    // // Update qtilde using bulk flow and ubar divergence: dq/dt = -q * div(ubar)
    // float step_x = - ubar_x_avg * timeStep / cellSize; // unitless (cells)
    // float step_y = - ubar_y_avg * timeStep / cellSize; // unitless (cells)
    // float2 samplePos = float2(id.x + step_x, id.y + step_y);
    // out1[CURR] = SampleCubicClamped2D(in4, samplePos) * exp(-div_ubar * timeStep);
    // out2[CURR] = SampleCubicClamped2D(in5, samplePos) * exp(-div_ubar * timeStep);

    // Limit flow to prevent negative water heights and enforce terrain boundaries
    out1[CURR] = LimitFlowRate(out1[CURR], in6[CURR], in6[RIGHT]); 
    out2[CURR] = LimitFlowRate(out2[CURR], in6[CURR], in6[UP]); 
    if (((ubar_x_avg >= 0.f) && (in6[CURR] <= minWaterHeight)) ||
        ((ubar_x_avg < 0.f)  && (in6[RIGHT] <= minWaterHeight)))
        out1[CURR] = 0.f;
    if (((ubar_y_avg >= 0.f) && (in6[CURR] <= minWaterHeight)) ||
        ((ubar_y_avg < 0.f)  && (in6[UP] <= minWaterHeight)))
        out2[CURR] = 0.f;        

    // Update htilde using ubar divergence
    // div_ubar  = (in0[CURR] - in0[LEFT] + in2[CURR] - in2[DOWN]) / cellSize;
    // div_ubar *= (div_ubar < 0.f) ? gammaTransport : 1; // dampen if converging to avoid breaking waves
    out0[CURR] = in7[CURR] * exp(-div_ubar * timeStep);
}

[numthreads(16, 16, 1)]
void CalcQAdvect(uint3 id : SV_DispatchThreadID) {
    // Inputs: in0 = ubarNew_x, in1 = ubarNew_y, in2 = htilde
    // Outputs: out0 = qAdvect_x, out1 = qAdvect_y
    // cubic reconstruction: 1/2 cell accounts for staggered grid, 0.5 * dt accounts for h and u at different 1/2 times
    if (id.x >= (uint)(gridSize) || id.y >= (uint)(gridSize)) return;
    // if (id.x <= 0 || id.x >= (uint)(gridSize-1) || id.y <= 0 || id.y >= (uint)(gridSize-1)) return;

    float step_x = 0.5f - in0[CURR] * (0.5f * timeStep) / cellSize;
    float step_y = 0.5f - in1[CURR] * (0.5f * timeStep) / cellSize;
    float2 samplePos = float2(id.x + step_x, id.y + step_y);
    float h_sample = SampleCubicClamped2D(in2, samplePos);
    out0[CURR] = in0[CURR] * h_sample;
    out1[CURR] = in1[CURR] * h_sample;
}

[numthreads(16, 16, 1)]
void IntegrateH(uint3 id : SV_DispatchThreadID) {
    // Inputs: in0 = qbar_x, in1 = qtilde_x, in2 = qAdvect_x, 
    //         in3 = qbar_y, in4 = qtilde_y, in5 = qAdvect_y, in6 = hPast, in7 = terrain
    // Outputs: out0 = h, out1 = q_x, out2 = q_y
    if (id.x >= (uint)(gridSize) || id.y >= (uint)(gridSize)) return;
    // if (id.x <= 0 || id.x >= (uint)(gridSize-1) || id.y <= 0 || id.y >= (uint)(gridSize-1)) return;

    // q = qbar + qtilde + qAdvect
    float q_x  = in0[CURR] + in1[CURR] + in2[CURR];
    float q_y  = in3[CURR] + in4[CURR] + in5[CURR];
    float q_xm = in0[LEFT] + in1[LEFT] + in2[LEFT];
    float q_ym = in3[DOWN] + in4[DOWN] + in5[DOWN];
    q_x  = LimitFlowRate(q_x,  in6[CURR], in6[RIGHT]);
    q_y  = LimitFlowRate(q_y,  in6[CURR], in6[UP]);
    q_xm = LimitFlowRate(q_xm, in6[LEFT], in6[CURR]);
    q_ym = LimitFlowRate(q_ym, in6[DOWN], in6[CURR]);
    if (StopFlowOnTerrainBoundary(in6, in7, CURR, false))
        q_x = 0.f;
    if (StopFlowOnTerrainBoundary(in6, in7, CURR, true ))
        q_y = 0.f;
    if (StopFlowOnTerrainBoundary(in6, in7, LEFT, false))
        q_xm = 0.f;
    if (StopFlowOnTerrainBoundary(in6, in7, DOWN, true ))
        q_ym = 0.f;

    // sponge layer: damp q near edges to absorb waves, simulate an open boundary
    uint thickness = asuint(floor(0.2f * gridSize)); // cells
    float alpha = 1.f;
    uint dist_x = thickness - min(thickness, (id.x < gridSize/2) ? id.x : gridSize-1-id.x);
    uint dist_y = thickness - min(thickness, (id.y < gridSize/2) ? id.y : gridSize-1-id.y);
    q_x  *= 1 - alpha * pow(dist_x / thickness, 2);
    q_xm *= 1 - alpha * pow(dist_x / thickness, 2);
    q_y  *= 1 - alpha * pow(dist_y / thickness, 2);
    q_ym *= 1 - alpha * pow(dist_y / thickness, 2);

    // update h with divergence of q
    float div_q = (q_x - q_xm + q_y - q_ym) / cellSize;
	out0[CURR] = max(0.f, in6[CURR] - timeStep * div_q);
    out1[CURR] = q_x - in2[CURR]; // qbar + qtilde, removing qAdvect
    out2[CURR] = q_y - in5[CURR];

    out1[CURR] = LimitFlowRate(out1[CURR], in6[CURR], in6[RIGHT]);
    out2[CURR] = LimitFlowRate(out2[CURR], in6[CURR], in6[UP]);
    if (StopFlowOnTerrainBoundary(in6, in7, CURR, false))
        out1[CURR] = 0.f;
    if (StopFlowOnTerrainBoundary(in6, in7, CURR, true))
        out2[CURR] = 0.f;
}


/////////// eWave ///////////
// (At bottom because requires different inputs/outputs for each kernel)
// These all have slight variants on the names so HLSL doesn't get mad for redefining a name used before

RWTexture2D<float2> hHat  : register(u1); // Complex
RWTexture2D<float2> qHat_x: register(u2); // Complex
RWTexture2D<float2> qHat_y: register(u3); // Complex
[numthreads(16, 16, 1)]
void TransferToFFT(uint3 id : SV_DispatchThreadID) {
    // Inputs: in0 = htilde, in1 = qtilde_x, in2 = qtilde_y
    // Outputs: out0 = htildeOld, out1 = hHat, out2 = qHat_x, out3 = qHat_y
    if (id.x >= (uint)(gridSize) || id.y >= (uint)(gridSize)) return;

    // Average htilde in time to get on same timestep as q, then prep variables for FFT
    float h_real = 0.5 * (in0[CURR] + out0[CURR]);
    out0[CURR]   = in0[CURR];
    hHat[CURR]   = float2(h_real, 0.0f);
    qHat_x[CURR] = float2(in1[CURR], 0.0f);
    qHat_y[CURR] = float2(in2[CURR], 0.0f);
}

Texture2D<float2> hhat   : register(t0);
Texture2D<float2> qhat_x : register(t1);
Texture2D<float2> qhat_y : register(t2);
RWTexture2DArray<float2> qhat_x_array: register(u0);
RWTexture2DArray<float2> qhat_y_array: register(u1);
[numthreads(16, 16, 1)]
void CalcEWave(uint3 id : SV_DispatchThreadID) {
    // Inputs: in0 = hHat, in1 = qHat_x, in2 = qHat_y, in3 = depth
    // Outputs: out0 = qHat_x_array, qHat_y_array
    if (id.x >= (uint)(gridSize) || id.y >= (uint)(gridSize)) return;

    // Early return for DC component 
    if (id.x == 0 && id.y == 0) {
        qhat_x_array[id] = qhat_x[CURR];
        qhat_y_array[id] = qhat_y[CURR];
        return;
    }

    ///// Wave Number ///////
    // Calculate the physical size of the grid and the frequency step (dK)
    float domainSize = (float)gridSize * cellSize;
    float dK = 2.0f * PI / domainSize; 
    float N_2 = (float)gridSize / 2.0f;
    // Calculate the physical 2D wavenumber vector components (signed, handling the Nyquist wrap-around)
    int freqX = (id.x < N_2) ? (int)id.x : (int)id.x - (int)gridSize;
    int freqY = (id.y < N_2) ? (int)id.y : (int)id.y - (int)gridSize;
    float kx = (float)freqX * dK; // spatial frequency: radians per meter
    float ky = (float)freqY * dK;
    float k = sqrt(kx * kx + ky * ky);

    // Unit vectors
    float kx_ = kx / k;
    float ky_ = ky / k;
    float kx2 = kx_ * kx_;
    float ky2 = ky_ * ky_;

    /////// Dispersion ///////
    // numerical dispersion correction
    float beta = sqrt((2.0 / (k * cellSize)) * sin(k * cellSize / 2.0)); // From their 1D code, but different from paper
    // float beta = sqrt((2.0 * k / cellSize) * sin(k * cellSize / 2.0)); // correct formula from paper, but not stable?
    // Angular frequency for dispersion relation
    float kh = k * in8[id.z];
    float tanh_kh = (kh > 20.f) ? 1.0f : tanh(kh);
    float omega = sqrt(GRAVITY * k * tanh_kh) / beta;
    float S = sin(omega * timeStep) * omega / (k * k);
    float C = cos(omega * timeStep);
    float Cx = 1 + (C-1) * kx2;
    float Cy = 1 + (C-1) * ky2;
    float Ck = (C-1) * kx_ * ky_;

    /////// Gradient of Shifted hHat ///////
    // Fourier gradient in X: dh/dx = hhat * (i * kx)
    float2 dhdx = ComplexMul(hhat[CURR], float2(0, kx));
    float2 dhdy = ComplexMul(hhat[CURR], float2(0, ky));
    // Phase shift in X: shift by dx/2 by multiplying by e^(i * shiftX) = cos(shiftX) + i*sin(shiftX)
    float shiftX = 0.5f * cellSize * kx;
    float shiftY = 0.5f * cellSize * ky;
    float2 e_mix = float2(cos(shiftX), sin(shiftX));
    float2 e_miy = float2(cos(shiftY), sin(shiftY));
    dhdx = ComplexMul(dhdx, e_mix);
    dhdy = ComplexMul(dhdy, e_miy);
    
    // Shift the q cross-terms to align with their target faces
    float theta_x = 0.5f * cellSize * (ky - kx);
    float theta_y = 0.5f * cellSize * (kx - ky);
    float2 shift_x = float2(cos(theta_x), sin(theta_x)); // e^{i * theta_x} 
    float2 shift_y = float2(cos(theta_y), sin(theta_y)); // e^{i * theta_y}
    float2 qx_shifted = ComplexMul(qhat_x[CURR], shift_x);
    float2 qy_shifted = ComplexMul(qhat_y[CURR], shift_y);
    // float2 qx_shifted = qhat_x[CURR];
    // float2 qy_shifted = qhat_y[CURR];

    /////// Update Q ///////
    // 1) Decompose q into parallel and perpendicular: q_|| = kx_*qx + ky_*qy, q_T = kx_*qy - ky_*qx
    // 2) Update in rotated basis: q_||+ = C*q_|| - S*dhdx (q_T unchanged bc Airy is irrotational)
    // 3) Recombine: q = q_T + q_||+
    float qx_r = Cx * qhat_x[CURR].x + Ck * qy_shifted.x - S * dhdx.x;
    float qx_i = Cx * qhat_x[CURR].y + Ck * qy_shifted.y - S * dhdx.y;
    float qy_r = Cy * qhat_y[CURR].x + Ck * qx_shifted.x - S * dhdy.x;
    float qy_i = Cy * qhat_y[CURR].y + Ck * qx_shifted.y - S * dhdy.y;
    // float qx_r = C * qhat_x[CURR].x - S * dhdx.x; // Naive 1D translation, but does just as well
    // float qx_i = C * qhat_x[CURR].y - S * dhdx.y;
    // float qy_r = C * qhat_y[CURR].x - S * dhdy.x;
    // float qy_i = C * qhat_y[CURR].y - S * dhdy.y;
    qhat_x_array[id] = float2(qx_r, qx_i);
    qhat_y_array[id] = float2(qy_r, qy_i);

    if ((id.x == N_2 && id.y == 0) || (id.x == 0 && id.y == N_2) || (id.x == N_2 && id.y == N_2)) {
        qhat_x_array[id] = float2(qhat_x_array[id].x, 0.0f);
        qhat_y_array[id] = float2(qhat_y_array[id].x, 0.0f);
    }
}

Texture2D<float>       hbar        : register(t0);
Texture2DArray<float2> qHat_x_array: register(t1);
Texture2DArray<float2> qHat_y_array: register(t2);
RWTexture2D<float>     qtilde_x    : register(u0);
RWTexture2D<float>     qtilde_y    : register(u1);
[numthreads(16, 16, 1)]
void InterpQ(uint3 id : SV_DispatchThreadID) {
    // Inputs: qtilde_x_array, qtilde_y_array, hbar, depth
    // Outputs: qtilde_x, qtilde_y
    if (id.x >= (uint)(gridSize) || id.y >= (uint)(gridSize)) return;

    float waterDepth_x = max(hbar[CURR], hbar[CURR + uint2(1, 0)]);
    float waterDepth_y = max(hbar[CURR], hbar[CURR + uint2(0, 1)]);
    int d1_x = 0;
    int d1_y = 0;
    for (int d = 0; d < depthNum; d++) {
        if (waterDepth_x >= in8[d]) d1_x = d;
        if (waterDepth_y >= in8[d]) d1_y = d;
    }
    int d2_x = min(depthNum - 1, d1_x + 1);
    int d2_y = min(depthNum - 1, d1_y + 1);
    float sx = 0.f;
    float sy = 0.f;
    if (d1_x != d2_x)
        sx = (in8[d2_x] - waterDepth_x) / (in8[d2_x] - in8[d1_x]);
    if (d1_y != d2_y)
        sy = (in8[d2_y] - waterDepth_y) / (in8[d2_y] - in8[d1_y]);
    qtilde_x[CURR] = sx * qHat_x_array[uint3(CURR, d1_x)].x + (1.f - sx) * qHat_x_array[uint3(CURR, d2_x)].x;
    qtilde_y[CURR] = sy * qHat_y_array[uint3(CURR, d1_y)].x + (1.f - sy) * qHat_y_array[uint3(CURR, d2_y)].x;
}