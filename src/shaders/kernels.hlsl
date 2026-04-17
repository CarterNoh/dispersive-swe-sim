// Constant buffer matching the C++ struct
cbuffer Constants : register(b0) {
    int gridSize;
    float cellSize;
    float deltaT;
    float diffusionPenalty;
    int diffusionIterations;
    int boundaryType;
    float buffer2;
    float buffer3;
};

// Texture Registers
Texture2D<float>   in0 : register(t0);
Texture2D<float>   in1 : register(t1);
Texture2D<float>   in2 : register(t2);
Texture2D<float>   in3 : register(t3);
Texture2D<float>   in4 : register(t4);
Texture2D<float>   in5 : register(t5);
Texture2D<float>   in6 : register(t6);
Texture2D<float>   in7 : register(t7);
RWTexture2D<float> out0: register(u0);
RWTexture2D<float> out1: register(u1);
RWTexture2D<float> out2: register(u2);
RWTexture2D<float> out3: register(u3);
RWTexture2D<float> out4: register(u4);
RWTexture2D<float> out5: register(u5);


//////////////////// Compute Shaders /////////////////////////

// Apply Boundary Conditions
[numthreads(16, 16, 1)]
void ApplyBoundaries(uint3 id : SV_DispatchThreadID){
    // Takes up to four fields
    uint x = id.x;
    uint y = id.y;
    if (x >= (uint)gridSize || y >= (uint)gridSize) return;
    bool left   = (x == 0);
    bool right  = (x == (uint)gridSize - 1);
    bool bottom = (y == 0);
    bool top    = (y == (uint)gridSize - 1);
    bool isEdge = left || right || top || bottom;
    if (!isEdge) return;

    // Pick a source interior cell depending on boundary side
    uint2 src;
    if (left && bottom)
        src = uint2(1, 1);
    else if (left && top)
        src = uint2(1, gridSize-2);
    else if (right && bottom)
        src = uint2(gridSize-2, 1);
    else if (right && top)
        src = uint2(gridSize-2, gridSize-2);
    else if (left)
        src = uint2(1, y);
    else if (right)
        src = uint2(gridSize-2, y);
    else if (bottom)
        src = uint2(x, 1);
    else
        src = uint2(x, gridSize-2);

    // Apply boundary condition
    if (boundaryType == 0) { // wall (copy neighbor)
        out0[id.xy] = in0[src]; 
        out1[id.xy] = in1[src];
        out2[id.xy] = in2[src];
        out3[id.xy] = in3[src];
    } 
    else if (boundaryType == 1) {// free (linear extrapolation)
        uint2 dir = uint2(
            (left ? 1 : right ? -1 : 0),
            (bottom ? 1 : top ? -1 : 0)
        );
        out0[id.xy] = 2.0f * in0[src] - in0[src + dir];
        out1[id.xy] = 2.0f * in1[src] - in1[src + dir];
        out2[id.xy] = 2.0f * in2[src] - in2[src + dir];
        out3[id.xy] = 2.0f * in3[src] - in3[src + dir];
    }
    else if (boundaryType == 2) { // zero (absorbing)
        out0[id.xy] = 0.0f;
        out1[id.xy] = 0.0f; 
        out2[id.xy] = 0.0f;
        out3[id.xy] = 0.0f;
    } 
}

// Calculate Diffusion Coefficients 
[numthreads(16, 16, 1)]
void CalcDiffusionCoeffs(uint3 id : SV_DispatchThreadID) {
    // Inputs: in0 = H
    // Outputs: out0 = alpha_H, out1 = alpha_Q_x, out2 = alpha_Q_y

    /*	NOTE: Their implementation uses an average height across neighbor cells (I assume to improve stability, 
	but reduces accuracy?). I am choosing to use only local cell values to align with eqn from paper, we'll 
	see how it does. 

	Their code for initial calculation:
	// Identify the correct height (sigma) to use for diffusivity calculation
	alpha_H[idx(x,y)] = 0.f;
	float maxGround = std::max(terrain[idx(x,y)], terrain[idx(x+1,y)], terrain[idx(x,y+1)]); // Why do this?
	float minWaterlevel = (H[idx(x,y)] + H[idx(x+1,y)] + H[idx(x,y+1)]) / 3.f; // Why average here?
	if ((h[idx(x,y)] > 0.f) && (h[idx(x+1,y)] > 0.f) && (h[idx(x,y+1)] > 0.f))
	{
		static const float sigma_max = 8.f;
		// they limit diffusion coefficient to between 0 and 1, maybe for stability?
		float sigma = std::min(sigma_max, std::max(0.f, minWaterlevel - maxGround));
		alpha_H[idx(x,y)] = sigma * sigma / (2*DELTA_T*DIFFUSION_ITERATIONS);
	} 
	*/

    if (id.x >= (uint)(gridSize - 1) || id.y >= (uint)(gridSize - 1)) return;    
    uint2 curr = id.xy;
    uint2 right = uint2(id.x + 1, id.y);
    uint2 up = uint2(id.x, id.y + 1);
    float max_alpha = (cellSize * cellSize) / (4.0f * deltaT);
    float denom = 2.0f * deltaT * diffusionIterations;
    // Alpha_H
    float grad_x = (in0[right] - in0[curr]) / cellSize;
    float grad_y = (in0[up] - in0[curr]) / cellSize;
    float penalty = -diffusionPenalty * (grad_x * grad_x + grad_y * grad_y);
    out0[curr] = min(max_alpha, in0[curr] * in0[curr] / denom) * exp(penalty);
    // Alpha_Q
    float avg_H_x = 0.5f * (in0[curr] + in0[right]);
    float avg_H_y = 0.5f * (in0[curr] + in0[up]);
    out1[curr] = min(max_alpha, avg_H_x * avg_H_x / denom) * exp(penalty);
    out2[curr] = min(max_alpha, avg_H_y * avg_H_y / denom) * exp(penalty);
}

// Diffusion Iteration
float calcGradient(Texture2D<float> f, Texture2D<float> a, uint2 curr) {
    uint2 left = curr - uint2(1,0);
    uint2 right = curr + uint2(1,0);
    uint2 up = curr + uint2(0,1);
    uint2 down = curr - uint2(0,1);
    float dF_x = a[curr] * (f[right] - f[curr]) - a[left] * (f[curr] - f[left]);
    float dF_y = a[curr] * (f[up] - f[curr]) - a[down] * (f[curr] - f[down]);
    float dFdT = (dF_x + dF_y) / (cellSize * cellSize);
    return dFdT;
}
[numthreads(16, 16, 1)]
void DiffusionStep(uint3 id : SV_DispatchThreadID) {
    // Inputs: in0 = terrain, in1 = HPast, in2 = QPast_x, in3 = QPast_y, 
    //         in4 = alpha_H, in5 = alpha_Q_x, in6 = alpha_Q_y
    // Outputs: out0 = H, out1 = Q_x, out3 = Q_y

    if (id.x < 1 || id.x >= (uint)(gridSize - 1) || id.y < 1 || id.y >= (uint)(gridSize - 1)) return;
    uint2 curr = id.xy;
    uint2 left = uint2(id.x - 1, id.y);
    uint2 right = uint2(id.x + 1, id.y);
    uint2 down = uint2(id.x, id.y - 1);
    uint2 up = uint2(id.x, id.y + 1);
    float newH = in1[curr] + deltaT * calcGradient(in1, in4, curr);
    out0[curr] = max(in0[curr], newH);
    out1[curr] = in2[curr] + deltaT * calcGradient(in2, in5, curr);
    out2[curr] = in3[curr] + deltaT * calcGradient(in3, in6, curr);
}




