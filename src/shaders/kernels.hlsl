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

    // Padding for 16-byte alignment 
    float buffer1; 
    float buffer2; 
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


//////////////////// HELPER FUNCTIONS /////////////////////////
float calcGradient(Texture2D<float> f, Texture2D<float> a, uint2 curr) {
    uint2 left = curr - uint2(1,0);
    uint2 right = curr + uint2(1,0);
    uint2 up = curr + uint2(0,1);
    uint2 down = curr - uint2(0,1);
    float dF_x = a[curr] * (f[right] - f[curr]) - a[left] * (f[curr] - f[left]);
    float dF_y = a[curr] * (f[up] - f[curr]) - a[down] * (f[curr] - f[down]);
    float dFdX = (dF_x + dF_y) / (cellSize * cellSize);
    return dFdX;
}

bool stopFlowOnTerrainBoundary(Texture2D<float> h, Texture2D<float> terrain, uint2 curr, bool isYDirection=false) {
    uint2 xplus = curr + uint2(1, 0);
    uint2 yplus = curr + uint2(0, 1);
    // test if the terrain boundary stops any flow across x+0.5
	if (!isYDirection) {
		if ((h[curr] <= minWaterHeight) && (terrain[curr] >= terrain[xplus] + h[xplus])) // positive q_x
			return true;
		if ((h[xplus] <= minWaterHeight) && (terrain[xplus] > terrain[curr] + h[curr])) // negative q_x
			return true;
		return false;
	}
	else {
		if ((h[curr] <= minWaterHeight) && (terrain[curr] >= terrain[yplus] + h[yplus])) // positive q_y
			return true;
		if ((h[yplus] <= minWaterHeight) && (terrain[yplus] > terrain[curr] + h[curr])) // negative q_y
			return true;
		return false;
	}
}

float Sim::SampleCubicClamped(std::vector<float>& dataField, float samplePos, uint2 curr, bool isYDirection=false) {
	// cubic interpolation with Catmull-Rom Spline https://en.wikipedia.org/wiki/Cubic_Hermite_spline#Interpolating_a_data_set
	int id_start = floor(samplePos) - 1;
	int id0 = max(0, min(id_start + 0, gridSize - 1));
	int id1 = max(0, min(id_start + 1, gridSize - 1));
	int id2 = max(0, min(id_start + 2, gridSize - 1));
	int id3 = max(0, min(id_start + 3, gridSize - 1));
	float fx = max(0.f, min(1.f, samplePos - floor(samplePos)));
	float x2 = fx * fx;
	float x3 = x2 * fx;
	float xcubicX = -x3 + 2.f * x2 - fx;
	float xcubicY =  3.f * x3 - 5.f * x2 + 2.f;
	float xcubicZ = -3.f * x3 + 4.f * x2 + fx;
	float xcubicW =  x3 - x2;
    if (isYDirection)
        float out = 0.5f * (xcubicX * dataField[id0 * gridSize + curr.x] + xcubicY * dataField[id1 * gridSize + curr.x] + 
                       xcubicZ * dataField[id2 * gridSize + curr.x] + xcubicW * dataField[id3 * gridSize + curr.x]);
        out = min(max(dataField[id1 * gridSize + curr.x], dataField[id2 * gridSize + curr.x]), out);  // value-limiting
	    out = max(min(dataField[id1 * gridSize + curr.x], dataField[id2 * gridSize + curr.x]), out);  // value-limiting
    else
        float out = 0.5f * (xcubicX * dataField[curr.y * gridSize + id0] + xcubicY * dataField[curr.y * gridSize + id1] + 
                       xcubicZ * dataField[curr.y * gridSize + id2] + xcubicW * dataField[curr.y * gridSize + id3]);
        out = min(max(dataField[curr.y * gridSize + id1], dataField[curr.y * gridSize + id2]), out);  // value-limiting
	    out = max(min(dataField[curr.y * gridSize + id1], dataField[curr.y * gridSize + id2]), out);  // value-limiting
	return out;
}

//////////////////// COMPUTE SHADERS /////////////////////////


// Apply Boundary Conditions
[numthreads(16, 16, 1)]
void ApplyBoundaries(uint3 id : SV_DispatchThreadID){
    uint x = id.x;
    uint y = id.y;
    if (x >= (uint)gridSize || y >= (uint)gridSize) return;
    bool left   = (x == 0);
    bool right  = (x == (uint)gridSize - 1);
    bool bottom = (y == 0);
    bool top    = (y == (uint)gridSize - 1);
    if (!(left || right || top || bottom)) return;

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
        // out3[id.xy] = in3[src];
    } 
    else if (boundaryType == 1) {// free (linear interpolation)
        uint2 dir = uint2(
            (left ? 1 : right ? -1 : 0),
            (bottom ? 1 : top ? -1 : 0)
        );
        out0[id.xy] = 2.0f * in0[src] - in0[src + dir];
        out1[id.xy] = 2.0f * in1[src] - in1[src + dir];
        out2[id.xy] = 2.0f * in2[src] - in2[src + dir];
        // out3[id.xy] = 2.0f * in3[src] - in3[src + dir];
    }
    else if (boundaryType == 2) { // zero (absorbing)
        out0[id.xy] = waterLevel;
        out1[id.xy] = waterLevel; 
        out2[id.xy] = waterLevel;
        // out3[id.xy] = waterLevel;
    } 
}

/////////// Decomposition ///////////

// Decomposition: Initialize 
[numthreads(16, 16, 1)]
void InitDecomp(uint3 id : SV_DispatchThreadID) {
    // Inputs: in0 = h, in1 = q_x, in2 = q_y, in3 = terrain
    // Outputs: out0 = H, out1 = Q_x, out2 = Q_y
    out0[id.xy] = in3[id.xy] + in0[id.xy];
    out1[id.xy] = in1[id.xy];
    out2[id.xy] = in2[id.xy];
}

// Decomposition: Calculate Diffusion Coefficients 
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

// Decomposition: Diffusion Iteration
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

// Decomposition: Decompose Fields
[numthreads(16, 16, 1)]
void DecomposeFields(uint3 id : SV_DispatchThreadID) {
    // Inputs: in0 = H, in1 = Q_x, in2 = Q_y
    //         in3 = h, in4 = q_x, in5 = q_y, in6 = terrain
    // Outputs: out0 = hbar, out1 = qbar_x, out2 = qbar_y, 
    //          out3 = htilde, out4 = qtilde_x, out5 = qtilde_y
    if (id.x < 1 || id.x >= (uint)(gridSize - 1) || id.y < 1 || id.y >= (uint)(gridSize - 1)) return;
    current = id.xy;

    out0[curr] = std::max(0.f, H[curr] - in0[curr]);
	out1[curr] = Q_x[curr];
	out2[curr] = Q_y[curr];
	out3[curr] = h[curr] - out0[curr];
	out4[curr] = q_x[curr] - out1[curr];
	out5[curr] = q_y[curr] - out2[curr];

    if (stopFlowOnTerrainBoundary(in3, in6, curr, false)) { // stop flow in x direction
        qbar_x[curr] = 0.f;
        qtilde_x[curr] = 0.f;
    }
    if (stopFlowOnTerrainBoundary(in3, in6, curr, true)) { // stop flow in y direction
        qbar_y[curr] = 0.f;
        qtilde_y[curr] = 0.f;
    }
}


/////////// eWave ///////////




/////////// SWE ///////////

// [numthreads(16, 16, 1)]
// void CalcDiffusionCoeffs(uint3 id : SV_DispatchThreadID) {}

// [numthreads(16, 16, 1)]
// void CalcDiffusionCoeffs(uint3 id : SV_DispatchThreadID) {}

// [numthreads(16, 16, 1)]
// void CalcDiffusionCoeffs(uint3 id : SV_DispatchThreadID) {}


/////////// Transport ///////////

[numthreads(16, 16, 1)]
void UpdateTilde(uint3 id : SV_DispatchThreadID) {
    // Inputs: in0 = ubarNew_x, in1 = ubar_x, in2 = ubarNew_y, in3 = ubar_y, 
    //         in4 = qtildePast_x, in5 = qtildePast_y, in6 = h
    // Outputs: out0 = qtilde_x, out1 = qtilde_y, out2 = htilde

    // Adjust qtilde to account for bulk flow
    float bulkVelocity_x = 0.5f * (in0[id.xy] + in1[id.xy]); // ubar is on same timestep as h, need to get back to timestep of q 
    float bulkVelocity_y = 0.5f * (in2[id.xy] + in3[id.xy]);
    float step_x = - bulkVelocity_x * timeStep / cellSize; // unitless (cells)
    float step_y = - bulkVelocity_y * timeStep / cellSize;

    // ubar divergence using central difference: dq/dt = -q * div(ubar)
    float ubar_x_m1 = 0.5f * (in0[uint2(id.x-1, id.y)] + in1[uint2(id.x-1, id.y)]);
    float ubar_x_p1 = 0.5f * (in0[uint2(id.x+1, id.y)] + in1[uint2(id.x+1, id.y)]);
    float ubar_y_m1 = 0.5f * (in2[uint2(id.x, id.y-1)] + in3[uint2(id.x, id.y-1)]);
    float ubar_y_p1 = 0.5f * (in2[uint2(id.x, id.y+1)] + in3[uint2(id.x, id.y+1)]);
    float div_ubar_x = (ubar_x_p1 - ubar_x_m1) / (2.f * cellSize);
    float div_ubar_y = (ubar_y_p1 - ubar_y_m1) / (2.f * cellSize);
    float div_ubar = div_ubar_x + div_ubar_y;
    if (div_ubar < 0.f) // dampen if converging to avoid breaking waves
        div_ubar *= gammaTransport;

    // Update qtilde using bulk flow and ubar divergence
    if (((bulkVelocity_x >= 0.f) && (in6[id.xy] < minWaterHeight)) ||
        ((bulkVelocity_x < 0.f) && (in6[uint2(id.x + 1, id.y)] < minWaterHeight)))
        out0[id.xy] = 0.f;
    else
        out0[id.xy] = SampleCubicClamped(id.x + step_x, in4, id.xy, false) * exp(-div_ubar * timeStep);
    if (((bulkVelocity_y >= 0.f) && (in6[id.xy] < minWaterHeight)) ||
        ((bulkVelocity_y < 0.f) && (in6[uint2(id.x, id.y + 1)] < minWaterHeight)))
        out1[id.xy] = 0.f;
    else
        out1[id.xy] = SampleCubicClamped(id.y + step_y, in5, id.xy, true) * exp(-div_ubar * timeStep);

    // Update htilde using ubar divergence
    div_ubar_x = (in0[id.xy] - in0[uint2(id.x-1, id.y)]) / cellSize;
    div_ubar_y = (in2[id.xy] - in2[uint2(id.x, id.y-1)]) / cellSize;
    div_ubar = div_ubar_x + div_ubar_y;
    if (div_ubar < 0.f) // dampen if converging to avoid breaking waves
        div_ubar *= GAMMA_TRANSPORT;
    htilde[id.xy] *= exp(-div_ubar * timeStep);
}

[numthreads(16, 16, 1)]
void UpdateHtilde(uint3 id : SV_DispatchThreadID) {

}




