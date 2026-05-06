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
    float depth[5];
    float depthNum;

    // // Padding for 16-byte alignment 
    // float buffer;  
};

#define GRAVITY 9.80665
#define PI 3.14159265358979323846f

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
float CalcGradient(Texture2D<float> f, Texture2D<float> a, uint2 curr) {
    uint2 left = curr - uint2(1,0);
    uint2 right = curr + uint2(1,0);
    uint2 up = curr + uint2(0,1);
    uint2 down = curr - uint2(0,1);
    float dF_x = a[curr] * (f[right] - f[curr]) - a[left] * (f[curr] - f[left]);
    float dF_y = a[curr] * (f[up] - f[curr]) - a[down] * (f[curr] - f[down]);
    float dFdX = (dF_x + dF_y) / (cellSize * cellSize);
    return dFdX;
}

bool StopFlowOnTerrainBoundary(Texture2D<float> h, Texture2D<float> terrain, uint2 curr, bool isYDirection=false) {
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

float LimitVelocity(float velocity_in) {
	if (velocity_in >= 0.f)
		return min(velocity_in, cflCondition * cellSize / timeStep);   // 0.25 since other neighbors might take from this source cell as well
	else
		return max(velocity_in, -cflCondition * cellSize / timeStep);
}

float LimitFlowRate(float flow_rate_in, float waterDepth_left, float waterDepth_right) {
	if (flow_rate_in >= 0.f)
		return min(flow_rate_in, cflCondition * waterDepth_left * cellSize / timeStep);
	else
		return max(flow_rate_in, -cflCondition * waterDepth_right * cellSize / timeStep);
}

float SampleCubicClamped(Texture2D<float> dataField, float samplePos, uint2 curr, bool isYDirection=false) {
	// cubic interpolation with Catmull-Rom Spline https://en.wikipedia.org/wiki/Cubic_Hermite_spline#Interpolating_a_data_set
	uint id_start = floor(samplePos) - 1;
	uint id0 = max(0, min(id_start + 0, (uint)gridSize - 1));
	uint id1 = max(0, min(id_start + 1, (uint)gridSize - 1));
	uint id2 = max(0, min(id_start + 2, (uint)gridSize - 1));
	uint id3 = max(0, min(id_start + 3, (uint)gridSize - 1));
	float fx = max(0.f, min(1.f, samplePos - floor(samplePos)));
	float x2 = fx * fx;
	float x3 = x2 * fx;
	float X = -x3 + 2.f * x2 - fx;
	float Y =  3.f * x3 - 5.f * x2 + 2.f;
	float Z = -3.f * x3 + 4.f * x2 + fx;
	float W =  x3 - x2;

    float result; 
    if (isYDirection) {
        float data0 = dataField[uint2(curr.x, id0)];
        float data1 = dataField[uint2(curr.x, id1)];
        float data2 = dataField[uint2(curr.x, id2)];
        float data3 = dataField[uint2(curr.x, id3)];
        result = 0.5f * (X * data0 + Y * data1 + Z * data2 + W * data3);
        result = min(max(data1, data2), result);  // value-limiting
	    result = max(min(data1, data2), result);  // value-limiting
    }
    else {
        float data0 = dataField[uint2(id0, curr.y)];
        float data1 = dataField[uint2(id1, curr.y)];
        float data2 = dataField[uint2(id2, curr.y)];
        float data3 = dataField[uint2(id3, curr.y)];
        result = 0.5f * (X * data0 + Y * data1 + Z * data2 + W * data3);
        result = min(max(data1, data2), result);  // value-limiting
	    result = max(min(data1, data2), result);  // value-limiting
    }
	return result;
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
        out0[id.xy] = out0[src];
        out1[id.xy] = out1[src];
        out2[id.xy] = out2[src];
        out3[id.xy] = out3[src];
    }
    else if (boundaryType == 1) {// free (linear interpolation)
        uint2 dir = uint2(
            (left ? 1 : right ? -1 : 0),
            (bottom ? 1 : top ? -1 : 0)
        );
        out0[id.xy] = 2.0f * out0[src] - out0[src + dir];
        out1[id.xy] = 2.0f * out1[src] - out1[src + dir];
        out2[id.xy] = 2.0f * out2[src] - out2[src + dir];
        out3[id.xy] = 2.0f * out3[src] - out3[src + dir];
    }
}

/////////// Decomposition ///////////

[numthreads(16, 16, 1)]
void InitDecomp(uint3 id : SV_DispatchThreadID) {
    // Inputs: in0 = h, in1 = q_x, in2 = q_y, in3 = terrain
    // Outputs: out0 = H, out1 = Q_x, out2 = Q_y
    out0[id.xy] = in0[id.xy] + in3[id.xy];
    out1[id.xy] = in1[id.xy];
    out2[id.xy] = in2[id.xy];
}

[numthreads(16, 16, 1)]
void CalcDiffusionCoeffs(uint3 id : SV_DispatchThreadID) {
    // Inputs: in0 = H, in1 = terrain
    // Outputs: out0 = alpha_H, out1 = alpha_Q_x, out2 = alpha_Q_y
    if (id.x < 1 || id.x >= (uint)(gridSize - 1) || id.y < 1 || id.y >= (uint)(gridSize - 1)) return;

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
 
    uint2 curr = id.xy;
    uint2 right = curr + uint2(1, 0);
    uint2 up = curr + uint2(0, 1);

    // float max_ground = max(in1[curr], in1[right], in1[up]);
    // float min_water = min(in0[curr], in0[right], in0[up]);
    float max_alpha = (cellSize * cellSize) / (4.0f * deltaT);
    float denom = 2.0f * deltaT * diffusionIterations;
    // Alpha_H
    float grad_x = (in0[right] - in0[curr]) / cellSize;
    float grad_y = (in0[up] - in0[curr]) / cellSize;
    float penalty = -diffusionPenalty * (grad_x * grad_x + grad_y * grad_y);
    // out0[curr] = max(0.f, min(max_alpha, in0[curr] * in0[curr] / denom) * exp(penalty));
    // Alpha_Q
    float avg_H_x = 0.5f * (in0[curr] + in0[right]);
    float avg_H_y = 0.5f * (in0[curr] + in0[up]);
    // out1[curr] = max(0.f, min(max_alpha, avg_H_x * avg_H_x / denom) * exp(penalty));
    // out2[curr] = max(0.f, min(max_alpha, avg_H_y * avg_H_y / denom) * exp(penalty));
    out0[curr] = min(max_alpha, in0[curr] * in0[curr] / denom) * exp(penalty);
    out1[curr] = min(max_alpha, avg_H_x * avg_H_x / denom) * exp(penalty);
    out2[curr] = min(max_alpha, avg_H_y * avg_H_y / denom) * exp(penalty);
}

[numthreads(16, 16, 1)]
void DiffusionStep(uint3 id : SV_DispatchThreadID) {
    // Inputs: in0 = terrain, in1 = HPast, in2 = QPast_x, in3 = QPast_y, 
    //         in4 = alpha_H, in5 = alpha_Q_x, in6 = alpha_Q_y
    // Outputs: out0 = H, out1 = Q_x, out3 = Q_y
    if (id.x < 1 || id.x >= (uint)(gridSize - 1) || id.y < 1 || id.y >= (uint)(gridSize - 1)) return;

    float newH = in1[id.xy] + deltaT * CalcGradient(in1, in4, id.xy);
    out0[id.xy] = max(in0[id.xy], newH);
    out1[id.xy] = in2[id.xy] + deltaT * CalcGradient(in2, in5, id.xy);
    out2[id.xy] = in3[id.xy] + deltaT * CalcGradient(in3, in6, id.xy);
}

[numthreads(16, 16, 1)]
void DecomposeFields(uint3 id : SV_DispatchThreadID) {
    // Inputs: in0 = H, in1 = Q_x, in2 = Q_y
    //         in3 = h, in4 = q_x, in5 = q_y, in6 = terrain
    // Outputs: out0 = hbar, out1 = qbar_x, out2 = qbar_y, 
    //          out3 = htilde, out4 = qtilde_x, out5 = qtilde_y
    if (id.x >= (uint)(gridSize - 1) || id.y >= (uint)(gridSize - 1)) return;

    float hbar = max(0.f, in0[id.xy] - in6[id.xy]);
    float qbar_x = in1[id.xy];
    float qbar_y = in2[id.xy];
    float htilde = in3[id.xy] - hbar;
    float qtilde_x = in4[id.xy] - qbar_x;
    float qtilde_y = in5[id.xy] - qbar_y;

    if (StopFlowOnTerrainBoundary(in3, in6, id.xy, false)) { // stop flow in x direction
        qbar_x = 0.f;
        qtilde_x = 0.f;
    }
    if (StopFlowOnTerrainBoundary(in3, in6, id.xy, true)) { // stop flow in y direction
        qbar_y = 0.f;
        qtilde_y = 0.f;
    }

    out0[id.xy] = hbar;
    out1[id.xy] = qbar_x;
    out2[id.xy] = qbar_y;
	out3[id.xy] = htilde;
    out4[id.xy] = qtilde_x;
    out5[id.xy] = qtilde_y;
}


/////////// SWE ///////////

[numthreads(16, 16, 1)]
void CalcUbar(uint3 id : SV_DispatchThreadID) {
    // Inputs: in0 = qbar_x, in1 = qbar_y, in2 = hbarOld
    // Outputs: out0 = ubar_x, out1 = ubar_y
    if (id.x < 1 || id.x >= (uint)(gridSize - 1) || id.y < 1 || id.y >= (uint)(gridSize - 1)) return;
    
    float ubar_x = in0[id.xy];
    float ubar_y = in1[id.xy];

    // First-Order Up-Winding
    if (ubar_x >= 0.f || id.x == gridSize-1)
        ubar_x /= max(minWaterHeight, in2[id.xy]);
    else
        ubar_x /= max(minWaterHeight, in2[id.xy + uint2(1,0)]);
    if (ubar_y >= 0.f || id.y == gridSize-1)
        ubar_y /= max(minWaterHeight, in2[id.xy]);
    else
        ubar_y /= max(minWaterHeight, in2[id.xy + uint2(0,1)]);

    // Enforcing CFL condition for later surface waves advection
    out0[id.xy] = LimitVelocity(ubar_x);  
    out1[id.xy] = LimitVelocity(ubar_y); 
}

[numthreads(16, 16, 1)]
void CalcSWE(uint3 id : SV_DispatchThreadID) {
    // Inputs: in0 = ubar_x, in1 = ubar_y, in2 = hbar, in3 = H
    // Outputs: out0 = ubarNew_x, out1 = ubarNew_y, out2 = qbar_x, out3 = qbar_y
    if (id.x < 1 || id.x >= (uint)(gridSize - 1) || id.y < 1 || id.y >= (uint)(gridSize - 1)) return;

    uint2 curr = id.xy;
    uint2 right = curr + uint2(1,0);
    uint2 left = curr - uint2(1,0);
    uint2 down = curr - uint2(0,1);
    uint2 up = curr + uint2(0,1);
    // Compute intermediate values needed for du/dt calculations
    // Need: q_x/y_ij, q_x_i-0.5_j, q_y_i_j-0.5, 
    // 	     h_i+0.5_j, h_i_j+0.5, h_ij=hbar[x,y], h_i+1_j=hbar[x+1,y], h_i_j+1=hbar[x,y+1],
    //       u_x_i+0.5_j=ubar_x[x,y], u_x_i-0.5_j=ubar_x[x-1,y], 
    //       u_y_i_j+0.5=ubar_y[x,y], u_y_i_j-0.5=ubar_y[x,y-1], 
    //       u_x_i+0.5_j-1=ubar_x[x,y-1], u_y_i-1_j+0.5=ubar_y[x-1,y]
    // Note: The 1D implementation uses far more intermeidate values that the paper, and I don't know why. 
    // The commented sections below are the conversion of their 1D code with all its complexity, and the 
    // uncommented sections are my simplified version that tries to stay more true to the paper. We'll see 
    // if it causes stability/accuracy issues.

    ////// X DIRECTION //////
    // Use upwinding to evaluate q_(i-0.5,j), q_(i+0.5,j), q_(i+1.5,j) 
    // for x direction to get q_x_(i,j), q_x_(i+1,j)
    float q_x_m05 = in0[left];
    if (q_x_m05 >= 0.f)
        q_x_m05 *= in2[left];
    else
        q_x_m05 *= in2[curr];
    float q_x_p05 = in0[curr];  
    if (q_x_p05 >= 0.f)
        q_x_p05 *= in2[curr];
    else
        q_x_p05 *= in2[right];
    float q_x_0 = 0.5f * (q_x_m05 + q_x_p05);
    
    // Calculate h_(i+0.5,j) and h_(i-0.5,j) using upwinding
    float h_x_p05 = 0.f;
    if (in0[curr] >= 0.f)
        h_x_p05 = in2[curr];
    else
        h_x_p05 = in2[right];
    

    // float h_x_p05 = (in2[curr] + in2[right]) / 2.f; // averaging, for some reason the old paper used this
    // float q_x_p15 = in0[right];  //q_(i+1.5,j) = hfr at position x
    // if (q_x_p15 >= 0.f)
    // 	q_x_p15 *= in2[right];
    // else
    // 	q_x_p15 *= in2[uint2(min(id.x + 2, gridSize-1), id.y)];
    // float q_x_p1 = 0.5f * (q_x_p05 + q_x_p15);
    // // Calculate corresponding values for u_x_(i,j) using upwinding
    // // (why do we use upwinding here instead of averaging like q?)
    // float u_star_x_0 = 0.f;
    // if (q_x_0 >= 0.f)
    // 	u_star_x_0 = in0[left];
    // else
    // 	u_star_x_0 = in0[curr];
    // float u_star_x_p1 = 0.f;
    // if (q_x_p1 > 0.f)
    // 	u_star_x_p1 = in0[curr];
    // else
    // 	u_star_x_p1 = in0[right];
    

    /////// Y DIRECTION //////
    // Use upwinding to evaluate q_(i,j-0.5), q_(i,j+0.5), q_(i,j+1.5) 
    // for y direction to get q_y_(i,j), q_y_(i,j+1)
    float q_y_m05 = in1[down];
    if (q_y_m05 >= 0.f)
        q_y_m05 *= in2[down];
    else
        q_y_m05 *= in2[curr];
    float q_y_p05 = in1[curr];  
    if (q_y_p05 >= 0.f)
        q_y_p05 *= in2[curr];
    else
        q_y_p05 *= in2[up];
    float q_y_0 = 0.5f * (q_y_m05 + q_y_p05);

    // Calculate h_(i,j+0.5) and h_(i,j-0.5) using upwinding
    float h_y_p05 = 0.f;
    if (in1[curr] >= 0.f)
        h_y_p05 = in2[curr];
    else
        h_y_p05 = in2[up];

    // float h_y_p05 = (in2[curr] + in2[right]) / 2.f; // averaging, for some reason the old paper used this
    // float q_y_p15 = in1[up];  //q_(i+1.5,j) = hfr at position x
    // if (q_y_p15 >= 0.f)
    // 	q_y_p15 *= in2[up];
    // else
    // 	q_y_p15 *= in2[uint2(id.x, min(id.y + 2, gridSize-1))];
    // float q_y_p1 = 0.5f * (q_y_p05 + q_y_p15);
    // // Calculate corresponding values for u_x_(i,j) using upwinding
    // // (why do we use upwinding here instead of averaging like q?)
    // float u_star_y_0 = 0.f;
    // if (q_y_0 >= 0.f)
    // 	u_star_y_0 = in0[down];
    // else
    // 	u_star_y_0 = in0[curr];
    // float u_star_y_p1 = 0.f;
    // if (q_y_p1 > 0.f)
    // 	u_star_y_p1 = in0[curr];
    // else
    // 	u_star_y_p1 = in0[up];


    // Compute dux_dt and duy_dt
    float dux_dt = - (1/cellSize) * ((q_x_0/h_x_p05) * (in0[curr] - in0[left]) + (q_y_m05/h_x_p05) * (in0[curr] - in0[down]));
    float duy_dt = - (1/cellSize) * ((q_y_0/h_y_p05) * (in1[curr] - in1[down]) + (q_x_m05/h_y_p05) * (in1[curr] - in1[left]));
    // float dux_dt = - (1/(cellSize*h_x_p05)) * ((q_x_p1 * u_star_x_p1 - q_x_0 * u_star_x_0) - in0[curr] * (q_x_p1 - q_x_0));
    // float duy_dt = - (1/(cellSize*h_y_p05)) * ((q_y_p1 * u_star_y_p1 - q_y_0 * u_star_y_0) - in1[curr] * (q_y_p1 - q_y_0));
    dux_dt -= (1/cellSize) * GRAVITY * (in3[right] - in3[curr]);
    duy_dt -= (1/cellSize) * GRAVITY * (in3[up] - in3[curr]);
    float ubarNew_x = LimitVelocity(in0[curr] + timeStep * dux_dt);  // Enforcing CFL condition
    float ubarNew_y = LimitVelocity(in1[curr] + timeStep * duy_dt);  // Enforcing CFL condition
    
    out0[curr] = ubarNew_x;
    out1[curr] = ubarNew_y;
    if (ubarNew_x >= 0.f || id.x == gridSize-1)
        out2[curr] = ubarNew_x * in2[curr];
    else
        out2[curr] = ubarNew_x * in2[right];
    if (ubarNew_y >= 0.f || id.y == gridSize-1)
        out3[curr] = ubarNew_y * in2[curr];
    else
        out3[curr] = ubarNew_y * in2[up];
}


/////////// Transport ///////////

[numthreads(16, 16, 1)]
void UpdateTilde(uint3 id : SV_DispatchThreadID) {
    // Inputs: in0 = ubarNew_x, in1 = ubar_x, in2 = ubarNew_y, in3 = ubar_y, 
    //         in4 = qtildePast_x, in5 = qtildePast_y, in6 = h, in7 = htilde
    // Outputs: out0 = qtilde_x, out1 = qtilde_y, out2 = htilde
    if (id.x < 1 || id.x >= (uint)(gridSize - 1) || id.y < 1 || id.y >= (uint)(gridSize - 1)) return;

    // Adjust qtilde to account for bulk flow
    float ubar_x_avg = 0.5f * (in0[id.xy] + in1[id.xy]); // ubar is on same timestep as h, need to get back to timestep of q 
    float ubar_y_avg = 0.5f * (in2[id.xy] + in3[id.xy]);
    float step_x = - ubar_x_avg * timeStep / cellSize; // unitless (cells)
    float step_y = - ubar_y_avg * timeStep / cellSize;

    // ubar divergence using central difference
    float ubar_x_m1 = 0.5f * (in0[id.xy - uint2(1, 0)] + in1[id.xy - uint2(1, 0)]);
    float ubar_x_p1 = 0.5f * (in0[id.xy + uint2(1, 0)] + in1[id.xy + uint2(1, 0)]);
    float ubar_y_m1 = 0.5f * (in2[id.xy - uint2(0, 1)] + in3[id.xy - uint2(0, 1)]);
    float ubar_y_p1 = 0.5f * (in2[id.xy + uint2(0, 1)] + in3[id.xy + uint2(0, 1)]);
    float div_ubar_x = (ubar_x_p1 - ubar_x_m1) / (2.f * cellSize);
    float div_ubar_y = (ubar_y_p1 - ubar_y_m1) / (2.f * cellSize);
    float div_ubar = div_ubar_x + div_ubar_y;
    if (div_ubar < 0.f) // dampen if converging to avoid breaking waves
        div_ubar *= gammaTransport;

    // Update qtilde using bulk flow and ubar divergence: dq/dt = -q * div(ubar)
    if (((ubar_x_avg >= 0.f) && (in6[id.xy] < minWaterHeight)) ||
        ((ubar_x_avg < 0.f) && (in6[id.xy + uint2(1, 0)] < minWaterHeight)))
        out0[id.xy] = 0.f;
    else
        out0[id.xy] = SampleCubicClamped(in4, id.x + step_x, id.xy, false) * exp(-div_ubar * timeStep);
        // out0[id.xy] = LimitFlowRate(out0[id.xy], in6[id.xy], in6[id.xy + uint2(1, 0)]); 
    if (((ubar_y_avg >= 0.f) && (in6[id.xy] < minWaterHeight)) ||
        ((ubar_y_avg < 0.f) && (in6[id.xy + uint2(0, 1)] < minWaterHeight)))
        out1[id.xy] = 0.f;
    else
        out1[id.xy] = SampleCubicClamped(in5, id.y + step_y, id.xy, true) * exp(-div_ubar * timeStep);
        // out1[id.xy] = LimitFlowRate(out1[id.xy], in6[id.xy], in6[id.xy + uint2(0, 1)]); 

    // Update htilde using ubar divergence
    div_ubar_x = (in0[id.xy] - in0[id.xy - uint2(1, 0)]) / cellSize;
    div_ubar_y = (in2[id.xy] - in2[id.xy - uint2(0, 1)]) / cellSize;
    div_ubar = div_ubar_x + div_ubar_y;
    if (div_ubar < 0.f) // dampen if converging to avoid breaking waves
        div_ubar *= gammaTransport;
    out2[id.xy] = in7[id.xy] * exp(-div_ubar * timeStep);
}

[numthreads(16, 16, 1)]
void CalcQAdvect(uint3 id : SV_DispatchThreadID) {
    // Inputs: in0 = ubarNew_x, in1 = ubarNew_y, in2 = htilde
    // Outputs: out0 = qAdvect_x, out1 = qAdvect_y
    // cubic reconstruction: 1/2 cell accounts for staggered grid, 0.5 * dt accounts for h and u at different 1/2 times
    float step_x = 0.5f - in0[id.xy] * (0.5f * timeStep) / cellSize;
    float step_y = 0.5f - in1[id.xy] * (0.5f * timeStep) / cellSize; 
    out0[id.xy] = in0[id.xy] * SampleCubicClamped(in2, id.x + step_x, id.xy, false);  
    out1[id.xy] = in1[id.xy] * SampleCubicClamped(in2, id.y + step_y, id.xy, true);  
}

[numthreads(16, 16, 1)]
void IntegrateH(uint3 id : SV_DispatchThreadID) {
    // Inputs: in0 = qbar_x, in1 = qtilde_x, in2 = qAdvect_x, 
    //         in3 = qbar_y, in4 = qtilde_y, in5 = qAdvect_y, in6 = hPast, in7 = terrain
    // Outputs: out0 = h, out1 = q_x, out2 = q_y
    if (id.x < 1 || id.x >= (uint)(gridSize - 1) || id.y < 1 || id.y >= (uint)(gridSize - 1)) return;

    uint2 curr = id.xy;
    uint2 curr_xm = id.xy - uint2(1, 0);
    uint2 curr_xp = id.xy + uint2(1, 0);
    uint2 curr_ym = id.xy - uint2(0, 1);
    uint2 curr_yp = id.xy + uint2(0, 1);

    // q = qbar + qtilde + qAdvect
    float q_x = in0[curr] + in1[curr] + in2[curr];
    float q_y = in3[curr] + in4[curr] + in5[curr];
    float q_xm = in0[curr_xm] + in1[curr_xm] + in2[curr_xm];
    float q_ym = in3[curr_ym] + in4[curr_ym] + in5[curr_ym];

    q_x  = LimitFlowRate(q_x, in6[curr], in6[curr_xp]);
    q_y  = LimitFlowRate(q_y, in6[curr], in6[curr_yp]);
    q_xm = LimitFlowRate(q_xm, in6[curr_xm], in6[curr]);
    q_ym = LimitFlowRate(q_ym, in6[curr_ym], in6[curr]);

    if (StopFlowOnTerrainBoundary(in6, in7, curr, false)) // || (id.x == gridSize - 2) // Volume Preserving (wall boundary)
        q_x = 0.f;
    if (StopFlowOnTerrainBoundary(in6, in7, curr, true )) // || (id.y == gridSize - 2) // Volume Preserving (wall boundary)
        q_y = 0.f;
    if (StopFlowOnTerrainBoundary(in6, in7, curr_xm, false))
        q_xm = 0.f;
    if (StopFlowOnTerrainBoundary(in6, in7, curr_ym, true )) 
        q_ym = 0.f;

    float div_q = (q_x - q_xm + q_y - q_ym) / cellSize;
	out0[curr] = max(0.f, in6[curr] - timeStep * div_q);
    out1[curr] = q_x;// - in2[curr]; // qbar + qtilde, removing qAdvect
    out2[curr] = q_y;// - in5[curr];
}


/////////// eWave ///////////
// (At bottom because requires different inputs/outputs for each kernel)
// These all have slight variants on teh names so HLSL doesn't get mad for redefining a name used before

// Average htilde in time to get on same timestep as q, then prep variables for FFT
RWTexture2D<float2> hHat  : register(u1); // Complex
RWTexture2D<float2> qHat_x: register(u2); // Complex
RWTexture2D<float2> qHat_y: register(u3); // Complex
[numthreads(16, 16, 1)]
void TransferToFFT(uint3 id : SV_DispatchThreadID) {
    // Inputs: in0 = htilde, in1 = qtilde_x, in2 = qtilde_y
    // Outputs: out0 = htildeOld, out1 = hHat, out2 = qHat_x, out3 = qHat_y
    float h_real = 0.5 * (in0[id.xy] + out0[id.xy]);
    out0  [id.xy] = in0[id.xy];
    hHat  [id.xy] = float2(h_real, 0.0f);
    qHat_x[id.xy] = float2(in1[id.xy], 0.0f);
    qHat_y[id.xy] = float2(in2[id.xy], 0.0f);
}

// FFT all three (see fft.hlsl)

// Update hHat to be at the right sample time and calculate wavenumber
RWTexture2D<float2> hhat   : register(u0);
RWTexture2D<float2> wavenum: register(u1);
[numthreads(16, 16, 1)]
void CalcHHat(uint3 id : SV_DispatchThreadID) {
    // Inputs: none (hhat updated in place)
    // Outputs: out0 = hHat, out1 = wavenum

    // Calculate k (wave number) from grid position
    float pos_x = id.x * 2 / gridSize; // range [0, 2), normalized by gridsize/2
    float pos_y = id.y * 2 / gridSize;
    float kx = 1 - abs(1 - pos_x);
    float ky = 1 - abs(1 - pos_y);
    float k = sqrt(kx*kx + ky*ky) * PI / cellSize;
    float kS = k;  // signed k
    if (id.x > (float)(gridSize) / 2.f)
        kS = -k;

    // Fourier gradient: dhhat/dx = hhat * -ik
    float real = hhat[id.xy].x;
    float imag = hhat[id.xy].y;
    hhat[id.xy].x = -kS * imag;
    hhat[id.xy].y =  kS * real;

    // Phase shift to translate hhat to cell boundaries: multiply hhat * e^-i*shift
    real = hhat[id.xy].x;
    imag = hhat[id.xy].y;
    float shift = 0.5 * cellSize * kS;
    hhat[id.xy] = float2(cos(shift) * real - sin(shift) * imag, // complex multiplication
                         sin(shift) * real + cos(shift) * imag);
    wavenum[id.xy] = k;
}

// eWave update: calculate resulting qHat at varying depths
Texture2D<float>  waveNum: register(t0);
Texture2D<float2> htHat  : register(t1);
Texture2D<float2> qtHat_x: register(t2);
Texture2D<float2> qtHat_y: register(t3);  
RWTexture2DArray<float2> qhat_x_array: register(u0);
RWTexture2DArray<float2> qhat_y_array: register(u1);
[numthreads(16, 16, 1)]
void CalcQHat(uint3 id : SV_DispatchThreadID) {
    // Inputs: in0 = wavenum, in1 = hHat, in2 = qHat_x, in3 = qHat_y
    // Outputs: out0 = qHat_x_array, qHat_y_array
    float k = waveNum[id.xy];
    float kNonZero = max(0.001, k);
    float beta = sqrt(2.0 / (kNonZero * cellSize) * sin(k * cellSize / 2.0)); // from 2D, I think this is wrong? 
    // float beta = sqrt(2.0 * k / cellSize * sin(k * cellSize / 2.0)); // grid dispersion correction
    float omega = sqrt(GRAVITY * k * tanh(k * depth[id.z])) / beta;
    float sin_term = (omega / (kNonZero * kNonZero)) * sin(omega * timeStep);
    float cos_term = cos(omega * timeStep);
    qhat_x_array[id] = float2(qtHat_x[id.xy].x * cos_term - htHat[id.xy].x * sin_term,
                              qtHat_x[id.xy].y * cos_term - htHat[id.xy].y * sin_term);
    qhat_y_array[id] = float2(qtHat_y[id.xy].x * cos_term - htHat[id.xy].x * sin_term,
                              qtHat_y[id.xy].y * cos_term - htHat[id.xy].y * sin_term);
}

// Inverse FFT the resulting qtilde arrays

// // Convert complex output to real
Texture2DArray<float2> qHat_x_array : register(t0);
Texture2DArray<float2> qHat_y_array : register(t1);
RWTexture2DArray<float> qtilde_x_array: register(u0);
RWTexture2DArray<float> qtilde_y_array: register(u1);
[numthreads(16, 16, 1)]
void CopyFromFFT(uint3 id : SV_DispatchThreadID) {
    qtilde_x_array[id] = qHat_x_array[id].x;
    qtilde_y_array[id] = qHat_y_array[id].x;
}

// Interpolate qtilde between depths
Texture2D<float>       hbar        : register(t0);
Texture2DArray<float2> Qhat_x_array: register(t1);
Texture2DArray<float2> Qhat_y_array: register(t2);
RWTexture2D<float>     qtilde_x    : register(u0);
RWTexture2D<float>     qtilde_y    : register(u1);
[numthreads(16, 16, 1)]
void InterpQ(uint3 id : SV_DispatchThreadID) {
    // Inputs: qtilde_x_array, qtilde_y_array, hbar
    // Outputs: qtilde_x, qtilde_y
    if (id.x < 0 || id.x >= (uint)(gridSize - 1) || id.y < 0 || id.y >= (uint)(gridSize - 1)) return;

    float waterDepth_x = max(hbar[id.xy], hbar[id.xy + uint2(1, 0)]);
    float waterDepth_y = max(hbar[id.xy], hbar[id.xy + uint2(0, 1)]);
    int d1_x = 0;
    int d1_y = 0;
    for (int d = 0; d < depthNum; d++)
        if (waterDepth_x >= depth[d])
            d1_x = d;
        if (waterDepth_y >= depth[d])
            d1_y = d;
    int d2_x = min(depthNum - 1, d1_x + 1);
    int d2_y = min(depthNum - 1, d1_y + 1);
    float sx = 0.f;
    float sy = 0.f;
    if (d1_x != d2_x)
        sx = (depth[d2_x] - waterDepth_x) / (depth[d2_x] - depth[d1_x]);
    if (d1_y != d2_y)
        sy = (depth[d2_y] - waterDepth_y) / (depth[d2_y] - depth[d1_y]);

    qtilde_x[id.xy] = sx * Qhat_x_array[uint3(id.x, id.y, d1_x)].x + (1.f - sx) * Qhat_x_array[uint3(id.x, id.y, d2_x)].x;
    qtilde_y[id.xy] = sy * Qhat_y_array[uint3(id.x, id.y, d1_y)].x + (1.f - sy) * Qhat_y_array[uint3(id.x, id.y, d2_y)].x;
}