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

float LimitVelocity(float velocity_in)
{
	if (velocity_in >= 0.f)
		return min(velocity_in, cflCondition * cellSize / timeStep);   // 0.25 since other neighbors might take from this source cell as well
	else
		return max(velocity_in, -cflCondition * cellSize / timeStep);
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
	float xcubicX = -x3 + 2.f * x2 - fx;
	float xcubicY =  3.f * x3 - 5.f * x2 + 2.f;
	float xcubicZ = -3.f * x3 + 4.f * x2 + fx;
	float xcubicW =  x3 - x2;

    float result; 
    if (isYDirection) {
        float data0 = dataField[uint2(curr.x, id0)];
        float data1 = dataField[uint2(curr.x, id1)];
        float data2 = dataField[uint2(curr.x, id2)];
        float data3 = dataField[uint2(curr.x, id3)];
        result = 0.5f * (xcubicX * data0 + xcubicY * data1 + xcubicZ * data2 + xcubicW * data3);
        result = min(max(data1, data2), result);  // value-limiting
	    result = max(min(data1, data2), result);  // value-limiting
    }
    else {
        float data0 = dataField[uint2(id0, curr.y)];
        float data1 = dataField[uint2(id1, curr.y)];
        float data2 = dataField[uint2(id2, curr.y)];
        float data3 = dataField[uint2(id3, curr.y)];
        result = 0.5f * (xcubicX * data0 + xcubicY * data1 + xcubicZ * data2 + xcubicW * data3);
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
        out0[id.xy] = 0.f;
        out1[id.xy] = 0.f;
        out2[id.xy] = 0.f;
        // out3[id.xy] = 0.f;
    } 
}

/////////// Decomposition ///////////

// Decomposition: Initialize 
[numthreads(16, 16, 1)]
void InitDecomp(uint3 id : SV_DispatchThreadID) {
    // Inputs: in0 = h, in1 = q_x, in2 = q_y, in3 = terrain
    // Outputs: out0 = H, out1 = Q_x, out2 = Q_y
    out0[id.xy] = in0[id.xy] + in3[id.xy];
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
    float hbar = max(0.f, in0[id.xy] - in6[id.xy]);
    out0[id.xy] = hbar;
    out1[id.xy] = in1[id.xy];
    out2[id.xy] = in2[id.xy];
	out3[id.xy] = in3[id.xy] - hbar;
    out4[id.xy] = in4[id.xy] - in1[id.xy];
    out5[id.xy] = in5[id.xy] - in2[id.xy];
}

[numthreads(16, 16, 1)]
void StopFlow(uint3 id : SV_DispatchThreadID) {
    // Inputs: in0 = h, in1 = terrain
    // Outputs: out0 = qbar_x, out1 = qbar_y, out2 = qtilde_x, out3 = qtilde_y
    if (id.x >= (uint)(gridSize - 1) || id.y >= (uint)(gridSize - 1)) return;
    if (stopFlowOnTerrainBoundary(in0, in1, id.xy, false)) { // stop flow in x direction
        out0[id.xy] = 0.f;
        out2[id.xy] = 0.f;
    }
    if (stopFlowOnTerrainBoundary(in0, in1, id.xy, true)) { // stop flow in y direction
        out1[id.xy] = 0.f;
        out3[id.xy] = 0.f;
    }
}


/////////// eWave ///////////




/////////// SWE ///////////

[numthreads(16, 16, 1)]
void Ubar(uint3 id : SV_DispatchThreadID) {
    // Inputs: in0 = qbar_x, in1 = qbar_y, in2 = hbarOld
    // Outputs: out0 = ubar_x, out1 = ubar_y
    float ubar_x = in0[id.xy];
    float ubar_y = in1[id.xy];

    // First-Order Up-Winding
    // SOMEDAY: Try interpolating h or using higher-order upwinding for better accuracy?
    // Technically u = q / H, not h?? Different derivations differ here
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
void SWE(uint3 id : SV_DispatchThreadID) {
    // Inputs: in0 = ubar_x, in1 = ubar_y, in2 = hbar, in3 = H
    // Outputs: out0 = ubarNew_x, out1 = ubarNew_y, out2 = qbar_x, out3 = qbar_y

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
    float q_x_m05 = in0[id.xy - uint2(1, 0)];
    if (q_x_m05 >= 0.f)
        q_x_m05 *= in2[id.xy - uint2(1, 0)];
    else
        q_x_m05 *= in2[id.xy];
    float q_x_p05 = in0[id.xy];  
    if (q_x_p05 >= 0.f)
        q_x_p05 *= in2[id.xy];
    else
        q_x_p05 *= in2[id.xy + uint2(1, 0)];
    float q_x_0 = 0.5f * (q_x_m05 + q_x_p05);
    // float q_x_p15 = ubar_x[idx(x+1,y)];  //q_(i+1.5,j) = hfr at position x
    // if (q_x_p15 >= 0.f)
    // 	q_x_p15 *= hbar[idx(x+1,y)];
    // else
    // 	q_x_p15 *= hbar[std::min(idx(x+1,y) + 1, GRIDSIZE*GRIDSIZE-1)];
    // float q_x_p1 = 0.5f * (q_x_p05 + q_x_p15);
    // Calculate corresponding vaules for u_x_(i,j) using upwinding
    // (why do we use upwinding here instead of averaging like q?)
    // float u_star_x_0 = 0.f;
    // if (q_x_0 >= 0.f)
    // 	u_star_x_0 = ubar_x[idx(x-1,y)];
    // else
    // 	u_star_x_0 = ubar_x[idx(x,y)];
    // float u_star_x_p1 = 0.f;
    // if (q_x_p1 > 0.f)
    // 	u_star_x_p1 = ubar_x[idx(x,y)];
    // else
    // 	u_star_x_p1 = ubar_x[idx(x+1,y)];
    // Calculate h_(i+0.5,j) and h_(i-0.5,j) using upwinding
    // float h_avg_x_p05 = (hbar[idx(x,y)] + hbar[idx(x+1,y)]) / 2.f; // averaging, for some reason the old paper used this
    float h_x_p05 = 0.f;
    if (in0[id.xy] >= 0.f)
        h_x_p05 = in2[id.xy];
    else
        h_x_p05 = in2[id.xy + uint2(1, 0)];

    /////// Y DIRECTION //////
    // Use upwinding to evaluate q_(i,j-0.5), q_(i,j+0.5), q_(i,j+1.5) 
    // for y direction to get q_y_(i,j), q_y_(i,j+1)
    float q_y_m05 = in1[id.xy - uint2(0, 1)];
    if (q_y_m05 >= 0.f)
        q_y_m05 *= in2[id.xy - uint2(0, 1)];
    else
        q_y_m05 *= in2[id.xy];
    float q_y_p05 = in1[id.xy];  
    if (q_y_p05 >= 0.f)
        q_y_p05 *= in2[id.xy];
    else
        q_y_p05 *= in2[id.xy + uint2(0, 1)];
    float q_y_0 = 0.5f * (q_y_m05 + q_y_p05);
    // Calculate h_(i,j+0.5) and h_(i,j-0.5) using upwinding
    float h_y_p05 = 0.f;
    if (in1[id.xy] >= 0.f)
        h_y_p05 = in2[id.xy];
    else
        h_y_p05 = in2[id.xy + uint2(0, 1)];

    // Compute dux_dt and duy_dt
    // X DIRECTION
    float gravity = 9.8066f;
    float dux_dt = - (1/cellSize) * ((q_x_0/h_x_p05) * (in0[id.xy] - in0[id.xy - uint2(1, 0)]) + (q_y_m05/h_x_p05) * (in0[id.xy] - in0[id.xy - uint2(0, 1)]) + gravity * (in3[id.xy + uint2(1, 0)] - in3[id.xy]));
    out0[id.xy] = LimitVelocity(in0[id.xy] + timeStep * dux_dt);  // Enforcing CFL condition
    // Y DIRECTION
    float duy_dt = - (1/cellSize) * ((q_y_0/h_y_p05) * (in1[id.xy] - in1[id.xy - uint2(0, 1)]) + (q_x_m05/h_y_p05) * (in1[id.xy] - in1[id.xy - uint2(1, 0)]) + gravity * (in3[id.xy + uint2(0, 1)] - in3[id.xy]));
    out1[id.xy] = LimitVelocity(in1[id.xy] + timeStep * duy_dt);  // Enforcing CFL condition


    if (ubarNew_x[idx(x,y)] >= 0.f || x == GRIDSIZE-1)
        qbar_x[idx(x,y)] = ubarNew_x[idx(x,y)] * hbar[idx(x,y)];
    else
        qbar_x[idx(x,y)] = ubarNew_x[idx(x,y)] * hbar[idx(x+1,y)];
    if (ubarNew_y[idx(x,y)] >= 0.f || y == GRIDSIZE-1)
        qbar_y[idx(x,y)] = ubarNew_y[idx(x,y)] * hbar[idx(x,y)];
    else
        qbar_y[idx(x,y)] = ubarNew_y[idx(x,y)] * hbar[idx(x,y+1)];
}


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
        out0[id.xy] = SampleCubicClamped(in4, id.x + step_x, id.xy, false) * exp(-div_ubar * timeStep);
    if (((bulkVelocity_y >= 0.f) && (in6[id.xy] < minWaterHeight)) ||
        ((bulkVelocity_y < 0.f) && (in6[uint2(id.x, id.y + 1)] < minWaterHeight)))
        out1[id.xy] = 0.f;
    else
        out1[id.xy] = SampleCubicClamped(in5, id.y + step_y, id.xy, true) * exp(-div_ubar * timeStep);

    // Update htilde using ubar divergence
    div_ubar_x = (in0[id.xy] - in0[uint2(id.x-1, id.y)]) / cellSize;
    div_ubar_y = (in2[id.xy] - in2[uint2(id.x, id.y-1)]) / cellSize;
    div_ubar = div_ubar_x + div_ubar_y;
    if (div_ubar < 0.f) // dampen if converging to avoid breaking waves
        div_ubar *= gammaTransport;
    out2[id.xy] *= exp(-div_ubar * timeStep);
}

[numthreads(16, 16, 1)]
void UpdateQAdvect(uint3 id : SV_DispatchThreadID) {

}




