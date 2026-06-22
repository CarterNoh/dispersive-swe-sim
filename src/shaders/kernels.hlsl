// Shaders executing steps of algorithm outlined in "Generalizing Shallow Water Simulations with Dispersive Surface Waves"

cbuffer Constants : register(b0) {
    float time;
    // Sim Params
    int gridSizeX; 
    int gridSizeY; 
    float cellSize;
    float timeStep;
    int boundaryType;
    float minWaterHeight;
    float surfaceTension;
    float density;
    // Decomposition Params
    int diffusionIterations;
    float deltaT;
    float diffusionPenalty;
    // SWE & Transport Params
    float slopeLimit;
    float cflCondition;
    float gammaTransport;
    // eWave Params
    int depthNum;
    // FFT Wave Params
    float fetch;
    float windSpeed;
    float windAngle;
    float swell;
    float swellAngle;
    float choppiness;
    float filterSmall;
    float filterBig;
    float filterWidth;
    float filterMin;
    float depthCutoff;
    int paddedGridSizeX;
    int paddedGridSizeY;
    float simConstantPadding[3];
};

#define GRAVITY 9.80665
#define PI 3.14159265358979323846f

// Directional helpers using signed int2 relative to SV_DispatchThreadID (id)
#define CURR  (int2(id.xy))
#define RIGHT (int2(id.xy) + int2(1, 0))
#define LEFT  (int2(id.xy) + int2(-1, 0))
#define UP    (int2(id.xy) + int2(0, 1))
#define DOWN  (int2(id.xy) + int2(0, -1))
#define RIGHT2 (int2(id.xy) + int2(2, 0))
#define UP2    (int2(id.xy) + int2(0, 2))
#define RDOWN (int2(id.xy) + int2(1, -1))
#define ULEFT (int2(id.xy) + int2(-1, 1))

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
float GetExtrapolatedValue(Texture2D<float> tex, int2 coord) {
    int tx = coord.x;
    int ty = coord.y;
    
    int x0 = clamp(tx, 0, gridSizeX - 1);
    int y0 = clamp(ty, 0, gridSizeY - 1);
    
    if (boundaryType == 0) { // wall (copy neighbor / clamp coordinate)
        return tex[uint2(x0, y0)];
    }
    // else if (boundaryType == 1) { // free (constant-slope linear extrapolation)
        int x1 = (tx < 0) ? 1 : (tx >= gridSizeX) ? gridSizeX - 2 : x0;
        int y1 = (ty < 0) ? 1 : (ty >= gridSizeY) ? gridSizeY - 2 : y0;
        
        float wx0 = (tx < 0 || tx >= gridSizeX) ? 2.0f : 1.0f;
        float wx1 = (tx < 0 || tx >= gridSizeX) ? -1.0f : 0.0f;
        
        float wy0 = (ty < 0 || ty >= gridSizeY) ? 2.0f : 1.0f;
        float wy1 = (ty < 0 || ty >= gridSizeY) ? -1.0f : 0.0f;
        
        float v00 = tex[uint2(x0, y0)];
        float v10 = tex[uint2(x1, y0)];
        float v01 = tex[uint2(x0, y1)];
        float v11 = tex[uint2(x1, y1)];
    
        return wy0 * (wx0 * v00 + wx1 * v10) + wy1 * (wx0 * v01 + wx1 * v11);
    // }
}

float CalcDiffusion(Texture2D<float> f, Texture2D<float> a, int2 curr) {
    float f_right = GetExtrapolatedValue(f, curr + int2(1, 0));
    float f_left  = GetExtrapolatedValue(f, curr + int2(-1, 0));
    float f_up    = GetExtrapolatedValue(f, curr + int2(0, 1));
    float f_down  = GetExtrapolatedValue(f, curr + int2(0, -1));

    float a_left  = GetExtrapolatedValue(a, curr + int2(-1, 0));
    float a_down  = GetExtrapolatedValue(a, curr + int2(0, -1));

    float dF_x = a[curr] * (f_right - f[curr]) - a_left * (f[curr] - f_left);
    float dF_y = a[curr] * (f_up - f[curr]) - a_down * (f[curr] - f_down);
    float dFdX = (dF_x + dF_y) / (cellSize * cellSize);
    return dFdX;
}

bool StopFlowOnTerrainBoundary(Texture2D<float> h, Texture2D<float> terrain, int2 curr, bool isYDirection=false) {
    int2 check_curr = curr;
    float h_curr;
    float terrain_curr;
    if (boundaryType == 0) { // wall (clamp coordinate)
        check_curr = clamp(curr, int2(0, 0), int2(gridSizeX - 1, gridSizeY - 1));
        h_curr = h[check_curr];
        terrain_curr = terrain[check_curr];
    } else {
        h_curr = GetExtrapolatedValue(h, check_curr);
        terrain_curr = GetExtrapolatedValue(terrain, check_curr);
    }
    
    if (!isYDirection) {
        float h_right = GetExtrapolatedValue(h, check_curr + int2(1, 0));
        float terrain_right = GetExtrapolatedValue(terrain, check_curr + int2(1, 0));
        if ((h_curr <= minWaterHeight) && (terrain_curr >= terrain_right + h_right)) // positive q_x
            return true;
        if ((h_right <= minWaterHeight) && (terrain_right > terrain_curr + h_curr)) // negative q_x
            return true;
        return false;
    }
    else {
        float h_up = GetExtrapolatedValue(h, check_curr + int2(0, 1));
        float terrain_up = GetExtrapolatedValue(terrain, check_curr + int2(0, 1));
        if ((h_curr <= minWaterHeight) && (terrain_curr >= terrain_up + h_up)) // positive q_y
            return true;
        if ((h_up <= minWaterHeight) && (terrain_up > terrain_curr + h_curr)) // negative q_y
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
    int2 maxGrid = int2(gridSizeX - 1, gridSizeY - 1);
    int2 id0 = clamp(id_start + 0, int2(0,0), maxGrid);
    int2 id1 = clamp(id_start + 1, int2(0,0), maxGrid);
    int2 id2 = clamp(id_start + 2, int2(0,0), maxGrid);
    int2 id3 = clamp(id_start + 3, int2(0,0), maxGrid);

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

float SafeTanh(float x) {
    if (x > 20.f)       return 1.f;
    else if (x < -20.f) return -1.f;
    else                return tanh(x);
}

//////////////////// COMPUTE SHADERS /////////////////////////

// Apply Boundary Conditions
[numthreads(16, 16, 1)]
void ApplyBoundaries(uint3 id : SV_DispatchThreadID){
    uint x = id.x;
    uint y = id.y;
    if (x >= (uint)(gridSizeX) || y >= (uint)(gridSizeY)) return;

    bool left   = (x == 0);
    bool right  = (x == (uint)gridSizeX - 1);
    bool bottom = (y == 0);
    bool top    = (y == (uint)gridSizeY - 1);
    if (!(left || right || top || bottom)) return;

    // Pick a source interior cell depending on boundary side
    uint2 src;
    if (left && bottom)
        src = uint2(1, 1);
    else if (left && top)
        src = uint2(1, gridSizeY-2);
    else if (right && bottom)
        src = uint2(gridSizeX-2, 1);
    else if (right && top)
        src = uint2(gridSizeX-2, gridSizeY-2);
    else if (left)
        src = uint2(1, y);
    else if (right)
        src = uint2(gridSizeX-2, y);
    else if (bottom)
        src = uint2(x, 1);
    else // top
        src = uint2(x, gridSizeY-2);

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
    if (id.x >= (uint)(gridSizeX) || id.y >= (uint)(gridSizeY)) return;
    out0[id.xy] = in0[id.xy] + in3[id.xy];
    out1[id.xy] = in1[id.xy];
    out2[id.xy] = in2[id.xy];
}

[numthreads(16, 16, 1)]
void CalcDiffusionCoeffs(uint3 id : SV_DispatchThreadID) {
    // Inputs: in0 = H, in1 = terrain
    // Outputs: out0 = alpha_H, out1 = alpha_Q_x, out2 = alpha_Q_y
    if (id.x >= (uint)(gridSizeX) || id.y >= (uint)(gridSizeY)) return;

    // Use constant-slope linear extrapolation at boundaries
    float h = in0[CURR] - in1[CURR];

    float hr = GetExtrapolatedValue(in0, RIGHT) - GetExtrapolatedValue(in1, RIGHT);
    float hu = GetExtrapolatedValue(in0, UP) - GetExtrapolatedValue(in1, UP);

    float grad_x = (GetExtrapolatedValue(in0, RIGHT) - in0[CURR]) / cellSize;
    float grad_y = (GetExtrapolatedValue(in0, UP) - in0[CURR]) / cellSize;

    hr = (h + hr) / 2;
    hu = (h + hu) / 2;
    float denom = 2.0f * deltaT * diffusionIterations;
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
    if (id.x >= (uint)(gridSizeX) || id.y >= (uint)(gridSizeY)) return;

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
    if (id.x >= (uint)(gridSizeX) || id.y >= (uint)(gridSizeY)) return;
    
    float hbar = max(0.f, in0[CURR] - in6[CURR]);
    out0[CURR] = hbar;
    out1[CURR] = in1[CURR];
    out2[CURR] = in2[CURR];
    out3[CURR] = in3[CURR] - hbar;
    out4[CURR] = in4[CURR] - in1[CURR];
    out5[CURR] = in5[CURR] - in2[CURR];

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
    if (id.x >= (uint)(gridSizeX) || id.y >= (uint)(gridSizeY)) return;

    // First-Order Up-Winding
    float hx = (in0[CURR] >= 0.f) ? in2[CURR] : GetExtrapolatedValue(in2, RIGHT);
    float hy = (in1[CURR] >= 0.f) ? in2[CURR] : GetExtrapolatedValue(in2, UP);
    float ubar_x = in0[CURR] / max(minWaterHeight, hx);
    float ubar_y = in1[CURR] / max(minWaterHeight, hy);

    // Enforcing CFL condition for later surface waves advection
    out0[CURR] = LimitVelocity(ubar_x);  
    out1[CURR] = LimitVelocity(ubar_y); 
}

[numthreads(16, 16, 1)]
void CalcSWE(uint3 id : SV_DispatchThreadID) {
    // Inputs: in0 = ubar_x, in1 = ubar_y, in2 = hbar, in3 = H, in4 = delHfft_x, delHfft_y
    // Outputs: out0 = ubarNew_x, out1 = ubarNew_y, out2 = qbar_x, out3 = qbar_y
    if (id.x >= (uint)(gridSizeX) || id.y >= (uint)(gridSizeY)) return;

    // Extrapolated neighbor values at boundaries using constant-slope linear extrapolation
    float in0_left  = GetExtrapolatedValue(in0, LEFT);
    float in0_right = GetExtrapolatedValue(in0, RIGHT);
    
    float in1_down  = GetExtrapolatedValue(in1, DOWN);
    float in1_up    = GetExtrapolatedValue(in1, UP);

    float in2_left  = GetExtrapolatedValue(in2, LEFT);
    float in2_right = GetExtrapolatedValue(in2, RIGHT);
    float in2_r2    = GetExtrapolatedValue(in2, RIGHT2);
    float in2_down  = GetExtrapolatedValue(in2, DOWN);
    float in2_up    = GetExtrapolatedValue(in2, UP);
    float in2_u2    = GetExtrapolatedValue(in2, UP2);

    float in3_right = GetExtrapolatedValue(in3, RIGHT);
    float in3_up    = GetExtrapolatedValue(in3, UP);

    // Compute intermediate values needed for du/dt calculations
    // Need: q_x/y_ij, q_x_i-0.5_j, q_y_i_j-0.5, 
    // 	     h_i+0.5_j, h_i_j+0.5, h_ij=hbar[x,y], h_i+1_j=hbar[x+1,y], h_i_j+1=hbar[x,y+1],
    //       u_x_i+0.5_j=ubar_x[x,y], u_x_i-0.5_j=ubar_x[x-1,y], 
    //       u_y_i_j+0.5=ubar_y[x,y], u_y_i_j-0.5=ubar_y[x,y-1], 
    //       u_x_i+0.5_j-1=ubar_x[x,y-1], u_y_i-1_j+0.5=ubar_y[x-1,y]

    // Use upwinding to evaluate q_(i-0.5,j), q_(i+0.5,j), q_(i+1.5,j) 
    // for x direction to get q_x_(i,j), q_x_(i+1,j)
    float q_x_m05 = in0_left * ((in0_left >= 0.f) ? in2_left : in2[CURR]);
    float q_x_p05 = in0[CURR] * ((in0[CURR] >= 0.f) ? in2[CURR] : in2_right);
    float q_x_p15 = in0_right* ((in0_right>= 0.f) ? in2_right: in2_r2);
    float q_x_0 = 0.5f * (q_x_m05 + q_x_p05);
    float q_x_1 = 0.5f * (q_x_p05 + q_x_p15);
    
    float q_y_m05 = in1_down * ((in1_down >= 0.f) ? in2_down : in2[CURR]);
    float q_y_p05 = in1[CURR] * ((in1[CURR] >= 0.f) ? in2[CURR] : in2_up);
    float q_y_p15 = in1_up   * ((in1_up   >= 0.f) ? in2_up   : in2_u2);
    float q_y_0 = 0.5f * (q_y_m05 + q_y_p05);
    float q_y_1 = 0.5f * (q_y_p05 + q_y_p15);
    
    // Calculate h_(i+0.5,j)
    // float h_x_p05 = (in0[curr] >= 0.f) ? in2[curr] : in2[right]; // upwinding
    // float h_y_p05 = (in1[curr] >= 0.f) ? in2[curr] : in2[up];
    float h_x_p05 = (in2[CURR] + in2_right) / 2.f;  
    float h_y_p05 = (in2[CURR] + in2_up) / 2.f;

    // Calculate corresponding values for u_x_(i,j) using upwinding
    float u_x_0 = (q_x_0 >= 0.f) ? in0_left : in0[CURR];
    float u_x_1 = (q_x_1 >= 0.f) ? in0[CURR] : in0_right; 
    float u_y_0 = (q_y_0 >= 0.f) ? in1_down : in1[CURR];
    float u_y_1 = (q_y_1 >= 0.f) ? in1[CURR] : in1_up;       

    // Compute dux_dt and duy_dt 
    // float dux_dt = - (1/cellSize) * ((q_x_0/h_x_p05) * (in0[curr] - in0[left]) + (q_y_m05/h_x_p05) * (in0[curr] - in0[down]));
    // float duy_dt = - (1/cellSize) * ((q_y_0/h_y_p05) * (in1[curr] - in1[down]) + (q_x_m05/h_y_p05) * (in1[curr] - in1[left]));
    float dux_dt = - (1/(cellSize * h_x_p05)) * ((q_x_1 * u_x_1 - q_x_0 * u_x_0) - in0[CURR] * (q_x_1 - q_x_0));
    float duy_dt = - (1/(cellSize * h_y_p05)) * ((q_y_1 * u_y_1 - q_y_0 * u_y_0) - in1[CURR] * (q_y_1 - q_y_0));
    
    // Incorporate gravity force and limit steep waves
    float gradh_x = (in3_right - in3[CURR]) / cellSize;
    float gradh_y = (in3_up - in3[CURR]) / cellSize;
    if (abs(gradh_x) > slopeLimit) gradh_x = sign(gradh_x) * slopeLimit; // When wave gets too steep, it "crashes"
    if (abs(gradh_y) > slopeLimit) gradh_y = sign(gradh_y) * slopeLimit; // https://arxiv.org/html/2503.03009v1
    dux_dt -= GRAVITY * gradh_x; // gravitational force
    duy_dt -= GRAVITY * gradh_y;

    // Calculate FFT wave forcing
    float depth_weight = SafeTanh(in2[CURR] / depthCutoff); // scaling term to reduce FFT waves in shallow water
    dux_dt += depth_weight * GRAVITY * in4[CURR]; // FFT wave pressure gradient 
    dux_dt += depth_weight * GRAVITY * in5[CURR];

    // Integrate u, calculate q
    float ubarNew_x = LimitVelocity(in0[CURR] + timeStep * dux_dt);
    float ubarNew_y = LimitVelocity(in1[CURR] + timeStep * duy_dt);
    out0[CURR] = ubarNew_x;
    out1[CURR] = ubarNew_y;
    out2[CURR] = ubarNew_x * ((ubarNew_x >= 0.f) ? in2[CURR] : in2_right);
    out3[CURR] = ubarNew_y * ((ubarNew_y >= 0.f) ? in2[CURR] : in2_up);
}


/////////// Transport & Combine ///////////

[numthreads(16, 16, 1)]
void UpdateTilde(uint3 id : SV_DispatchThreadID) {
    // Inputs: in0 = ubarNew_x, in1 = ubar_x, in2 = ubarNew_y, in3 = ubar_y, 
    //         in4 = qtildePast_x, in5 = qtildePast_y, in6 = h, in7 = htilde
    // Outputs: out0 = htilde, out1 = qtilde_x, out2 = qtilde_y
    if (id.x >= (uint)(gridSizeX) || id.y >= (uint)(gridSizeY)) return;

    // integrate at midpoint of timestep 
    float ubar_x_avg = 0.5f * (in0[CURR] + in1[CURR]); 
    float ubar_y_avg = 0.5f * (in2[CURR] + in3[CURR]);
    float ubar_x_left  = 0.5f * (GetExtrapolatedValue(in0, LEFT) + GetExtrapolatedValue(in1, LEFT));
    float ubar_x_right = 0.5f * (GetExtrapolatedValue(in0, RIGHT) + GetExtrapolatedValue(in1, RIGHT));
    float ubar_y_right = 0.5f * (GetExtrapolatedValue(in2, RIGHT) + GetExtrapolatedValue(in3, RIGHT));
    float ubar_y_down  = 0.5f * (GetExtrapolatedValue(in2, DOWN) + GetExtrapolatedValue(in3, DOWN));
    float ubar_x_up  = 0.5f * (GetExtrapolatedValue(in0, UP) + GetExtrapolatedValue(in1, UP));
    float ubar_y_up  = 0.5f * (GetExtrapolatedValue(in2, UP) + GetExtrapolatedValue(in3, UP));
    float ubar_y_rdown = 0.5f * (GetExtrapolatedValue(in2, RDOWN) + GetExtrapolatedValue(in3, RDOWN));
    float ubar_x_uleft = 0.5f * (GetExtrapolatedValue(in0, ULEFT) + GetExtrapolatedValue(in1, ULEFT));

    // ubar divergence
    float div_ubar  = (ubar_x_avg   - ubar_x_left  + ubar_y_avg   - ubar_y_down)  / cellSize;  // at cell center
    float div_right = (ubar_x_right - ubar_x_avg   + ubar_y_right - ubar_y_rdown) / cellSize; // at center of right cell
    float div_up    = (ubar_x_up    - ubar_x_uleft + ubar_y_up    - ubar_y_avg)   / cellSize;  // at center of up cell
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

    // Limit flow to prevent negative water heights and enforce terrain boundaries
    float h_right = GetExtrapolatedValue(in6, RIGHT);
    float h_up    = GetExtrapolatedValue(in6, UP);
    out1[CURR] = LimitFlowRate(out1[CURR], in6[CURR], h_right); 
    out2[CURR] = LimitFlowRate(out2[CURR], in6[CURR], h_up); 
    if (((ubar_x_avg >= 0.f) && (in6[CURR] <= minWaterHeight)) ||
        ((ubar_x_avg < 0.f)  && (h_right <= minWaterHeight)))
        out1[CURR] = 0.f;
    if (((ubar_y_avg >= 0.f) && (in6[CURR] <= minWaterHeight)) ||
        ((ubar_y_avg < 0.f)  && (h_up <= minWaterHeight)))
        out2[CURR] = 0.f;   
      
    // Update htilde using ubar divergence (not at middle of timestep)
    float in0_left = GetExtrapolatedValue(in0, LEFT);
    float in2_down = GetExtrapolatedValue(in2, DOWN);
    div_ubar  = (in0[CURR] - in0_left + in2[CURR] - in2_down) / cellSize;
    div_ubar *= (div_ubar < 0.f) ? gammaTransport : 1; // dampen if converging to avoid breaking waves
    out0[CURR] = in7[CURR] * exp(-div_ubar * timeStep);
}

[numthreads(16, 16, 1)]
void CalcQAdvect(uint3 id : SV_DispatchThreadID) {
    // Inputs: in0 = ubarNew_x, in1 = ubarNew_y, in2 = htilde
    // Outputs: out0 = qAdvect_x, out1 = qAdvect_y
    // cubic reconstruction: 1/2 cell accounts for staggered grid, 0.5 * dt accounts for h and u at different 1/2 times
    if (id.x >= (uint)(gridSizeX) || id.y >= (uint)(gridSizeY)) return;
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
    if (id.x >= (uint)(gridSizeX) || id.y >= (uint)(gridSizeY)) return;

    // q = qbar + qtilde + qAdvect
    float q_x  = in0[CURR] + in1[CURR] + in2[CURR];
    float q_y  = in3[CURR] + in4[CURR] + in5[CURR];

    float qbar_x_left   = GetExtrapolatedValue(in0, LEFT);
    float qtilde_x_left = GetExtrapolatedValue(in1, LEFT);
    float qadv_x_left   = GetExtrapolatedValue(in2, LEFT);
    
    float qbar_y_down   = GetExtrapolatedValue(in3, DOWN);
    float qtilde_y_down = GetExtrapolatedValue(in4, DOWN);
    float qadv_y_down   = GetExtrapolatedValue(in5, DOWN);
    
    float q_xm = qbar_x_left + qtilde_x_left + qadv_x_left;
    float q_ym = qbar_y_down + qtilde_y_down + qadv_y_down;

    float hPast_right = GetExtrapolatedValue(in6, RIGHT);
    float hPast_up    = GetExtrapolatedValue(in6, UP);
    float hPast_left  = GetExtrapolatedValue(in6, LEFT);
    float hPast_down  = GetExtrapolatedValue(in6, DOWN);

    q_x  = LimitFlowRate(q_x,  in6[CURR], hPast_right);
    q_y  = LimitFlowRate(q_y,  in6[CURR], hPast_up);
    q_xm = LimitFlowRate(q_xm, hPast_left, in6[CURR]);
    q_ym = LimitFlowRate(q_ym, hPast_down, in6[CURR]);

    if (StopFlowOnTerrainBoundary(in6, in7, CURR, false))
        q_x = 0.f;
    if (StopFlowOnTerrainBoundary(in6, in7, CURR, true ))
        q_y = 0.f;
    if (StopFlowOnTerrainBoundary(in6, in7, LEFT, false))
        q_xm = 0.f;
    if (StopFlowOnTerrainBoundary(in6, in7, DOWN, true ))
        q_ym = 0.f;

    float div_q = (q_x - q_xm + q_y - q_ym) / cellSize;
	out0[CURR] = max(0.f, in6[CURR] - timeStep * div_q);
    out1[CURR] = q_x - in2[CURR]; // qbar + qtilde, removing qAdvect
    out2[CURR] = q_y - in5[CURR];
    out1[CURR] = LimitFlowRate(out1[CURR], in6[CURR], hPast_right);
    out2[CURR] = LimitFlowRate(out2[CURR], in6[CURR], hPast_up);
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
    if (id.x >= (uint)(paddedGridSizeX) || id.y >= (uint)(paddedGridSizeY)) return;

    if (id.x >= (uint)(gridSizeX) || id.y >= (uint)(gridSizeY)) {
        hHat[id.xy]   = float2(0.0f, 0.0f);
        qHat_x[id.xy] = float2(0.0f, 0.0f);
        qHat_y[id.xy] = float2(0.0f, 0.0f);
        return;
    }

    // Average htilde in time to get on same timestep as q, then prep variables for FFT
    float h_real = 0.5 * (in0[id.xy] + out0[id.xy]);
    out0[id.xy]   = in0[id.xy];
    hHat[id.xy]   = float2(h_real, 0.0f);
    qHat_x[id.xy] = float2(in1[id.xy], 0.0f);
    qHat_y[id.xy] = float2(in2[id.xy], 0.0f);
}

Texture2D<float2> hhat   : register(t0);
Texture2D<float2> qhat_x : register(t1);
Texture2D<float2> qhat_y : register(t2);
Texture2D<float2> qhatFFT_x : register(t3);
Texture2D<float2> qhatFFT_y : register(t4);
RWTexture2DArray<float2> qhat_x_array: register(u0);
RWTexture2DArray<float2> qhat_y_array: register(u1);
[numthreads(16, 16, 1)]
void CalcEWave(uint3 id : SV_DispatchThreadID) {
    // Inputs: in0 = hHat, in1 = qHat_x, in2 = qHat_y, in3 = depth
    // Outputs: out0 = qHat_x_array, qHat_y_array
    if (id.x >= (uint)(paddedGridSizeX) || id.y >= (uint)(paddedGridSizeY)) return;

    // Early return for DC component 
    if (id.x == 0 && id.y == 0) {
        qhat_x_array[id] = qhat_x[id.xy];
        qhat_y_array[id] = qhat_y[id.xy];
        return;
    }

    ///// Wave Number ///////
    // Calculate the physical size of the grid and the frequency step (dK) in each dimension
    float domainSizeX = (float)paddedGridSizeX * cellSize;
    float domainSizeY = (float)paddedGridSizeY * cellSize;
    float dKx = 2.0f * PI / domainSizeX; 
    float dKy = 2.0f * PI / domainSizeY;
    float NXdiv2 = (float)paddedGridSizeX / 2.0f;
    float NYdiv2 = (float)paddedGridSizeY / 2.0f;
    // Calculate the physical 2D wavenumber vector components (signed, handling the Nyquist wrap-around)
    int freqX = (id.x < NXdiv2) ? (int)id.x : (int)id.x - (int)paddedGridSizeX;
    int freqY = (id.y < NYdiv2) ? (int)id.y : (int)id.y - (int)paddedGridSizeY;
    float kx = (float)freqX * dKx; // spatial frequency: radians per meter
    float ky = (float)freqY * dKy;
    float k = sqrt(kx * kx + ky * ky);
    if (k < 1e-6) {
        qhat_x_array[id] = qhat_x[id.xy];
        qhat_y_array[id] = qhat_y[id.xy];
        return;
    }
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
    float omega = sqrt(GRAVITY * k * SafeTanh(k * in8[id.z])) / beta;
    float S = sin(omega * timeStep) * omega / (k * k);
    float C = cos(omega * timeStep);
    float Cx = 1 + (C-1) * kx2;
    float Cy = 1 + (C-1) * ky2;
    float Ck = (C-1) * kx_ * ky_;

    /////// Gradient of Shifted hHat ///////
    // Fourier gradient in X: dh/dx = hhat * (i * kx)
    float2 dhdx = ComplexMul(hhat[id.xy], float2(0, kx));
    float2 dhdy = ComplexMul(hhat[id.xy], float2(0, ky));
    // Phase shift in X: shift by dx/2 by multiplying by e^(i * shiftX) = cos(shiftX) + i*sin(shiftX)
    float shiftX = 0.5f * cellSize * kx;
    float shiftY = 0.5f * cellSize * ky;
    float2 e_ix = float2(cos(shiftX), sin(shiftX));
    float2 e_iy = float2(cos(shiftY), sin(shiftY));
    dhdx = ComplexMul(dhdx, e_ix);
    dhdy = ComplexMul(dhdy, e_iy);
    
    // Shift the q cross-terms to align with their target faces
    float theta_x = 0.5f * cellSize * (ky - kx);
    float theta_y = 0.5f * cellSize * (kx - ky);
    float2 shift_x = float2(cos(theta_x), sin(theta_x)); // e^{i * theta_x} 
    float2 shift_y = float2(cos(theta_y), sin(theta_y)); // e^{i * theta_y}
    float2 qx_shifted = ComplexMul(qhat_x[id.xy], shift_x);
    float2 qy_shifted = ComplexMul(qhat_y[id.xy], shift_y);
    // float2 qx_shifted = qhat_x[id.xy];
    // float2 qy_shifted = qhat_y[id.xy];

    /////// Update Q ///////
    // 1) Decompose q into parallel and perpendicular: q_|| = kx_*qx + ky_*qy, q_T = kx_*qy - ky_*qx
    // 2) Update in rotated basis: q_||+ = C*q_|| - S*dhdx (q_T unchanged bc Airy is irrotational)
    // 3) Recombine: q = q_T + q_||+
    float qx_r = Cx * qhat_x[id.xy].x + Ck * qy_shifted.x - S * dhdx.x;
    float qx_i = Cx * qhat_x[id.xy].y + Ck * qy_shifted.y - S * dhdx.y;
    float qy_r = Cy * qhat_y[id.xy].x + Ck * qx_shifted.x - S * dhdy.x;
    float qy_i = Cy * qhat_y[id.xy].y + Ck * qx_shifted.y - S * dhdy.y;
    // float qx_r = C * qhat_x[id.xy].x - S * dhdx.x; // Naive 1D translation, but does just as well
    // float qx_i = C * qhat_x[id.xy].y - S * dhdx.y;
    // float qy_r = C * qhat_y[id.xy].x - S * dhdy.x;
    // float qy_i = C * qhat_y[id.xy].y - S * dhdy.y;
    qhat_x_array[id] = float2(qx_r, qx_i);
    qhat_y_array[id] = float2(qy_r, qy_i);

    if ((id.x == NXdiv2 && id.y == 0) || (id.x == 0 && id.y == NYdiv2) || (id.x == NXdiv2 && id.y == NYdiv2)) {
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
    if (id.x >= (uint)(gridSizeX) || id.y >= (uint)(gridSizeY)) return;
    uint2 right = uint2(min(id.x + 1, gridSizeX - 1), id.y);
    uint2 up    = uint2(id.x, min(id.y + 1, gridSizeY - 1));

    float waterDepth_x = max(hbar[id.xy], hbar[right]);
    float waterDepth_y = max(hbar[id.xy], hbar[up]);
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
    qtilde_x[id.xy] = sx * qHat_x_array[uint3(id.x, id.y, d1_x)].x + (1.f - sx) * qHat_x_array[uint3(id.x, id.y, d2_x)].x;
    qtilde_y[id.xy] = sy * qHat_y_array[uint3(id.x, id.y, d1_y)].x + (1.f - sy) * qHat_y_array[uint3(id.x, id.y, d2_y)].x;
}