// Shaders executing steps of algorithm outlined in "Generalizing Shallow Water Simulations with Dispersive Surface Waves"

cbuffer Constants : register(b0) {
    float time;
    // Sim Params
    int gridSizeX; 
    int gridSizeY; 
    float cellSize;
    float timeStep;
    float minWaterHeight;
    float surfaceTension;
    float density;
    // Decomposition Params
    int diffusionIterations;
    int maxDiffusionCells;
    float diffusionPenalty;
    // SWE & Transport Params
    float slopeLimit;
    float cflCondition;
    float gammaTransport;
    int spongeThickness;
    float laplacianDamping;
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
    float maxSafeDepth;
    float simConstantPadding[1];
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
Texture2D<float>   in8 : register(t8);
Texture2D<float>   in9 : register(t9);
Texture2D<float>   in10 : register(t10);
Texture2D<float>   in11 : register(t11);
StructuredBuffer<float> depth : register(t12);
RWTexture2D<float> out0: register(u0);
RWTexture2D<float> out1: register(u1);
RWTexture2D<float> out2: register(u2);
RWTexture2D<float> out3: register(u3);
RWTexture2D<float> out4: register(u4);
RWTexture2D<float> out5: register(u5);


//////////////////// HELPER FUNCTIONS /////////////////////////
float LimitVelocity(float velocity_in, float cflFactor) {
	if (velocity_in >= 0.f)
		return min(velocity_in, cflFactor);
	else
		return max(velocity_in, -cflFactor);
}

float LimitFlowRate(float flow_rate_in, float waterDepth_left, float waterDepth_right, float cflFactor) {
	if (flow_rate_in >= 0.f)
		return min(flow_rate_in, cflFactor * waterDepth_left);
	else
		return max(flow_rate_in, -cflFactor * waterDepth_right);
}

float LimitFlowRate(float flow_rate_in, float waterDepth_left, float terrain_left, float waterDepth_right, float terrain_right, float cflFactor) {
	if (flow_rate_in >= 0.f) {
		if (waterDepth_left <= minWaterHeight || (terrain_left + waterDepth_left <= terrain_right))
			return 0.f;
		return min(flow_rate_in, cflFactor * waterDepth_left);
	} else {
		if (waterDepth_right <= minWaterHeight || (terrain_right + waterDepth_right <= terrain_left))
			return 0.f;
		return max(flow_rate_in, -cflFactor * waterDepth_right);
	}
}

float SampleCubicClamped2D(Texture2D<float> dataField, float2 samplePos) {
    // 1. Calculate integer base and fractional coordinates
    int2 id_start = (int2)floor(samplePos) - 1;
    float2 frac = samplePos - floor(samplePos);
    
    // 2. Safely clamp all 4x4 grid indices to the domain boundaries
    int2 maxGrid = int2(gridSizeX - 1, gridSizeY - 1);
    int2 id0 = clamp(id_start,     int2(0,0), maxGrid);
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

    // 4. Load the 16 values explicitly (reusing registers)
    float v00 = dataField[int2(id0.x, id0.y)];
    float v10 = dataField[int2(id1.x, id0.y)];
    float v20 = dataField[int2(id2.x, id0.y)];
    float v30 = dataField[int2(id3.x, id0.y)];

    float v01 = dataField[int2(id0.x, id1.y)];
    float v11 = dataField[int2(id1.x, id1.y)]; // c00
    float v21 = dataField[int2(id2.x, id1.y)]; // c10
    float v31 = dataField[int2(id3.x, id1.y)];

    float v02 = dataField[int2(id0.x, id2.y)];
    float v12 = dataField[int2(id1.x, id2.y)]; // c01
    float v22 = dataField[int2(id2.x, id2.y)]; // c11
    float v32 = dataField[int2(id3.x, id2.y)];

    float v03 = dataField[int2(id0.x, id3.y)];
    float v13 = dataField[int2(id1.x, id3.y)];
    float v23 = dataField[int2(id2.x, id3.y)];
    float v33 = dataField[int2(id3.x, id3.y)];

    // 5. Interpolate along X for each row
    float row0 = wX.x * v00 + wX.y * v10 + wX.z * v20 + wX.w * v30;
    float row1 = wX.x * v01 + wX.y * v11 + wX.z * v21 + wX.w * v31;
    float row2 = wX.x * v02 + wX.y * v12 + wX.z * v22 + wX.w * v32;
    float row3 = wX.x * v03 + wX.y * v13 + wX.z * v23 + wX.w * v33;

    // 6. Interpolate the 4 row results along Y
    float result = wY.x * row0 + wY.y * row1 + wY.z * row2 + wY.w * row3;

    // 7. Monotonicity clamp (reusing the already-loaded values)
    float minVal = min(min(v11, v21), min(v12, v22));
    float maxVal = max(max(v11, v21), max(v12, v22));

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

float SolveJacobi(float u_orig, float sigma, float u_E, float u_W, float u_N, float u_S, float C_E, float C_W, float C_N, float C_S) {
    return (u_orig + sigma * (C_E * u_E + C_W * u_W + C_N * u_N + C_S * u_S)) / (1.0f + sigma * (C_E + C_W + C_N + C_S));
}

float SpongeDamping(uint2 id, bool isYDir = false) {
    // If the cell is in the interior (outside the sponge boundary), return 1.0 immediately
    if (id.x >= (uint)spongeThickness && id.x < (uint)(gridSizeX - spongeThickness) &&
        id.y >= (uint)spongeThickness && id.y < (uint)(gridSizeY - spongeThickness)) {
        return 1.f;
    }

    float x = (float)id.x;
    float y = (float)id.y;
    float distX = min((float)spongeThickness, (x < (float)gridSizeX / 2.f) ? x : (float)gridSizeX - x);
    float distY = min((float)spongeThickness, (y < (float)gridSizeY / 2.f) ? y : (float)gridSizeY - y);
    
    float invSpongeThickness = 1.f / (float)spongeThickness;
    float ratioX = ((float)spongeThickness - distX) * invSpongeThickness;
    float ratioY = ((float)spongeThickness - distY) * invSpongeThickness;
    
    return (1.f - ratioX * ratioX) * (1.f - ratioY * ratioY);
}

//////////////////// COMPUTE SHADERS /////////////////////////

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

    int2 curr  = int2(id.xy);
    int2 maxGrid = int2(gridSizeX - 1, gridSizeY - 1);
    int2 right = clamp(curr + int2(1, 0), int2(0, 0), maxGrid);
    int2 up    = clamp(curr + int2(0, 1), int2(0, 0), maxGrid);

    float h = clamp(in0[curr] - in1[curr], 0.f, maxSafeDepth);
    float hr = clamp(in0[right] - in1[right], 0.f, maxSafeDepth);
    float hu = clamp(in0[up] - in1[up], 0.f, maxSafeDepth);
    hr = (h + hr) * 0.5f;
    hu = (h + hu) * 0.5f;
    float invCellSize = 1.0f / cellSize;
    float grad_x = (in0[right] - in0[curr]) * invCellSize;
    float grad_y = (in0[up] - in0[curr]) * invCellSize;
    float penalty = exp(- diffusionPenalty * (grad_x * grad_x + grad_y * grad_y));
    
    float max_alpha = 0.5 * (maxDiffusionCells * cellSize) * (maxDiffusionCells * cellSize);
    out0[curr]   = clamp(0.5 * h  * h * penalty,  0.0f, max_alpha);
    out1[curr] = clamp(0.5 * hr * hr * penalty, 0.0f, max_alpha);
    out2[curr] = clamp(0.5 * hu * hu * penalty, 0.0f, max_alpha);
}

[numthreads(16, 16, 1)]
void DiffusionStep(uint3 id : SV_DispatchThreadID) {
    // Inputs: in0 = terrain, in1 = HOrig, in2 = QOrig_x, in3 = QOrig_y,
    //         in4 = HPast, in5 = QPast_x, in6 = QPast_y,
    //         in7 = alpha_H, in8 = alpha_Q_x, in9 = alpha_Q_y,
    // Outputs: out0 = H, out1 = Q_x, out2 = Q_y
    if (id.x >= (uint)(gridSizeX) || id.y >= (uint)(gridSizeY)) return;

    int2 curr = int2(id.xy);
    int2 maxGrid = int2(gridSizeX - 1, gridSizeY - 1);
    int2 right = clamp(curr + int2(1, 0), int2(0, 0), maxGrid);
    int2 left  = clamp(curr + int2(-1, 0), int2(0, 0), maxGrid);
    int2 up    = clamp(curr + int2(0, 1), int2(0, 0), maxGrid);
    int2 down  = clamp(curr + int2(0, -1), int2(0, 0), maxGrid);
    int2 ru    = clamp(curr + int2(1, 1), int2(0, 0), maxGrid);
    int2 rd    = clamp(curr + int2(1, -1), int2(0, 0), maxGrid);
    int2 lu    = clamp(curr + int2(-1, 1), int2(0, 0), maxGrid);

    float aH_curr  = in7[curr];
    float aH_right = in7[right];
    float aH_left  = in7[left];
    float aH_up    = in7[up];
    float aH_down  = in7[down];
    float a_topright = 0.25f * (aH_curr + aH_up   + aH_right + in7[ru]);
    float a_botright = 0.25f * (aH_curr + aH_down + aH_right + in7[rd]);
    float a_topleft  = 0.25f * (aH_curr + aH_up   + aH_left  + in7[lu]);

    float invCellSizeSq = 1.0f / (cellSize * cellSize);   
    float newH = SolveJacobi(in1[curr], invCellSizeSq, in4[right], in4[left], in4[up], in4[down], in8[curr], in8[left], in9[curr], in9[down]);
    out0[curr] = max(in0[curr], newH);
    out1[curr] = SolveJacobi(in2[curr], invCellSizeSq, in5[right], in5[left], in5[up], in5[down], aH_right,   aH_curr,   a_topright, a_botright);
    out2[curr] = SolveJacobi(in3[curr], invCellSizeSq, in6[right], in6[left], in6[up], in6[down], a_topright, a_topleft, aH_up,      aH_curr);
}

[numthreads(16, 16, 1)]
void DecomposeFields(uint3 id : SV_DispatchThreadID) {
    // Inputs: in0 = H, in1 = Q_x, in2 = Q_y
    //         in3 = h, in4 = q_x, in5 = q_y, in6 = terrain
    // Outputs: out0 = hbar, out1 = qbar_x, out2 = qbar_y, 
    //          out3 = htilde, out4 = qtilde_x, out5 = qtilde_y
    if (id.x >= (uint)(gridSizeX) || id.y >= (uint)(gridSizeY)) return;
    
    int2 curr = int2(id.xy);
    int2 maxGrid = int2(gridSizeX - 1, gridSizeY - 1);
    
    float t_curr = in6[curr];
    float hbar = clamp(in0[curr] - t_curr, 0.f, maxSafeDepth);
    out0[curr] = hbar;
    out1[curr] = in1[curr];
    out2[curr] = in2[curr];
    out3[curr] = in3[curr] - hbar;
    out4[curr] = in4[curr] - in1[curr];
    out5[curr] = in5[curr] - in2[curr];

    float h_curr = in3[curr];
    int2 right = clamp(curr + int2(1, 0), int2(0, 0), maxGrid);
    float h_right = in3[right];
    float t_right = in6[right];

    float cflFactor = cflCondition * cellSize / timeStep;
    out1[curr] = LimitFlowRate(out1[curr], h_curr, t_curr, clamp(in0[right] - t_right, 0.f, maxSafeDepth), t_right, cflFactor);
    out4[curr] = LimitFlowRate(out4[curr], h_curr, t_curr, h_right, t_right, cflFactor);

    int2 up = clamp(curr + int2(0, 1), int2(0, 0), maxGrid);
    float h_up = in3[up];
    float t_up = in6[up];

    out2[curr] = LimitFlowRate(out2[curr], h_curr, t_curr, clamp(in0[up] - t_up, 0.f, maxSafeDepth), t_up, cflFactor);
    out5[curr] = LimitFlowRate(out5[curr], h_curr, t_curr, h_up, t_up, cflFactor);
}


/////////// SWE ///////////

[numthreads(16, 16, 1)]
void CalcUbar(uint3 id : SV_DispatchThreadID) {
    // Inputs: in0 = qbar_x, in1 = qbar_y, in2 = hbarOld
    // Outputs: out0 = ubar_x, out1 = ubar_y
    if (id.x >= (uint)(gridSizeX) || id.y >= (uint)(gridSizeY)) return;

    int2 curr = int2(id.xy);
    int2 maxGrid = int2(gridSizeX - 1, gridSizeY - 1);
    int2 right = clamp(curr + int2(1, 0), int2(0, 0), maxGrid);
    int2 up    = clamp(curr + int2(0, 1), int2(0, 0), maxGrid);

    // First-Order Up-Winding with regularized division near zero depth
    float hx = (in0[curr] >= 0.f || curr.x == gridSizeX-1) ? in2[curr] : in2[right];
    float hy = (in1[curr] >= 0.f || curr.y == gridSizeY-1) ? in2[curr] : in2[up];
    hx = min(hx, maxSafeDepth);
    hy = min(hy, maxSafeDepth);
    float ubar_x = (hx <= minWaterHeight) ? 0.f : (in0[curr] * hx / (hx * hx + minWaterHeight));
    float ubar_y = (hy <= minWaterHeight) ? 0.f : (in1[curr] * hy / (hy * hy + minWaterHeight));

    // Enforcing CFL condition for later surface waves advection
    float cflFactor = cflCondition * cellSize / timeStep;
    out0[curr] = LimitVelocity(ubar_x, cflFactor);  
    out1[curr] = LimitVelocity(ubar_y, cflFactor); 
}

[numthreads(16, 16, 1)]
void CalcSWE(uint3 id : SV_DispatchThreadID) {
    // Inputs: in0 = ubar_x, in1 = ubar_y, in2 = hbar, in3 = H, in4 = delHfft_x, delHfft_y
    // Outputs: out0 = ubarNew_x, out1 = ubarNew_y, out2 = qbar_x, out3 = qbar_y
    if (id.x >= (uint)(gridSizeX) || id.y >= (uint)(gridSizeY)) return;

    int2 curr = int2(id.xy);
    int2 maxGrid = int2(gridSizeX - 1, gridSizeY - 1);
    int2 right = clamp(curr + int2(1, 0), int2(0, 0), maxGrid);
    int2 left  = clamp(curr + int2(-1, 0), int2(0, 0), maxGrid);
    int2 up    = clamp(curr + int2(0, 1), int2(0, 0), maxGrid);
    int2 down  = clamp(curr + int2(0, -1), int2(0, 0), maxGrid);
    int2 r2    = clamp(curr + int2(2, 0), int2(0, 0), maxGrid);
    int2 u2    = clamp(curr + int2(0, 2), int2(0, 0), maxGrid);

    // Use upwinding to evaluate q_(i-0.5,j), q_(i+0.5,j), q_(i+1.5,j) 
    // for x direction to get q_x_(i,j), q_x_(i+1,j)
    float q_x_m05 = in0[left] * ((in0[left] >= 0.f) ? in2[left] : in2[curr]);
    float q_x_p05 = in0[curr] * ((in0[curr] >= 0.f) ? in2[curr] : in2[right]);
    float q_x_p15 = in0[right]* ((in0[right]>= 0.f) ? in2[right]: in2[r2]);
    float q_x_0 = 0.5f * (q_x_m05 + q_x_p05);
    float q_x_1 = 0.5f * (q_x_p05 + q_x_p15);
    float q_y_m05 = in1[down] * ((in1[down] >= 0.f) ? in2[down] : in2[curr]);
    float q_y_p05 = in1[curr] * ((in1[curr] >= 0.f) ? in2[curr] : in2[up]);
    float q_y_p15 = in1[up]   * ((in1[up]   >= 0.f) ? in2[up]   : in2[u2]);
    float q_y_0 = 0.5f * (q_y_m05 + q_y_p05);
    float q_y_1 = 0.5f * (q_y_p05 + q_y_p15);
    
    // Calculate h_(i+0.5,j)
    float h_x_p05 = clamp((in2[curr] + in2[right]) * 0.5f, 0.f, maxSafeDepth);  
    float h_y_p05 = clamp((in2[curr] + in2[up]) * 0.5f, 0.f, maxSafeDepth);

    // Calculate corresponding values for u_x_(i,j) using upwinding
    float u_x_0 = (q_x_0 >= 0.f) ? in0[left] : in0[curr];
    float u_x_1 = (q_x_1 >= 0.f) ? in0[curr] : in0[right]; 
    float u_y_0 = (q_y_0 >= 0.f) ? in1[down] : in1[curr];
    float u_y_1 = (q_y_1 >= 0.f) ? in1[curr] : in1[up];       

    // Compute dux_dt and duy_dt 
    float invCellSize = 1.0f / cellSize;
    float invH_x = h_x_p05 / (h_x_p05 * h_x_p05 + minWaterHeight); // adding minWaterHeight to avoid division by zero
    float invH_y = h_y_p05 / (h_y_p05 * h_y_p05 + minWaterHeight);
    float dux_dt = (h_x_p05 <= minWaterHeight) ? 0.f : (- (invCellSize * invH_x) * ((q_x_1 * u_x_1 - q_x_0 * u_x_0) - in0[curr] * (q_x_1 - q_x_0)));
    float duy_dt = (h_y_p05 <= minWaterHeight) ? 0.f : (- (invCellSize * invH_y) * ((q_y_1 * u_y_1 - q_y_0 * u_y_0) - in1[curr] * (q_y_1 - q_y_0)));
    
    // Calculate wave gradient
    float gradh_x = (in3[right] - in3[curr]) * invCellSize;
    float gradh_y = (in3[up] - in3[curr]) * invCellSize;

    // Calculate FFT wave forcing
    float invDepthCutoff = 1.0f / depthCutoff;
    float depthWeight = SafeTanh(in2[curr] * invDepthCutoff); // scaling term to reduce FFT waves in shallow water
    gradh_x += depthWeight * in4[curr]; // FFT wave pressure gradient 
    gradh_y += depthWeight * in5[curr];

    // Limit steep waves: When wave gets too steep, it "crashes"
    gradh_x = clamp(gradh_x, -slopeLimit, slopeLimit);
    gradh_y = clamp(gradh_y, -slopeLimit, slopeLimit);

    // Incorporate gravity force and bottom friction
    float bottomFriction = 0.05f; // subtle bottom drag to allow water to settle to rest
    dux_dt -= GRAVITY * gradh_x + bottomFriction * in0[curr];
    duy_dt -= GRAVITY * gradh_y + bottomFriction * in1[curr];

    // Integrate u, calculate q
    float cflFactor = cflCondition * cellSize / timeStep;
    float ubarNew_x = LimitVelocity(in0[curr] + timeStep * dux_dt, cflFactor);
    float ubarNew_y = LimitVelocity(in1[curr] + timeStep * duy_dt, cflFactor);

    float h_x_face = (ubarNew_x >= 0.f) ? in2[curr] : in2[right];
    float h_y_face = (ubarNew_y >= 0.f) ? in2[curr] : in2[up];
    if (h_x_face <= minWaterHeight) ubarNew_x = 0.f;
    if (h_y_face <= minWaterHeight) ubarNew_y = 0.f;

    out0[curr] = ubarNew_x;
    out1[curr] = ubarNew_y;

    float t_curr = in3[curr] - in2[curr];
    float t_right = in3[right] - in2[right];
    float t_up = in3[up] - in2[up];

    out2[curr] = LimitFlowRate(ubarNew_x * h_x_face, in2[curr], t_curr, in2[right], t_right, cflFactor);
    out3[curr] = LimitFlowRate(ubarNew_y * h_y_face, in2[curr], t_curr, in2[up], t_up, cflFactor);
}


/////////// Transport & Combine ///////////

[numthreads(16, 16, 1)]
void UpdateTilde(uint3 id : SV_DispatchThreadID) {
    // Inputs: in0 = ubarNew_x, in1 = ubar_x, in2 = ubarNew_y, in3 = ubar_y, 
    //         in4 = qtildePast_x, in5 = qtildePast_y, in6 = h, in7 = htilde
    // Outputs: out0 = htilde, out1 = qtilde_x, out2 = qtilde_y
    if (id.x >= (uint)(gridSizeX) || id.y >= (uint)(gridSizeY)) return;

    int2 curr = int2(id.xy);
    int2 maxGrid = int2(gridSizeX - 1, gridSizeY - 1);
    int2 right = clamp(curr + int2(1, 0), int2(0, 0), maxGrid);
    int2 left  = clamp(curr + int2(-1, 0), int2(0, 0), maxGrid);
    int2 up    = clamp(curr + int2(0, 1), int2(0, 0), maxGrid);
    int2 down  = clamp(curr + int2(0, -1), int2(0, 0), maxGrid);
    int2 rdown = clamp(curr + int2(1, -1), int2(0, 0), maxGrid);
    int2 uleft = clamp(curr + int2(-1, 1), int2(0, 0), maxGrid);

    // integrate at midpoint of timestep 
    float ubar_x_avg = 0.5f * (in0[curr] + in1[curr]); 
    float ubar_y_avg = 0.5f * (in2[curr] + in3[curr]);
    float ubar_x_left  = 0.5f * (in0[left] + in1[left]);
    float ubar_x_right  = 0.5f * (in0[right] + in1[right]);
    float ubar_y_right  = 0.5f * (in2[right] + in3[right]);
    float ubar_y_down  = 0.5f * (in2[down] + in3[down]);
    float ubar_x_up  = 0.5f * (in0[up] + in1[up]);
    float ubar_y_up  = 0.5f * (in2[up] + in3[up]);
    float ubar_y_rdown = 0.5f * (in2[rdown] + in3[rdown]);
    float ubar_x_uleft = 0.5f * (in0[uleft] + in1[uleft]);

    // ubar divergence
    float invCellSize = 1.0f / cellSize;
    float div_ubar  = (ubar_x_avg   - ubar_x_left  + ubar_y_avg   - ubar_y_down)  * invCellSize;  // at cell center
    float div_right = (ubar_x_right - ubar_x_avg   + ubar_y_right - ubar_y_rdown) * invCellSize; // at center of right cell
    float div_up    = (ubar_x_up    - ubar_x_uleft + ubar_y_up    - ubar_y_avg)   * invCellSize;  // at center of up cell
    float div_ux = 0.5f * (div_ubar + div_right); // at right boundary
    float div_uy = 0.5f * (div_ubar + div_up); // at up boundary

    // Clamp divergence to prevent exponential growth during flow convergence over shallow bathymetry
    float maxDivergence = 1.0f / max(1e-4f, timeStep);
    div_ubar = clamp(div_ubar, -maxDivergence, maxDivergence);
    div_ux   = clamp(div_ux,   -maxDivergence, maxDivergence);
    div_uy   = clamp(div_uy,   -maxDivergence, maxDivergence);

    div_ubar *= (div_ubar < 0.f) ? gammaTransport : 1; // Dampen if converging to avoid breaking waves
    div_ux   *= (div_ux   < 0.f) ? gammaTransport : 1;
    div_uy   *= (div_uy   < 0.f) ? gammaTransport : 1;

    // Bound exponential amplification factor to prevent runaway growth
    float exp_ux = exp(clamp(-div_ux * timeStep, -5.0f, 1.0f));
    float exp_uy = exp(clamp(-div_uy * timeStep, -5.0f, 1.0f));

    // Update qtilde using bulk flow and ubar divergence: dq/dt = -q * div(ubar)
    float timeStepOverCellSize = timeStep * invCellSize;
    float step_x = - ubar_x_avg * timeStepOverCellSize; // unitless (cells)
    float step_y = - ubar_y_avg * timeStepOverCellSize; // unitless (cells)
    float2 samplePos = float2(id.x + step_x, id.y + step_y);
    out1[curr] = SampleCubicClamped2D(in4, samplePos) * exp_ux;
    out2[curr] = SampleCubicClamped2D(in5, samplePos) * exp_uy;

    // Limit flow to prevent negative water heights and enforce terrain boundaries
    float cflFactor = cflCondition * cellSize / timeStep;
    out1[curr] = LimitFlowRate(out1[curr], in6[curr], in6[right], cflFactor); 
    out2[curr] = LimitFlowRate(out2[curr], in6[curr], in6[up], cflFactor);
    // float t_curr = terrain[curr];
    // float t_right = terrain[right];
    // float t_up = terrain[up];
    // qtildeOut_x[curr] = LimitFlowRate(qtildeOut_x[curr], hIn[curr], t_curr, hIn[right], t_right, cflFactor); 
    // qtildeOut_y[curr] = LimitFlowRate(qtildeOut_y[curr], hIn[curr], t_curr, hIn[up], t_up, cflFactor); 
      
    // Update htilde using ubar divergence (not at middle of timestep)
    div_ubar  = (in0[curr] - in0[left] + in2[curr] - in2[down]) * invCellSize;
    div_ubar  = clamp(div_ubar, -maxDivergence, maxDivergence);
    div_ubar *= (div_ubar < 0.f) ? gammaTransport : 1; // dampen if converging to avoid breaking waves
    float exp_h = exp(clamp(-div_ubar * timeStep, -5.0f, 1.0f));
    out0[curr] = in7[curr] * exp_h;
}

[numthreads(16, 16, 1)]
void CalcQAdvect(uint3 id : SV_DispatchThreadID) {
    // Inputs: in0 = ubarNew_x, in1 = ubarNew_y, in2 = htilde
    // Outputs: out0 = qAdvect_x, out1 = qAdvect_y
    // cubic reconstruction: 1/2 cell accounts for staggered grid, 0.5 * dt accounts for h and u at different 1/2 times
    if (id.x >= (uint)(gridSizeX) || id.y >= (uint)(gridSizeY)) return;
    float halfTimeStepOverCellSize = (0.5f * timeStep) / cellSize;
    float step_x = 0.5f - in0[id.xy] * halfTimeStepOverCellSize;
    float step_y = 0.5f - in1[id.xy] * halfTimeStepOverCellSize;
    float2 samplePos = float2(id.x + step_x, id.y + step_y);
    float h_sample = SampleCubicClamped2D(in2, samplePos);
    out0[id.xy] = in0[id.xy] * h_sample;
    out1[id.xy] = in1[id.xy] * h_sample;
}

[numthreads(16, 16, 1)]
void IntegrateH(uint3 id : SV_DispatchThreadID) {
    // Inputs: in0 = qbar_x, in1 = qtilde_x, in2 = qAdvect_x, 
    //         in3 = qbar_y, in4 = qtilde_y, in5 = qAdvect_y, 
    //         in6 = hPast, in7 = terrain
    // Outputs: out0 = h, out1 = q_x, out2 = q_y
    if (id.x >= (uint)(gridSizeX) || id.y >= (uint)(gridSizeY)) return;

    int2 curr  = int2(id.xy);
    int2 maxGrid = int2(gridSizeX - 1, gridSizeY - 1);
    int2 right = clamp(curr + int2(1, 0), int2(0, 0), maxGrid);
    int2 left  = clamp(curr + int2(-1, 0), int2(0, 0), maxGrid);
    int2 up    = clamp(curr + int2(0, 1), int2(0, 0), maxGrid);
    int2 down  = clamp(curr + int2(0, -1), int2(0, 0), maxGrid);

    // q = qbar + qtilde + qAdvect
    float q_x  = in0[curr] + in1[curr] + in2[curr];
    float q_y  = in3[curr] + in4[curr] + in5[curr];
    float q_xm = in0[left] + in1[left] + in2[left];
    float q_ym = in3[down] + in4[down] + in5[down];

    // sponge layer: damp q near edges to absorb waves, simulate an open boundary
    float damping = SpongeDamping(id.xy);
    q_x  *= damping;
    q_xm *= damping;
    q_y  *= damping;
    q_ym *= damping;

    float h_curr = in6[curr];
    float t_curr = in7[curr];
    float h_right = in6[right];
    float t_right = in7[right];
    float h_left = in6[left];
    float t_left = in7[left];
    float h_up = in6[up];
    float t_up = in7[up];
    float h_down = in6[down];
    float t_down = in7[down];

    float cflFactor = cflCondition * cellSize / timeStep;
    q_x  = LimitFlowRate(q_x,  h_curr, t_curr, h_right, t_right, cflFactor);
    q_y  = LimitFlowRate(q_y,  h_curr, t_curr, h_up,    t_up,    cflFactor);
    q_xm = LimitFlowRate(q_xm, h_left, t_left, h_curr,  t_curr,  cflFactor);
    q_ym = LimitFlowRate(q_ym, h_down, t_down, h_curr,  t_curr,  cflFactor);

    float invCellSize = 1.0f / cellSize;
    float div_q = (q_x - q_xm + q_y - q_ym) * invCellSize;

    // Wetting-Aware Laplacian to reduce spikes and unstable grid-scale ripples
    bool isFullyWet = (h_curr > minWaterHeight && h_left > minWaterHeight && 
                    h_right > minWaterHeight && h_up > minWaterHeight && h_down > minWaterHeight);
    if (isFullyWet) {
        float laplacian_h = (h_right + h_left + h_up + h_down - 4.0f * h_curr);
        out0[curr] = clamp(h_curr - timeStep * div_q + laplacianDamping * laplacian_h, 0.f, maxSafeDepth);
    } else {
        out0[curr] = clamp(h_curr - timeStep * div_q, 0.f, maxSafeDepth);
    }

    out1[curr] = q_x - in2[curr]; // qbar + qtilde, removing qAdvect
    out2[curr] = q_y - in5[curr];
    out1[curr] = LimitFlowRate(out1[curr], h_curr, t_curr, h_right, t_right, cflFactor);
    out2[curr] = LimitFlowRate(out2[curr], h_curr, t_curr, h_up,    t_up,    cflFactor);
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
    float k2 = kx * kx + ky * ky;
    if (k2 < 1e-12) {
        qhat_x_array[id] = qhat_x[id.xy];
        qhat_y_array[id] = qhat_y[id.xy];
        return;
    }
    float invK = rsqrt(k2);
    float k = k2 * invK;
    // Unit vectors
    float kx_ = kx * invK;
    float ky_ = ky * invK;
    float kx2 = kx_ * kx_;
    float ky2 = ky_ * ky_;

    /////// Dispersion ///////
    // numerical dispersion correction
    float beta = sqrt((2.0 / (k * cellSize)) * sin(k * cellSize / 2.0)); // From their 1D code, but different from paper
    // float beta = sqrt((2.0 * k / cellSize) * sin(k * cellSize / 2.0)); // correct formula from paper, but not stable?
    // Angular frequency for dispersion relation
    float omega = sqrt(GRAVITY * k * SafeTanh(k * min(depth[id.z], maxSafeDepth))) / beta;
    float S = sin(omega * timeStep) * omega / k2;
    float C = cos(omega * timeStep);
    float Cx = C * kx2 + ky2;
    float Cy = C * ky2 + kx2;
    float Ck = (C-1) * kx_ * ky_;

    /////// Gradient of Shifted hHat ///////
    // Fourier gradient along wave direction: dh/dx = hhat * (i * k)
    float2 dh = ComplexMul(hhat[id.xy], float2(0, k));
    // Phase shift in X: shift by dx/2 by multiplying by e^(i * shiftX) = cos(shiftX) + i*sin(shiftX)
    float theta_x = -0.5f * cellSize * kx;
    float theta_y = -0.5f * cellSize * ky;
    float2 shift_x = float2(cos(theta_x), sin(theta_x));
    float2 shift_y = float2(cos(theta_y), sin(theta_y));
    float2 dhdx = ComplexMul(dh, shift_x);
    float2 dhdy = ComplexMul(dh, shift_y);

    // Shift the q cross-terms to align with their target faces
    theta_x = 0.5f * cellSize * (kx - ky);
    theta_y = 0.5f * cellSize * (ky - kx);
    shift_x = float2(cos(theta_x), sin(theta_x));
    shift_y = float2(cos(theta_y), sin(theta_y));
    float2 qx_shifted = ComplexMul(qhat_x[id.xy], shift_x); // now on top face
    float2 qy_shifted = ComplexMul(qhat_y[id.xy], shift_y); // now on right face
    // float2 qx_shifted = qhat_x[id.xy];
    // float2 qy_shifted = qhat_y[id.xy];

    /////// Update Q ///////
    // Shift to cell centers, do the following, shift back to faces
    // 1) Decompose q into parallel and perpendicular: q_|| = kx_*qx + ky_*qy, q_T = kx_*qy - ky_*qx
    // 2) Update in rotated basis: q_||+ = C*q_|| - S*dhdx (q_T unchanged bc Airy is irrotational)
    // 3) Recombine: qx = kx_*q_||+ - ky_*q_T, qy = ky_*q_||+ + kx_*q_T
    qhat_x_array[id] = Cx * qhat_x[id.xy] + Ck * qy_shifted - kx_ * S * dhdx;
    qhat_y_array[id] = Cy * qhat_y[id.xy] + Ck * qx_shifted - ky_ * S * dhdy;
    // qhat_x_array[id] = C * qhat_x[id.xy] - S * dhdx; // Naive 1D translation, but does just as well?
    // qhat_y_array[id] = C * qhat_y[id.xy] - S * dhdy;

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

    float waterDepth_x = min(max(hbar[id.xy], hbar[right]), maxSafeDepth);
    float waterDepth_y = min(max(hbar[id.xy], hbar[up]), maxSafeDepth);
    int d1_x = 0;
    int d1_y = 0;
    for (int d = 0; d < depthNum; d++) {
        if (waterDepth_x >= depth[d]) d1_x = d;
        if (waterDepth_y >= depth[d]) d1_y = d;
    }
    int d2_x = min(depthNum - 1, d1_x + 1);
    int d2_y = min(depthNum - 1, d1_y + 1);
    float sx = 0.f;
    float sy = 0.f;

    // Depth Interpolation (handle shallowest case separately)
    if (waterDepth_x < depth[0]) {
        sx = waterDepth_x / depth[0];
        qtilde_x[id.xy] = sx * qHat_x_array[uint3(id.x, id.y, 0)].x;
    } else {
        if (d1_x != d2_x)
            sx = (depth[d2_x] - waterDepth_x) / (depth[d2_x] - depth[d1_x]);
        qtilde_x[id.xy] = sx * qHat_x_array[uint3(id.x, id.y, d1_x)].x + (1.f - sx) * qHat_x_array[uint3(id.x, id.y, d2_x)].x;
    }
    if (waterDepth_y < depth[0]) {
        sy = waterDepth_y / depth[0];
        qtilde_y[id.xy] = sy * qHat_y_array[uint3(id.x, id.y, 0)].x;
    } else {
        if (d1_y != d2_y)
            sy = (depth[d2_y] - waterDepth_y) / (depth[d2_y] - depth[d1_y]);
        qtilde_y[id.xy] = sy * qHat_y_array[uint3(id.x, id.y, d1_y)].x + (1.f - sy) * qHat_y_array[uint3(id.x, id.y, d2_y)].x;
    }
}