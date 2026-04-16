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
RWTexture2D<float> terrain   : register(u0);
RWTexture2D<float> H         : register(u1);
RWTexture2D<float> Q_x       : register(u2);
RWTexture2D<float> Q_y       : register(u3);
RWTexture2D<float> HPast     : register(u4);
RWTexture2D<float> QPast_x   : register(u5);
RWTexture2D<float> QPast_y   : register(u6);
RWTexture2D<float> alpha_H   : register(u7);
RWTexture2D<float> alpha_Q_x : register(u8);
RWTexture2D<float> alpha_Q_y : register(u9);

// Helper Functions
float calcGradient(RWTexture2D<float> field, RWTexture2D<float> alpha, uint2 pos) {
    float curr = field[pos];
    float left = field[pos - uint2(1,0)];
    float right = field[pos + uint2(1,0)];
    float up = field[pos + uint2(0,1)];
    float down = field[pos - uint2(0,1)];
    float dF_x = alpha[pos] * (right - curr) - left * (curr - left);
    float dF_y = alpha[pos] * (up - curr) - down * (curr - down);
    float dFdT = (dF_x + dF_y) / (cellSize * cellSize);
    return dFdT;
}

// Calculate Diffusion Coefficients 
[numthreads(16, 16, 1)]
void CalcDiffusionCoeffs(uint3 id : SV_DispatchThreadID) {
    if (id.x >= (uint)(gridSize - 1) || id.y >= (uint)(gridSize - 1)) return;

    uint2 curr = id.xy;
    uint2 right = uint2(id.x + 1, id.y);
    uint2 up = uint2(id.x, id.y + 1);

    float max_alpha = (cellSize * cellSize) / (4.0f * deltaT);
    float denom = 2.0f * deltaT * diffusionIterations;
    
    // Alpha_H
    float grad_x = (H[right] - H[curr]) / cellSize;
    float grad_y = (H[up] - H[curr]) / cellSize;
    float penalty = -diffusionPenalty * (grad_x * grad_x + grad_y * grad_y);
    alpha_H[curr] = min(max_alpha, H[curr] * H[curr] / denom) * exp(penalty);

    // Alpha_Q
    float avg_H_x = 0.5f * (H[curr] + H[right]);
    float avg_H_y = 0.5f * (H[curr] + H[up]);
    alpha_Q_x[curr] = min(max_alpha, avg_H_x * avg_H_x / denom) * exp(penalty);
    alpha_Q_y[curr] = min(max_alpha, avg_H_y * avg_H_y / denom) * exp(penalty);
}

// Diffusion Iteration
[numthreads(16, 16, 1)]
void DiffusionStep(uint3 id : SV_DispatchThreadID) {
    if (id.x < 1 || id.x >= (uint)(gridSize - 1) || id.y < 1 || id.y >= (uint)(gridSize - 1)) return;

    uint2 curr = id.xy;
    uint2 left = uint2(id.x - 1, id.y);
    uint2 right = uint2(id.x + 1, id.y);
    uint2 down = uint2(id.x, id.y - 1);
    uint2 up = uint2(id.x, id.y + 1);

    float dH_x = alpha_H[curr] * (HPast[right] - HPast[curr]) - alpha_H[left] * (HPast[curr] - HPast[left]);
    float dH_y = alpha_H[curr] * (HPast[up] - HPast[curr]) - alpha_H[down] * (HPast[curr] - HPast[down]);
    float dHdT = (dH_x + dH_y) / (cellSize * cellSize);
    
    float newH = HPast[curr] + deltaT * calcGradient(HPast, alpha_H, id.xy);
    H[curr] = max(terrain[curr], newH);
    Q_x[curr] = QPast_x[curr] + deltaT * calcGradient(QPast_x, alpha_Q_x, id.xy);
    Q_y[curr] = QPast_y[curr] + deltaT * calcGradient(QPast_y, alpha_Q_y, id.xy);
}

[numthreads(16, 16, 1)]
void ApplyBoundary(uint3 id : SV_DispatchThreadID) {

}

