/*
Implementation of FFT waves. Based on Cougar's FFTWave implementation in HoloOcean
and the paper "Empirical Directional Wave Spectra for Computer Graphics". 
 */


cbuffer Constants : register(b0) {
    // Sim params
    int gridSize; 
    float cellSize;
    float timeStep;
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
    float surfaceTension;
    float density;

    // Padding for 16-byte alignment 
    // float buffer;
};

cbuffer FFTWaveConstants : register(b1) {
    float depth; 
    float windSpeed;
    float windAngle;
    float windDirection[2];
    float period;
    float amplitude;
    float windDirectionality;
    float windTighten;
    float choppiness;

    float surfaceTension;
    float density;
}

#define G 9.80665
#define PI 3.14159265358979323846f

// Texture Registers
Texture2D<float>   in0 : register(t0);
Texture2D<float>   in1 : register(t1);
Texture2D<float>   in2 : register(t2);
RWTexture2D<float4> out0: register(u0);
RWTexture2D<float> out1: register(u1);
RWTexture2D<float> out2: register(u2);

float2 Omega(k) {
    float kh = k * depth;
    float tanh_kh = (kh > 20.f) ? 1.f : tanh(kh);
    float sech_kh = (kh > 20.f) ? 0.f : sech(kh);
    float capillaryTerm = pow(k, 3) * surfaceTension / density;
    float omega = sqrt((G * k + capillaryTerm) * tanh_kh);
    float dwdk = (depth * (G * k + capillaryTerm) * sech_kh * sech_kh + omega * omega) / (2 * omega);
    return float2(omega, dwdk);
}

[numthreads(16, 16, 1)]
void PopulateSpectrum(uint3 id : SV_DispatchThreadID) {
    // Inputs: in0 = h, in1 = q_x, in2 = q_y, in3 = terrain
    // Outputs: out0 = H, out1 = Q_x, out2 = Q_y
    if (id.x >= (uint)(gridSize) || id.y >= (uint)(gridSize)) return;

    if (id.x == 0 && id.y == 0) {
        out0[id.xy] = float4(0.f, 0.f, 0.f, 0.f);
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

    ///// Frequency & Wave Angle /////
    float2 omegaData = Omega(k); //(omega, dwdk)
    float omega = omegaData.x; 
    float theta = atan2(ky_, kx_) - windAngle;
    theta = (theta + PI) % (2 * PI) - PI; // wrap to [-pi, pi]      // theta = atan2(sin(theta), cos(theta)); 

    ///// Get Spectrum /////
    float S = Spectrum(omega, theta);
    // Convert S(w,theta) to S(kx, ky)
    float dwdk = omegaData.y;
    S *= dwdk / k;

    ///// Convert to amplitude & phase /////


    // // Dampen waves outside user-defined threshold
    // S *= exp(-(k * k) * cutoff[0]);
    // S *= k < cutoff[1] ? 0 : 1;

    // // Dampen waves not aligned with wind
    // float2 windDir = float2(windDirection[0], windDirection[1]);
    // float windAlignment = kx_ * windDirection[0] + ky_ * windDirection[1];
    // float windFactor = pow(abs(windAlignment), windTighten);
    // float2 windFactor2 = float2(windFactor * ( windAlignment > 0 ? 1 : (1 - windDirectionality)),
    //                             windFactor * (-windAlignment > 0 ? 1 : (1 - windDirectionality)));
    // float2 result = float2(S * windFactor2[0], S * windFactor2[1]);
    // result = sqrt(result / 2.0f);

    // // Generate random amplitude and phase
    // float idx = y * gridSize + x;
    // float r1 = random(idx) * 2 * PI;
    // float r2 = random(idx) * 2 * PI;
    // float r3 = random(idx) * 2 * PI;
    // float r4 = random(idx) * 2 * PI;
    // float4 randomMagnitude = float4(r1, r1, r2, r2);
    // float4 randomPhase = float4 (sin(r3), cos(r3), sin(r4), cos(r4));

    // // Get complex amplitude & phase: (pos.real, pos.imag, neg.real, neg.imag)
    // float4 posAndNegSpectrum = float4(result.x, result.x, result.y, result.y);
    // posAndNegSpectrum *= randomAmplitudes; // scale by random seed
    // posAndNegSpectrum *= float4(1.0f, 1.0f, 1.0f, -1.0f); // get conjugate of negative spectrum by flipping imaginary part
    
    // out0[id.xy] = posAndNegSpectrum;
}

