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
    float fetch;
    float windSpeed;
    float windAngle;
    float windDirectionality;
    float windTighten;
    float choppiness;
    float swell;
    float swellAngle;
    // float period;
    float amplitude;


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

float SafeTanh(x) {
    if (x > 20.f) {
        return 1.f;
    } else if (x < -20.f) {
        return -1.f;
    } else {
        return tanh(x);
    }
}

float Gamma(float s) {
    // Euler Gamma function (complex version of a factorial, or something, idk)
    if s <= 0:
        return 0.f;
    
    // Reflection formula for s < 0.5 to improve accuracy
    if s < 0.5:
        return PI / (sin(PI * s) * Gamma(1 - s))

    // Lanczos coefficients (g=7, n=9)
    // Ultra janky because HLSL doesn't have arrays
    float g = 7;
    s -= 0.5;
    float x = 0.99999999999980993;
    for (int i = 1; i < 9; i++) {
        if (i == 1) float p = 676.5203681218851;
        if (i == 2) float p = -1259.1392167224028;
        if (i == 3) float p = 771.32342877765313;
        if (i == 4) float p = -176.61502916214059;
        if (i == 5) float p = 12.507343278686905;
        if (i == 6) float p = -0.13857109526572012;
        if (i == 7) float p = 9.9843695780195716e-6;
        if (i == 8) float p = 1.5056327351493116e-7;
        x += p / (s + i);
    }
    t = s + g;
    return sqrt(2 * PI) * pow(t, s) * exp(-t) * x;
}

float LonguetHiggins(float theta, float s) {
    float Q = (pow(2, 2 * s - 1) / PI) * pow(Gamma(s + 1), 2) / Gamma(2 * s + 1);
    return Q * pow(abs(cos(theta / 2.f)), 2.f * s);
}

float Integrate(float a, float b, int n) {
    // Directional Spreading
    // numerically integrate the product of D * D_swell from -pi to pi to normalize
    float nf = (float)n;
    float step = (b - a) / nf;
    float sum = 0;
    for (int i = 0; 1 < n; i++) {
        float thetaI = a + i * step;
        sum += TotalDirectional(thetaI, thetaI, w, wp);
    }
    return step * sum;
    // float avg = (TotalDirectional(a, a, w, wp) + TotalDirectional(b, b, w, wp)) / 2.f;
    // return step * (avg + sum);
}

float2 Dispersion(float k) {
    // // Deep Water Dispersion
    // float omega = sqrt(G * k);
    // float dwdk = G / (2 * omega)
    
    // // Depth-Limited Dispersion
    // float kh = k * depth;
    // float tanh_kh = SafeTanh(kh)
    // float omega = sqrt(G * k * tanh_kh);
    // float dwdk = G * (tanh_kh + kh / (cosh(kh) * cosh(kh))) / (2 * omega);

    // Capillary Dispersion
    float kh = k * depth;
    float tanh_kh = (kh > 20.f) ? 1.f : tanh(kh);
    float term = G * k + pow(k, 3) * surfaceTension / density;
    float omega = sqrt((term) * tanh_kh);
    float dwdk = term * (depth * pow(1.f / cosh(kh), 2) + tanh_kh) / (2 * omega);

    return float2(omega, dwdk);
}

float2 WaveSpectra(float w) {
    // // Pierson-Moskowitz
    // float a = 8.1e-3;
    // float b = 0.74;
    // float wp = G / (1.026 * windSpeed);
    // float S = (a * G * G / pow(w, 5)) * exp(-b * pow((wp / w), 4));
    
    // JONSWAP
    float F = fetch * 1000 // convert km to m
    float FBar = (G * F) / pow(windSpeed, 2); // dimensionless fetch
    float wp = 7 * PI * pow(FBar, 1.f/3.f);
    float a = 0.076f * pow(FBar, -0.22);
    float gamma = 3.3f;
    float s = (w > wp) ? 0.09f : 0.07f;
    float r = exp(- pow(w - wp, 2) / (2 * pow(s * wp, 2)));
    float b = 5.f / 4.f;
    float S = (a * G * G / pow(w, 5)) * exp(-b * pow((wp / w), 4)) * pow(gamma, r)

    // Texel MARSEN ARSLOE (TMA)
    float wh = w * sqrt(depth / G);
    float z = 1.8f * (wh - 1.125f);
    float tanh_z = (z > 20) ? 1 : (z < -20) ? -1 : tanh(z);
    float Phi = 0.5f + 0.5f * tanh(z); // tanh approx. of kitaigorodskii depth attenuation
    S *= Phi;
    
    return float2(S, wp);
}

float DirectionalSpreading(float theta, float w, float wp) {
    float D = 0.f;

    // Positive Cosine-Squared
    if ((theta <= PI / 2) || (theta >= -PI / 2)) {
        D = (2.f / PI) * pow(cos(theta), 2);
    }

    // Longuet-Higgins
    // // Mitsuyasu
    // float sp = 11.5 * pow(wp * windSpeed / G, -2.5f);
    // float exponent = (w > wp) ? -2.5f : 5.f;
    // Hasselmann
    float sp = 6.97f;
    float exponent = 4.06f;
    if (w > wp) {
        exponent = -2.33f - 1.45f * ((wp * windSpeed / G) - 1.17f);
        sp = 9.77f;
    }
    float s = sp * pow(w / wp, exponent);
    D = LonguetHiggins(theta, s);

    // Donelan-Banner
    float w_wp = w / wp;
    float eps = -0.4 + 0.8393 * exp(-0.567 * log(pow(w_wp, 2)));
    float beta = pow(10, eps);
    if (w_wp < 0.95) {
        beta = 2.61 * pow(w_wp, 1.3);
    } else if (w_wp < 1.6) {
        beta = 2.28 * pow(w_wp, -1.3);
    }
    D = beta * pow(1.f / cosh(beta * theta), 2) / (2 * SafeTanh(beta * PI));
    
    return D 
}

float Swell(float theta, float w, float wp) {
    ///// Swell ///// 
    float safeSwell = clamp(swell, 0, 1);
    float shape = 16.1f * SafeTanh(wp / w) * pow(safeSwell, 2);
    return = LonguetHiggins(theta, shape)
}

float TotalDirectional(float th, float th_s, float w, float wp) {
    return DirectionalSpreading(th, w, wp) * Swell(th_s, w, wp);
}

float GetSpectrum(float w, float theta, float theta_swell) {
    // Wave Spectra
    float2 spectraData = WaveSpectra(w);
    float S = spectraData.x;
    float wp = spectraData.y;

    // Directional Spreading & Swell
    float D = TotalDirectional(theta, theta_swell, w, wp) / Integrate(-PI, PI, 36);

    return S * D;
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
    float2 omegaData = Dispersion(k); //(omega, dwdk)
    float omega = omegaData.x; 
    float dwdk = omegaData.y;

    float thetaPos      = atan2( ky_,  kx_) - windAngle;
    float thetaNeg      = atan2(-ky_, -kx_) - windAngle;
    float thetaSwellPos = atan2( ky_,  kx_) - swellAngle;
    float thetaSwellNeg = atan2(-ky_, -kx_) - swellAngle;
    thetaPos = (thetaPos + PI) % (2 * PI) - PI; // wrap to [-pi, pi]      // theta = atan2(sin(theta), cos(theta)); 
    thetaNeg = (thetaNeg + PI) % (2 * PI) - PI;
    thetaSwellPos = (thetaSwellPos + PI) % (2 * PI) - PI;
    thetaSwellNeg = (thetaSwellNeg + PI) % (2 * PI) - PI;

    ///// Get Spectrum /////
    float SPos = GetSpectrum(omega, thetaPos, thetaSwellPos);
    float SNeg = GetSpectrum(omega, thetaNeg, thetaSwellNeg);
    // Convert S(w,theta) to S(kx, ky)
    SPos *= dwdk / k;
    SNeg *= dwdk / k;

    ///// Amplitude /////
    float2 ampPos = GetRandomAmp() * sqrt(2 * SPos * dk * dk);
    float2 ampPos = GetRandomAmp() * sqrt(2 * SNeg * dk * dk);
    float2 phasePos = GetRandomPhase();
    float2 phasePos = GetRandomPhase();

    // Filter amplitudes outside of user-defined thresholds
    // // Dampen waves outside user-defined thresholds
    // S *= exp(-(k * k) * cutoff[0]);
    // S *= k < cutoff[1] ? 0 : 1;


    // need a complex number for positive and negative waves, uniformly sampled for phase and gaussian sampled for amplitude


    
    
    // // Dampen waves not aligned with wind
    // float windAlignment = kx_ * cos(windAngle) + ky_ * sin(windAngle);
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

