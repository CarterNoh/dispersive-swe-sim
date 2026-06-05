/*
Implementation of FFT waves. Based on Cougar's FFTWave implementation in HoloOcean
and the paper "Empirical Directional Wave Spectra for Computer Graphics". 
 */


cbuffer Constants : register(b0) {
    float time;

    // Sim params
    int gridSize; 
    float cellSize;
    float timeStep;

    // eWave params
    int depthNum;
    float surfaceTension;
    float density;

    // Padding for 16-byte alignment 
    float buffer;

    float depth; 
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

    float buffer2;
};

#define G 9.80665f
#define PI 3.14159265358979323846f


///////////////// Random Number Generation /////////////////
uint Seed(float kx, float ky, uint seed) {
    static const uint p1 = 73856093u;
    static const uint p2 = 19349663u;
    static const uint p3 = 83492791u;
    uint hashX = (asuint(kx) ^ p1) * p2;
    uint hashY = (asuint(ky) ^ p2) * p3;
    uint hashTime = (seed ^ p3) * p1;
    return hashX ^ hashY ^ hashTime;
}

float Rand(inout uint state) {
    // Convert a uint to a uniform float in the range (0.0, 1.0]
    uint word = ((state >> ((state >> 28u) + 4u)) ^ state) * 277803737u;
    state = state * 747796405u + 2891336453u; // Advance state for the next call
    uint hash = (word >> 22u) ^ word;
    return max(float(hash) / 4294967295.0f, 1e-7f); 
}

float2 Rand2(inout uint state) {
    return float2(Rand(state), Rand(state));
}

float2 Randn2(inout uint state) {
    // Generate two standard normal floats (mean 0, variance 1) using Box-Muller
    float2 unif = Rand2(state); // generate two random uniform numbers
    float r = sqrt(-2.0f * log(unif.x));
    float theta = 2.0f * PI * unif.y;
    return float2(r * cos(theta), r * sin(theta));
}

///////////////// Math Functions /////////////////
float SafeTanh(float x) {
    if (x > 20.f) {
        return 1.f;
    } else if (x < -20.f) {
        return -1.f;
    } else {
        return tanh(x);
    }
}

float Gamma(float z) {
    float x = z;
    float reflectOffset = 0.0;
    bool reflected = false;

    static const float GAMMA_SQRT_2PI = 2.5066282746310002;
    // Lanczos coefficients (g=7, n=9)
    static const float GAMMA_LANCZOS_G = 7.0;
    static const float GAMMA_C[9] = {
        0.99999999999980993,
        676.5203681218851,
        -1259.1392167224028,
        771.32342877765313,
        -176.61502916214059,
        12.507343278686905,
        -0.13857109526572012,
        9.9843695780195716e-6,
        1.5056327351493116e-7
    };

    [branch]
    if (x < 0.5) {
        reflected = true;
        x = 1.0 - x;
    }

    x -= 1.0;

    float a = GAMMA_C[0]
            + GAMMA_C[1] / (x + 1.0)
            + GAMMA_C[2] / (x + 2.0)
            + GAMMA_C[3] / (x + 3.0)
            + GAMMA_C[4] / (x + 4.0)
            + GAMMA_C[5] / (x + 5.0)
            + GAMMA_C[6] / (x + 6.0)
            + GAMMA_C[7] / (x + 7.0)
            + GAMMA_C[8] / (x + 8.0);

    float t = x + GAMMA_LANCZOS_G + 0.5;

    float logG = 0.5 * log(2.0 * PI)
               + (x + 0.5) * log(t) - t
               + log(a);

    [branch]
    if (reflected) {
        // log(Gamma(z)) = log(pi) - log(|sin(pi*z)|) - log(Gamma(1-z))
        logG = log(PI) - log(abs(sin(PI * z))) - logG;
    }

    return exp(logG);
}

float Smoothstep(float x) {
    // Cubic function that is nearly linear between 0 and 1, but 
    // smoothly transitions with a slope of zero at each end
    if (x <= 0.f) return 0.f;
    if (x >= 1.f) return 1.f;
    return x * x * (3.f - 2.f * x);
}

float LonguetHiggins(float theta, float s) {
    float Q = (pow(2, 2 * s - 1) / PI) * pow(Gamma(s + 1), 2) / Gamma(2 * s + 1);
    return Q * pow(abs(cos(theta / 2.f)), 2.f * s);
}

float2 ComplexMul(float2 a, float2 b) {
    return float2(a.x * b.x - a.y * b.y,
                  a.x * b.y + a.y * b.x);
}

///////////////// Spectrum Functions /////////////////
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
    float tanh_kh = SafeTanh(kh);
    float term = G * k + pow(k, 3) * surfaceTension / density;
    float omega = sqrt((term) * tanh_kh);
    float dwdk = term * (depth * pow(1.f / cosh(kh), 2) + tanh_kh) / (2 * omega);

    return float2(omega, dwdk);
}

float2 WaveSpectra(float w) {
    // Pierson-Moskowitz
    float a = 8.1 * pow(10, -3);
    float b = 0.74;
    float wp = G / (1.026 * windSpeed);
    float S = (a * G * G / pow(w, 5)) * exp(-b * pow((wp / w), 4));
    
    // // JONSWAP
    // float F = fetch * 1000; // convert km to m
    // float FBar = (G * F) / pow(windSpeed, 2); // dimensionless fetch
    // float wp = 7 * PI * pow(abs(FBar), 1.f/3.f);
    // float a = 0.076f * pow(abs(FBar), -0.22);
    // float gamma = 3.3f;
    // float s = (w > wp) ? 0.09f : 0.07f;
    // float r = exp(- pow(w - wp, 2) / (2 * pow(s * wp, 2)));
    // float b = 5.f / 4.f;
    // float S = (a * G * G / pow(w, 5)) * exp(-b * pow((wp / w), 4)) * pow(gamma, r);

    // // Texel MARSEN ARSLOE (TMA)
    // float wh = w * sqrt(depth / G);
    // float z = 1.8f * (wh - 1.125f);
    // float tanh_z = (z > 20) ? 1 : (z < -20) ? -1 : tanh(z);
    // float Phi = 0.5f + 0.5f * tanh_z; // tanh approx. of kitaigorodskii depth attenuation
    // S *= Phi;
    
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
    
    return D;
}

float Swell(float theta, float w, float wp) {
    ///// Swell ///// 
    float safeSwell = clamp(swell, 0, 1);
    float shape = 16.1f * SafeTanh(wp / w) * pow(safeSwell, 2);
    return LonguetHiggins(theta, shape);
}

float TotalDirectional(float th, float th_s, float w, float wp) {
    return DirectionalSpreading(th, w, wp) * Swell(th_s, w, wp);
}

float IntegrateDirectional(float a, float b, int n, float w, float wp) {
    // Directional Spreading
    // numerically integrate the product of D * D_swell from -pi to pi to normalize
    float nf = (float)n;
    float step = (b - a) / nf;
    float sum = 0;
    for (int i = 0; i < n; i++) {
        float thetaI = a + i * step;
        sum += TotalDirectional(thetaI, thetaI, w, wp);
    }
    return step * sum;
    // float avg = (TotalDirectional(a, a, w, wp) + TotalDirectional(b, b, w, wp)) / 2.f;
    // return step * (avg + sum);
}

float GetSpectrum(float w, float theta, float theta_swell) {
    // Wave Spectra
    float2 spectraData = WaveSpectra(w);
    float S = spectraData.x;
    float wp = spectraData.y;

    // Directional Spreading & Swell
    float D = TotalDirectional(theta, theta_swell, w, wp) / IntegrateDirectional(-PI, PI, 36, w, wp);

    return S * D;
}

float Filter(float k, int filterInvert) {
    // Smooth invertible band-pass filter
    float wavelength = 2 * PI / k; // convert wavenumber to wavelength
    float fracSmall = (wavelength - (filterSmall - filterWidth)) / filterWidth;
    float fracBig = (wavelength - filterBig) / filterWidth;
    float t = Smoothstep(fracSmall) - Smoothstep(fracBig);
    float f = clamp(filterMin + (1.f - filterMin) * t, 0.f, 1.f);
    if (filterInvert == 1) f = 1.f - f;
    return f;
}


///////////////// Shaders /////////////////
RWTexture2D<float>  omegaOut: register(u0);
RWTexture2D<float2> HPosOut: register(u1);
RWTexture2D<float2> HNegOut: register(u2);
[numthreads(16, 16, 1)]
void PopulateSpectrum(uint3 id : SV_DispatchThreadID) {
    // Outputs: omegaOut = omega, HPosOut = HPos, HNegOut = HNeg (complex amplitudes)
    if (id.x >= (uint)(gridSize) || id.y >= (uint)(gridSize)) return;

    if (id.x == 0 && id.y == 0) {
        omegaOut[id.xy] = 0.f;
        HPosOut[id.xy] = float2(0.f, 0.f);
        HNegOut[id.xy] = float2(0.f, 0.f);
        return;
    }

    ///// Wave Number ///////
    // Calculate the physical size of the grid and the frequency step (dK)
    float domainSize = (float)gridSize * cellSize;
    float dK = 2.0f * PI / domainSize; 
    float N_2 = (float)gridSize / 2.0f;
    // Calculate the physical 2D wavenumber vector components (signed, handling the Nyquist wrap-around)
    int freqX = ((float)id.x < N_2) ? (int)id.x : (int)id.x - (int)gridSize;
    int freqY = ((float)id.y < N_2) ? (int)id.y : (int)id.y - (int)gridSize;
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

    float thetaPos      = atan2( ky_, -kx_) - windAngle;
    float thetaNeg      = atan2(-ky_,  kx_) - windAngle;
    float thetaSwellPos = atan2( ky_, -kx_) - swellAngle;
    float thetaSwellNeg = atan2(-ky_,  kx_) - swellAngle;
    thetaPos = (thetaPos + PI) % (2 * PI) - PI; // wrap to [-pi, pi]      // theta = atan2(sin(theta), cos(theta)); 
    thetaNeg = (thetaNeg + PI) % (2 * PI) - PI;
    thetaSwellPos = (thetaSwellPos + PI) % (2 * PI) - PI;
    thetaSwellNeg = (thetaSwellNeg + PI) % (2 * PI) - PI;

    ///// Spectrum /////
    float SPos = GetSpectrum(omega, thetaPos, thetaSwellPos);
    float SNeg = GetSpectrum(omega, thetaNeg, thetaSwellNeg);
    SPos *= dwdk / k; // Convert S(w,theta) to S(kx, ky)
    SNeg *= dwdk / k;

    ///// Convert to Amplitude & Phase /////
    uint seed = Seed(kx, ky, 0);
    float2 unif = Rand2(seed);
    float2 norm = Randn2(seed); 
    float ampPos = norm.x * dK * sqrt(2 * SPos);
    float ampNeg = norm.y * dK * sqrt(2 * SNeg);
    float phasePos = unif.x * 2 * PI;
    float phaseNeg = unif.y * 2 * PI;

    // Filter amplitudes of wavelengths outside of thresholds
    float filter = Filter(k, 0);
    ampPos *= filter;
    ampNeg *= filter;

    // Scale by gridSize^2, since iFFT will normalize by this
    ampPos *= gridSize * gridSize;
    ampNeg *= gridSize * gridSize;

    // Store results
    omegaOut[id.xy] = omega;
    HPosOut[id.xy] = ampPos * float2(cos(phasePos), -sin(phasePos));
    HNegOut[id.xy] = ampNeg * float2(cos(phaseNeg), -sin(phaseNeg));
}

Texture2D<float>    omegaIn : register(t0);
Texture2D<float2>   HPosIn : register(t1);
Texture2D<float2>   HNegIn : register(t2);
RWTexture2D<float2> HPropOut: register(u0);
RWTexture2D<float2> DxPropOut: register(u1);
RWTexture2D<float2> DyPropOut: register(u2);
[numthreads(16, 16, 1)]
void PropagateWaves(uint3 id : SV_DispatchThreadID) {
    // Inputs: omegaIn = omega, HPosIn = HPos, HNegIn = HNeg
    // Outputs: HPropOut = HProp, DxPropOut = DxProp, DyPropOut = DyProp

    // Calculate Wavevector
    float domainSize = (float)gridSize * cellSize;
    float dK = 2.0f * PI / domainSize; 
    float N_2 = (float)gridSize / 2.0f;
    int freqX = ((float)id.x < N_2) ? (int)id.x : (int)id.x - (int)gridSize;
    int freqY = ((float)id.y < N_2) ? (int)id.y : (int)id.y - (int)gridSize;
    float kx = (float)freqX * dK; // spatial frequency: radians per meter
    float ky = (float)freqY * dK;
    float k = sqrt(kx * kx + ky * ky);
    float kx_ = kx / k;
    float ky_ = ky / k;

    // Propagate Height H
    float w = omegaIn[id.xy];
    // w = floor() * w0; //modify w to be multiple of w0?
    float S = sin(w * time);
    float C = cos(w * time);
    float2 fwd = float2(C, -S); // why is forward the negative one? Try switching if doesn't look right
    float2 bkwd = float2(C, S);
    float2 HProp = ComplexMul(HPosIn[id.xy], fwd) + ComplexMul(HNegIn[id.xy], bkwd);
    HPropOut[id.xy] = HProp;
    
    // Calculate Horizontal Displacement Dx, Dy
    DxPropOut[id.xy] = ComplexMul(HProp, float2(0.f, -kx_ * choppiness));
    DyPropOut[id.xy] = ComplexMul(HProp, float2(0.f, -ky_ * choppiness));
}

Texture2D<float2>  in0 : register(t0);
RWTexture2D<float> out0: register(u0);
[numthreads(16, 16, 1)]
void ComplexToReal(uint3 id : SV_DispatchThreadID) {
    out0[id.xy] = in0[id.xy].x; // real part of complex number
}