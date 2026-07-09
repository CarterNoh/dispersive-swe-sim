#include "DispersiveSWEShaders.h"

IMPLEMENT_GLOBAL_SHADER_PARAMETER_STRUCT(FSimConstants, "SimConstants");

IMPLEMENT_GLOBAL_SHADER(FInitializeWaterCS,     "/Plugin/DispersiveSWESim/kernels.usf",  "InitializeWater",    SF_Compute);

IMPLEMENT_GLOBAL_SHADER(FInitDecompCS,          "/Plugin/DispersiveSWESim/kernels.usf",  "InitDecomp",             SF_Compute);
IMPLEMENT_GLOBAL_SHADER(FRecomputeHCS,          "/Plugin/DispersiveSWESim/kernels.usf",  "RecomputeH",             SF_Compute);
IMPLEMENT_GLOBAL_SHADER(FCalcDiffusionCoeffsCS, "/Plugin/DispersiveSWESim/kernels.usf",  "CalcDiffusionCoeffs",    SF_Compute);
IMPLEMENT_GLOBAL_SHADER(FDiffuseTerrainCS,      "/Plugin/DispersiveSWESim/kernels.usf",  "DiffuseTerrain",         SF_Compute);
IMPLEMENT_GLOBAL_SHADER(FDiffusionStepCS,       "/Plugin/DispersiveSWESim/kernels.usf",  "DiffusionStep",          SF_Compute);
IMPLEMENT_GLOBAL_SHADER(FDecomposeFieldsCS,     "/Plugin/DispersiveSWESim/kernels.usf",  "DecomposeFields",        SF_Compute);
IMPLEMENT_GLOBAL_SHADER(FCalcUbarCS,            "/Plugin/DispersiveSWESim/kernels.usf",  "CalcUbar",               SF_Compute);
IMPLEMENT_GLOBAL_SHADER(FCalcSWECS,             "/Plugin/DispersiveSWESim/kernels.usf",  "CalcSWE",                SF_Compute);
IMPLEMENT_GLOBAL_SHADER(FUpdateTildeCS,         "/Plugin/DispersiveSWESim/kernels.usf",  "UpdateTilde",            SF_Compute);
IMPLEMENT_GLOBAL_SHADER(FCalcQAdvectCS,         "/Plugin/DispersiveSWESim/kernels.usf",  "CalcQAdvect",            SF_Compute);
IMPLEMENT_GLOBAL_SHADER(FIntegrateHCS,          "/Plugin/DispersiveSWESim/kernels.usf",  "IntegrateH",             SF_Compute);
IMPLEMENT_GLOBAL_SHADER(FTransferToFFTCS,       "/Plugin/DispersiveSWESim/kernels.usf",  "TransferToFFT",          SF_Compute);
IMPLEMENT_GLOBAL_SHADER(FCalcEWaveCS,           "/Plugin/DispersiveSWESim/kernels.usf",  "CalcEWave",              SF_Compute);
IMPLEMENT_GLOBAL_SHADER(FInterpQCS,             "/Plugin/DispersiveSWESim/kernels.usf",  "InterpQ",                SF_Compute);
IMPLEMENT_GLOBAL_SHADER(FScaleCopyDisplacementCS,    "/Plugin/DispersiveSWESim/kernels.usf",  "ScaleCopyDisplacement",    SF_Compute);
IMPLEMENT_GLOBAL_SHADER(FCalcSurfaceNormalAndFoamCS, "/Plugin/DispersiveSWESim/kernels.usf",  "CalcSurfaceNormalAndFoam", SF_Compute);
IMPLEMENT_GLOBAL_SHADER(FCalcRoughnessLUTCS,    "/Plugin/DispersiveSWESim/kernels.usf",  "CalcRoughnessLUT",       SF_Compute);

IMPLEMENT_GLOBAL_SHADER(FPopulateSpectrumCS,    "/Plugin/DispersiveSWESim/fftwaves.usf", "PopulateSpectrum",       SF_Compute);
IMPLEMENT_GLOBAL_SHADER(FPropagateWavesCS,      "/Plugin/DispersiveSWESim/fftwaves.usf", "PropagateWaves",         SF_Compute);
IMPLEMENT_GLOBAL_SHADER(FInterpCS,              "/Plugin/DispersiveSWESim/fftwaves.usf", "Interp",                 SF_Compute);

IMPLEMENT_GLOBAL_SHADER(FFFTKernel1DCS,         "/Plugin/DispersiveSWESim/fft.usf",      "FFTKernel_1D",           SF_Compute);



