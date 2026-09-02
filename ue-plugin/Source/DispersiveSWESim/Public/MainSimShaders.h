#pragma once

#include "CoreMinimal.h"
#include "GlobalShader.h"
#include "ShaderParameterStruct.h"
#include "RenderGraphResources.h"
#include "RenderGraphUtils.h"

// Uniform buffer containing constants for SWE simulation kernels
BEGIN_GLOBAL_SHADER_PARAMETER_STRUCT(FSimConstants, )
	SHADER_PARAMETER(int32, gridSizeX)
	SHADER_PARAMETER(int32, gridSizeY)
	SHADER_PARAMETER(float, cellSize)
	SHADER_PARAMETER(float, timeStep)
	SHADER_PARAMETER(float, minWaterHeight)
	SHADER_PARAMETER(float, surfaceTension)
	SHADER_PARAMETER(float, density)
	SHADER_PARAMETER(int32, diffusionIterations)
	SHADER_PARAMETER(int32, maxDiffusionCells)
	SHADER_PARAMETER(float, diffusionPenalty)
	SHADER_PARAMETER(float, slopeLimit)
	SHADER_PARAMETER(float, cflCondition)
	SHADER_PARAMETER(float, gammaTransport)
	SHADER_PARAMETER(int32, spongeThickness)
	SHADER_PARAMETER(float, laplacianDamping)
	SHADER_PARAMETER(int32, depthNum)
	SHADER_PARAMETER(float, depthCutoff)
	SHADER_PARAMETER(int32, paddedGridSizeX)
	SHADER_PARAMETER(int32, paddedGridSizeY)
	SHADER_PARAMETER(float, maxSafeDepth)
	SHADER_PARAMETER_ARRAY(FVector4f, depthLevels, [4])
END_GLOBAL_SHADER_PARAMETER_STRUCT()

// ============================================================================
// Shader Classes from kernels.usf
// ============================================================================

class FInitializeWaterCS : public FGlobalShader
{
public:
	DECLARE_GLOBAL_SHADER(FInitializeWaterCS);
	SHADER_USE_PARAMETER_STRUCT(FInitializeWaterCS, FGlobalShader);

	BEGIN_SHADER_PARAMETER_STRUCT(FParameters, )
		SHADER_PARAMETER_STRUCT_REF(FSimConstants, SimConstants)
		SHADER_PARAMETER(float, WaterLevel)
		SHADER_PARAMETER(float, TerrainCaptureCameraZ)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, terrain)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, hOut)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, H_Out)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, terrainOut)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, terrainOutCM)
	END_SHADER_PARAMETER_STRUCT()
};

class FInitDecompCS : public FGlobalShader
{
public:
	DECLARE_GLOBAL_SHADER(FInitDecompCS);
	SHADER_USE_PARAMETER_STRUCT(FInitDecompCS, FGlobalShader);

	BEGIN_SHADER_PARAMETER_STRUCT(FParameters, )
		SHADER_PARAMETER_STRUCT_REF(FSimConstants, SimConstants)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, hIn)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, qIn_x)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, qIn_y)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, terrain)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, H_Out)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, Q_Out_x)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, Q_Out_y)
	END_SHADER_PARAMETER_STRUCT()
};

class FCalcDiffusionCoeffsCS : public FGlobalShader
{
public:
	DECLARE_GLOBAL_SHADER(FCalcDiffusionCoeffsCS);
	SHADER_USE_PARAMETER_STRUCT(FCalcDiffusionCoeffsCS, FGlobalShader);

	BEGIN_SHADER_PARAMETER_STRUCT(FParameters, )
		SHADER_PARAMETER_STRUCT_REF(FSimConstants, SimConstants)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, H_In)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, terrain)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, alpha_HOut)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, alpha_QOut_x)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, alpha_QOut_y)
	END_SHADER_PARAMETER_STRUCT()
};

class FDiffusionStepCS : public FGlobalShader
{
public:
	DECLARE_GLOBAL_SHADER(FDiffusionStepCS);
	SHADER_USE_PARAMETER_STRUCT(FDiffusionStepCS, FGlobalShader);

	BEGIN_SHADER_PARAMETER_STRUCT(FParameters, )
		SHADER_PARAMETER_STRUCT_REF(FSimConstants, SimConstants)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, terrain)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, H_Orig)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, Q_Orig_x)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, Q_Orig_y)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, H_Past)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, Q_Past_x)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, Q_Past_y)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, alpha_HIn)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, alpha_QIn_x)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, alpha_QIn_y)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, H_Out)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, Q_Out_x)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, Q_Out_y)
	END_SHADER_PARAMETER_STRUCT()
};

class FDecomposeFieldsCS : public FGlobalShader
{
public:
	DECLARE_GLOBAL_SHADER(FDecomposeFieldsCS);
	SHADER_USE_PARAMETER_STRUCT(FDecomposeFieldsCS, FGlobalShader);

	BEGIN_SHADER_PARAMETER_STRUCT(FParameters, )
		SHADER_PARAMETER_STRUCT_REF(FSimConstants, SimConstants)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, H_In)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, Q_In_x)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, Q_In_y)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, hIn)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, qIn_x)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, qIn_y)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, terrain)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, hbarOut)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, qbarOut_x)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, qbarOut_y)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, htildeOut)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, qtildeOut_x)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, qtildeOut_y)
	END_SHADER_PARAMETER_STRUCT()
};

class FRecomputeHCS : public FGlobalShader
{
public:
	DECLARE_GLOBAL_SHADER(FRecomputeHCS);
	SHADER_USE_PARAMETER_STRUCT(FRecomputeHCS, FGlobalShader);

	BEGIN_SHADER_PARAMETER_STRUCT(FParameters, )
		SHADER_PARAMETER_STRUCT_REF(FSimConstants, SimConstants)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, hIn)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, terrain)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, H_Out)
	END_SHADER_PARAMETER_STRUCT()
};

class FRecomputeHAvgCS : public FGlobalShader
{
public:
	DECLARE_GLOBAL_SHADER(FRecomputeHAvgCS);
	SHADER_USE_PARAMETER_STRUCT(FRecomputeHAvgCS, FGlobalShader);

	BEGIN_SHADER_PARAMETER_STRUCT(FParameters, )
		SHADER_PARAMETER_STRUCT_REF(FSimConstants, SimConstants)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, hIn)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, hPast)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, terrain)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, H_Out)
	END_SHADER_PARAMETER_STRUCT()
};

class FTransferToFFTCS : public FGlobalShader
{
public:
	DECLARE_GLOBAL_SHADER(FTransferToFFTCS);
	SHADER_USE_PARAMETER_STRUCT(FTransferToFFTCS, FGlobalShader);

	BEGIN_SHADER_PARAMETER_STRUCT(FParameters, )
		SHADER_PARAMETER_STRUCT_REF(FSimConstants, SimConstants)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, htildeIn)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, htildeOldIn)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, qtildeIn_x)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, qtildeIn_y)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, htildeOldNext)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float2>, hHat)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float2>, qHat_x)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float2>, qHat_y)
	END_SHADER_PARAMETER_STRUCT()
};

class FCalcEWaveCS : public FGlobalShader
{
public:
	DECLARE_GLOBAL_SHADER(FCalcEWaveCS);
	SHADER_USE_PARAMETER_STRUCT(FCalcEWaveCS, FGlobalShader);

	BEGIN_SHADER_PARAMETER_STRUCT(FParameters, )
		SHADER_PARAMETER_STRUCT_REF(FSimConstants, SimConstants)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float2>, hhat)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float2>, qhat_x)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float2>, qhat_y)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2DArray<float2>, Flow_x)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2DArray<float2>, Flow_y)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2DArray<float2>, qhat_x_array)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2DArray<float2>, qhat_y_array)
	END_SHADER_PARAMETER_STRUCT()
};

class FInterpQCS : public FGlobalShader
{
public:
	DECLARE_GLOBAL_SHADER(FInterpQCS);
	SHADER_USE_PARAMETER_STRUCT(FInterpQCS, FGlobalShader);

	BEGIN_SHADER_PARAMETER_STRUCT(FParameters, )
		SHADER_PARAMETER_STRUCT_REF(FSimConstants, SimConstants)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, hbarIn)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2DArray<float2>, qHat_x_array)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2DArray<float2>, qHat_y_array)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, qtildeOut_x)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, qtildeOut_y)
	END_SHADER_PARAMETER_STRUCT()
};

class FCalcUbarCS : public FGlobalShader
{
public:
	DECLARE_GLOBAL_SHADER(FCalcUbarCS);
	SHADER_USE_PARAMETER_STRUCT(FCalcUbarCS, FGlobalShader);

	BEGIN_SHADER_PARAMETER_STRUCT(FParameters, )
		SHADER_PARAMETER_STRUCT_REF(FSimConstants, SimConstants)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, qbarIn_x)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, qbarIn_y)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, hbarIn)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, ubarOut_x)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, ubarOut_y)
	END_SHADER_PARAMETER_STRUCT()
};

class FCalcSWECS : public FGlobalShader
{
public:
	DECLARE_GLOBAL_SHADER(FCalcSWECS);
	SHADER_USE_PARAMETER_STRUCT(FCalcSWECS, FGlobalShader);

	BEGIN_SHADER_PARAMETER_STRUCT(FParameters, )
		SHADER_PARAMETER_STRUCT_REF(FSimConstants, SimConstants)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, ubarIn_x)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, ubarIn_y)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, hbarIn)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, H_In)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, delH_x)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, delH_y)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, terrain)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, ubarNewOut_x)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, ubarNewOut_y)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, qbarOut_x)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, qbarOut_y)
	END_SHADER_PARAMETER_STRUCT()
};

class FUpdateTildeCS : public FGlobalShader
{
public:
	DECLARE_GLOBAL_SHADER(FUpdateTildeCS);
	SHADER_USE_PARAMETER_STRUCT(FUpdateTildeCS, FGlobalShader);

	BEGIN_SHADER_PARAMETER_STRUCT(FParameters, )
		SHADER_PARAMETER_STRUCT_REF(FSimConstants, SimConstants)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, ubarNewIn_x)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, ubarIn_x)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, ubarNewIn_y)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, ubarIn_y)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, qtildePast_x)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, qtildePast_y)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, hIn)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, htildePast)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, htildeOut)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, qtildeOut_x)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, qtildeOut_y)
	END_SHADER_PARAMETER_STRUCT()
};

class FCalcQAdvectCS : public FGlobalShader
{
public:
	DECLARE_GLOBAL_SHADER(FCalcQAdvectCS);
	SHADER_USE_PARAMETER_STRUCT(FCalcQAdvectCS, FGlobalShader);

	BEGIN_SHADER_PARAMETER_STRUCT(FParameters, )
		SHADER_PARAMETER_STRUCT_REF(FSimConstants, SimConstants)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, ubarNewIn_x)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, ubarNewIn_y)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, htildeIn)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, qAdvectOut_x)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, qAdvectOut_y)
	END_SHADER_PARAMETER_STRUCT()
};

class FIntegrateHCS : public FGlobalShader
{
public:
	DECLARE_GLOBAL_SHADER(FIntegrateHCS);
	SHADER_USE_PARAMETER_STRUCT(FIntegrateHCS, FGlobalShader);

	BEGIN_SHADER_PARAMETER_STRUCT(FParameters, )
		SHADER_PARAMETER_STRUCT_REF(FSimConstants, SimConstants)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, qbarIn_x)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, qtildeIn_x)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, qAdvectIn_x)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, qbarIn_y)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, qtildeIn_y)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, qAdvectIn_y)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, hPast)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, terrain)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, qOut_x)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, qOut_y)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, hOut)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, hdot_out)
	END_SHADER_PARAMETER_STRUCT()
};


// ============================================================================
// RDG Helper Structures & Functions
// ============================================================================

// --- Water Initialization ---
struct FInitializeWaterInputs
{
	FRDGTextureRef TerrainInput = nullptr;
	float WaterLevel = 0.0f;
	float TerrainCaptureCameraZ = 0.0f;
	int32 GridSizeX = 0;
	int32 GridSizeY = 0;
};

struct FInitializeWaterOutputs
{
	FRDGTextureRef hOut = nullptr;
	FRDGTextureRef H_Out = nullptr;
	FRDGTextureRef terrainOut = nullptr;
	FRDGTextureRef terrainOutCM = nullptr;
	FRDGTextureRef hbarOut = nullptr;
	FRDGTextureRef hbarOldOut = nullptr;
	FRDGTextureRef FoamOut = nullptr;      // Cleared to 0 if valid
	FRDGTextureRef RoughnessOut = nullptr; // Cleared to 0 if valid
};

void AddInitializeWaterPass(
	FRDGBuilder& GraphBuilder,
	TUniformBufferRef<FSimConstants> ConstantBuffer,
	const FInitializeWaterInputs& Inputs,
	const FInitializeWaterOutputs& Outputs);

// --- Field Decomposition & Jacobi Diffusion ---
struct FDecompositionInputs
{
	FRDGTextureRef hIn = nullptr;
	FRDGTextureRef qIn_x = nullptr;
	FRDGTextureRef qIn_y = nullptr;
	FRDGTextureRef terrain = nullptr;
	FRDGTextureRef HOrig = nullptr;
	FRDGTextureRef QOrig_x = nullptr;
	FRDGTextureRef QOrig_y = nullptr;
	FRDGTextureRef alpha_H = nullptr;
	FRDGTextureRef alpha_Q_x = nullptr;
	FRDGTextureRef alpha_Q_y = nullptr;
	int32 GridSizeX = 0;
	int32 GridSizeY = 0;
	int32 DiffusionIterations = 128;
};

struct FDecompositionOutputs
{
	FRDGTextureRef H_SrcDst = nullptr;      // Ping-pongs with HPast
	FRDGTextureRef HPast_SrcDst = nullptr;
	FRDGTextureRef Qx_SrcDst = nullptr;     // Ping-pongs with QPast_x
	FRDGTextureRef QPast_x_SrcDst = nullptr;
	FRDGTextureRef Qy_SrcDst = nullptr;     // Ping-pongs with QPast_y
	FRDGTextureRef QPast_y_SrcDst = nullptr;
	FRDGTextureRef hbarOut = nullptr;
	FRDGTextureRef qbarOut_x = nullptr;
	FRDGTextureRef qbarOut_y = nullptr;
	FRDGTextureRef htildeOut = nullptr;
	FRDGTextureRef qtildeOut_x = nullptr;
	FRDGTextureRef qtildeOut_y = nullptr;
};

void AddDecompositionPasses(
	FRDGBuilder& GraphBuilder,
	TUniformBufferRef<FSimConstants> ConstantBuffer,
	const FDecompositionInputs& Inputs,
	FDecompositionOutputs& Outputs);

// --- eWave Dispersion Solver ---
struct FEWaveInputs
{
	FRDGTextureRef htildeIn = nullptr;
	FRDGTextureRef htildeOldIn = nullptr;
	FRDGTextureRef qtildeIn_x = nullptr;
	FRDGTextureRef qtildeIn_y = nullptr;
	FRDGTextureRef Flow_x = nullptr;
	FRDGTextureRef Flow_y = nullptr;
	FRDGTextureRef hbarIn = nullptr;
	FRDGTextureRef hHat = nullptr;
	FRDGTextureRef qHat_x = nullptr;
	FRDGTextureRef qHat_y = nullptr;
	FRDGTextureRef qHat_x_array = nullptr;
	FRDGTextureRef qHat_y_array = nullptr;
	int32 PaddedSizeX = 0;
	int32 PaddedSizeY = 0;
	int32 GridSizeX = 0;
	int32 GridSizeY = 0;
	int32 DepthNum = 0;
};

struct FEWaveOutputs
{
	FRDGTextureRef htildeOldNext = nullptr;
	FRDGTextureRef qtildeOut_x = nullptr;
	FRDGTextureRef qtildeOut_y = nullptr;
};

void AddEWavePasses(
	FRDGBuilder& GraphBuilder,
	TUniformBufferRef<FSimConstants> ConstantBuffer,
	const FEWaveInputs& Inputs,
	const FEWaveOutputs& Outputs);

// --- Bulk Flow (SWE Velocity & Momentum) ---
struct FSWEBulkInputs
{
	FRDGTextureRef qbarIn_x = nullptr;
	FRDGTextureRef qbarIn_y = nullptr;
	FRDGTextureRef hbarIn = nullptr;
	FRDGTextureRef hbarOldIn = nullptr;
	FRDGTextureRef H_In = nullptr;
	FRDGTextureRef delH_x = nullptr;
	FRDGTextureRef delH_y = nullptr;
	FRDGTextureRef terrain = nullptr;
	int32 GridSizeX = 0;
	int32 GridSizeY = 0;
};

struct FSWEBulkOutputs
{
	FRDGTextureRef ubarOut_x = nullptr;
	FRDGTextureRef ubarOut_y = nullptr;
	FRDGTextureRef ubarNewOut_x = nullptr;
	FRDGTextureRef ubarNewOut_y = nullptr;
	FRDGTextureRef qbarOut_x = nullptr;
	FRDGTextureRef qbarOut_y = nullptr;
};

void AddSWEBulkPasses(
	FRDGBuilder& GraphBuilder,
	TUniformBufferRef<FSimConstants> ConstantBuffer,
	const FSWEBulkInputs& Inputs,
	const FSWEBulkOutputs& Outputs);

// --- Advective Transport & Height Integration ---
struct FTransportAndIntegrateInputs
{
	FRDGTextureRef ubarNewIn_x = nullptr;
	FRDGTextureRef ubarIn_x = nullptr;
	FRDGTextureRef ubarNewIn_y = nullptr;
	FRDGTextureRef ubarIn_y = nullptr;
	FRDGTextureRef qtildePast_x = nullptr;
	FRDGTextureRef qtildePast_y = nullptr;
	FRDGTextureRef htildePast = nullptr;
	FRDGTextureRef hPast = nullptr;
	FRDGTextureRef terrain = nullptr;
	FRDGTextureRef qbarIn_x = nullptr;
	FRDGTextureRef qbarIn_y = nullptr;
	int32 GridSizeX = 0;
	int32 GridSizeY = 0;
};

struct FTransportAndIntegrateOutputs
{
	FRDGTextureRef htildeOut = nullptr;
	FRDGTextureRef qtildeOut_x = nullptr;
	FRDGTextureRef qtildeOut_y = nullptr;
	FRDGTextureRef qAdvectOut_x = nullptr;
	FRDGTextureRef qAdvectOut_y = nullptr;
	FRDGTextureRef hOut = nullptr;
	FRDGTextureRef qOut_x = nullptr;
	FRDGTextureRef qOut_y = nullptr;
	FRDGTextureRef hdot_out = nullptr;
	FRDGTextureRef H_Out = nullptr;
	FRDGTextureRef HAvg_Out = nullptr;
};

void AddTransportAndIntegratePasses(
	FRDGBuilder& GraphBuilder,
	TUniformBufferRef<FSimConstants> ConstantBuffer,
	const FTransportAndIntegrateInputs& Inputs,
	const FTransportAndIntegrateOutputs& Outputs);
