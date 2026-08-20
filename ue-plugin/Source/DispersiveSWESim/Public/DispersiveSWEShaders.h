#pragma once

#include "CoreMinimal.h"
#include "ShaderParameters.h"
#include "ShaderParameterUtils.h"
#include "ShaderParameterStruct.h"
#include "GlobalShader.h"

// Uniform buffer containing all simulation constants
BEGIN_GLOBAL_SHADER_PARAMETER_STRUCT(FSimConstants, )
	SHADER_PARAMETER(float, time)
	SHADER_PARAMETER(int32, gridSizeX)
	SHADER_PARAMETER(int32, gridSizeY)
	SHADER_PARAMETER(float, cellSize)
	SHADER_PARAMETER(float, timeStep)
	SHADER_PARAMETER(float, minWaterHeight)
	SHADER_PARAMETER(float, surfaceTension)
	SHADER_PARAMETER(float, density)
	SHADER_PARAMETER(int32, diffusionIterations)
	SHADER_PARAMETER(float, diffusionTime)
	SHADER_PARAMETER(float, diffusionPenalty)
	SHADER_PARAMETER(float, slopeLimit)
	SHADER_PARAMETER(float, cflCondition)
	SHADER_PARAMETER(float, gammaTransport)
	SHADER_PARAMETER(int32, spongeThickness)
	SHADER_PARAMETER(float, laplacianDamping)
	SHADER_PARAMETER(int32, depthNum)
	SHADER_PARAMETER(float, fetch)
	SHADER_PARAMETER(float, windSpeed)
	SHADER_PARAMETER(float, windAngle)
	SHADER_PARAMETER(float, swell)
	SHADER_PARAMETER(float, swellAngle)
	SHADER_PARAMETER(float, choppiness)
	SHADER_PARAMETER(float, filterSmall)
	SHADER_PARAMETER(float, filterBig)
	SHADER_PARAMETER(float, filterWidth)
	SHADER_PARAMETER(float, filterMin)
	SHADER_PARAMETER(float, depthCutoff)
	SHADER_PARAMETER(int32, paddedGridSizeX)
	SHADER_PARAMETER(int32, paddedGridSizeY)
	SHADER_PARAMETER(float, maxSafeDepth)
	SHADER_PARAMETER(float, foamThreshold)
	SHADER_PARAMETER(float, foamMultiplier)
	SHADER_PARAMETER(float, foamFade)
	SHADER_PARAMETER(float, foamBlur)
	SHADER_PARAMETER_ARRAY(FVector4f, depthLevels, [4])
END_GLOBAL_SHADER_PARAMETER_STRUCT()

// Base class for our compute shaders to reduce duplicate code
class FDispersiveSWEComputeShader : public FGlobalShader
{
public:
	FDispersiveSWEComputeShader() {}
	FDispersiveSWEComputeShader(const ShaderMetaType::CompiledShaderInitializerType& Initializer)
		: FGlobalShader(Initializer) {}
};

// --- Shaders from kernels.usf ---

class FInitializeWaterCS : public FDispersiveSWEComputeShader
{
public:
	DECLARE_GLOBAL_SHADER(FInitializeWaterCS);
	SHADER_USE_PARAMETER_STRUCT(FInitializeWaterCS, FDispersiveSWEComputeShader);

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

class FInitDecompCS : public FDispersiveSWEComputeShader
{
public:
	DECLARE_GLOBAL_SHADER(FInitDecompCS);
	SHADER_USE_PARAMETER_STRUCT(FInitDecompCS, FDispersiveSWEComputeShader);

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

class FRecomputeHCS : public FDispersiveSWEComputeShader
{
public:
	DECLARE_GLOBAL_SHADER(FRecomputeHCS);
	SHADER_USE_PARAMETER_STRUCT(FRecomputeHCS, FDispersiveSWEComputeShader);

	BEGIN_SHADER_PARAMETER_STRUCT(FParameters, )
		SHADER_PARAMETER_STRUCT_REF(FSimConstants, SimConstants)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, hIn)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, terrain)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, H_Out)
	END_SHADER_PARAMETER_STRUCT()
};

class FCalcDiffusionCoeffsCS : public FDispersiveSWEComputeShader
{
public:
	DECLARE_GLOBAL_SHADER(FCalcDiffusionCoeffsCS);
	SHADER_USE_PARAMETER_STRUCT(FCalcDiffusionCoeffsCS, FDispersiveSWEComputeShader);

	BEGIN_SHADER_PARAMETER_STRUCT(FParameters, )
		SHADER_PARAMETER_STRUCT_REF(FSimConstants, SimConstants)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, H_In)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, terrain)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, alpha_HOut)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, alpha_QOut_x)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, alpha_QOut_y)
	END_SHADER_PARAMETER_STRUCT()
};

class FDiffusionStepCS : public FDispersiveSWEComputeShader
{
public:
	DECLARE_GLOBAL_SHADER(FDiffusionStepCS);
	SHADER_USE_PARAMETER_STRUCT(FDiffusionStepCS, FDispersiveSWEComputeShader);

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

class FDecomposeFieldsCS : public FDispersiveSWEComputeShader
{
public:
	DECLARE_GLOBAL_SHADER(FDecomposeFieldsCS);
	SHADER_USE_PARAMETER_STRUCT(FDecomposeFieldsCS, FDispersiveSWEComputeShader);

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

class FTransferToFFTCS : public FDispersiveSWEComputeShader
{
public:
	DECLARE_GLOBAL_SHADER(FTransferToFFTCS);
	SHADER_USE_PARAMETER_STRUCT(FTransferToFFTCS, FDispersiveSWEComputeShader);

	BEGIN_SHADER_PARAMETER_STRUCT(FParameters, )
		SHADER_PARAMETER_STRUCT_REF(FSimConstants, SimConstants)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, htildeIn)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, qtildeIn_x)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, qtildeIn_y)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, htildeOldCopy)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, htildeOldOut)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float2>, hHat)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float2>, qHat_x)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float2>, qHat_y)
	END_SHADER_PARAMETER_STRUCT()
};

class FCalcEWaveCS : public FDispersiveSWEComputeShader
{
public:
	DECLARE_GLOBAL_SHADER(FCalcEWaveCS);
	SHADER_USE_PARAMETER_STRUCT(FCalcEWaveCS, FDispersiveSWEComputeShader);

	BEGIN_SHADER_PARAMETER_STRUCT(FParameters, )
		SHADER_PARAMETER_STRUCT_REF(FSimConstants, SimConstants)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float2>, hhat)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float2>, qhat_x)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float2>, qhat_y)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2DArray<float2>, qhat_x_array)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2DArray<float2>, qhat_y_array)
	END_SHADER_PARAMETER_STRUCT()
};

class FInterpQCS : public FDispersiveSWEComputeShader
{
public:
	DECLARE_GLOBAL_SHADER(FInterpQCS);
	SHADER_USE_PARAMETER_STRUCT(FInterpQCS, FDispersiveSWEComputeShader);

	BEGIN_SHADER_PARAMETER_STRUCT(FParameters, )
		SHADER_PARAMETER_STRUCT_REF(FSimConstants, SimConstants)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, hbarIn)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2DArray<float2>, qHat_x_array)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2DArray<float2>, qHat_y_array)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, qtildeOut_x)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, qtildeOut_y)
	END_SHADER_PARAMETER_STRUCT()
};

class FCalcUbarCS : public FDispersiveSWEComputeShader
{
public:
	DECLARE_GLOBAL_SHADER(FCalcUbarCS);
	SHADER_USE_PARAMETER_STRUCT(FCalcUbarCS, FDispersiveSWEComputeShader);

	BEGIN_SHADER_PARAMETER_STRUCT(FParameters, )
		SHADER_PARAMETER_STRUCT_REF(FSimConstants, SimConstants)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, qbarIn_x)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, qbarIn_y)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, hbarIn)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, ubarOut_x)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, ubarOut_y)
	END_SHADER_PARAMETER_STRUCT()
};

class FCalcSWECS : public FDispersiveSWEComputeShader
{
public:
	DECLARE_GLOBAL_SHADER(FCalcSWECS);
	SHADER_USE_PARAMETER_STRUCT(FCalcSWECS, FDispersiveSWEComputeShader);

	BEGIN_SHADER_PARAMETER_STRUCT(FParameters, )
		SHADER_PARAMETER_STRUCT_REF(FSimConstants, SimConstants)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, ubarIn_x)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, ubarIn_y)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, hbarIn)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, H_In)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, delH_x)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, delH_y)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, ubarNewOut_x)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, ubarNewOut_y)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, qbarOut_x)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, qbarOut_y)
	END_SHADER_PARAMETER_STRUCT()
};

class FUpdateTildeCS : public FDispersiveSWEComputeShader
{
public:
	DECLARE_GLOBAL_SHADER(FUpdateTildeCS);
	SHADER_USE_PARAMETER_STRUCT(FUpdateTildeCS, FDispersiveSWEComputeShader);

	BEGIN_SHADER_PARAMETER_STRUCT(FParameters, )
		SHADER_PARAMETER_STRUCT_REF(FSimConstants, SimConstants)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, ubarNewIn_x)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, ubarIn_x)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, ubarNewIn_y)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, ubarIn_y)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, qtildePast_x)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, qtildePast_y)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, hIn)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, htildeCopy)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, terrain)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, htildeOut)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, qtildeOut_x)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, qtildeOut_y)
	END_SHADER_PARAMETER_STRUCT()
};

class FCalcQAdvectCS : public FDispersiveSWEComputeShader
{
public:
	DECLARE_GLOBAL_SHADER(FCalcQAdvectCS);
	SHADER_USE_PARAMETER_STRUCT(FCalcQAdvectCS, FDispersiveSWEComputeShader);

	BEGIN_SHADER_PARAMETER_STRUCT(FParameters, )
		SHADER_PARAMETER_STRUCT_REF(FSimConstants, SimConstants)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, ubarNewIn_x)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, ubarNewIn_y)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, htildeIn)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, qAdvectOut_x)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, qAdvectOut_y)
	END_SHADER_PARAMETER_STRUCT()
};

class FIntegrateHCS : public FDispersiveSWEComputeShader
{
public:
	DECLARE_GLOBAL_SHADER(FIntegrateHCS);
	SHADER_USE_PARAMETER_STRUCT(FIntegrateHCS, FDispersiveSWEComputeShader);

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
	END_SHADER_PARAMETER_STRUCT()
};

class FScaleCopyDisplacementCS : public FDispersiveSWEComputeShader
{
public:
	DECLARE_GLOBAL_SHADER(FScaleCopyDisplacementCS);
	SHADER_USE_PARAMETER_STRUCT(FScaleCopyDisplacementCS, FDispersiveSWEComputeShader);

	BEGIN_SHADER_PARAMETER_STRUCT(FParameters, )
		SHADER_PARAMETER_STRUCT_REF(FSimConstants, SimConstants)
		SHADER_PARAMETER(float, ScaleFactor)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, inDispX)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, inDispY)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, inHeight)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float4>, outDisp4)
	END_SHADER_PARAMETER_STRUCT()
};

class FCalcSurfaceNormalAndFoamCS : public FDispersiveSWEComputeShader
{
public:
	DECLARE_GLOBAL_SHADER(FCalcSurfaceNormalAndFoamCS);
	SHADER_USE_PARAMETER_STRUCT(FCalcSurfaceNormalAndFoamCS, FDispersiveSWEComputeShader);

	BEGIN_SHADER_PARAMETER_STRUCT(FParameters, )
		SHADER_PARAMETER_STRUCT_REF(FSimConstants, SimConstants)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, inDispX)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, inDispY)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, inHeight)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, inPreviousFoam)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float4>, outNormal)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, outFoam)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, outJacobianDet)
	END_SHADER_PARAMETER_STRUCT()
};

class FCalcRoughnessLUTCS : public FDispersiveSWEComputeShader
{
public:
	DECLARE_GLOBAL_SHADER(FCalcRoughnessLUTCS);
	SHADER_USE_PARAMETER_STRUCT(FCalcRoughnessLUTCS, FDispersiveSWEComputeShader);

	BEGIN_SHADER_PARAMETER_STRUCT(FParameters, )
		SHADER_PARAMETER_STRUCT_REF(FSimConstants, SimConstants)
		SHADER_PARAMETER(float, IntegrationSamples)
		SHADER_PARAMETER(float, RoughnessPower)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float4>, inNormal)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, inPreviousRoughness)
		SHADER_PARAMETER_SAMPLER(SamplerState, BilinearWrapSampler)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, outRoughness)
	END_SHADER_PARAMETER_STRUCT()
};


// --- Shaders from fftwaves.usf ---

class FPopulateSpectrumCS : public FDispersiveSWEComputeShader
{
public:
	DECLARE_GLOBAL_SHADER(FPopulateSpectrumCS);
	SHADER_USE_PARAMETER_STRUCT(FPopulateSpectrumCS, FDispersiveSWEComputeShader);

	BEGIN_SHADER_PARAMETER_STRUCT(FParameters, )
		SHADER_PARAMETER_STRUCT_REF(FSimConstants, SimConstants)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2DArray<float2>, HPosOut)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2DArray<float2>, HNegOut)
	END_SHADER_PARAMETER_STRUCT()
};

class FPropagateWavesCS : public FDispersiveSWEComputeShader
{
public:
	DECLARE_GLOBAL_SHADER(FPropagateWavesCS);
	SHADER_USE_PARAMETER_STRUCT(FPropagateWavesCS, FDispersiveSWEComputeShader);

	BEGIN_SHADER_PARAMETER_STRUCT(FParameters, )
		SHADER_PARAMETER_STRUCT_REF(FSimConstants, SimConstants)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2DArray<float2>, HPosIn)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2DArray<float2>, HNegIn)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2DArray<float2>, DelHxOut)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2DArray<float2>, DelHyOut)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2DArray<float2>, DispXOut)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2DArray<float2>, DispYOut)
	END_SHADER_PARAMETER_STRUCT()
};

class FInterpCS : public FDispersiveSWEComputeShader
{
public:
	DECLARE_GLOBAL_SHADER(FInterpCS);
	SHADER_USE_PARAMETER_STRUCT(FInterpCS, FDispersiveSWEComputeShader);

	BEGIN_SHADER_PARAMETER_STRUCT(FParameters, )
		SHADER_PARAMETER_STRUCT_REF(FSimConstants, SimConstants)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2DArray<float2>, HxIn)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2DArray<float2>, HyIn)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2DArray<float2>, DxIn)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2DArray<float2>, DyIn)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, hbarIn)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, HxOut)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, HyOut)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, DxOut)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, DyOut)
	END_SHADER_PARAMETER_STRUCT()
};


// --- Shaders from fft.usf ---

class FFFTKernel1DCS : public FGlobalShader
{
public:
	DECLARE_GLOBAL_SHADER(FFFTKernel1DCS);
	SHADER_USE_PARAMETER_STRUCT(FFFTKernel1DCS, FGlobalShader);

	// Sparse permutations supporting different power-of-two FFT sizes
	class FFFTSizeDim : SHADER_PERMUTATION_SPARSE_INT("FFT_SIZE", 32, 64, 128, 256, 512, 1024, 2048);
	class FIsArrayDim : SHADER_PERMUTATION_BOOL("IS_ARRAY");

	using FPermutationDomain = TShaderPermutationDomain<FFFTSizeDim, FIsArrayDim>;

	BEGIN_SHADER_PARAMETER_STRUCT(FParameters, )
		SHADER_PARAMETER(int32, cb_Nx)
		SHADER_PARAMETER(int32, cb_Ny)
		SHADER_PARAMETER(int32, cb_BitsX)
		SHADER_PARAMETER(int32, cb_BitsY)
		SHADER_PARAMETER(int32, cb_Inverse)
		SHADER_PARAMETER(int32, cb_IsRow)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float2>, fft)
	END_SHADER_PARAMETER_STRUCT()

	static bool ShouldCompilePermutation(const FGlobalShaderPermutationParameters& Parameters)
	{
		return true;
	}

	static void ModifyCompilationEnvironment(const FGlobalShaderPermutationParameters& Parameters, FShaderCompilerEnvironment& OutEnvironment)
	{
		FGlobalShader::ModifyCompilationEnvironment(Parameters, OutEnvironment);
		
		FPermutationDomain PermutationVector(Parameters.PermutationId);
		OutEnvironment.SetDefine(TEXT("FFT_SIZE"), PermutationVector.Get<FFFTSizeDim>());
		OutEnvironment.SetDefine(TEXT("IS_ARRAY"), PermutationVector.Get<FIsArrayDim>() ? 1 : 0);
	}
};

