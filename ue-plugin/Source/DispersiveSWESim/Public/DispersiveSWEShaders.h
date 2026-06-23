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
	SHADER_PARAMETER(int32, spongeThickness)
	SHADER_PARAMETER(float, minWaterHeight)
	SHADER_PARAMETER(float, surfaceTension)
	SHADER_PARAMETER(float, density)
	SHADER_PARAMETER(int32, diffusionIterations)
	SHADER_PARAMETER(float, deltaT)
	SHADER_PARAMETER(float, diffusionPenalty)
	SHADER_PARAMETER(float, slopeLimit)
	SHADER_PARAMETER(float, cflCondition)
	SHADER_PARAMETER(float, gammaTransport)
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
END_GLOBAL_SHADER_PARAMETER_STRUCT()

// Base class for our compute shaders to reduce duplicate code
class FDispersiveSWEComputeShader : public FGlobalShader
{
public:
	FDispersiveSWEComputeShader() {}
	FDispersiveSWEComputeShader(const ShaderMetaType::CompiledShaderInitializerType& Initializer)
		: FGlobalShader(Initializer)
	{}
};

// --- Shaders from kernels.usf ---

class FInitDecompCS : public FDispersiveSWEComputeShader
{
public:
	DECLARE_GLOBAL_SHADER(FInitDecompCS);
	SHADER_USE_PARAMETER_STRUCT(FInitDecompCS, FDispersiveSWEComputeShader);

	BEGIN_SHADER_PARAMETER_STRUCT(FParameters, )
		SHADER_PARAMETER_STRUCT_REF(FSimConstants, SimConstants)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, in0)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, in1)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, in2)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, in3)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, out0)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, out1)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, out2)
	END_SHADER_PARAMETER_STRUCT()
};

class FCalcDiffusionCoeffsCS : public FDispersiveSWEComputeShader
{
public:
	DECLARE_GLOBAL_SHADER(FCalcDiffusionCoeffsCS);
	SHADER_USE_PARAMETER_STRUCT(FCalcDiffusionCoeffsCS, FDispersiveSWEComputeShader);

	BEGIN_SHADER_PARAMETER_STRUCT(FParameters, )
		SHADER_PARAMETER_STRUCT_REF(FSimConstants, SimConstants)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, in0)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, in1)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, out0)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, out1)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, out2)
	END_SHADER_PARAMETER_STRUCT()
};

class FDiffusionStepCS : public FDispersiveSWEComputeShader
{
public:
	DECLARE_GLOBAL_SHADER(FDiffusionStepCS);
	SHADER_USE_PARAMETER_STRUCT(FDiffusionStepCS, FDispersiveSWEComputeShader);

	BEGIN_SHADER_PARAMETER_STRUCT(FParameters, )
		SHADER_PARAMETER_STRUCT_REF(FSimConstants, SimConstants)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, in0)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, in1)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, in2)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, in3)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, in4)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, in5)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, in6)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, out0)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, out1)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, out2)
	END_SHADER_PARAMETER_STRUCT()
};

class FDecomposeFieldsCS : public FDispersiveSWEComputeShader
{
public:
	DECLARE_GLOBAL_SHADER(FDecomposeFieldsCS);
	SHADER_USE_PARAMETER_STRUCT(FDecomposeFieldsCS, FDispersiveSWEComputeShader);

	BEGIN_SHADER_PARAMETER_STRUCT(FParameters, )
		SHADER_PARAMETER_STRUCT_REF(FSimConstants, SimConstants)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, in0)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, in1)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, in2)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, in3)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, in4)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, in5)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, in6)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, out0)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, out1)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, out2)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, out3)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, out4)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, out5)
	END_SHADER_PARAMETER_STRUCT()
};

class FCalcUbarCS : public FDispersiveSWEComputeShader
{
public:
	DECLARE_GLOBAL_SHADER(FCalcUbarCS);
	SHADER_USE_PARAMETER_STRUCT(FCalcUbarCS, FDispersiveSWEComputeShader);

	BEGIN_SHADER_PARAMETER_STRUCT(FParameters, )
		SHADER_PARAMETER_STRUCT_REF(FSimConstants, SimConstants)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, in0)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, in1)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, in2)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, out0)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, out1)
	END_SHADER_PARAMETER_STRUCT()
};

class FCalcSWECS : public FDispersiveSWEComputeShader
{
public:
	DECLARE_GLOBAL_SHADER(FCalcSWECS);
	SHADER_USE_PARAMETER_STRUCT(FCalcSWECS, FDispersiveSWEComputeShader);

	BEGIN_SHADER_PARAMETER_STRUCT(FParameters, )
		SHADER_PARAMETER_STRUCT_REF(FSimConstants, SimConstants)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, in0)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, in1)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, in2)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, in3)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, in4)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, in5)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, out0)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, out1)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, out2)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, out3)
	END_SHADER_PARAMETER_STRUCT()
};

class FUpdateTildeCS : public FDispersiveSWEComputeShader
{
public:
	DECLARE_GLOBAL_SHADER(FUpdateTildeCS);
	SHADER_USE_PARAMETER_STRUCT(FUpdateTildeCS, FDispersiveSWEComputeShader);

	BEGIN_SHADER_PARAMETER_STRUCT(FParameters, )
		SHADER_PARAMETER_STRUCT_REF(FSimConstants, SimConstants)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, in0)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, in1)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, in2)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, in3)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, in4)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, in5)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, in6)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, in7)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, out0)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, out1)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, out2)
	END_SHADER_PARAMETER_STRUCT()
};

class FCalcQAdvectCS : public FDispersiveSWEComputeShader
{
public:
	DECLARE_GLOBAL_SHADER(FCalcQAdvectCS);
	SHADER_USE_PARAMETER_STRUCT(FCalcQAdvectCS, FDispersiveSWEComputeShader);

	BEGIN_SHADER_PARAMETER_STRUCT(FParameters, )
		SHADER_PARAMETER_STRUCT_REF(FSimConstants, SimConstants)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, in0)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, in1)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, in2)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, out0)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, out1)
	END_SHADER_PARAMETER_STRUCT()
};

class FIntegrateHCS : public FDispersiveSWEComputeShader
{
public:
	DECLARE_GLOBAL_SHADER(FIntegrateHCS);
	SHADER_USE_PARAMETER_STRUCT(FIntegrateHCS, FDispersiveSWEComputeShader);

	BEGIN_SHADER_PARAMETER_STRUCT(FParameters, )
		SHADER_PARAMETER_STRUCT_REF(FSimConstants, SimConstants)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, in0)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, in1)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, in2)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, in3)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, in4)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, in5)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, in6)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, in7)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, out0)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, out1)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, out2)
	END_SHADER_PARAMETER_STRUCT()
};

class FTransferToFFTCS : public FDispersiveSWEComputeShader
{
public:
	DECLARE_GLOBAL_SHADER(FTransferToFFTCS);
	SHADER_USE_PARAMETER_STRUCT(FTransferToFFTCS, FDispersiveSWEComputeShader);

	BEGIN_SHADER_PARAMETER_STRUCT(FParameters, )
		SHADER_PARAMETER_STRUCT_REF(FSimConstants, SimConstants)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, in0)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, in1)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, in2)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, out0)
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
		SHADER_PARAMETER_RDG_BUFFER_SRV(StructuredBuffer<float>, in8) // depth buffer
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
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, hbar)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2DArray<float2>, qHat_x_array)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2DArray<float2>, qHat_y_array)
		SHADER_PARAMETER_RDG_BUFFER_SRV(StructuredBuffer<float>, in8)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, qtilde_x)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, qtilde_y)
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
		SHADER_PARAMETER_RDG_BUFFER_SRV(StructuredBuffer<float>, depth)
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
		SHADER_PARAMETER_RDG_BUFFER_SRV(StructuredBuffer<float>, depth)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2DArray<float2>, HPropOut)
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
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2DArray<float2>, HIn)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2DArray<float2>, HxIn)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2DArray<float2>, HyIn)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2DArray<float2>, DxIn)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2DArray<float2>, DyIn)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, hbar)
		SHADER_PARAMETER_RDG_BUFFER_SRV(StructuredBuffer<float>, depth)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, HOut)
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
	class FFFTSizeDim : SHADER_PERMUTATION_SPARSE_INT("FFT_SIZE", 32, 64, 128, 256, 512, 1024);
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

class FInitializeWaterHeightCS : public FDispersiveSWEComputeShader
{
public:
	DECLARE_GLOBAL_SHADER(FInitializeWaterHeightCS);
	SHADER_USE_PARAMETER_STRUCT(FInitializeWaterHeightCS, FDispersiveSWEComputeShader);

	BEGIN_SHADER_PARAMETER_STRUCT(FParameters, )
		SHADER_PARAMETER_STRUCT_REF(FSimConstants, SimConstants)
		SHADER_PARAMETER(float, WaterLevel)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, in3)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, out0)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, out1)
	END_SHADER_PARAMETER_STRUCT()
};

class FScaleCopyTextureCS : public FDispersiveSWEComputeShader
{
public:
	DECLARE_GLOBAL_SHADER(FScaleCopyTextureCS);
	SHADER_USE_PARAMETER_STRUCT(FScaleCopyTextureCS, FDispersiveSWEComputeShader);

	BEGIN_SHADER_PARAMETER_STRUCT(FParameters, )
		SHADER_PARAMETER_STRUCT_REF(FSimConstants, SimConstants)
		SHADER_PARAMETER(float, ScaleFactor)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, in0)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, out0)
	END_SHADER_PARAMETER_STRUCT()
};
