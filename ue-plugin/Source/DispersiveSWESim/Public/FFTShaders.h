#pragma once

#include "CoreMinimal.h"
#include "GlobalShader.h"
#include "ShaderParameterStruct.h"
#include "RenderGraphResources.h"
#include "RenderGraphUtils.h"

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

/**
 * Executes a 2D Fast Fourier Transform (Row pass followed by Column pass) using RDG.
 * Supports both single 2D textures and 2D texture arrays.
 */
void Add2DFFTPasses(
	FRDGBuilder& GraphBuilder,
	FRDGTextureRef TargetTexture,
	int32 SizeX,
	int32 SizeY,
	bool bInverse,
	int32 NumLayers = 1);
