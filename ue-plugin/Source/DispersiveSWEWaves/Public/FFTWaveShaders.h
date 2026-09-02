#pragma once

#include "CoreMinimal.h"
#include "GlobalShader.h"
#include "ShaderParameterStruct.h"
#include "RenderGraphResources.h"
#include "RenderGraphUtils.h"

// Uniform buffer containing constants for FFT spectral wave generation and propagation
BEGIN_GLOBAL_SHADER_PARAMETER_STRUCT(FFFTWaveConstants, )
	SHADER_PARAMETER(float, time)
	SHADER_PARAMETER(int32, gridSizeX)
	SHADER_PARAMETER(int32, gridSizeY)
	SHADER_PARAMETER(float, cellSize)
	SHADER_PARAMETER(int32, paddedGridSizeX)
	SHADER_PARAMETER(int32, paddedGridSizeY)
	SHADER_PARAMETER(float, minWaterHeight)
	SHADER_PARAMETER(float, maxSafeDepth)
	SHADER_PARAMETER(int32, depthNum)
	SHADER_PARAMETER(float, surfaceTension)
	SHADER_PARAMETER(float, density)
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
	SHADER_PARAMETER_ARRAY(FVector4f, depthLevels, [4])
END_GLOBAL_SHADER_PARAMETER_STRUCT()

// --- Shaders from fftwaves.usf ---

class FPopulateSpectrumCS : public FGlobalShader {
public:
	DECLARE_GLOBAL_SHADER(FPopulateSpectrumCS);
	SHADER_USE_PARAMETER_STRUCT(FPopulateSpectrumCS, FGlobalShader);

	BEGIN_SHADER_PARAMETER_STRUCT(FParameters, )
		SHADER_PARAMETER_STRUCT_REF(FFFTWaveConstants, FFTWaveConstants)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2DArray<float2>, HPosOut)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2DArray<float2>, HNegOut)
	END_SHADER_PARAMETER_STRUCT()
};

class FPropagateWavesCS : public FGlobalShader {
public:
	DECLARE_GLOBAL_SHADER(FPropagateWavesCS);
	SHADER_USE_PARAMETER_STRUCT(FPropagateWavesCS, FGlobalShader);

	BEGIN_SHADER_PARAMETER_STRUCT(FParameters, )
		SHADER_PARAMETER_STRUCT_REF(FFFTWaveConstants, FFTWaveConstants)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2DArray<float2>, HPosIn)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2DArray<float2>, HNegIn)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2DArray<float2>, DispXOut)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2DArray<float2>, DispYOut)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2DArray<float2>, DelHXOut)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2DArray<float2>, DelHYOut)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2DArray<float2>, FlowXOut)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2DArray<float2>, FlowYOut)
	END_SHADER_PARAMETER_STRUCT()
};

class FInterpCS : public FGlobalShader {
public:
	DECLARE_GLOBAL_SHADER(FInterpCS);
	SHADER_USE_PARAMETER_STRUCT(FInterpCS, FGlobalShader);

	BEGIN_SHADER_PARAMETER_STRUCT(FParameters, )
		SHADER_PARAMETER_STRUCT_REF(FFFTWaveConstants, FFTWaveConstants)
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

// --- RDG Helper Structures & Functions ---

struct FPropagateFFTWavesInputs {
	FRDGTextureRef HPosIn = nullptr;
	FRDGTextureRef HNegIn = nullptr;
	FRDGTextureRef hbarIn = nullptr;

	int32 PaddedSizeX = 0;
	int32 PaddedSizeY = 0;
	int32 GridSizeX = 0;
	int32 GridSizeY = 0;
	int32 DepthNum = 0;
};

struct FPropagateFFTWavesOutputs {
	// Spatial 2D outputs (after iFFT & depth interpolation)
	FRDGTextureRef disp_x_Out = nullptr;
	FRDGTextureRef disp_y_Out = nullptr;
	FRDGTextureRef delH_x_Out = nullptr;
	FRDGTextureRef delH_y_Out = nullptr;

	// Spectral outputs for eWave coupling (2D complex array)
	FRDGTextureRef Flow_x_Out = nullptr;
	FRDGTextureRef Flow_y_Out = nullptr;
};

/**
 * Initializes spectrum HPos and HNeg textures.
 */
void AddPopulateSpectrumPass(
	FRDGBuilder& GraphBuilder,
	TUniformBufferRef<FFFTWaveConstants> ConstantBuffer,
	FRDGTextureRef HPosOut,
	FRDGTextureRef HNegOut,
	int32 PaddedSizeX,
	int32 PaddedSizeY,
	int32 DepthLevelsCount);

/**
 * Propagates spectral waves, runs inverse FFTs on arrays, and interpolates between depth levels.
 */
void AddPropagateFFTWavesPasses(
	FRDGBuilder& GraphBuilder,
	TUniformBufferRef<FFFTWaveConstants> ConstantBuffer,
	const FPropagateFFTWavesInputs& Inputs,
	const FPropagateFFTWavesOutputs& Outputs);
