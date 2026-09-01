#pragma once

#include "CoreMinimal.h"
#include "GlobalShader.h"
#include "ShaderParameterStruct.h"
#include "RenderGraphResources.h"
#include "RenderGraphUtils.h"

// Uniform buffer containing constants for visual export (displacement, foam, normal, roughness)
BEGIN_GLOBAL_SHADER_PARAMETER_STRUCT(FExportConstants, )
	SHADER_PARAMETER(float, time)
	SHADER_PARAMETER(int32, gridSizeX)
	SHADER_PARAMETER(int32, gridSizeY)
	SHADER_PARAMETER(float, cellSize)
	SHADER_PARAMETER(float, timeStep)
	SHADER_PARAMETER(float, foamThreshold)
	SHADER_PARAMETER(float, foamMultiplier)
	SHADER_PARAMETER(float, foamFade)
	SHADER_PARAMETER(float, foamBlur)
END_GLOBAL_SHADER_PARAMETER_STRUCT()

// --- Shaders from export.usf ---

class FScaleCopyDisplacementCS : public FGlobalShader
{
public:
	DECLARE_GLOBAL_SHADER(FScaleCopyDisplacementCS);
	SHADER_USE_PARAMETER_STRUCT(FScaleCopyDisplacementCS, FGlobalShader);

	BEGIN_SHADER_PARAMETER_STRUCT(FParameters, )
		SHADER_PARAMETER_STRUCT_REF(FExportConstants, ExportConstants)
		SHADER_PARAMETER(float, ScaleFactor)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, inDispX)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, inDispY)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, inHeight)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float4>, outDisp4)
	END_SHADER_PARAMETER_STRUCT()
};

class FCalcSurfaceNormalAndFoamCS : public FGlobalShader
{
public:
	DECLARE_GLOBAL_SHADER(FCalcSurfaceNormalAndFoamCS);
	SHADER_USE_PARAMETER_STRUCT(FCalcSurfaceNormalAndFoamCS, FGlobalShader);

	BEGIN_SHADER_PARAMETER_STRUCT(FParameters, )
		SHADER_PARAMETER_STRUCT_REF(FExportConstants, ExportConstants)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, inDispX)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, inDispY)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, inHeight)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, inPreviousFoam)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float4>, outNormal)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, outFoam)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, outJacobianDet)
	END_SHADER_PARAMETER_STRUCT()
};

class FCalcRoughnessLUTCS : public FGlobalShader
{
public:
	DECLARE_GLOBAL_SHADER(FCalcRoughnessLUTCS);
	SHADER_USE_PARAMETER_STRUCT(FCalcRoughnessLUTCS, FGlobalShader);

	BEGIN_SHADER_PARAMETER_STRUCT(FParameters, )
		SHADER_PARAMETER_STRUCT_REF(FExportConstants, ExportConstants)
		SHADER_PARAMETER(float, IntegrationSamples)
		SHADER_PARAMETER(float, RoughnessPower)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float4>, inNormal)
		SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, inPreviousRoughness)
		SHADER_PARAMETER_SAMPLER(SamplerState, BilinearWrapSampler)
		SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, outRoughness)
	END_SHADER_PARAMETER_STRUCT()
};

// --- RDG Helper Structures & Functions ---

struct FVisualExportInputs
{
	FRDGTextureRef inDispX = nullptr;
	FRDGTextureRef inDispY = nullptr;
	FRDGTextureRef inHeight = nullptr;
	FRDGTextureRef inPreviousFoam = nullptr;
	FRDGTextureRef inPreviousRoughness = nullptr;

	FRDGTextureRef ExportDispDest = nullptr;
	FRDGTextureRef ExportNormalDest = nullptr;
	FRDGTextureRef ExportFoamDest = nullptr;
	FRDGTextureRef ExportJacobianDest = nullptr;
	FRDGTextureRef ExportRoughnessDest = nullptr;

	float ScaleFactor = 100.0f; // m to cm
	float IntegrationSamples = 100.0f;
	float RoughnessPower = 1.0f;
	int32 GridSizeX = 0;
	int32 GridSizeY = 0;
};

struct FVisualExportOutputs
{
	FRDGTextureRef outNewFoam = nullptr;
	FRDGTextureRef outNewRoughness = nullptr;
};

/**
 * Executes all visual export passes using FExportConstants.
 */
void AddVisualExportPasses(
	FRDGBuilder& GraphBuilder,
	TUniformBufferRef<FExportConstants> ConstantBuffer,
	const FVisualExportInputs& Inputs,
	const FVisualExportOutputs& Outputs);
