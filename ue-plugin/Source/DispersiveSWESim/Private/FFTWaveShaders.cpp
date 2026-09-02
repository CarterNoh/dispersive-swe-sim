#include "FFTWaveShaders.h"
#include "FFTShaders.h"
#include "RenderGraphBuilder.h"

IMPLEMENT_GLOBAL_SHADER_PARAMETER_STRUCT(FFFTWaveConstants, "FFTWaveConstants");

IMPLEMENT_GLOBAL_SHADER(FPopulateSpectrumCS, "/Plugin/DispersiveSWESim/fftwaves.usf", "PopulateSpectrum", SF_Compute);
IMPLEMENT_GLOBAL_SHADER(FPropagateWavesCS,   "/Plugin/DispersiveSWESim/fftwaves.usf", "PropagateWaves",   SF_Compute);
IMPLEMENT_GLOBAL_SHADER(FInterpCS,           "/Plugin/DispersiveSWESim/fftwaves.usf", "Interp",           SF_Compute);

void AddPopulateSpectrumPass(
	FRDGBuilder& GraphBuilder,
	TUniformBufferRef<FFFTWaveConstants> ConstantBuffer,
	FRDGTextureRef HPosOut,
	FRDGTextureRef HNegOut,
	int32 PaddedSizeX,
	int32 PaddedSizeY,
	int32 DepthLevelsCount)
{
	FGlobalShaderMap* ShaderMap = GetGlobalShaderMap(GMaxRHIFeatureLevel);
	TShaderMapRef<FPopulateSpectrumCS> PopulateSpectrumCS(ShaderMap);
	FPopulateSpectrumCS::FParameters* PassParams = GraphBuilder.AllocParameters<FPopulateSpectrumCS::FParameters>();
	PassParams->FFTWaveConstants = ConstantBuffer;
	PassParams->HPosOut = GraphBuilder.CreateUAV(HPosOut);
	PassParams->HNegOut = GraphBuilder.CreateUAV(HNegOut);

	FComputeShaderUtils::AddPass(
		GraphBuilder,
		RDG_EVENT_NAME("SWE_PopulateSpectrum"),
		ERDGPassFlags::Compute,
		PopulateSpectrumCS,
		PassParams,
		FIntVector(
			FMath::DivideAndRoundUp(PaddedSizeX, 16),
			FMath::DivideAndRoundUp(PaddedSizeY, 16),
			DepthLevelsCount
		)
	);
}

void AddPropagateFFTWavesPasses(
	FRDGBuilder& GraphBuilder,
	TUniformBufferRef<FFFTWaveConstants> ConstantBuffer,
	const FPropagateFFTWavesInputs& Inputs,
	const FPropagateFFTWavesOutputs& Outputs)
{
	FGlobalShaderMap* ShaderMap = GetGlobalShaderMap(GMaxRHIFeatureLevel);

	FIntVector ComplexArrayGroups(
		FMath::DivideAndRoundUp(Inputs.PaddedSizeX, 16),
		FMath::DivideAndRoundUp(Inputs.PaddedSizeY, 16),
		Inputs.DepthNum
	);
	FIntVector GridGroups(
		FMath::DivideAndRoundUp(Inputs.GridSizeX, 16),
		FMath::DivideAndRoundUp(Inputs.GridSizeY, 16),
		1
	);

	// Propagate spectral waves
	{
		TShaderMapRef<FPropagateWavesCS> Shader(ShaderMap);
		FPropagateWavesCS::FParameters* Params = GraphBuilder.AllocParameters<FPropagateWavesCS::FParameters>();
		Params->FFTWaveConstants = ConstantBuffer;
		Params->HPosIn = GraphBuilder.CreateSRV(Inputs.HPosIn);
		Params->HNegIn = GraphBuilder.CreateSRV(Inputs.HNegIn);
		Params->DispXOut = GraphBuilder.CreateUAV(Inputs.Disp_x_Array);
		Params->DispYOut = GraphBuilder.CreateUAV(Inputs.Disp_y_Array);
		Params->DelHXOut = GraphBuilder.CreateUAV(Inputs.DelH_x_Array);
		Params->DelHYOut = GraphBuilder.CreateUAV(Inputs.DelH_y_Array);
		Params->FlowXOut = GraphBuilder.CreateUAV(Inputs.Flow_x_Array);
		Params->FlowYOut = GraphBuilder.CreateUAV(Inputs.Flow_y_Array);

		FComputeShaderUtils::AddPass(
			GraphBuilder,
			RDG_EVENT_NAME("SWE_FFTWaves_Propagate"),
			ERDGPassFlags::Compute,
			Shader,
			Params,
			ComplexArrayGroups
		);
	}

	// Inverse FFTs on Disp and DelH 2D Texture Arrays
	Add2DFFTPasses(GraphBuilder, Inputs.Disp_x_Array, Inputs.PaddedSizeX, Inputs.PaddedSizeY, true, Inputs.DepthNum);
	Add2DFFTPasses(GraphBuilder, Inputs.Disp_y_Array, Inputs.PaddedSizeX, Inputs.PaddedSizeY, true, Inputs.DepthNum);
	Add2DFFTPasses(GraphBuilder, Inputs.DelH_x_Array, Inputs.PaddedSizeX, Inputs.PaddedSizeY, true, Inputs.DepthNum);
	Add2DFFTPasses(GraphBuilder, Inputs.DelH_y_Array, Inputs.PaddedSizeX, Inputs.PaddedSizeY, true, Inputs.DepthNum);

	// Interpolate wind wave outputs between depths
	{
		TShaderMapRef<FInterpCS> Shader(ShaderMap);
		FInterpCS::FParameters* Params = GraphBuilder.AllocParameters<FInterpCS::FParameters>();
		Params->FFTWaveConstants = ConstantBuffer;
		Params->HxIn = GraphBuilder.CreateSRV(Inputs.DelH_x_Array);
		Params->HyIn = GraphBuilder.CreateSRV(Inputs.DelH_y_Array);
		Params->DxIn = GraphBuilder.CreateSRV(Inputs.Disp_x_Array);
		Params->DyIn = GraphBuilder.CreateSRV(Inputs.Disp_y_Array);
		Params->hbarIn = GraphBuilder.CreateSRV(Inputs.hbarIn);
		Params->HxOut = GraphBuilder.CreateUAV(Outputs.delH_x_Out);
		Params->HyOut = GraphBuilder.CreateUAV(Outputs.delH_y_Out);
		Params->DxOut = GraphBuilder.CreateUAV(Outputs.disp_x_Out);
		Params->DyOut = GraphBuilder.CreateUAV(Outputs.disp_y_Out);

		FComputeShaderUtils::AddPass(
			GraphBuilder,
			RDG_EVENT_NAME("SWE_FFTWaves_Interp"),
			ERDGPassFlags::Compute,
			Shader,
			Params,
			GridGroups
		);
	}
}
