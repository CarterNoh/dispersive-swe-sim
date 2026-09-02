#include "FFTWaveShaders.h"
#include "FFTShaders.h"
#include "RenderGraphBuilder.h"

IMPLEMENT_GLOBAL_SHADER_PARAMETER_STRUCT(FFFTWaveConstants, "FFTWaveConstants");

IMPLEMENT_GLOBAL_SHADER(FPopulateSpectrumCS, "/Plugin/DispersiveSWEWaves/fftwaves.usf", "PopulateSpectrum", SF_Compute);
IMPLEMENT_GLOBAL_SHADER(FPropagateWavesCS,   "/Plugin/DispersiveSWEWaves/fftwaves.usf", "PropagateWaves",   SF_Compute);
IMPLEMENT_GLOBAL_SHADER(FInterpCS,           "/Plugin/DispersiveSWEWaves/fftwaves.usf", "Interp",           SF_Compute);

void AddPopulateSpectrumPass(
	FRDGBuilder& GraphBuilder,
	TUniformBufferRef<FFFTWaveConstants> ConstantBuffer,
	FRDGTextureRef HPosOut,
	FRDGTextureRef HNegOut,
	int32 PaddedSizeX,
	int32 PaddedSizeY,
	int32 DepthLevelsCount) {
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
	const FPropagateFFTWavesOutputs& Outputs) {
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

	// Transient complex arrays for FFT wave workspace
	FRDGTextureDesc ComplexArrayDesc = FRDGTextureDesc::Create2DArray(
		FIntPoint(Inputs.PaddedSizeX, Inputs.PaddedSizeY),
		PF_G32R32F,
		FClearValueBinding::None,
		TexCreate_ShaderResource | TexCreate_UAV,
		Inputs.DepthNum
	);

	FRDGTextureRef Disp_x_Array = GraphBuilder.CreateTexture(ComplexArrayDesc, TEXT("Disp_x_Array"));
	FRDGTextureRef Disp_y_Array = GraphBuilder.CreateTexture(ComplexArrayDesc, TEXT("Disp_y_Array"));
	FRDGTextureRef DelH_x_Array = GraphBuilder.CreateTexture(ComplexArrayDesc, TEXT("DelH_x_Array"));
	FRDGTextureRef DelH_y_Array = GraphBuilder.CreateTexture(ComplexArrayDesc, TEXT("DelH_y_Array"));

	// Propagate spectral waves
	TShaderMapRef<FPropagateWavesCS> PropagateShader(ShaderMap);
	FPropagateWavesCS::FParameters* PropagateParams = GraphBuilder.AllocParameters<FPropagateWavesCS::FParameters>();
	PropagateParams->FFTWaveConstants = ConstantBuffer;
	PropagateParams->HPosIn = GraphBuilder.CreateSRV(Inputs.HPosIn);
	PropagateParams->HNegIn = GraphBuilder.CreateSRV(Inputs.HNegIn);
	PropagateParams->DispXOut = GraphBuilder.CreateUAV(Disp_x_Array);
	PropagateParams->DispYOut = GraphBuilder.CreateUAV(Disp_y_Array);
	PropagateParams->DelHXOut = GraphBuilder.CreateUAV(DelH_x_Array);
	PropagateParams->DelHYOut = GraphBuilder.CreateUAV(DelH_y_Array);
	PropagateParams->FlowXOut = GraphBuilder.CreateUAV(Outputs.Flow_x_Out);
	PropagateParams->FlowYOut = GraphBuilder.CreateUAV(Outputs.Flow_y_Out);

	FComputeShaderUtils::AddPass(
		GraphBuilder,
		RDG_EVENT_NAME("SWE_FFTWaves_Propagate"),
		ERDGPassFlags::Compute,
		PropagateShader,
		PropagateParams,
		ComplexArrayGroups
	);

	// Inverse FFTs on Disp and DelH 2D Texture Arrays
	Add2DFFTPasses(GraphBuilder, Disp_x_Array, Inputs.PaddedSizeX, Inputs.PaddedSizeY, true, Inputs.DepthNum);
	Add2DFFTPasses(GraphBuilder, Disp_y_Array, Inputs.PaddedSizeX, Inputs.PaddedSizeY, true, Inputs.DepthNum);
	Add2DFFTPasses(GraphBuilder, DelH_x_Array, Inputs.PaddedSizeX, Inputs.PaddedSizeY, true, Inputs.DepthNum);
	Add2DFFTPasses(GraphBuilder, DelH_y_Array, Inputs.PaddedSizeX, Inputs.PaddedSizeY, true, Inputs.DepthNum);

	// Interpolate wind wave outputs between depths
	TShaderMapRef<FInterpCS> InterpShader(ShaderMap);
	FInterpCS::FParameters* InterpParams = GraphBuilder.AllocParameters<FInterpCS::FParameters>();
	InterpParams->FFTWaveConstants = ConstantBuffer;
	InterpParams->HxIn = GraphBuilder.CreateSRV(DelH_x_Array);
	InterpParams->HyIn = GraphBuilder.CreateSRV(DelH_y_Array);
	InterpParams->DxIn = GraphBuilder.CreateSRV(Disp_x_Array);
	InterpParams->DyIn = GraphBuilder.CreateSRV(Disp_y_Array);
	InterpParams->hbarIn = GraphBuilder.CreateSRV(Inputs.hbarIn);
	InterpParams->HxOut = GraphBuilder.CreateUAV(Outputs.delH_x_Out);
	InterpParams->HyOut = GraphBuilder.CreateUAV(Outputs.delH_y_Out);
	InterpParams->DxOut = GraphBuilder.CreateUAV(Outputs.disp_x_Out);
	InterpParams->DyOut = GraphBuilder.CreateUAV(Outputs.disp_y_Out);

	FComputeShaderUtils::AddPass(
		GraphBuilder,
		RDG_EVENT_NAME("SWE_FFTWaves_Interp"),
		ERDGPassFlags::Compute,
		InterpShader,
		InterpParams,
		GridGroups
	);
}
