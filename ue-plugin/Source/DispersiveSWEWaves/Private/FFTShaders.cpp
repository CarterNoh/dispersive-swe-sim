#include "FFTShaders.h"
#include "RenderGraphBuilder.h"

IMPLEMENT_GLOBAL_SHADER(FFFTKernel1DCS, "/Plugin/DispersiveSWEWaves/fft.usf", "FFTKernel_1D", SF_Compute);

void Add2DFFTPasses(
	FRDGBuilder& GraphBuilder,
	FRDGTextureRef TargetTexture,
	int32 SizeX,
	int32 SizeY,
	bool bInverse,
	int32 NumLayers) {
	uint32 BitsX = 0, BitsY = 0;
	int32 TmpX = SizeX;
	while (TmpX >>= 1) ++BitsX;
	int32 TmpY = SizeY;
	while (TmpY >>= 1) ++BitsY;

	bool bIsArray = NumLayers > 1;

	// 1. Row Pass: Transform each row
	FFFTKernel1DCS::FPermutationDomain RowPermutationVector;
	RowPermutationVector.Set<FFFTKernel1DCS::FFFTSizeDim>(SizeX);
	RowPermutationVector.Set<FFFTKernel1DCS::FIsArrayDim>(bIsArray);

	TShaderMapRef<FFFTKernel1DCS> RowFFTShader(GetGlobalShaderMap(GMaxRHIFeatureLevel), RowPermutationVector);

	FFFTKernel1DCS::FParameters* RowPassParams = GraphBuilder.AllocParameters<FFFTKernel1DCS::FParameters>();
	RowPassParams->cb_Nx = SizeX;
	RowPassParams->cb_Ny = SizeY;
	RowPassParams->cb_BitsX = BitsX;
	RowPassParams->cb_BitsY = BitsY;
	RowPassParams->cb_Inverse = bInverse ? 1 : 0;
	RowPassParams->cb_IsRow = 1;
	RowPassParams->fft = GraphBuilder.CreateUAV(TargetTexture);

	FComputeShaderUtils::AddPass(
		GraphBuilder,
		RDG_EVENT_NAME("FFT_RowPass"),
		ERDGPassFlags::Compute,
		RowFFTShader,
		RowPassParams,
		FIntVector(1, SizeY, NumLayers)
	);

	// 2. Column Pass: Transform each column
	FFFTKernel1DCS::FPermutationDomain ColPermutationVector;
	ColPermutationVector.Set<FFFTKernel1DCS::FFFTSizeDim>(SizeY);
	ColPermutationVector.Set<FFFTKernel1DCS::FIsArrayDim>(bIsArray);

	TShaderMapRef<FFFTKernel1DCS> ColFFTShader(GetGlobalShaderMap(GMaxRHIFeatureLevel), ColPermutationVector);

	FFFTKernel1DCS::FParameters* ColPassParams = GraphBuilder.AllocParameters<FFFTKernel1DCS::FParameters>();
	ColPassParams->cb_Nx = SizeX;
	ColPassParams->cb_Ny = SizeY;
	ColPassParams->cb_BitsX = BitsX;
	ColPassParams->cb_BitsY = BitsY;
	ColPassParams->cb_Inverse = bInverse ? 1 : 0;
	ColPassParams->cb_IsRow = 0;
	ColPassParams->fft = GraphBuilder.CreateUAV(TargetTexture);

	FComputeShaderUtils::AddPass(
		GraphBuilder,
		RDG_EVENT_NAME("FFT_ColPass"),
		ERDGPassFlags::Compute,
		ColFFTShader,
		ColPassParams,
		FIntVector(1, SizeX, NumLayers)
	);
}
