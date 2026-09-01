#include "ExportShaders.h"
#include "RenderGraphBuilder.h"

IMPLEMENT_GLOBAL_SHADER_PARAMETER_STRUCT(FExportConstants, "ExportConstants");

IMPLEMENT_GLOBAL_SHADER(FScaleCopyDisplacementCS,    "/Plugin/DispersiveSWESim/export.usf", "ScaleCopyDisplacement",    SF_Compute);
IMPLEMENT_GLOBAL_SHADER(FCalcSurfaceNormalAndFoamCS, "/Plugin/DispersiveSWESim/export.usf", "CalcSurfaceNormalAndFoam", SF_Compute);
IMPLEMENT_GLOBAL_SHADER(FCalcRoughnessLUTCS,         "/Plugin/DispersiveSWESim/export.usf", "CalcRoughnessLUT",         SF_Compute);

void AddVisualExportPasses(
	FRDGBuilder& GraphBuilder,
	TUniformBufferRef<FExportConstants> ConstantBuffer,
	const FVisualExportInputs& Inputs,
	const FVisualExportOutputs& Outputs)
{
	FGlobalShaderMap* ShaderMap = GetGlobalShaderMap(GMaxRHIFeatureLevel);
	FIntVector GridGroups(
		FMath::DivideAndRoundUp(Inputs.GridSizeX, 16),
		FMath::DivideAndRoundUp(Inputs.GridSizeY, 16),
		1
	);

	// 1. Export Displacement (Pack X, Y, Z/Height into single FloatRGBA target)
	if (Inputs.ExportDispDest && Inputs.inHeight)
	{
		TShaderMapRef<FScaleCopyDisplacementCS> ScaleCopyCS(ShaderMap);
		FScaleCopyDisplacementCS::FParameters* PassParams = GraphBuilder.AllocParameters<FScaleCopyDisplacementCS::FParameters>();
		PassParams->ExportConstants = ConstantBuffer;
		PassParams->ScaleFactor = Inputs.ScaleFactor;
		PassParams->inDispX = GraphBuilder.CreateSRV(Inputs.inDispX);
		PassParams->inDispY = GraphBuilder.CreateSRV(Inputs.inDispY);
		PassParams->inHeight = GraphBuilder.CreateSRV(Inputs.inHeight);
		PassParams->outDisp4 = GraphBuilder.CreateUAV(Inputs.ExportDispDest);

		FComputeShaderUtils::AddPass(
			GraphBuilder,
			RDG_EVENT_NAME("SWE_ExportDisp_Scale"),
			ERDGPassFlags::Compute,
			ScaleCopyCS,
			PassParams,
			GridGroups
		);
	}

	// 2. Export Surface Normal, Foam & Jacobian Determinant
	if (Inputs.ExportNormalDest && Inputs.ExportFoamDest && Inputs.ExportJacobianDest && Outputs.outNewFoam)
	{
		TShaderMapRef<FCalcSurfaceNormalAndFoamCS> NormalAndFoamCS(ShaderMap);
		FCalcSurfaceNormalAndFoamCS::FParameters* PassParams = GraphBuilder.AllocParameters<FCalcSurfaceNormalAndFoamCS::FParameters>();
		PassParams->ExportConstants = ConstantBuffer;
		PassParams->inDispX = GraphBuilder.CreateSRV(Inputs.inDispX);
		PassParams->inDispY = GraphBuilder.CreateSRV(Inputs.inDispY);
		PassParams->inHeight = GraphBuilder.CreateSRV(Inputs.inHeight);
		PassParams->inPreviousFoam = GraphBuilder.CreateSRV(Inputs.inPreviousFoam);
		PassParams->outNormal = GraphBuilder.CreateUAV(Inputs.ExportNormalDest);
		PassParams->outFoam = GraphBuilder.CreateUAV(Outputs.outNewFoam);
		PassParams->outJacobianDet = GraphBuilder.CreateUAV(Inputs.ExportJacobianDest);

		FComputeShaderUtils::AddPass(
			GraphBuilder,
			RDG_EVENT_NAME("SWE_CalcSurfaceNormalAndFoam"),
			ERDGPassFlags::Compute,
			NormalAndFoamCS,
			PassParams,
			GridGroups
		);

		// Copy updated foam to user-facing render target
		AddCopyTexturePass(GraphBuilder, Outputs.outNewFoam, Inputs.ExportFoamDest);
	}

	// 3. Calculate Roughness LUT
	if (Inputs.ExportRoughnessDest && Inputs.ExportNormalDest && Outputs.outNewRoughness && Inputs.inPreviousRoughness)
	{
		TShaderMapRef<FCalcRoughnessLUTCS> RoughnessCS(ShaderMap);
		FCalcRoughnessLUTCS::FParameters* PassParams = GraphBuilder.AllocParameters<FCalcRoughnessLUTCS::FParameters>();
		PassParams->ExportConstants = ConstantBuffer;
		PassParams->IntegrationSamples = Inputs.IntegrationSamples;
		PassParams->RoughnessPower = Inputs.RoughnessPower;
		PassParams->inNormal = GraphBuilder.CreateSRV(Inputs.ExportNormalDest);
		PassParams->inPreviousRoughness = GraphBuilder.CreateSRV(Inputs.inPreviousRoughness);
		PassParams->BilinearWrapSampler = TStaticSamplerState<SF_Bilinear, AM_Wrap, AM_Wrap, AM_Wrap>::GetRHI();
		PassParams->outRoughness = GraphBuilder.CreateUAV(Outputs.outNewRoughness);

		FComputeShaderUtils::AddPass(
			GraphBuilder,
			RDG_EVENT_NAME("SWE_CalcRoughnessLUT"),
			ERDGPassFlags::Compute,
			RoughnessCS,
			PassParams,
			FIntVector(FMath::DivideAndRoundUp(Inputs.GridSizeX, 16), 1, 1)
		);

		AddCopyTexturePass(GraphBuilder, Outputs.outNewRoughness, Inputs.ExportRoughnessDest);
	}
}
