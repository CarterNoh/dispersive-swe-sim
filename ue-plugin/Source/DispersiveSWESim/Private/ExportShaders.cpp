#include "ExportShaders.h"
#include "RenderGraphBuilder.h"

IMPLEMENT_GLOBAL_SHADER_PARAMETER_STRUCT(FExportConstants, "ExportConstants");

IMPLEMENT_GLOBAL_SHADER(FScaleCopyDisplacementCS,    "/Plugin/DispersiveSWESim/export.usf", "ScaleCopyDisplacement",    SF_Compute);
IMPLEMENT_GLOBAL_SHADER(FScaleCopyVelocityCS,        "/Plugin/DispersiveSWESim/export.usf", "ScaleCopyVelocity",        SF_Compute);
IMPLEMENT_GLOBAL_SHADER(FCalcAccelerationCS,        "/Plugin/DispersiveSWESim/export.usf", "CalcAcceleration",        SF_Compute);
IMPLEMENT_GLOBAL_SHADER(FCalcSurfaceNormalAndFoamCS, "/Plugin/DispersiveSWESim/export.usf", "CalcSurfaceNormalAndFoam", SF_Compute);
IMPLEMENT_GLOBAL_SHADER(FCalcRoughnessLUTCS,         "/Plugin/DispersiveSWESim/export.usf", "CalcRoughnessLUT",         SF_Compute);

void AddVisualExportPasses(
	FRDGBuilder& GraphBuilder,
	TUniformBufferRef<FExportConstants> ConstantBuffer,
	const FVisualExportInputs& Inputs,
	const FVisualExportOutputs& Outputs) {
	FGlobalShaderMap* ShaderMap = GetGlobalShaderMap(GMaxRHIFeatureLevel);
	FIntVector GridGroups(
		FMath::DivideAndRoundUp(Inputs.GridSizeX, 16),
		FMath::DivideAndRoundUp(Inputs.GridSizeY, 16),
		1
	);

	// Export Displacement (Pack X, Y, Z/Height into single FloatRGBA target)
	if (Outputs.ExportDispDest && Inputs.inHeight) {
		TShaderMapRef<FScaleCopyDisplacementCS> ScaleCopyCS(ShaderMap);
		FScaleCopyDisplacementCS::FParameters* PassParams = GraphBuilder.AllocParameters<FScaleCopyDisplacementCS::FParameters>();
		PassParams->ExportConstants = ConstantBuffer;
		PassParams->ScaleFactor = Inputs.ScaleFactor;
		PassParams->inDispX = GraphBuilder.CreateSRV(Inputs.inDispX);
		PassParams->inDispY = GraphBuilder.CreateSRV(Inputs.inDispY);
		PassParams->inHeight = GraphBuilder.CreateSRV(Inputs.inHeight);
		PassParams->outDisp4 = GraphBuilder.CreateUAV(Outputs.ExportDispDest);

		FComputeShaderUtils::AddPass(
			GraphBuilder,
			RDG_EVENT_NAME("SWE_ExportDisp_Scale"),
			ERDGPassFlags::Compute,
			ScaleCopyCS,
			PassParams,
			GridGroups
		);
	}

	// Export Fluid Velocity (Pack u_x, u_y, w into single FloatRGBA target in m/s)
	if (Outputs.ExportVelocityDest && Inputs.inQx && Inputs.inQy && Inputs.inHPast && Inputs.inHNew && Inputs.inHdot) {
		TShaderMapRef<FScaleCopyVelocityCS> VelocityCS(ShaderMap);
		FScaleCopyVelocityCS::FParameters* PassParams = GraphBuilder.AllocParameters<FScaleCopyVelocityCS::FParameters>();
		PassParams->ExportConstants = ConstantBuffer;
		PassParams->inQx = GraphBuilder.CreateSRV(Inputs.inQx);
		PassParams->inQy = GraphBuilder.CreateSRV(Inputs.inQy);
		PassParams->inHPast = GraphBuilder.CreateSRV(Inputs.inHPast);
		PassParams->inHNew = GraphBuilder.CreateSRV(Inputs.inHNew);
		PassParams->inHdot = GraphBuilder.CreateSRV(Inputs.inHdot);
		PassParams->outVel4 = GraphBuilder.CreateUAV(Outputs.ExportVelocityDest);

		FComputeShaderUtils::AddPass(
			GraphBuilder,
			RDG_EVENT_NAME("SWE_ExportVelocity"),
			ERDGPassFlags::Compute,
			VelocityCS,
			PassParams,
			GridGroups
		);
	}

	// Export Fluid Acceleration (Forward Euler from current and past velocity in m/s^2)
	if (Outputs.ExportAccelDest && Inputs.inVel && Inputs.inVelPast) {
		TShaderMapRef<FCalcAccelerationCS> AccelCS(ShaderMap);
		FCalcAccelerationCS::FParameters* PassParams = GraphBuilder.AllocParameters<FCalcAccelerationCS::FParameters>();
		PassParams->ExportConstants = ConstantBuffer;
		PassParams->inVel = GraphBuilder.CreateSRV(Inputs.inVel);
		PassParams->inVelPast = GraphBuilder.CreateSRV(Inputs.inVelPast);
		PassParams->outAccel4 = GraphBuilder.CreateUAV(Outputs.ExportAccelDest);

		FComputeShaderUtils::AddPass(
			GraphBuilder,
			RDG_EVENT_NAME("SWE_ExportAcceleration"),
			ERDGPassFlags::Compute,
			AccelCS,
			PassParams,
			GridGroups
		);
	}

	// Export Surface Normal, Foam & Jacobian Determinant
	if (Outputs.outNewFoam && Inputs.inPreviousFoam) {
		FRDGTextureRef NormalDest = Outputs.ExportNormalDest;
		if (!NormalDest) {
			NormalDest = GraphBuilder.CreateTexture(FRDGTextureDesc::Create2D(FIntPoint(Inputs.GridSizeX, Inputs.GridSizeY), PF_FloatRGBA, FClearValueBinding::None, TexCreate_ShaderResource | TexCreate_UAV), TEXT("NormalTemp"));
		}
		FRDGTextureRef JacobianDest = Outputs.ExportJacobianDest;
		if (!JacobianDest) {
			JacobianDest = GraphBuilder.CreateTexture(FRDGTextureDesc::Create2D(FIntPoint(Inputs.GridSizeX, Inputs.GridSizeY), PF_FloatRGBA, FClearValueBinding::None, TexCreate_ShaderResource | TexCreate_UAV), TEXT("JacobianTemp"));
		}

		TShaderMapRef<FCalcSurfaceNormalAndFoamCS> NormalAndFoamCS(ShaderMap);
		FCalcSurfaceNormalAndFoamCS::FParameters* PassParams = GraphBuilder.AllocParameters<FCalcSurfaceNormalAndFoamCS::FParameters>();
		PassParams->ExportConstants = ConstantBuffer;
		PassParams->inDispX = GraphBuilder.CreateSRV(Inputs.inDispX);
		PassParams->inDispY = GraphBuilder.CreateSRV(Inputs.inDispY);
		PassParams->inHeight = GraphBuilder.CreateSRV(Inputs.inHeight);
		PassParams->inPreviousFoam = GraphBuilder.CreateSRV(Inputs.inPreviousFoam);
		PassParams->outNormal = GraphBuilder.CreateUAV(NormalDest);
		PassParams->outFoam = GraphBuilder.CreateUAV(Outputs.outNewFoam);
		PassParams->outJacobianDet = GraphBuilder.CreateUAV(JacobianDest);

		FComputeShaderUtils::AddPass(
			GraphBuilder,
			RDG_EVENT_NAME("SWE_CalcSurfaceNormalAndFoam"),
			ERDGPassFlags::Compute,
			NormalAndFoamCS,
			PassParams,
			GridGroups
		);

		if (Outputs.ExportFoamDest) {
			AddCopyTexturePass(GraphBuilder, Outputs.outNewFoam, Outputs.ExportFoamDest);
		}
	}

	// Calculate Roughness LUT
	if (Outputs.outNewRoughness && Inputs.inPreviousRoughness && Outputs.ExportNormalDest) {
		TShaderMapRef<FCalcRoughnessLUTCS> RoughnessCS(ShaderMap);
		FCalcRoughnessLUTCS::FParameters* PassParams = GraphBuilder.AllocParameters<FCalcRoughnessLUTCS::FParameters>();
		PassParams->ExportConstants = ConstantBuffer;
		PassParams->IntegrationSamples = Inputs.IntegrationSamples;
		PassParams->RoughnessPower = Inputs.RoughnessPower;
		PassParams->inNormal = GraphBuilder.CreateSRV(Outputs.ExportNormalDest);
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

		if (Outputs.ExportRoughnessDest) {
			AddCopyTexturePass(GraphBuilder, Outputs.outNewRoughness, Outputs.ExportRoughnessDest);
		}
	}
}
