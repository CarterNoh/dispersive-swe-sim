#include "MainSimShaders.h"
#include "FFTShaders.h"
#include "RenderGraphBuilder.h"

IMPLEMENT_GLOBAL_SHADER_PARAMETER_STRUCT(FSimConstants, "SimConstants");

IMPLEMENT_GLOBAL_SHADER(FInitializeWaterCS,     "/Plugin/DispersiveSWESim/kernels.usf", "InitializeWater",     SF_Compute);
IMPLEMENT_GLOBAL_SHADER(FInitDecompCS,          "/Plugin/DispersiveSWESim/kernels.usf", "InitDecomp",          SF_Compute);
IMPLEMENT_GLOBAL_SHADER(FCalcDiffusionCoeffsCS, "/Plugin/DispersiveSWESim/kernels.usf", "CalcDiffusionCoeffs", SF_Compute);
IMPLEMENT_GLOBAL_SHADER(FDiffusionStepCS,       "/Plugin/DispersiveSWESim/kernels.usf", "DiffusionStep",       SF_Compute);
IMPLEMENT_GLOBAL_SHADER(FDecomposeFieldsCS,     "/Plugin/DispersiveSWESim/kernels.usf", "DecomposeFields",     SF_Compute);
IMPLEMENT_GLOBAL_SHADER(FRecomputeHCS,          "/Plugin/DispersiveSWESim/kernels.usf", "RecomputeH",          SF_Compute);
IMPLEMENT_GLOBAL_SHADER(FTransferToFFTCS,       "/Plugin/DispersiveSWESim/kernels.usf", "TransferToFFT",       SF_Compute);
IMPLEMENT_GLOBAL_SHADER(FCalcEWaveCS,           "/Plugin/DispersiveSWESim/kernels.usf", "CalcEWave",           SF_Compute);
IMPLEMENT_GLOBAL_SHADER(FInterpQCS,             "/Plugin/DispersiveSWESim/kernels.usf", "InterpQ",             SF_Compute);
IMPLEMENT_GLOBAL_SHADER(FCalcUbarCS,            "/Plugin/DispersiveSWESim/kernels.usf", "CalcUbar",            SF_Compute);
IMPLEMENT_GLOBAL_SHADER(FCalcSWECS,             "/Plugin/DispersiveSWESim/kernels.usf", "CalcSWE",             SF_Compute);
IMPLEMENT_GLOBAL_SHADER(FUpdateTildeCS,         "/Plugin/DispersiveSWESim/kernels.usf", "UpdateTilde",         SF_Compute);
IMPLEMENT_GLOBAL_SHADER(FCalcQAdvectCS,         "/Plugin/DispersiveSWESim/kernels.usf", "CalcQAdvect",         SF_Compute);
IMPLEMENT_GLOBAL_SHADER(FIntegrateHCS,          "/Plugin/DispersiveSWESim/kernels.usf", "IntegrateH",          SF_Compute);


void AddInitializeWaterPass(
	FRDGBuilder& GraphBuilder,
	TUniformBufferRef<FSimConstants> ConstantBuffer,
	const FInitializeWaterInputs& Inputs,
	const FInitializeWaterOutputs& Outputs)
{
	FGlobalShaderMap* ShaderMap = GetGlobalShaderMap(GMaxRHIFeatureLevel);
	FIntVector GridGroups(
		FMath::DivideAndRoundUp(Inputs.GridSizeX, 16),
		FMath::DivideAndRoundUp(Inputs.GridSizeY, 16),
		1
	);

	TShaderMapRef<FInitializeWaterCS> InitWaterHeightCS(ShaderMap);
	FInitializeWaterCS::FParameters* InitParams = GraphBuilder.AllocParameters<FInitializeWaterCS::FParameters>();
	InitParams->SimConstants = ConstantBuffer;
	InitParams->WaterLevel = Inputs.WaterLevel;
	InitParams->TerrainCaptureCameraZ = Inputs.TerrainCaptureCameraZ;
	InitParams->terrain = GraphBuilder.CreateSRV(Inputs.TerrainInput);
	InitParams->hOut = GraphBuilder.CreateUAV(Outputs.hOut);
	InitParams->H_Out = GraphBuilder.CreateUAV(Outputs.H_Out);
	InitParams->terrainOut = GraphBuilder.CreateUAV(Outputs.terrainOut);
	InitParams->terrainOutCM = GraphBuilder.CreateUAV(Outputs.terrainOutCM);

	FComputeShaderUtils::AddPass(
		GraphBuilder,
		RDG_EVENT_NAME("SWE_InitializeWater_GPU"),
		ERDGPassFlags::Compute,
		InitWaterHeightCS,
		InitParams,
		GridGroups
	);

	// Copy initial water height h to hbar and hbarOld
	if (Outputs.hbarOut)
	{
		AddCopyTexturePass(GraphBuilder, Outputs.hOut, Outputs.hbarOut);
	}
	if (Outputs.hbarOldOut)
	{
		AddCopyTexturePass(GraphBuilder, Outputs.hOut, Outputs.hbarOldOut);
	}

	// Initialize persistent foam and roughness targets to 0
	if (Outputs.FoamOut)
	{
		AddClearUAVPass(GraphBuilder, GraphBuilder.CreateUAV(Outputs.FoamOut), FLinearColor::Black);
	}
	if (Outputs.RoughnessOut)
	{
		AddClearUAVPass(GraphBuilder, GraphBuilder.CreateUAV(Outputs.RoughnessOut), FLinearColor::Black);
	}
}

void AddDecompositionPasses(
	FRDGBuilder& GraphBuilder,
	TUniformBufferRef<FSimConstants> ConstantBuffer,
	const FDecompositionInputs& Inputs,
	FDecompositionOutputs& Outputs)
{
	FGlobalShaderMap* ShaderMap = GetGlobalShaderMap(GMaxRHIFeatureLevel);
	FIntVector GridGroups(
		FMath::DivideAndRoundUp(Inputs.GridSizeX, 16),
		FMath::DivideAndRoundUp(Inputs.GridSizeY, 16),
		1
	);

	// 1. InitDecomp
	{
		TShaderMapRef<FInitDecompCS> Shader(ShaderMap);
		FInitDecompCS::FParameters* Params = GraphBuilder.AllocParameters<FInitDecompCS::FParameters>();
		Params->SimConstants = ConstantBuffer;
		Params->hIn = GraphBuilder.CreateSRV(Inputs.hIn);
		Params->qIn_x = GraphBuilder.CreateSRV(Inputs.qIn_x);
		Params->qIn_y = GraphBuilder.CreateSRV(Inputs.qIn_y);
		Params->terrain = GraphBuilder.CreateSRV(Inputs.terrain);
		Params->H_Out = GraphBuilder.CreateUAV(Outputs.H_SrcDst);
		Params->Q_Out_x = GraphBuilder.CreateUAV(Outputs.Qx_SrcDst);
		Params->Q_Out_y = GraphBuilder.CreateUAV(Outputs.Qy_SrcDst);

		FComputeShaderUtils::AddPass(
			GraphBuilder,
			RDG_EVENT_NAME("SWE_Decomp_Init"),
			ERDGPassFlags::Compute,
			Shader,
			Params,
			GridGroups
		);
	}

	// 2. Copy initial fields to Orig fields before Jacobi solver
	AddCopyTexturePass(GraphBuilder, Outputs.H_SrcDst, Inputs.HOrig);
	AddCopyTexturePass(GraphBuilder, Outputs.Qx_SrcDst, Inputs.QOrig_x);
	AddCopyTexturePass(GraphBuilder, Outputs.Qy_SrcDst, Inputs.QOrig_y);

	// 3. CalcDiffusionCoeffs
	{
		TShaderMapRef<FCalcDiffusionCoeffsCS> Shader(ShaderMap);
		FCalcDiffusionCoeffsCS::FParameters* Params = GraphBuilder.AllocParameters<FCalcDiffusionCoeffsCS::FParameters>();
		Params->SimConstants = ConstantBuffer;
		Params->H_In = GraphBuilder.CreateSRV(Outputs.H_SrcDst);
		Params->terrain = GraphBuilder.CreateSRV(Inputs.terrain);
		Params->alpha_HOut = GraphBuilder.CreateUAV(Inputs.alpha_H);
		Params->alpha_QOut_x = GraphBuilder.CreateUAV(Inputs.alpha_Q_x);
		Params->alpha_QOut_y = GraphBuilder.CreateUAV(Inputs.alpha_Q_y);

		FComputeShaderUtils::AddPass(
			GraphBuilder,
			RDG_EVENT_NAME("SWE_Decomp_Coeffs"),
			ERDGPassFlags::Compute,
			Shader,
			Params,
			GridGroups
		);
	}

	// 4. Diffusion loop - implicit Jacobi solver ping-ponging across iterations
	{
		FRDGTextureRef H_Src = Outputs.H_SrcDst;
		FRDGTextureRef H_Dst = Outputs.HPast_SrcDst;
		FRDGTextureRef Qx_Src = Outputs.Qx_SrcDst;
		FRDGTextureRef Qx_Dst = Outputs.QPast_x_SrcDst;
		FRDGTextureRef Qy_Src = Outputs.Qy_SrcDst;
		FRDGTextureRef Qy_Dst = Outputs.QPast_y_SrcDst;

		TShaderMapRef<FDiffusionStepCS> Shader(ShaderMap);
		for (int32 j = 0; j < Inputs.DiffusionIterations; j++)
		{
			FDiffusionStepCS::FParameters* Params = GraphBuilder.AllocParameters<FDiffusionStepCS::FParameters>();
			Params->SimConstants = ConstantBuffer;
			Params->terrain = GraphBuilder.CreateSRV(Inputs.terrain);
			Params->H_Orig = GraphBuilder.CreateSRV(Inputs.HOrig);
			Params->Q_Orig_x = GraphBuilder.CreateSRV(Inputs.QOrig_x);
			Params->Q_Orig_y = GraphBuilder.CreateSRV(Inputs.QOrig_y);
			Params->H_Past = GraphBuilder.CreateSRV(H_Src);
			Params->Q_Past_x = GraphBuilder.CreateSRV(Qx_Src);
			Params->Q_Past_y = GraphBuilder.CreateSRV(Qy_Src);
			Params->alpha_HIn = GraphBuilder.CreateSRV(Inputs.alpha_H);
			Params->alpha_QIn_x = GraphBuilder.CreateSRV(Inputs.alpha_Q_x);
			Params->alpha_QIn_y = GraphBuilder.CreateSRV(Inputs.alpha_Q_y);
			Params->H_Out = GraphBuilder.CreateUAV(H_Dst);
			Params->Q_Out_x = GraphBuilder.CreateUAV(Qx_Dst);
			Params->Q_Out_y = GraphBuilder.CreateUAV(Qy_Dst);

			FComputeShaderUtils::AddPass(
				GraphBuilder,
				RDG_EVENT_NAME("SWE_Decomp_Diffusion"),
				ERDGPassFlags::Compute,
				Shader,
				Params,
				GridGroups
			);

			Swap(H_Src, H_Dst);
			Swap(Qx_Src, Qx_Dst);
			Swap(Qy_Src, Qy_Dst);
		}

		Outputs.H_SrcDst = H_Src;
		Outputs.HPast_SrcDst = H_Dst;
		Outputs.Qx_SrcDst = Qx_Src;
		Outputs.QPast_x_SrcDst = Qx_Dst;
		Outputs.Qy_SrcDst = Qy_Src;
		Outputs.QPast_y_SrcDst = Qy_Dst;
	}

	// 5. DecomposeFields
	{
		TShaderMapRef<FDecomposeFieldsCS> Shader(ShaderMap);
		FDecomposeFieldsCS::FParameters* Params = GraphBuilder.AllocParameters<FDecomposeFieldsCS::FParameters>();
		Params->SimConstants = ConstantBuffer;
		Params->H_In = GraphBuilder.CreateSRV(Outputs.H_SrcDst);
		Params->Q_In_x = GraphBuilder.CreateSRV(Outputs.Qx_SrcDst);
		Params->Q_In_y = GraphBuilder.CreateSRV(Outputs.Qy_SrcDst);
		Params->hIn = GraphBuilder.CreateSRV(Inputs.hIn);
		Params->qIn_x = GraphBuilder.CreateSRV(Inputs.qIn_x);
		Params->qIn_y = GraphBuilder.CreateSRV(Inputs.qIn_y);
		Params->terrain = GraphBuilder.CreateSRV(Inputs.terrain);
		Params->hbarOut = GraphBuilder.CreateUAV(Outputs.hbarOut);
		Params->qbarOut_x = GraphBuilder.CreateUAV(Outputs.qbarOut_x);
		Params->qbarOut_y = GraphBuilder.CreateUAV(Outputs.qbarOut_y);
		Params->htildeOut = GraphBuilder.CreateUAV(Outputs.htildeOut);
		Params->qtildeOut_x = GraphBuilder.CreateUAV(Outputs.qtildeOut_x);
		Params->qtildeOut_y = GraphBuilder.CreateUAV(Outputs.qtildeOut_y);

		FComputeShaderUtils::AddPass(
			GraphBuilder,
			RDG_EVENT_NAME("SWE_Decomp_Final"),
			ERDGPassFlags::Compute,
			Shader,
			Params,
			GridGroups
		);
	}

	// 6. Recompute H
	{
		TShaderMapRef<FRecomputeHCS> Shader(ShaderMap);
		FRecomputeHCS::FParameters* Params = GraphBuilder.AllocParameters<FRecomputeHCS::FParameters>();
		Params->SimConstants = ConstantBuffer;
		Params->hIn = GraphBuilder.CreateSRV(Inputs.hIn);
		Params->terrain = GraphBuilder.CreateSRV(Inputs.terrain);
		Params->H_Out = GraphBuilder.CreateUAV(Outputs.H_SrcDst);

		FComputeShaderUtils::AddPass(
			GraphBuilder,
			RDG_EVENT_NAME("SWE_Decomp_ReH"),
			ERDGPassFlags::Compute,
			Shader,
			Params,
			GridGroups
		);
	}
}

void AddEWavePasses(
	FRDGBuilder& GraphBuilder,
	TUniformBufferRef<FSimConstants> ConstantBuffer,
	const FEWaveInputs& Inputs,
	const FEWaveOutputs& Outputs)
{
	FGlobalShaderMap* ShaderMap = GetGlobalShaderMap(GMaxRHIFeatureLevel);
	FIntVector PaddedGroups(
		FMath::DivideAndRoundUp(Inputs.PaddedSizeX, 16),
		FMath::DivideAndRoundUp(Inputs.PaddedSizeY, 16),
		1
	);
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

	// 1. Transfer variables to Fourier domain
	{
		TShaderMapRef<FTransferToFFTCS> Shader(ShaderMap);
		FTransferToFFTCS::FParameters* Params = GraphBuilder.AllocParameters<FTransferToFFTCS::FParameters>();
		Params->SimConstants = ConstantBuffer;
		Params->htildeIn = GraphBuilder.CreateSRV(Inputs.htildeIn);
		Params->htildeOldIn = GraphBuilder.CreateSRV(Inputs.htildeOldIn);
		Params->qtildeIn_x = GraphBuilder.CreateSRV(Inputs.qtildeIn_x);
		Params->qtildeIn_y = GraphBuilder.CreateSRV(Inputs.qtildeIn_y);
		Params->htildeOldNext = GraphBuilder.CreateUAV(Outputs.htildeOldNext);
		Params->hHat = GraphBuilder.CreateUAV(Inputs.hHat);
		Params->qHat_x = GraphBuilder.CreateUAV(Inputs.qHat_x);
		Params->qHat_y = GraphBuilder.CreateUAV(Inputs.qHat_y);

		FComputeShaderUtils::AddPass(
			GraphBuilder,
			RDG_EVENT_NAME("SWE_eWave_Transfer"),
			ERDGPassFlags::Compute,
			Shader,
			Params,
			PaddedGroups
		);
	}

	// 2. Forward 2D FFTs
	Add2DFFTPasses(GraphBuilder, Inputs.hHat, Inputs.PaddedSizeX, Inputs.PaddedSizeY, false, 1);
	Add2DFFTPasses(GraphBuilder, Inputs.qHat_x, Inputs.PaddedSizeX, Inputs.PaddedSizeY, false, 1);
	Add2DFFTPasses(GraphBuilder, Inputs.qHat_y, Inputs.PaddedSizeX, Inputs.PaddedSizeY, false, 1);

	// 3. Compute eWave dispersion updates across depth levels
	{
		TShaderMapRef<FCalcEWaveCS> Shader(ShaderMap);
		FCalcEWaveCS::FParameters* Params = GraphBuilder.AllocParameters<FCalcEWaveCS::FParameters>();
		Params->SimConstants = ConstantBuffer;
		Params->hhat = GraphBuilder.CreateSRV(Inputs.hHat);
		Params->qhat_x = GraphBuilder.CreateSRV(Inputs.qHat_x);
		Params->qhat_y = GraphBuilder.CreateSRV(Inputs.qHat_y);
		Params->Flow_x = GraphBuilder.CreateSRV(Inputs.Flow_x);
		Params->Flow_y = GraphBuilder.CreateSRV(Inputs.Flow_y);
		Params->qhat_x_array = GraphBuilder.CreateUAV(Inputs.qHat_x_array);
		Params->qhat_y_array = GraphBuilder.CreateUAV(Inputs.qHat_y_array);

		FComputeShaderUtils::AddPass(
			GraphBuilder,
			RDG_EVENT_NAME("SWE_eWave_Calc"),
			ERDGPassFlags::Compute,
			Shader,
			Params,
			ComplexArrayGroups
		);
	}

	// 4. Inverse 2D FFTs
	Add2DFFTPasses(GraphBuilder, Inputs.qHat_x_array, Inputs.PaddedSizeX, Inputs.PaddedSizeY, true, Inputs.DepthNum);
	Add2DFFTPasses(GraphBuilder, Inputs.qHat_y_array, Inputs.PaddedSizeX, Inputs.PaddedSizeY, true, Inputs.DepthNum);

	// 5. Interpolate between depths to get new qtilde
	{
		TShaderMapRef<FInterpQCS> Shader(ShaderMap);
		FInterpQCS::FParameters* Params = GraphBuilder.AllocParameters<FInterpQCS::FParameters>();
		Params->SimConstants = ConstantBuffer;
		Params->hbarIn = GraphBuilder.CreateSRV(Inputs.hbarIn);
		Params->qHat_x_array = GraphBuilder.CreateSRV(Inputs.qHat_x_array);
		Params->qHat_y_array = GraphBuilder.CreateSRV(Inputs.qHat_y_array);
		Params->qtildeOut_x = GraphBuilder.CreateUAV(Outputs.qtildeOut_x);
		Params->qtildeOut_y = GraphBuilder.CreateUAV(Outputs.qtildeOut_y);

		FComputeShaderUtils::AddPass(
			GraphBuilder,
			RDG_EVENT_NAME("SWE_eWave_InterpQ"),
			ERDGPassFlags::Compute,
			Shader,
			Params,
			GridGroups
		);
	}
}

void AddSWEBulkPasses(
	FRDGBuilder& GraphBuilder,
	TUniformBufferRef<FSimConstants> ConstantBuffer,
	const FSWEBulkInputs& Inputs,
	const FSWEBulkOutputs& Outputs)
{
	FGlobalShaderMap* ShaderMap = GetGlobalShaderMap(GMaxRHIFeatureLevel);
	FIntVector GridGroups(
		FMath::DivideAndRoundUp(Inputs.GridSizeX, 16),
		FMath::DivideAndRoundUp(Inputs.GridSizeY, 16),
		1
	);

	// 1. CalcUbar
	{
		TShaderMapRef<FCalcUbarCS> Shader(ShaderMap);
		FCalcUbarCS::FParameters* Params = GraphBuilder.AllocParameters<FCalcUbarCS::FParameters>();
		Params->SimConstants = ConstantBuffer;
		Params->qbarIn_x = GraphBuilder.CreateSRV(Inputs.qbarIn_x);
		Params->qbarIn_y = GraphBuilder.CreateSRV(Inputs.qbarIn_y);
		Params->hbarIn = GraphBuilder.CreateSRV(Inputs.hbarOldIn);
		Params->ubarOut_x = GraphBuilder.CreateUAV(Outputs.ubarOut_x);
		Params->ubarOut_y = GraphBuilder.CreateUAV(Outputs.ubarOut_y);

		FComputeShaderUtils::AddPass(
			GraphBuilder,
			RDG_EVENT_NAME("SWE_CalcUbar"),
			ERDGPassFlags::Compute,
			Shader,
			Params,
			GridGroups
		);
	}

	// 2. CalcSWE
	{
		TShaderMapRef<FCalcSWECS> Shader(ShaderMap);
		FCalcSWECS::FParameters* Params = GraphBuilder.AllocParameters<FCalcSWECS::FParameters>();
		Params->SimConstants = ConstantBuffer;
		Params->ubarIn_x = GraphBuilder.CreateSRV(Outputs.ubarOut_x);
		Params->ubarIn_y = GraphBuilder.CreateSRV(Outputs.ubarOut_y);
		Params->hbarIn = GraphBuilder.CreateSRV(Inputs.hbarIn);
		Params->H_In = GraphBuilder.CreateSRV(Inputs.H_In);
		Params->delH_x = GraphBuilder.CreateSRV(Inputs.delH_x);
		Params->delH_y = GraphBuilder.CreateSRV(Inputs.delH_y);
		Params->terrain = GraphBuilder.CreateSRV(Inputs.terrain);
		Params->ubarNewOut_x = GraphBuilder.CreateUAV(Outputs.ubarNewOut_x);
		Params->ubarNewOut_y = GraphBuilder.CreateUAV(Outputs.ubarNewOut_y);
		Params->qbarOut_x = GraphBuilder.CreateUAV(Outputs.qbarOut_x);
		Params->qbarOut_y = GraphBuilder.CreateUAV(Outputs.qbarOut_y);

		FComputeShaderUtils::AddPass(
			GraphBuilder,
			RDG_EVENT_NAME("SWE_CalcSWE"),
			ERDGPassFlags::Compute,
			Shader,
			Params,
			GridGroups
		);
	}
}

void AddTransportAndIntegratePasses(
	FRDGBuilder& GraphBuilder,
	TUniformBufferRef<FSimConstants> ConstantBuffer,
	const FTransportAndIntegrateInputs& Inputs,
	const FTransportAndIntegrateOutputs& Outputs)
{
	FGlobalShaderMap* ShaderMap = GetGlobalShaderMap(GMaxRHIFeatureLevel);
	FIntVector GridGroups(
		FMath::DivideAndRoundUp(Inputs.GridSizeX, 16),
		FMath::DivideAndRoundUp(Inputs.GridSizeY, 16),
		1
	);

	// 1. UpdateTilde (Advect wave height and flow rate)
	{
		TShaderMapRef<FUpdateTildeCS> Shader(ShaderMap);
		FUpdateTildeCS::FParameters* Params = GraphBuilder.AllocParameters<FUpdateTildeCS::FParameters>();
		Params->SimConstants = ConstantBuffer;
		Params->ubarNewIn_x = GraphBuilder.CreateSRV(Inputs.ubarNewIn_x);
		Params->ubarIn_x = GraphBuilder.CreateSRV(Inputs.ubarIn_x);
		Params->ubarNewIn_y = GraphBuilder.CreateSRV(Inputs.ubarNewIn_y);
		Params->ubarIn_y = GraphBuilder.CreateSRV(Inputs.ubarIn_y);
		Params->qtildePast_x = GraphBuilder.CreateSRV(Inputs.qtildePast_x);
		Params->qtildePast_y = GraphBuilder.CreateSRV(Inputs.qtildePast_y);
		Params->hIn = GraphBuilder.CreateSRV(Inputs.hPast);
		Params->htildePast = GraphBuilder.CreateSRV(Inputs.htildePast);
		Params->htildeOut = GraphBuilder.CreateUAV(Outputs.htildeOut);
		Params->qtildeOut_x = GraphBuilder.CreateUAV(Outputs.qtildeOut_x);
		Params->qtildeOut_y = GraphBuilder.CreateUAV(Outputs.qtildeOut_y);

		FComputeShaderUtils::AddPass(
			GraphBuilder,
			RDG_EVENT_NAME("SWE_UpdateTilde"),
			ERDGPassFlags::Compute,
			Shader,
			Params,
			GridGroups
		);
	}

	// 2. CalcQAdvect
	{
		TShaderMapRef<FCalcQAdvectCS> Shader(ShaderMap);
		FCalcQAdvectCS::FParameters* Params = GraphBuilder.AllocParameters<FCalcQAdvectCS::FParameters>();
		Params->SimConstants = ConstantBuffer;
		Params->ubarNewIn_x = GraphBuilder.CreateSRV(Inputs.ubarNewIn_x);
		Params->ubarNewIn_y = GraphBuilder.CreateSRV(Inputs.ubarNewIn_y);
		Params->htildeIn = GraphBuilder.CreateSRV(Outputs.htildeOut);
		Params->qAdvectOut_x = GraphBuilder.CreateUAV(Outputs.qAdvectOut_x);
		Params->qAdvectOut_y = GraphBuilder.CreateUAV(Outputs.qAdvectOut_y);

		FComputeShaderUtils::AddPass(
			GraphBuilder,
			RDG_EVENT_NAME("SWE_CalcQAdvect"),
			ERDGPassFlags::Compute,
			Shader,
			Params,
			GridGroups
		);
	}

	// 3. IntegrateH
	{
		TShaderMapRef<FIntegrateHCS> Shader(ShaderMap);
		FIntegrateHCS::FParameters* Params = GraphBuilder.AllocParameters<FIntegrateHCS::FParameters>();
		Params->SimConstants = ConstantBuffer;
		Params->qbarIn_x = GraphBuilder.CreateSRV(Inputs.qbarIn_x);
		Params->qtildeIn_x = GraphBuilder.CreateSRV(Outputs.qtildeOut_x);
		Params->qAdvectIn_x = GraphBuilder.CreateSRV(Outputs.qAdvectOut_x);
		Params->qbarIn_y = GraphBuilder.CreateSRV(Inputs.qbarIn_y);
		Params->qtildeIn_y = GraphBuilder.CreateSRV(Outputs.qtildeOut_y);
		Params->qAdvectIn_y = GraphBuilder.CreateSRV(Outputs.qAdvectOut_y);
		Params->hPast = GraphBuilder.CreateSRV(Inputs.hPast);
		Params->terrain = GraphBuilder.CreateSRV(Inputs.terrain);
		Params->hOut = GraphBuilder.CreateUAV(Outputs.hOut);
		Params->qOut_x = GraphBuilder.CreateUAV(Outputs.qOut_x);
		Params->qOut_y = GraphBuilder.CreateUAV(Outputs.qOut_y);

		FComputeShaderUtils::AddPass(
			GraphBuilder,
			RDG_EVENT_NAME("SWE_IntegrateH"),
			ERDGPassFlags::Compute,
			Shader,
			Params,
			GridGroups
		);
	}

	// 4. Recompute H (Total free surface elevation for rendering)
	if (Outputs.H_Out)
	{
		TShaderMapRef<FRecomputeHCS> Shader(ShaderMap);
		FRecomputeHCS::FParameters* Params = GraphBuilder.AllocParameters<FRecomputeHCS::FParameters>();
		Params->SimConstants = ConstantBuffer;
		Params->hIn = GraphBuilder.CreateSRV(Outputs.hOut);
		Params->terrain = GraphBuilder.CreateSRV(Inputs.terrain);
		Params->H_Out = GraphBuilder.CreateUAV(Outputs.H_Out);

		FComputeShaderUtils::AddPass(
			GraphBuilder,
			RDG_EVENT_NAME("SWE_Final_ReH"),
			ERDGPassFlags::Compute,
			Shader,
			Params,
			GridGroups
		);
	}
}
