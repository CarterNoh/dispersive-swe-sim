#include "MainSimShaders.h"
#include "FFTShaders.h"
#include "RenderGraphBuilder.h"

IMPLEMENT_GLOBAL_SHADER_PARAMETER_STRUCT(FSimConstants, "SimConstants");

IMPLEMENT_GLOBAL_SHADER(FInitializeWaterCS,     "/Plugin/DispersiveSWEWaves/dispersiveswe.usf", "InitializeWater",     SF_Compute);
IMPLEMENT_GLOBAL_SHADER(FInitDecompCS,          "/Plugin/DispersiveSWEWaves/dispersiveswe.usf", "InitDecomp",          SF_Compute);
IMPLEMENT_GLOBAL_SHADER(FCalcDiffusionCoeffsCS, "/Plugin/DispersiveSWEWaves/dispersiveswe.usf", "CalcDiffusionCoeffs", SF_Compute);
IMPLEMENT_GLOBAL_SHADER(FDiffusionStepCS,       "/Plugin/DispersiveSWEWaves/dispersiveswe.usf", "DiffusionStep",       SF_Compute);
IMPLEMENT_GLOBAL_SHADER(FDecomposeFieldsCS,     "/Plugin/DispersiveSWEWaves/dispersiveswe.usf", "DecomposeFields",     SF_Compute);
IMPLEMENT_GLOBAL_SHADER(FRecomputeHCS,          "/Plugin/DispersiveSWEWaves/dispersiveswe.usf", "RecomputeH",          SF_Compute);
IMPLEMENT_GLOBAL_SHADER(FRecomputeHAvgCS,       "/Plugin/DispersiveSWEWaves/dispersiveswe.usf", "RecomputeHAvg",       SF_Compute);
IMPLEMENT_GLOBAL_SHADER(FTransferToFFTCS,       "/Plugin/DispersiveSWEWaves/dispersiveswe.usf", "TransferToFFT",       SF_Compute);
IMPLEMENT_GLOBAL_SHADER(FCalcEWaveCS,           "/Plugin/DispersiveSWEWaves/dispersiveswe.usf", "CalcEWave",           SF_Compute);
IMPLEMENT_GLOBAL_SHADER(FInterpQCS,             "/Plugin/DispersiveSWEWaves/dispersiveswe.usf", "InterpQ",             SF_Compute);
IMPLEMENT_GLOBAL_SHADER(FCalcUbarCS,            "/Plugin/DispersiveSWEWaves/dispersiveswe.usf", "CalcUbar",            SF_Compute);
IMPLEMENT_GLOBAL_SHADER(FCalcSWECS,             "/Plugin/DispersiveSWEWaves/dispersiveswe.usf", "CalcSWE",             SF_Compute);
IMPLEMENT_GLOBAL_SHADER(FUpdateTildeCS,         "/Plugin/DispersiveSWEWaves/dispersiveswe.usf", "UpdateTilde",         SF_Compute);
IMPLEMENT_GLOBAL_SHADER(FCalcQAdvectCS,         "/Plugin/DispersiveSWEWaves/dispersiveswe.usf", "CalcQAdvect",         SF_Compute);
IMPLEMENT_GLOBAL_SHADER(FIntegrateHCS,          "/Plugin/DispersiveSWEWaves/dispersiveswe.usf", "IntegrateH",          SF_Compute);


void AddInitializeWaterPass(
	FRDGBuilder& GraphBuilder,
	TUniformBufferRef<FSimConstants> ConstantBuffer,
	const FInitializeWaterInputs& Inputs,
	const FInitializeWaterOutputs& Outputs) {
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
	if (Outputs.hbarOut) {
		AddCopyTexturePass(GraphBuilder, Outputs.hOut, Outputs.hbarOut);
	}
	if (Outputs.hbarOldOut) {
		AddCopyTexturePass(GraphBuilder, Outputs.hOut, Outputs.hbarOldOut);
	}

	// Initialize persistent foam and roughness targets to 0
	if (Outputs.FoamOut) {
		AddClearUAVPass(GraphBuilder, GraphBuilder.CreateUAV(Outputs.FoamOut), FLinearColor::Black);
	}
	if (Outputs.RoughnessOut) {
		AddClearUAVPass(GraphBuilder, GraphBuilder.CreateUAV(Outputs.RoughnessOut), FLinearColor::Black);
	}
}

void AddDecompositionPasses(
	FRDGBuilder& GraphBuilder,
	TUniformBufferRef<FSimConstants> ConstantBuffer,
	const FDecompositionInputs& Inputs,
	FDecompositionOutputs& Outputs) {
	FGlobalShaderMap* ShaderMap = GetGlobalShaderMap(GMaxRHIFeatureLevel);
	FIntVector GridGroups(
		FMath::DivideAndRoundUp(Inputs.GridSizeX, 16),
		FMath::DivideAndRoundUp(Inputs.GridSizeY, 16),
		1
	);

	// Transient textures for decomposition and Jacobi diffusion
	FRDGTextureDesc GridFloatDesc = FRDGTextureDesc::Create2D(
		FIntPoint(Inputs.GridSizeX, Inputs.GridSizeY),
		PF_R32_FLOAT,
		FClearValueBinding::None,
		TexCreate_ShaderResource | TexCreate_UAV
	);

	FRDGTextureRef HOrig = GraphBuilder.CreateTexture(GridFloatDesc, TEXT("HOrig"));
	FRDGTextureRef QOrig_x = GraphBuilder.CreateTexture(GridFloatDesc, TEXT("QOrig_x"));
	FRDGTextureRef QOrig_y = GraphBuilder.CreateTexture(GridFloatDesc, TEXT("QOrig_y"));
	FRDGTextureRef alpha_H = GraphBuilder.CreateTexture(GridFloatDesc, TEXT("alpha_H"));
	FRDGTextureRef alpha_Q_x = GraphBuilder.CreateTexture(GridFloatDesc, TEXT("alpha_Q_x"));
	FRDGTextureRef alpha_Q_y = GraphBuilder.CreateTexture(GridFloatDesc, TEXT("alpha_Q_y"));

	// Initialize Decomposition Fields
	TShaderMapRef<FInitDecompCS> InitDecompShader(ShaderMap);
	FInitDecompCS::FParameters* InitParams = GraphBuilder.AllocParameters<FInitDecompCS::FParameters>();
	InitParams->SimConstants = ConstantBuffer;
	InitParams->hIn = GraphBuilder.CreateSRV(Inputs.hIn);
	InitParams->qIn_x = GraphBuilder.CreateSRV(Inputs.qIn_x);
	InitParams->qIn_y = GraphBuilder.CreateSRV(Inputs.qIn_y);
	InitParams->terrain = GraphBuilder.CreateSRV(Inputs.terrain);
	InitParams->H_Out = GraphBuilder.CreateUAV(Outputs.H_SrcDst);
	InitParams->Q_Out_x = GraphBuilder.CreateUAV(Outputs.Qx_SrcDst);
	InitParams->Q_Out_y = GraphBuilder.CreateUAV(Outputs.Qy_SrcDst);

	FComputeShaderUtils::AddPass(
		GraphBuilder,
		RDG_EVENT_NAME("SWE_Decomp_Init"),
		ERDGPassFlags::Compute,
		InitDecompShader,
		InitParams,
		GridGroups
	);

	// Copy initial fields to Orig fields before Jacobi solver
	AddCopyTexturePass(GraphBuilder, Outputs.H_SrcDst, HOrig);
	AddCopyTexturePass(GraphBuilder, Outputs.Qx_SrcDst, QOrig_x);
	AddCopyTexturePass(GraphBuilder, Outputs.Qy_SrcDst, QOrig_y);

	// Calculate Diffusion Coefficients
	TShaderMapRef<FCalcDiffusionCoeffsCS> CoeffsShader(ShaderMap);
	FCalcDiffusionCoeffsCS::FParameters* CoeffsParams = GraphBuilder.AllocParameters<FCalcDiffusionCoeffsCS::FParameters>();
	CoeffsParams->SimConstants = ConstantBuffer;
	CoeffsParams->H_In = GraphBuilder.CreateSRV(Outputs.H_SrcDst);
	CoeffsParams->terrain = GraphBuilder.CreateSRV(Inputs.terrain);
	CoeffsParams->alpha_HOut = GraphBuilder.CreateUAV(alpha_H);
	CoeffsParams->alpha_QOut_x = GraphBuilder.CreateUAV(alpha_Q_x);
	CoeffsParams->alpha_QOut_y = GraphBuilder.CreateUAV(alpha_Q_y);

	FComputeShaderUtils::AddPass(
		GraphBuilder,
		RDG_EVENT_NAME("SWE_Decomp_Coeffs"),
		ERDGPassFlags::Compute,
		CoeffsShader,
		CoeffsParams,
		GridGroups
	);

	// Diffusion loop - implicit Jacobi solver ping-ponging across iterations
	FRDGTextureRef H_Src = Outputs.H_SrcDst;
	FRDGTextureRef H_Dst = Outputs.HPast_SrcDst;
	FRDGTextureRef Qx_Src = Outputs.Qx_SrcDst;
	FRDGTextureRef Qx_Dst = Outputs.QPast_x_SrcDst;
	FRDGTextureRef Qy_Src = Outputs.Qy_SrcDst;
	FRDGTextureRef Qy_Dst = Outputs.QPast_y_SrcDst;

	TShaderMapRef<FDiffusionStepCS> DiffusionShader(ShaderMap);
	for (int32 j = 0; j < Inputs.DiffusionIterations; j++) {
		FDiffusionStepCS::FParameters* Params = GraphBuilder.AllocParameters<FDiffusionStepCS::FParameters>();
		Params->SimConstants = ConstantBuffer;
		Params->terrain = GraphBuilder.CreateSRV(Inputs.terrain);
		Params->H_Orig = GraphBuilder.CreateSRV(HOrig);
		Params->Q_Orig_x = GraphBuilder.CreateSRV(QOrig_x);
		Params->Q_Orig_y = GraphBuilder.CreateSRV(QOrig_y);
		Params->H_Past = GraphBuilder.CreateSRV(H_Src);
		Params->Q_Past_x = GraphBuilder.CreateSRV(Qx_Src);
		Params->Q_Past_y = GraphBuilder.CreateSRV(Qy_Src);
		Params->alpha_HIn = GraphBuilder.CreateSRV(alpha_H);
		Params->alpha_QIn_x = GraphBuilder.CreateSRV(alpha_Q_x);
		Params->alpha_QIn_y = GraphBuilder.CreateSRV(alpha_Q_y);
		Params->H_Out = GraphBuilder.CreateUAV(H_Dst);
		Params->Q_Out_x = GraphBuilder.CreateUAV(Qx_Dst);
		Params->Q_Out_y = GraphBuilder.CreateUAV(Qy_Dst);

		FComputeShaderUtils::AddPass(
			GraphBuilder,
			RDG_EVENT_NAME("SWE_Decomp_Diffusion"),
			ERDGPassFlags::Compute,
			DiffusionShader,
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

	// Decompose Fields into low-frequency and high-frequency wave parts
	TShaderMapRef<FDecomposeFieldsCS> DecompShader(ShaderMap);
	FDecomposeFieldsCS::FParameters* DecompParams = GraphBuilder.AllocParameters<FDecomposeFieldsCS::FParameters>();
	DecompParams->SimConstants = ConstantBuffer;
	DecompParams->H_In = GraphBuilder.CreateSRV(Outputs.H_SrcDst);
	DecompParams->Q_In_x = GraphBuilder.CreateSRV(Outputs.Qx_SrcDst);
	DecompParams->Q_In_y = GraphBuilder.CreateSRV(Outputs.Qy_SrcDst);
	DecompParams->hIn = GraphBuilder.CreateSRV(Inputs.hIn);
	DecompParams->qIn_x = GraphBuilder.CreateSRV(Inputs.qIn_x);
	DecompParams->qIn_y = GraphBuilder.CreateSRV(Inputs.qIn_y);
	DecompParams->terrain = GraphBuilder.CreateSRV(Inputs.terrain);
	DecompParams->hbarOut = GraphBuilder.CreateUAV(Outputs.hbarOut);
	DecompParams->qbarOut_x = GraphBuilder.CreateUAV(Outputs.qbarOut_x);
	DecompParams->qbarOut_y = GraphBuilder.CreateUAV(Outputs.qbarOut_y);
	DecompParams->htildeOut = GraphBuilder.CreateUAV(Outputs.htildeOut);
	DecompParams->qtildeOut_x = GraphBuilder.CreateUAV(Outputs.qtildeOut_x);
	DecompParams->qtildeOut_y = GraphBuilder.CreateUAV(Outputs.qtildeOut_y);

	FComputeShaderUtils::AddPass(
		GraphBuilder,
		RDG_EVENT_NAME("SWE_Decomp_Final"),
		ERDGPassFlags::Compute,
		DecompShader,
		DecompParams,
		GridGroups
	);

	// Recompute H
	TShaderMapRef<FRecomputeHCS> RecomputeHShader(ShaderMap);
	FRecomputeHCS::FParameters* RecomputeHParams = GraphBuilder.AllocParameters<FRecomputeHCS::FParameters>();
	RecomputeHParams->SimConstants = ConstantBuffer;
	RecomputeHParams->hIn = GraphBuilder.CreateSRV(Inputs.hIn);
	RecomputeHParams->terrain = GraphBuilder.CreateSRV(Inputs.terrain);
	RecomputeHParams->H_Out = GraphBuilder.CreateUAV(Outputs.H_SrcDst);

	FComputeShaderUtils::AddPass(
		GraphBuilder,
		RDG_EVENT_NAME("SWE_Decomp_ReH"),
		ERDGPassFlags::Compute,
		RecomputeHShader,
		RecomputeHParams,
		GridGroups
	);
}

void AddEWavePasses(
	FRDGBuilder& GraphBuilder,
	TUniformBufferRef<FSimConstants> ConstantBuffer,
	const FEWaveInputs& Inputs,
	const FEWaveOutputs& Outputs) {
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

	// Transient complex Fourier buffers
	FRDGTextureDesc FloatPaddedDesc = FRDGTextureDesc::Create2D(
		FIntPoint(Inputs.PaddedSizeX, Inputs.PaddedSizeY),
		PF_G32R32F,
		FClearValueBinding::None,
		TexCreate_ShaderResource | TexCreate_UAV
	);

	FRDGTextureDesc ComplexArrayDesc = FRDGTextureDesc::Create2DArray(
		FIntPoint(Inputs.PaddedSizeX, Inputs.PaddedSizeY),
		PF_G32R32F,
		FClearValueBinding::None,
		TexCreate_ShaderResource | TexCreate_UAV,
		Inputs.DepthNum
	);

	FRDGTextureRef hHat = GraphBuilder.CreateTexture(FloatPaddedDesc, TEXT("hHat"));
	FRDGTextureRef qHat_x = GraphBuilder.CreateTexture(FloatPaddedDesc, TEXT("qHat_x"));
	FRDGTextureRef qHat_y = GraphBuilder.CreateTexture(FloatPaddedDesc, TEXT("qHat_y"));
	FRDGTextureRef qHat_x_array = GraphBuilder.CreateTexture(ComplexArrayDesc, TEXT("qHat_x_array"));
	FRDGTextureRef qHat_y_array = GraphBuilder.CreateTexture(ComplexArrayDesc, TEXT("qHat_y_array"));

	// Transfer variables to Fourier domain
	TShaderMapRef<FTransferToFFTCS> TransferShader(ShaderMap);
	FTransferToFFTCS::FParameters* TransferParams = GraphBuilder.AllocParameters<FTransferToFFTCS::FParameters>();
	TransferParams->SimConstants = ConstantBuffer;
	TransferParams->htildeIn = GraphBuilder.CreateSRV(Inputs.htildeIn);
	TransferParams->htildeOldIn = GraphBuilder.CreateSRV(Inputs.htildeOldIn);
	TransferParams->qtildeIn_x = GraphBuilder.CreateSRV(Inputs.qtildeIn_x);
	TransferParams->qtildeIn_y = GraphBuilder.CreateSRV(Inputs.qtildeIn_y);
	TransferParams->htildeOldNext = GraphBuilder.CreateUAV(Outputs.htildeOldNext);
	TransferParams->hHat = GraphBuilder.CreateUAV(hHat);
	TransferParams->qHat_x = GraphBuilder.CreateUAV(qHat_x);
	TransferParams->qHat_y = GraphBuilder.CreateUAV(qHat_y);

	FComputeShaderUtils::AddPass(
		GraphBuilder,
		RDG_EVENT_NAME("SWE_eWave_Transfer"),
		ERDGPassFlags::Compute,
		TransferShader,
		TransferParams,
		PaddedGroups
	);

	// Forward 2D FFTs
	Add2DFFTPasses(GraphBuilder, hHat, Inputs.PaddedSizeX, Inputs.PaddedSizeY, false, 1);
	Add2DFFTPasses(GraphBuilder, qHat_x, Inputs.PaddedSizeX, Inputs.PaddedSizeY, false, 1);
	Add2DFFTPasses(GraphBuilder, qHat_y, Inputs.PaddedSizeX, Inputs.PaddedSizeY, false, 1);

	// Compute eWave dispersion updates across depth levels
	TShaderMapRef<FCalcEWaveCS> CalcEWaveShader(ShaderMap);
	FCalcEWaveCS::FParameters* CalcEWaveParams = GraphBuilder.AllocParameters<FCalcEWaveCS::FParameters>();
	CalcEWaveParams->SimConstants = ConstantBuffer;
	CalcEWaveParams->hhat = GraphBuilder.CreateSRV(hHat);
	CalcEWaveParams->qhat_x = GraphBuilder.CreateSRV(qHat_x);
	CalcEWaveParams->qhat_y = GraphBuilder.CreateSRV(qHat_y);
	CalcEWaveParams->Flow_x = GraphBuilder.CreateSRV(Inputs.Flow_x);
	CalcEWaveParams->Flow_y = GraphBuilder.CreateSRV(Inputs.Flow_y);
	CalcEWaveParams->qhat_x_array = GraphBuilder.CreateUAV(qHat_x_array);
	CalcEWaveParams->qhat_y_array = GraphBuilder.CreateUAV(qHat_y_array);

	FComputeShaderUtils::AddPass(
		GraphBuilder,
		RDG_EVENT_NAME("SWE_eWave_Calc"),
		ERDGPassFlags::Compute,
		CalcEWaveShader,
		CalcEWaveParams,
		ComplexArrayGroups
	);

	// Inverse 2D FFTs
	Add2DFFTPasses(GraphBuilder, qHat_x_array, Inputs.PaddedSizeX, Inputs.PaddedSizeY, true, Inputs.DepthNum);
	Add2DFFTPasses(GraphBuilder, qHat_y_array, Inputs.PaddedSizeX, Inputs.PaddedSizeY, true, Inputs.DepthNum);

	// Interpolate between depths to get new qtilde
	TShaderMapRef<FInterpQCS> InterpQShader(ShaderMap);
	FInterpQCS::FParameters* InterpQParams = GraphBuilder.AllocParameters<FInterpQCS::FParameters>();
	InterpQParams->SimConstants = ConstantBuffer;
	InterpQParams->hbarIn = GraphBuilder.CreateSRV(Inputs.hbarIn);
	InterpQParams->qHat_x_array = GraphBuilder.CreateSRV(qHat_x_array);
	InterpQParams->qHat_y_array = GraphBuilder.CreateSRV(qHat_y_array);
	InterpQParams->qtildeOut_x = GraphBuilder.CreateUAV(Outputs.qtildeOut_x);
	InterpQParams->qtildeOut_y = GraphBuilder.CreateUAV(Outputs.qtildeOut_y);

	FComputeShaderUtils::AddPass(
		GraphBuilder,
		RDG_EVENT_NAME("SWE_eWave_InterpQ"),
		ERDGPassFlags::Compute,
		InterpQShader,
		InterpQParams,
		GridGroups
	);
}

void AddSWEBulkPasses(
	FRDGBuilder& GraphBuilder,
	TUniformBufferRef<FSimConstants> ConstantBuffer,
	const FSWEBulkInputs& Inputs,
	const FSWEBulkOutputs& Outputs) {
	FGlobalShaderMap* ShaderMap = GetGlobalShaderMap(GMaxRHIFeatureLevel);
	FIntVector GridGroups(
		FMath::DivideAndRoundUp(Inputs.GridSizeX, 16),
		FMath::DivideAndRoundUp(Inputs.GridSizeY, 16),
		1
	);

	// Calculate SWE Velocity (ubar) from momentum (qbar)
	TShaderMapRef<FCalcUbarCS> CalcUbarShader(ShaderMap);
	FCalcUbarCS::FParameters* CalcUbarParams = GraphBuilder.AllocParameters<FCalcUbarCS::FParameters>();
	CalcUbarParams->SimConstants = ConstantBuffer;
	CalcUbarParams->qbarIn_x = GraphBuilder.CreateSRV(Inputs.qbarIn_x);
	CalcUbarParams->qbarIn_y = GraphBuilder.CreateSRV(Inputs.qbarIn_y);
	CalcUbarParams->hbarIn = GraphBuilder.CreateSRV(Inputs.hbarOldIn);
	CalcUbarParams->ubarOut_x = GraphBuilder.CreateUAV(Outputs.ubarOut_x);
	CalcUbarParams->ubarOut_y = GraphBuilder.CreateUAV(Outputs.ubarOut_y);

	FComputeShaderUtils::AddPass(
		GraphBuilder,
		RDG_EVENT_NAME("SWE_CalcUbar"),
		ERDGPassFlags::Compute,
		CalcUbarShader,
		CalcUbarParams,
		GridGroups
	);

	// Solve SWE Non-Linear Wave Equations
	TShaderMapRef<FCalcSWECS> CalcSWEEShader(ShaderMap);
	FCalcSWECS::FParameters* CalcSWEParams = GraphBuilder.AllocParameters<FCalcSWECS::FParameters>();
	CalcSWEParams->SimConstants = ConstantBuffer;
	CalcSWEParams->ubarIn_x = GraphBuilder.CreateSRV(Outputs.ubarOut_x);
	CalcSWEParams->ubarIn_y = GraphBuilder.CreateSRV(Outputs.ubarOut_y);
	CalcSWEParams->hbarIn = GraphBuilder.CreateSRV(Inputs.hbarIn);
	CalcSWEParams->H_In = GraphBuilder.CreateSRV(Inputs.H_In);
	CalcSWEParams->delH_x = GraphBuilder.CreateSRV(Inputs.delH_x);
	CalcSWEParams->delH_y = GraphBuilder.CreateSRV(Inputs.delH_y);
	CalcSWEParams->terrain = GraphBuilder.CreateSRV(Inputs.terrain);
	CalcSWEParams->ubarNewOut_x = GraphBuilder.CreateUAV(Outputs.ubarNewOut_x);
	CalcSWEParams->ubarNewOut_y = GraphBuilder.CreateUAV(Outputs.ubarNewOut_y);
	CalcSWEParams->qbarOut_x = GraphBuilder.CreateUAV(Outputs.qbarOut_x);
	CalcSWEParams->qbarOut_y = GraphBuilder.CreateUAV(Outputs.qbarOut_y);

	FComputeShaderUtils::AddPass(
		GraphBuilder,
		RDG_EVENT_NAME("SWE_CalcSWE"),
		ERDGPassFlags::Compute,
		CalcSWEEShader,
		CalcSWEParams,
		GridGroups
	);
}

void AddTransportAndIntegratePasses(
	FRDGBuilder& GraphBuilder,
	TUniformBufferRef<FSimConstants> ConstantBuffer,
	const FTransportAndIntegrateInputs& Inputs,
	const FTransportAndIntegrateOutputs& Outputs) {
	FGlobalShaderMap* ShaderMap = GetGlobalShaderMap(GMaxRHIFeatureLevel);
	FIntVector GridGroups(
		FMath::DivideAndRoundUp(Inputs.GridSizeX, 16),
		FMath::DivideAndRoundUp(Inputs.GridSizeY, 16),
		1
	);

	// UpdateTilde (Advect wave height and flow rate)
	TShaderMapRef<FUpdateTildeCS> UpdateTildeShader(ShaderMap);
	FUpdateTildeCS::FParameters* UpdateTildeParams = GraphBuilder.AllocParameters<FUpdateTildeCS::FParameters>();
	UpdateTildeParams->SimConstants = ConstantBuffer;
	UpdateTildeParams->ubarNewIn_x = GraphBuilder.CreateSRV(Inputs.ubarNewIn_x);
	UpdateTildeParams->ubarIn_x = GraphBuilder.CreateSRV(Inputs.ubarIn_x);
	UpdateTildeParams->ubarNewIn_y = GraphBuilder.CreateSRV(Inputs.ubarNewIn_y);
	UpdateTildeParams->ubarIn_y = GraphBuilder.CreateSRV(Inputs.ubarIn_y);
	UpdateTildeParams->qtildePast_x = GraphBuilder.CreateSRV(Inputs.qtildePast_x);
	UpdateTildeParams->qtildePast_y = GraphBuilder.CreateSRV(Inputs.qtildePast_y);
	UpdateTildeParams->hIn = GraphBuilder.CreateSRV(Inputs.hPast);
	UpdateTildeParams->htildePast = GraphBuilder.CreateSRV(Inputs.htildePast);
	UpdateTildeParams->htildeOut = GraphBuilder.CreateUAV(Outputs.htildeOut);
	UpdateTildeParams->qtildeOut_x = GraphBuilder.CreateUAV(Outputs.qtildeOut_x);
	UpdateTildeParams->qtildeOut_y = GraphBuilder.CreateUAV(Outputs.qtildeOut_y);

	FComputeShaderUtils::AddPass(
		GraphBuilder,
		RDG_EVENT_NAME("SWE_UpdateTilde"),
		ERDGPassFlags::Compute,
		UpdateTildeShader,
		UpdateTildeParams,
		GridGroups
	);

	// Transient advection intermediate flux
	FRDGTextureDesc GridFloatDesc = FRDGTextureDesc::Create2D(
		FIntPoint(Inputs.GridSizeX, Inputs.GridSizeY),
		PF_R32_FLOAT,
		FClearValueBinding::None,
		TexCreate_ShaderResource | TexCreate_UAV
	);

	FRDGTextureRef qAdvectOut_x = GraphBuilder.CreateTexture(GridFloatDesc, TEXT("qAdvect_x"));
	FRDGTextureRef qAdvectOut_y = GraphBuilder.CreateTexture(GridFloatDesc, TEXT("qAdvect_y"));

	// CalcQAdvect (Compute advective flux)
	TShaderMapRef<FCalcQAdvectCS> CalcQAdvectShader(ShaderMap);
	FCalcQAdvectCS::FParameters* CalcQAdvectParams = GraphBuilder.AllocParameters<FCalcQAdvectCS::FParameters>();
	CalcQAdvectParams->SimConstants = ConstantBuffer;
	CalcQAdvectParams->ubarNewIn_x = GraphBuilder.CreateSRV(Inputs.ubarNewIn_x);
	CalcQAdvectParams->ubarNewIn_y = GraphBuilder.CreateSRV(Inputs.ubarNewIn_y);
	CalcQAdvectParams->htildeIn = GraphBuilder.CreateSRV(Outputs.htildeOut);
	CalcQAdvectParams->qAdvectOut_x = GraphBuilder.CreateUAV(qAdvectOut_x);
	CalcQAdvectParams->qAdvectOut_y = GraphBuilder.CreateUAV(qAdvectOut_y);

	FComputeShaderUtils::AddPass(
		GraphBuilder,
		RDG_EVENT_NAME("SWE_CalcQAdvect"),
		ERDGPassFlags::Compute,
		CalcQAdvectShader,
		CalcQAdvectParams,
		GridGroups
	);

	// IntegrateH (Time-integrate water depth)
	TShaderMapRef<FIntegrateHCS> IntegrateHShader(ShaderMap);
	FIntegrateHCS::FParameters* IntegrateHParams = GraphBuilder.AllocParameters<FIntegrateHCS::FParameters>();
	IntegrateHParams->SimConstants = ConstantBuffer;
	IntegrateHParams->qbarIn_x = GraphBuilder.CreateSRV(Inputs.qbarIn_x);
	IntegrateHParams->qtildeIn_x = GraphBuilder.CreateSRV(Outputs.qtildeOut_x);
	IntegrateHParams->qAdvectIn_x = GraphBuilder.CreateSRV(qAdvectOut_x);
	IntegrateHParams->qbarIn_y = GraphBuilder.CreateSRV(Inputs.qbarIn_y);
	IntegrateHParams->qtildeIn_y = GraphBuilder.CreateSRV(Outputs.qtildeOut_y);
	IntegrateHParams->qAdvectIn_y = GraphBuilder.CreateSRV(qAdvectOut_y);
	IntegrateHParams->hPast = GraphBuilder.CreateSRV(Inputs.hPast);
	IntegrateHParams->terrain = GraphBuilder.CreateSRV(Inputs.terrain);
	IntegrateHParams->hOut = GraphBuilder.CreateUAV(Outputs.hOut);
	IntegrateHParams->qOut_x = GraphBuilder.CreateUAV(Outputs.qOut_x);
	IntegrateHParams->qOut_y = GraphBuilder.CreateUAV(Outputs.qOut_y);
	IntegrateHParams->hdot_out = GraphBuilder.CreateUAV(Outputs.hdot_out);

	FComputeShaderUtils::AddPass(
		GraphBuilder,
		RDG_EVENT_NAME("SWE_IntegrateH"),
		ERDGPassFlags::Compute,
		IntegrateHShader,
		IntegrateHParams,
		GridGroups
	);

	// 4. Recompute H (Total free surface elevation for next tick)
	if (Outputs.H_Out) {
		TShaderMapRef<FRecomputeHCS> RecomputeHShader(ShaderMap);
		FRecomputeHCS::FParameters* RecomputeHParams = GraphBuilder.AllocParameters<FRecomputeHCS::FParameters>();
		RecomputeHParams->SimConstants = ConstantBuffer;
		RecomputeHParams->hIn = GraphBuilder.CreateSRV(Outputs.hOut);
		RecomputeHParams->terrain = GraphBuilder.CreateSRV(Inputs.terrain);
		RecomputeHParams->H_Out = GraphBuilder.CreateUAV(Outputs.H_Out);

		FComputeShaderUtils::AddPass(
			GraphBuilder,
			RDG_EVENT_NAME("SWE_Final_ReH"),
			ERDGPassFlags::Compute,
			RecomputeHShader,
			RecomputeHParams,
			GridGroups
		);
	}

	// 5. Recompute HAvg (Time-averaged free surface elevation for visual export & CPU readback)
	if (Outputs.HAvg_Out) {
		TShaderMapRef<FRecomputeHAvgCS> RecomputeHAvgShader(ShaderMap);
		FRecomputeHAvgCS::FParameters* RecomputeHAvgParams = GraphBuilder.AllocParameters<FRecomputeHAvgCS::FParameters>();
		RecomputeHAvgParams->SimConstants = ConstantBuffer;
		RecomputeHAvgParams->hIn = GraphBuilder.CreateSRV(Outputs.hOut);
		RecomputeHAvgParams->hPast = GraphBuilder.CreateSRV(Inputs.hPast);
		RecomputeHAvgParams->terrain = GraphBuilder.CreateSRV(Inputs.terrain);
		RecomputeHAvgParams->H_Out = GraphBuilder.CreateUAV(Outputs.HAvg_Out);

		FComputeShaderUtils::AddPass(
			GraphBuilder,
			RDG_EVENT_NAME("SWE_Final_ReHAvg"),
			ERDGPassFlags::Compute,
			RecomputeHAvgShader,
			RecomputeHAvgParams,
			GridGroups
		);
	}
}
