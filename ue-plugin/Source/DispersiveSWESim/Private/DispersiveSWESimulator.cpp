#include "DispersiveSWESimulator.h"
#include "DispersiveSWEShaders.h"
#include "JsonObjectConverter.h"
#include "Misc/FileHelper.h"
#include "Misc/Paths.h"
#include "HAL/FileManager.h"
#include "Serialization/JsonReader.h"
#include "Serialization/JsonSerializer.h"
#include "RenderGraphBuilder.h"
#include "RenderGraphUtils.h"
#include "RenderTargetPool.h"
#include "RHI.h"
#include "RHIResources.h"
#include "RHIGPUReadback.h"
#include "RHIStaticStates.h"
#include "ShaderParameterUtils.h"
#include "ShaderCompilerCore.h"
#include "GlobalShader.h"
#include "CommonRenderResources.h"
#include "TextureResource.h"

// DECLARE_GPU_STAT_NAMED(DispersiveSWESim, TEXT("Dispersive SWE Simulation"));

UDispersiveSWESimulator::UDispersiveSWESimulator() 
{
	PrimaryComponentTick.bCanEverTick = true;
	PrimaryComponentTick.bStartWithTickEnabled = true;
}

void UDispersiveSWESimulator::BeginPlay() 
{
	UE_LOG(LogTemp, Warning, TEXT("UDispersiveSWESimulator::BeginPlay() called"));
	Super::BeginPlay();
	InitializeSimulation();
}

int32 NextPowerOf2(int32 n) 
{
	if (n <= 0) return 1;
	int32 p = 1;
	while (p < n) p <<= 1;
	return p;
}

void UDispersiveSWESimulator::InitializeSimulation() 
{
	UE_LOG(LogTemp, Warning, TEXT("InitializeSimulation() running: GridSizeX=%d, GridSizeY=%d, CellSize=%f, CapturedWorldWidth=%f"), GridSizeX, GridSizeY, CellSize, CapturedWorldWidth);

	// // Disable RDG resource aliasing and extend lifetimes to ensure pooled persistent render targets remain stable
	// auto SafeSetCVar = [](const TCHAR* Name, int32 Value) {
	// 	if (IConsoleVariable* CVar = IConsoleManager::Get().FindConsoleVariable(Name))
	// 	{
	// 		CVar->Set(Value, ECVF_SetByCode);
	// 	}
	// };

	// SafeSetCVar(TEXT("r.RDG.TransientAllocator"), 0);
	// SafeSetCVar(TEXT("r.RDG.ResourceAliasing"), 0);
	// SafeSetCVar(TEXT("r.RDG.Debug.ExtendResourceLifetimes"), 1);
	// SafeSetCVar(TEXT("r.RDG.Debug.ResourceClobber"), 1);
	// SafeSetCVar(TEXT("r.RDG.ClobberResources"), 1);
	// SafeSetCVar(TEXT("r.RDG.ParallelExecute"), 0);

	// Load configuration from JSON if path is provided
	if (!JsonConfigFilePath.IsEmpty()) {
		LoadParametersFromJson(JsonConfigFilePath);
	}

	// Automatically calculate CellSize based on CapturedWorldWidth and GridSizeX if enabled
	if (bAutoCalculateCellSize && GridSizeX > 0) {
		CellSize = CapturedWorldWidth / (float)GridSizeX;
	} 

	// Automatically calculate CalculatedMaxSafeDepth based on MaxSafeDepth
	float CellSizeMeters = CellSize * 0.01f;
	if (MaxSafeDepth > 0.0f) {
		CalculatedMaxSafeDepth = MaxSafeDepth * 0.01f;
	} else {
		float a = 182.80027907467993f;
		float b = 0.045464332332812774f;
		float c = -0.14717654147795045f;
		CalculatedMaxSafeDepth = a * CellSizeMeters * CellSizeMeters + b * CellSizeMeters + c;
		CalculatedMaxSafeDepth *= StabilitySafetyFactor;
		if (CalculatedMaxSafeDepth < 0.0f) CalculatedMaxSafeDepth = 0.0f;
	}

	PaddedSizeX = NextPowerOf2(GridSizeX);
	PaddedSizeY = NextPowerOf2(GridSizeY);
	SimulationTime = 0.0f;

	ENQUEUE_RENDER_COMMAND(InitializeSWESimulation)(
		[this](FRHICommandListImmediate& RHICmdList) {
			UE_LOG(LogTemp, Warning, TEXT("InitializeSWESimulation RenderCommand queueing AllocatePersistentTargets and SetupInitialStates"));
			AllocatePersistentTargets(RHICmdList);
			SetupInitialStates(RHICmdList);
			bInitialized = true;
			UE_LOG(LogTemp, Warning, TEXT("InitializeSWESimulation RenderCommand finished. bInitialized=true"));
		}
	);
}

void UDispersiveSWESimulator::AllocatePersistentTargets(FRHICommandListImmediate& RHICmdList) 
{
	FPooledRenderTargetDesc Desc = FPooledRenderTargetDesc::Create2DDesc(
		FIntPoint(GridSizeX, GridSizeY),
		PF_R32_FLOAT,
		FClearValueBinding::None,
		TexCreate_None,
		TexCreate_ShaderResource | TexCreate_UAV,
		false
	);

	GRenderTargetPool.FindFreeElement(RHICmdList, Desc, TexTerrain, TEXT("Terrain"));
	GRenderTargetPool.FindFreeElement(RHICmdList, Desc, TexH, TEXT("H"));
	GRenderTargetPool.FindFreeElement(RHICmdList, Desc, TexQ_x, TEXT("Q_x"));
	GRenderTargetPool.FindFreeElement(RHICmdList, Desc, TexQ_y, TEXT("Q_y"));
	GRenderTargetPool.FindFreeElement(RHICmdList, Desc, Texh, TEXT("h"));
	GRenderTargetPool.FindFreeElement(RHICmdList, Desc, Texq_x, TEXT("q_x"));
	GRenderTargetPool.FindFreeElement(RHICmdList, Desc, Texq_y, TEXT("q_y"));
	GRenderTargetPool.FindFreeElement(RHICmdList, Desc, TexHOrig, TEXT("HOrig"));
	GRenderTargetPool.FindFreeElement(RHICmdList, Desc, TexQOrig_x, TEXT("QOrig_x"));
	GRenderTargetPool.FindFreeElement(RHICmdList, Desc, TexQOrig_y, TEXT("QOrig_y"));
	GRenderTargetPool.FindFreeElement(RHICmdList, Desc, TexHPast, TEXT("HPast"));
	GRenderTargetPool.FindFreeElement(RHICmdList, Desc, TexQPast_x, TEXT("QPast_x"));
	GRenderTargetPool.FindFreeElement(RHICmdList, Desc, TexQPast_y, TEXT("QPast_y"));
	GRenderTargetPool.FindFreeElement(RHICmdList, Desc, TexAlpha_H, TEXT("alpha_H"));
	GRenderTargetPool.FindFreeElement(RHICmdList, Desc, TexAlpha_Q_x, TEXT("alpha_Q_x"));
	GRenderTargetPool.FindFreeElement(RHICmdList, Desc, TexAlpha_Q_y, TEXT("alpha_Q_y"));
	GRenderTargetPool.FindFreeElement(RHICmdList, Desc, Texhbar, TEXT("hbar"));
	GRenderTargetPool.FindFreeElement(RHICmdList, Desc, TexhbarOld, TEXT("hbarOld"));
	GRenderTargetPool.FindFreeElement(RHICmdList, Desc, Texqbar_x, TEXT("qbar_x"));
	GRenderTargetPool.FindFreeElement(RHICmdList, Desc, Texqbar_y, TEXT("qbar_y"));
	GRenderTargetPool.FindFreeElement(RHICmdList, Desc, Texhtilde, TEXT("htilde"));
	GRenderTargetPool.FindFreeElement(RHICmdList, Desc, TexhtildePast, TEXT("htildePast"));
	GRenderTargetPool.FindFreeElement(RHICmdList, Desc, TexhtildeOld, TEXT("htildeOld"));
	GRenderTargetPool.FindFreeElement(RHICmdList, Desc, TexhtildeOldNext, TEXT("htildeOldNext"));
	GRenderTargetPool.FindFreeElement(RHICmdList, Desc, Texqtilde_x, TEXT("qtilde_x"));
	GRenderTargetPool.FindFreeElement(RHICmdList, Desc, Texqtilde_y, TEXT("qtilde_y"));
	GRenderTargetPool.FindFreeElement(RHICmdList, Desc, Texubar_x, TEXT("ubar_x"));
	GRenderTargetPool.FindFreeElement(RHICmdList, Desc, Texubar_y, TEXT("ubar_y"));
	GRenderTargetPool.FindFreeElement(RHICmdList, Desc, TexubarNew_x, TEXT("ubarNew_x"));
	GRenderTargetPool.FindFreeElement(RHICmdList, Desc, TexubarNew_y, TEXT("ubarNew_y"));
	GRenderTargetPool.FindFreeElement(RHICmdList, Desc, TexqtildePast_x, TEXT("qtildePast_x"));
	GRenderTargetPool.FindFreeElement(RHICmdList, Desc, TexqtildePast_y, TEXT("qtildePast_y"));
	GRenderTargetPool.FindFreeElement(RHICmdList, Desc, TexqAdvect_x, TEXT("qAdvect_x"));
	GRenderTargetPool.FindFreeElement(RHICmdList, Desc, TexqAdvect_y, TEXT("qAdvect_y"));
	GRenderTargetPool.FindFreeElement(RHICmdList, Desc, TexhPast, TEXT("hPast"));
	GRenderTargetPool.FindFreeElement(RHICmdList, Desc, TexdelH_x, TEXT("delH_x"));
	GRenderTargetPool.FindFreeElement(RHICmdList, Desc, TexdelH_y, TEXT("delH_y"));
	GRenderTargetPool.FindFreeElement(RHICmdList, Desc, Texdisp_x, TEXT("disp_x"));
	GRenderTargetPool.FindFreeElement(RHICmdList, Desc, Texdisp_y, TEXT("disp_y"));
	GRenderTargetPool.FindFreeElement(RHICmdList, Desc, TexTerrainExportDummy, TEXT("TerrainExportDummy"));

	FPooledRenderTargetDesc FoamDesc = FPooledRenderTargetDesc::Create2DDesc(
		FIntPoint(GridSizeX, GridSizeY),
		PF_FloatRGBA,
		FClearValueBinding::None,
		TexCreate_None,
		TexCreate_ShaderResource | TexCreate_UAV,
		false
	);
	GRenderTargetPool.FindFreeElement(RHICmdList, FoamDesc, TexFoam, TEXT("FoamState"));
	GRenderTargetPool.FindFreeElement(RHICmdList, FoamDesc, TexNewFoam, TEXT("NewFoamState"));

	FPooledRenderTargetDesc RoughnessDesc = FPooledRenderTargetDesc::Create2DDesc(
		FIntPoint(GridSizeX, 1),
		PF_FloatRGBA,
		FClearValueBinding::None,
		TexCreate_None,
		TexCreate_ShaderResource | TexCreate_UAV,
		false
	);
	GRenderTargetPool.FindFreeElement(RHICmdList, RoughnessDesc, TexRoughness, TEXT("RoughnessState"));
	GRenderTargetPool.FindFreeElement(RHICmdList, RoughnessDesc, TexNewRoughness, TEXT("NewRoughnessState"));

	// Complex/array targets use padded sizes and PF_G32R32F format
	FPooledRenderTargetDesc ComplexArrayDesc = FPooledRenderTargetDesc::Create2DArrayDesc(
		FIntPoint(PaddedSizeX, PaddedSizeY),
		PF_G32R32F,
		FClearValueBinding::None,
		TexCreate_None,
		TexCreate_ShaderResource | TexCreate_UAV,
		false,
		DepthLevels.Num()
	);

	GRenderTargetPool.FindFreeElement(RHICmdList, ComplexArrayDesc, TexHPos, TEXT("HPos"));
	GRenderTargetPool.FindFreeElement(RHICmdList, ComplexArrayDesc, TexHNeg, TEXT("HNeg"));
	GRenderTargetPool.FindFreeElement(RHICmdList, ComplexArrayDesc, TexDisp_x, TEXT("Disp_x"));
	GRenderTargetPool.FindFreeElement(RHICmdList, ComplexArrayDesc, TexDisp_y, TEXT("Disp_y"));
	GRenderTargetPool.FindFreeElement(RHICmdList, ComplexArrayDesc, TexDelH_x, TEXT("DelH_x"));
	GRenderTargetPool.FindFreeElement(RHICmdList, ComplexArrayDesc, TexDelH_y, TEXT("DelH_y"));
	GRenderTargetPool.FindFreeElement(RHICmdList, ComplexArrayDesc, TexFlow_x, TEXT("Flow_x"));
	GRenderTargetPool.FindFreeElement(RHICmdList, ComplexArrayDesc, TexFlow_y, TEXT("Flow_y"));
	GRenderTargetPool.FindFreeElement(RHICmdList, ComplexArrayDesc, TexqHat_x_array, TEXT("qHat_x_array"));
	GRenderTargetPool.FindFreeElement(RHICmdList, ComplexArrayDesc, TexqHat_y_array, TEXT("qHat_y_array"));

	FPooledRenderTargetDesc ComplexPaddedDesc = FPooledRenderTargetDesc::Create2DDesc(
		FIntPoint(PaddedSizeX, PaddedSizeY),
		PF_G32R32F,
		FClearValueBinding::None,
		TexCreate_None,
		TexCreate_ShaderResource | TexCreate_UAV,
		false
	);
	GRenderTargetPool.FindFreeElement(RHICmdList, ComplexPaddedDesc, TexhHat, TEXT("hHat"));
	GRenderTargetPool.FindFreeElement(RHICmdList, ComplexPaddedDesc, TexqHat_x, TEXT("qHat_x"));
	GRenderTargetPool.FindFreeElement(RHICmdList, ComplexPaddedDesc, TexqHat_y, TEXT("qHat_y"));

	// Initialize staging textures for double-buffered async readback
	FRHITextureCreateDesc Desc0 = FRHITextureCreateDesc::Create2D(TEXT("Staging0"), GridSizeX, GridSizeY, PF_R32_FLOAT)
		.SetFlags(TexCreate_CPUReadback)
		.SetNumMips(1);
	StagingTextures[0] = RHICreateTexture(Desc0);

	FRHITextureCreateDesc Desc1 = FRHITextureCreateDesc::Create2D(TEXT("Staging1"), GridSizeX, GridSizeY, PF_R32_FLOAT)
		.SetFlags(TexCreate_CPUReadback)
		.SetNumMips(1);
	StagingTextures[1] = RHICreateTexture(Desc1);

	CPUHeightData[0].Reset();
	CPUHeightData[0].SetNumZeroed(GridSizeX * GridSizeY);
	CPUHeightData[1].Reset();
	CPUHeightData[1].SetNumZeroed(GridSizeX * GridSizeY);
	ActiveCPUBufferIndex.store(0);
	StagingWriteIndex = 0;
	StagingReadIndex = 1;
}

void UDispersiveSWESimulator::AssignConstants(FSimConstants& Constants) const
{
	Constants.time = SimulationTime;
	Constants.gridSizeX = GridSizeX;
	Constants.gridSizeY = GridSizeY;
	Constants.cellSize = CellSize * 0.01f; // Convert cm to meters
	Constants.timeStep = TimeStep;
	Constants.minWaterHeight = MinWaterHeight * 0.01f; // Convert cm to meters
	Constants.maxSafeDepth = CalculatedMaxSafeDepth;
	Constants.surfaceTension = SurfaceTension;
	Constants.density = Density;
	Constants.diffusionIterations = DiffusionIterations;
	Constants.maxDiffusionCells = MaxDiffusionCells;
	Constants.diffusionPenalty = DiffusionPenalty;
	Constants.slopeLimit = SlopeLimit;
	Constants.cflCondition = CFLCondition;
	Constants.gammaTransport = GammaTransport;
	Constants.spongeThickness = SpongeThickness;
	Constants.laplacianDamping = LaplacianDamping;
	Constants.depthNum = DepthLevels.Num();
	Constants.fetch = Fetch;
	Constants.windSpeed = WindSpeed;
	Constants.windAngle = WindAngle;
	Constants.swell = Swell;
	Constants.swellAngle = SwellAngle;
	Constants.choppiness = Choppiness;
	Constants.filterSmall = FilterSmall;
	Constants.filterBig = FilterBig;
	Constants.filterWidth = FilterWidth;
	Constants.filterMin = FilterMin;
	Constants.depthCutoff = DepthCutoff;
	Constants.paddedGridSizeX = PaddedSizeX;
	Constants.paddedGridSizeY = PaddedSizeY;
	Constants.foamThreshold = FoamThreshold;
	Constants.foamMultiplier = FoamMultiplier;
	Constants.foamFade = FoamFade;
	Constants.foamBlur = FoamBlur;
	for (int32 i = 0; i < 16; ++i) {
		float Val = i < DepthLevels.Num() ? DepthLevels[i] : 0.0f;
		Constants.depthLevels[i / 4][i % 4] = Val;
	}
}

void UDispersiveSWESimulator::SetupInitialStates(FRHICommandListImmediate& RHICmdList) 
{

	// Setup the common uniform buffer
	FSimConstants CPUConstants = {};
	AssignConstants(CPUConstants);

	TUniformBufferRef<FSimConstants> ConstantBuffer = CreateUniformBufferImmediate(CPUConstants, EUniformBufferUsage::UniformBuffer_SingleFrame);

	// Clear all persistent simulation state textures to zero to avoid uninitialized GPU garbage (matching Sim2D.cpp)
	{
		FRDGBuilder GraphBuilder(RHICmdList);

		// 1. 2D Scalar simulation fields
		if (TexQ_x.IsValid()) AddClearUAVPass(GraphBuilder, GraphBuilder.CreateUAV(GraphBuilder.RegisterExternalTexture(TexQ_x)), FLinearColor::Black);
		if (TexQ_y.IsValid()) AddClearUAVPass(GraphBuilder, GraphBuilder.CreateUAV(GraphBuilder.RegisterExternalTexture(TexQ_y)), FLinearColor::Black);
		if (Texq_x.IsValid()) AddClearUAVPass(GraphBuilder, GraphBuilder.CreateUAV(GraphBuilder.RegisterExternalTexture(Texq_x)), FLinearColor::Black);
		if (Texq_y.IsValid()) AddClearUAVPass(GraphBuilder, GraphBuilder.CreateUAV(GraphBuilder.RegisterExternalTexture(Texq_y)), FLinearColor::Black);
		if (TexHOrig.IsValid()) AddClearUAVPass(GraphBuilder, GraphBuilder.CreateUAV(GraphBuilder.RegisterExternalTexture(TexHOrig)), FLinearColor::Black);
		if (TexQOrig_x.IsValid()) AddClearUAVPass(GraphBuilder, GraphBuilder.CreateUAV(GraphBuilder.RegisterExternalTexture(TexQOrig_x)), FLinearColor::Black);
		if (TexQOrig_y.IsValid()) AddClearUAVPass(GraphBuilder, GraphBuilder.CreateUAV(GraphBuilder.RegisterExternalTexture(TexQOrig_y)), FLinearColor::Black);
		if (TexHPast.IsValid()) AddClearUAVPass(GraphBuilder, GraphBuilder.CreateUAV(GraphBuilder.RegisterExternalTexture(TexHPast)), FLinearColor::Black);
		if (TexQPast_x.IsValid()) AddClearUAVPass(GraphBuilder, GraphBuilder.CreateUAV(GraphBuilder.RegisterExternalTexture(TexQPast_x)), FLinearColor::Black);
		if (TexQPast_y.IsValid()) AddClearUAVPass(GraphBuilder, GraphBuilder.CreateUAV(GraphBuilder.RegisterExternalTexture(TexQPast_y)), FLinearColor::Black);
		if (TexAlpha_H.IsValid()) AddClearUAVPass(GraphBuilder, GraphBuilder.CreateUAV(GraphBuilder.RegisterExternalTexture(TexAlpha_H)), FLinearColor::Black);
		if (TexAlpha_Q_x.IsValid()) AddClearUAVPass(GraphBuilder, GraphBuilder.CreateUAV(GraphBuilder.RegisterExternalTexture(TexAlpha_Q_x)), FLinearColor::Black);
		if (TexAlpha_Q_y.IsValid()) AddClearUAVPass(GraphBuilder, GraphBuilder.CreateUAV(GraphBuilder.RegisterExternalTexture(TexAlpha_Q_y)), FLinearColor::Black);
		if (Texqbar_x.IsValid()) AddClearUAVPass(GraphBuilder, GraphBuilder.CreateUAV(GraphBuilder.RegisterExternalTexture(Texqbar_x)), FLinearColor::Black);
		if (Texqbar_y.IsValid()) AddClearUAVPass(GraphBuilder, GraphBuilder.CreateUAV(GraphBuilder.RegisterExternalTexture(Texqbar_y)), FLinearColor::Black);
		if (Texhtilde.IsValid()) AddClearUAVPass(GraphBuilder, GraphBuilder.CreateUAV(GraphBuilder.RegisterExternalTexture(Texhtilde)), FLinearColor::Black);
		if (TexhtildePast.IsValid()) AddClearUAVPass(GraphBuilder, GraphBuilder.CreateUAV(GraphBuilder.RegisterExternalTexture(TexhtildePast)), FLinearColor::Black);
		if (TexhtildeOld.IsValid()) AddClearUAVPass(GraphBuilder, GraphBuilder.CreateUAV(GraphBuilder.RegisterExternalTexture(TexhtildeOld)), FLinearColor::Black);
		if (TexhtildeOldNext.IsValid()) AddClearUAVPass(GraphBuilder, GraphBuilder.CreateUAV(GraphBuilder.RegisterExternalTexture(TexhtildeOldNext)), FLinearColor::Black);
		if (Texqtilde_x.IsValid()) AddClearUAVPass(GraphBuilder, GraphBuilder.CreateUAV(GraphBuilder.RegisterExternalTexture(Texqtilde_x)), FLinearColor::Black);
		if (Texqtilde_y.IsValid()) AddClearUAVPass(GraphBuilder, GraphBuilder.CreateUAV(GraphBuilder.RegisterExternalTexture(Texqtilde_y)), FLinearColor::Black);
		if (TexqtildePast_x.IsValid()) AddClearUAVPass(GraphBuilder, GraphBuilder.CreateUAV(GraphBuilder.RegisterExternalTexture(TexqtildePast_x)), FLinearColor::Black);
		if (TexqtildePast_y.IsValid()) AddClearUAVPass(GraphBuilder, GraphBuilder.CreateUAV(GraphBuilder.RegisterExternalTexture(TexqtildePast_y)), FLinearColor::Black);
		if (Texubar_x.IsValid()) AddClearUAVPass(GraphBuilder, GraphBuilder.CreateUAV(GraphBuilder.RegisterExternalTexture(Texubar_x)), FLinearColor::Black);
		if (Texubar_y.IsValid()) AddClearUAVPass(GraphBuilder, GraphBuilder.CreateUAV(GraphBuilder.RegisterExternalTexture(Texubar_y)), FLinearColor::Black);
		if (TexubarNew_x.IsValid()) AddClearUAVPass(GraphBuilder, GraphBuilder.CreateUAV(GraphBuilder.RegisterExternalTexture(TexubarNew_x)), FLinearColor::Black);
		if (TexubarNew_y.IsValid()) AddClearUAVPass(GraphBuilder, GraphBuilder.CreateUAV(GraphBuilder.RegisterExternalTexture(TexubarNew_y)), FLinearColor::Black);
		if (TexqAdvect_x.IsValid()) AddClearUAVPass(GraphBuilder, GraphBuilder.CreateUAV(GraphBuilder.RegisterExternalTexture(TexqAdvect_x)), FLinearColor::Black);
		if (TexqAdvect_y.IsValid()) AddClearUAVPass(GraphBuilder, GraphBuilder.CreateUAV(GraphBuilder.RegisterExternalTexture(TexqAdvect_y)), FLinearColor::Black);
		if (TexhPast.IsValid()) AddClearUAVPass(GraphBuilder, GraphBuilder.CreateUAV(GraphBuilder.RegisterExternalTexture(TexhPast)), FLinearColor::Black);
		if (TexdelH_x.IsValid()) AddClearUAVPass(GraphBuilder, GraphBuilder.CreateUAV(GraphBuilder.RegisterExternalTexture(TexdelH_x)), FLinearColor::Black);
		if (TexdelH_y.IsValid()) AddClearUAVPass(GraphBuilder, GraphBuilder.CreateUAV(GraphBuilder.RegisterExternalTexture(TexdelH_y)), FLinearColor::Black);
		if (Texdisp_x.IsValid()) AddClearUAVPass(GraphBuilder, GraphBuilder.CreateUAV(GraphBuilder.RegisterExternalTexture(Texdisp_x)), FLinearColor::Black);
		if (Texdisp_y.IsValid()) AddClearUAVPass(GraphBuilder, GraphBuilder.CreateUAV(GraphBuilder.RegisterExternalTexture(Texdisp_y)), FLinearColor::Black);
		if (TexTerrainExportDummy.IsValid()) AddClearUAVPass(GraphBuilder, GraphBuilder.CreateUAV(GraphBuilder.RegisterExternalTexture(TexTerrainExportDummy)), FLinearColor::Black);

		// 2. 2D Complex float2 fields
		if (TexhHat.IsValid()) AddClearUAVPass(GraphBuilder, GraphBuilder.CreateUAV(GraphBuilder.RegisterExternalTexture(TexhHat)), FLinearColor::Black);
		if (TexqHat_x.IsValid()) AddClearUAVPass(GraphBuilder, GraphBuilder.CreateUAV(GraphBuilder.RegisterExternalTexture(TexqHat_x)), FLinearColor::Black);
		if (TexqHat_y.IsValid()) AddClearUAVPass(GraphBuilder, GraphBuilder.CreateUAV(GraphBuilder.RegisterExternalTexture(TexqHat_y)), FLinearColor::Black);

		// 3. 2D Complex float2 array fields
		if (TexDisp_x.IsValid()) AddClearUAVPass(GraphBuilder, GraphBuilder.CreateUAV(GraphBuilder.RegisterExternalTexture(TexDisp_x)), FLinearColor::Black);
		if (TexDisp_y.IsValid()) AddClearUAVPass(GraphBuilder, GraphBuilder.CreateUAV(GraphBuilder.RegisterExternalTexture(TexDisp_y)), FLinearColor::Black);
		if (TexDelH_x.IsValid()) AddClearUAVPass(GraphBuilder, GraphBuilder.CreateUAV(GraphBuilder.RegisterExternalTexture(TexDelH_x)), FLinearColor::Black);
		if (TexDelH_y.IsValid()) AddClearUAVPass(GraphBuilder, GraphBuilder.CreateUAV(GraphBuilder.RegisterExternalTexture(TexDelH_y)), FLinearColor::Black);
		if (TexFlow_x.IsValid()) AddClearUAVPass(GraphBuilder, GraphBuilder.CreateUAV(GraphBuilder.RegisterExternalTexture(TexFlow_x)), FLinearColor::Black);
		if (TexFlow_y.IsValid()) AddClearUAVPass(GraphBuilder, GraphBuilder.CreateUAV(GraphBuilder.RegisterExternalTexture(TexFlow_y)), FLinearColor::Black);
		if (TexqHat_x_array.IsValid()) AddClearUAVPass(GraphBuilder, GraphBuilder.CreateUAV(GraphBuilder.RegisterExternalTexture(TexqHat_x_array)), FLinearColor::Black);
		if (TexqHat_y_array.IsValid()) AddClearUAVPass(GraphBuilder, GraphBuilder.CreateUAV(GraphBuilder.RegisterExternalTexture(TexqHat_y_array)), FLinearColor::Black);
		if (TexHPos.IsValid()) AddClearUAVPass(GraphBuilder, GraphBuilder.CreateUAV(GraphBuilder.RegisterExternalTexture(TexHPos)), FLinearColor::Black);
		if (TexHNeg.IsValid()) AddClearUAVPass(GraphBuilder, GraphBuilder.CreateUAV(GraphBuilder.RegisterExternalTexture(TexHNeg)), FLinearColor::Black);

		// 4. Foam & Roughness
		if (TexFoam.IsValid()) AddClearUAVPass(GraphBuilder, GraphBuilder.CreateUAV(GraphBuilder.RegisterExternalTexture(TexFoam)), FLinearColor::Black);
		if (TexNewFoam.IsValid()) AddClearUAVPass(GraphBuilder, GraphBuilder.CreateUAV(GraphBuilder.RegisterExternalTexture(TexNewFoam)), FLinearColor::Black);
		if (TexRoughness.IsValid()) AddClearUAVPass(GraphBuilder, GraphBuilder.CreateUAV(GraphBuilder.RegisterExternalTexture(TexRoughness)), FLinearColor::Black);
		if (TexNewRoughness.IsValid()) AddClearUAVPass(GraphBuilder, GraphBuilder.CreateUAV(GraphBuilder.RegisterExternalTexture(TexNewRoughness)), FLinearColor::Black);

		GraphBuilder.Execute();
	}

	// Initialize Terrain and Water Height fields
	if (TerrainHeightInputRT && 
		TerrainHeightInputRT->GetRenderTargetResource() && 
		TerrainHeightInputRT->GetRenderTargetResource()->GetTexture2DRHI()) 
	{
		FRDGBuilder GraphBuilder(RHICmdList);

		// Import the live terrain height input render target
		FRDGTextureRef TerrainInput_RDG = GraphBuilder.RegisterExternalTexture(
			CreateRenderTarget(TerrainHeightInputRT->GetRenderTargetResource()->GetTexture2DRHI(), TEXT("TerrainInput"))
		);

		// Import the persistent states
		FRDGTextureRef Terrain_RDG = GraphBuilder.RegisterExternalTexture(TexTerrain);
		FRDGTextureRef H_RDG = GraphBuilder.RegisterExternalTexture(TexH);
		FRDGTextureRef h_RDG = GraphBuilder.RegisterExternalTexture(Texh);
		FRDGTextureRef hbar_RDG = GraphBuilder.RegisterExternalTexture(Texhbar);
		FRDGTextureRef hbarOld_RDG = GraphBuilder.RegisterExternalTexture(TexhbarOld);

		// Run the compute shader to initialize water height on GPU
		TShaderMapRef<FInitializeWaterCS> InitWaterHeightCS(GetGlobalShaderMap(GMaxRHIFeatureLevel));
		FInitializeWaterCS::FParameters* InitParams = GraphBuilder.AllocParameters<FInitializeWaterCS::FParameters>();
		InitParams->SimConstants = ConstantBuffer;
		InitParams->WaterLevel = WaterLevel * 0.01f; // Convert cm to meters
		InitParams->TerrainCaptureCameraZ = TerrainCaptureCameraZ;
		InitParams->terrain = GraphBuilder.CreateSRV(TerrainInput_RDG);
		InitParams->hOut = GraphBuilder.CreateUAV(h_RDG);
		InitParams->H_Out = GraphBuilder.CreateUAV(H_RDG);
		InitParams->terrainOut = GraphBuilder.CreateUAV(Terrain_RDG);

		FRDGTextureRef TerrainExportRDG;
		if (TerrainRT && TerrainRT->GetRenderTargetResource())
		{
			TerrainExportRDG = GraphBuilder.RegisterExternalTexture(
				CreateRenderTarget(TerrainRT->GetRenderTargetResource()->GetTexture2DRHI(), TEXT("TerrainExport"))
			);
		}
		else
		{
			TerrainExportRDG = GraphBuilder.RegisterExternalTexture(TexTerrainExportDummy);
		}
		InitParams->terrainOutCM = GraphBuilder.CreateUAV(TerrainExportRDG);

		FComputeShaderUtils::AddPass(
			GraphBuilder,
			RDG_EVENT_NAME("SWE_InitializeWater_GPU"),
			ERDGPassFlags::Compute,
			InitWaterHeightCS,
			InitParams,
			FIntVector(FMath::DivideAndRoundUp(GridSizeX, 16), FMath::DivideAndRoundUp(GridSizeY, 16), 1)
		);

		// Copy the computed water height h to hbar and hbarOld
		AddCopyTexturePass(GraphBuilder, h_RDG, hbar_RDG);
		AddCopyTexturePass(GraphBuilder, h_RDG, hbarOld_RDG);

		// Initialize persistent foam target to 0
		if (TexFoam.IsValid()) {
			FRDGTextureRef Foam_RDG = GraphBuilder.RegisterExternalTexture(TexFoam);
			AddClearUAVPass(GraphBuilder, GraphBuilder.CreateUAV(Foam_RDG), FLinearColor::Black);
		}

		// Initialize persistent roughness target to 0 (smooth water by default)
		if (TexRoughness.IsValid()) {
			FRDGTextureRef Roughness_RDG = GraphBuilder.RegisterExternalTexture(TexRoughness);
			AddClearUAVPass(GraphBuilder, GraphBuilder.CreateUAV(Roughness_RDG), FLinearColor::Black);
		}

		GraphBuilder.Execute();
	}
	else
	{
		TArray<float> TerrainData;
		TerrainData.SetNumZeroed(GridSizeX * GridSizeY);
		TArray<float> WaterData;
		WaterData.SetNumZeroed(GridSizeX * GridSizeY);
		TArray<float> HData;
		HData.SetNumZeroed(GridSizeX * GridSizeY);

		float TerrainHeight = -13.0f;
		float TerrainScale = 20.0f;
		float WaterLevelLocal = WaterLevel * 0.01f; // Convert cm to meters

		for (int32 y = 0; y < GridSizeY; y++)
		{
			for (int32 x = 0; x < GridSizeX; x++)
			{
				float xf = (float)x / (GridSizeX - 1);
				float yf = (float)y / (GridSizeY - 1);
				int32 i = y * GridSizeX + x;

				// Dunes beach terrain
				float dunes = 0.05f * FMath::Sin(20.f * yf);
				TerrainData[i] = TerrainHeight + TerrainScale * (xf * (1.0f + dunes));

				// Flat initial water level
				float WaterSurface = WaterLevelLocal;
				WaterData[i] = FMath::Max(0.0f, WaterSurface - TerrainData[i]);
				HData[i] = TerrainData[i] + WaterData[i];
			}
		}

		FUpdateTextureRegion2D Region(0, 0, 0, 0, GridSizeX, GridSizeY);

		FRHITexture2D* TerrainRHI = static_cast<FRHITexture2D*>(TexTerrain->GetRHI());
		RHIUpdateTexture2D(TerrainRHI, 0, Region, GridSizeX * sizeof(float), (uint8*)TerrainData.GetData());

		FRHITexture2D* HRHI = static_cast<FRHITexture2D*>(TexH->GetRHI());
		RHIUpdateTexture2D(HRHI, 0, Region, GridSizeX * sizeof(float), (uint8*)HData.GetData());

		FRHITexture2D* hRHI = static_cast<FRHITexture2D*>(Texh->GetRHI());
		RHIUpdateTexture2D(hRHI, 0, Region, GridSizeX * sizeof(float), (uint8*)WaterData.GetData());

		FRHITexture2D* hbarRHI = static_cast<FRHITexture2D*>(Texhbar->GetRHI());
		RHIUpdateTexture2D(hbarRHI, 0, Region, GridSizeX * sizeof(float), (uint8*)WaterData.GetData());

		FRHITexture2D* hbarOldRHI = static_cast<FRHITexture2D*>(TexhbarOld->GetRHI());
		RHIUpdateTexture2D(hbarOldRHI, 0, Region, GridSizeX * sizeof(float), (uint8*)WaterData.GetData());
	}

	// Initialize the Wave Spectrum
	{
		FRDGBuilder GraphBuilder(RHICmdList);

		FRDGTextureRef HPos_RDG = GraphBuilder.RegisterExternalTexture(TexHPos);
		FRDGTextureRef HNeg_RDG = GraphBuilder.RegisterExternalTexture(TexHNeg);

		TShaderMapRef<FPopulateSpectrumCS> PopulateSpectrumCS(GetGlobalShaderMap(GMaxRHIFeatureLevel));
		FPopulateSpectrumCS::FParameters* PassParams = GraphBuilder.AllocParameters<FPopulateSpectrumCS::FParameters>();
		PassParams->SimConstants = ConstantBuffer;
		PassParams->HPosOut = GraphBuilder.CreateUAV(HPos_RDG);
		PassParams->HNegOut = GraphBuilder.CreateUAV(HNeg_RDG);

		FComputeShaderUtils::AddPass(
			GraphBuilder,
			RDG_EVENT_NAME("SWE_PopulateSpectrum"),
			ERDGPassFlags::Compute,
			PopulateSpectrumCS,
			PassParams,
			FIntVector(FMath::DivideAndRoundUp(PaddedSizeX, 16), FMath::DivideAndRoundUp(PaddedSizeY, 16), DepthLevels.Num())
		);

		GraphBuilder.Execute();
	}

	// Export the captured terrain to a raw float32 file for standalone simulation testing
	if (TexTerrain.IsValid())
	{
		FRHITexture2D* TerrainRHI = static_cast<FRHITexture2D*>(TexTerrain->GetRHI());
		
		FRHIGPUTextureReadback* GPUReadback = new FRHIGPUTextureReadback(TEXT("TerrainExportReadback"));
		GPUReadback->EnqueueCopy(RHICmdList, TerrainRHI);

		// Force GPU to finish all work so we can read immediately (fine since this is one-time startup)
		RHICmdList.BlockUntilGPUIdle();

		if (GPUReadback->IsReady())
		{
			int32 OutWidth = 0;
			int32 OutHeight = 0;
			float* FloatData = static_cast<float*>(GPUReadback->Lock(OutWidth, &OutHeight));
			if (FloatData)
			{
				TArray<float> RawFloatData;
				RawFloatData.SetNumUninitialized(GridSizeX * GridSizeY);
				// Stride / Row pitch can sometimes be larger than GridSizeX, so copy row by row
				for (int32 y = 0; y < GridSizeY; ++y)
				{
					FMemory::Memcpy(&RawFloatData[y * GridSizeX], &FloatData[y * OutWidth], GridSizeX * sizeof(float));
				}
				
				GPUReadback->Unlock();

				FString ExportFilePath = FPaths::ProjectDir() / TEXT("terrain_captured.raw");
				FArchive* Ar = IFileManager::Get().CreateFileWriter(*ExportFilePath);
				if (Ar)
				{
					Ar->Serialize(RawFloatData.GetData(), RawFloatData.Num() * sizeof(float));
					delete Ar;
					UE_LOG(LogTemp, Warning, TEXT("DispersiveSWESim: Successfully exported captured terrain to %s"), *ExportFilePath);
				}
				else
				{
					UE_LOG(LogTemp, Error, TEXT("DispersiveSWESim: Failed to create export file for terrain!"));
				}
			}
		}
		delete GPUReadback;
	}
}

void UDispersiveSWESimulator::TickComponent(float DeltaTime, ELevelTick TickType, FActorComponentTickFunction* ThisTickFunction) 
{
	Super::TickComponent(DeltaTime, TickType, ThisTickFunction);

	if (!bInitialized) return;

	SimulationTime += TimeStep;

	FSimConstants Constants = {};
	AssignConstants(Constants);

	ENQUEUE_RENDER_COMMAND(ExecuteSWESimulation)(
		[this, Constants](FRHICommandListImmediate& RHICmdList) {
			ExecuteSimulation_RenderThread(RHICmdList, Constants);
		}
	);
}

void UDispersiveSWESimulator::ExecuteSimulation_RenderThread(
	FRHICommandListImmediate& RHICmdList,
	const FSimConstants& Constants) 
{
	FRDGBuilder GraphBuilder(RHICmdList);

	// Import persistent buffers (matching Sim2D.h)
	FRDGTextureRef Terrain_RDG = GraphBuilder.RegisterExternalTexture(TexTerrain);
	FRDGTextureRef H_RDG = GraphBuilder.RegisterExternalTexture(TexH);
	FRDGTextureRef Qx_RDG = GraphBuilder.RegisterExternalTexture(TexQ_x);
	FRDGTextureRef Qy_RDG = GraphBuilder.RegisterExternalTexture(TexQ_y);
	FRDGTextureRef h_RDG = GraphBuilder.RegisterExternalTexture(Texh);
	FRDGTextureRef qx_RDG = GraphBuilder.RegisterExternalTexture(Texq_x);
	FRDGTextureRef qy_RDG = GraphBuilder.RegisterExternalTexture(Texq_y);
	FRDGTextureRef HOrig_RDG = GraphBuilder.RegisterExternalTexture(TexHOrig);
	FRDGTextureRef QOrig_x_RDG = GraphBuilder.RegisterExternalTexture(TexQOrig_x);
	FRDGTextureRef QOrig_y_RDG = GraphBuilder.RegisterExternalTexture(TexQOrig_y);
	FRDGTextureRef HPast_RDG = GraphBuilder.RegisterExternalTexture(TexHPast);
	FRDGTextureRef QPast_x_RDG = GraphBuilder.RegisterExternalTexture(TexQPast_x);
	FRDGTextureRef QPast_y_RDG = GraphBuilder.RegisterExternalTexture(TexQPast_y);
	FRDGTextureRef alpha_H_RDG = GraphBuilder.RegisterExternalTexture(TexAlpha_H);
	FRDGTextureRef alpha_Q_x_RDG = GraphBuilder.RegisterExternalTexture(TexAlpha_Q_x);
	FRDGTextureRef alpha_Q_y_RDG = GraphBuilder.RegisterExternalTexture(TexAlpha_Q_y);
	FRDGTextureRef hbar_RDG = GraphBuilder.RegisterExternalTexture(Texhbar);
	FRDGTextureRef hbarOld_RDG = GraphBuilder.RegisterExternalTexture(TexhbarOld);
	FRDGTextureRef qbarx_RDG = GraphBuilder.RegisterExternalTexture(Texqbar_x);
	FRDGTextureRef qbary_RDG = GraphBuilder.RegisterExternalTexture(Texqbar_y);
	FRDGTextureRef htilde_RDG = GraphBuilder.RegisterExternalTexture(Texhtilde);
	FRDGTextureRef htildePast_RDG = GraphBuilder.RegisterExternalTexture(TexhtildePast);
	FRDGTextureRef htildeOld_RDG = GraphBuilder.RegisterExternalTexture(TexhtildeOld);
	FRDGTextureRef htildeOldNext_RDG = GraphBuilder.RegisterExternalTexture(TexhtildeOldNext);
	FRDGTextureRef qtildex_RDG = GraphBuilder.RegisterExternalTexture(Texqtilde_x);
	FRDGTextureRef qtildey_RDG = GraphBuilder.RegisterExternalTexture(Texqtilde_y);
	FRDGTextureRef ubarx_RDG = GraphBuilder.RegisterExternalTexture(Texubar_x);
	FRDGTextureRef ubary_RDG = GraphBuilder.RegisterExternalTexture(Texubar_y);
	FRDGTextureRef ubarNew_x_RDG = GraphBuilder.RegisterExternalTexture(TexubarNew_x);
	FRDGTextureRef ubarNew_y_RDG = GraphBuilder.RegisterExternalTexture(TexubarNew_y);
	FRDGTextureRef qtildePast_x_RDG = GraphBuilder.RegisterExternalTexture(TexqtildePast_x);
	FRDGTextureRef qtildePast_y_RDG = GraphBuilder.RegisterExternalTexture(TexqtildePast_y);
	FRDGTextureRef qAdvect_x_RDG = GraphBuilder.RegisterExternalTexture(TexqAdvect_x);
	FRDGTextureRef qAdvect_y_RDG = GraphBuilder.RegisterExternalTexture(TexqAdvect_y);
	FRDGTextureRef hPast_RDG = GraphBuilder.RegisterExternalTexture(TexhPast);
	FRDGTextureRef hHat_RDG = GraphBuilder.RegisterExternalTexture(TexhHat);
	FRDGTextureRef qHat_x_RDG = GraphBuilder.RegisterExternalTexture(TexqHat_x);
	FRDGTextureRef qHat_y_RDG = GraphBuilder.RegisterExternalTexture(TexqHat_y);
	FRDGTextureRef qHat_x_array_RDG = GraphBuilder.RegisterExternalTexture(TexqHat_x_array);
	FRDGTextureRef qHat_y_array_RDG = GraphBuilder.RegisterExternalTexture(TexqHat_y_array);

	// FFT Wave propagation outputs (DEPTH_NUM layers)
	FRDGTextureRef HPos_RDG = GraphBuilder.RegisterExternalTexture(TexHPos);
	FRDGTextureRef HNeg_RDG = GraphBuilder.RegisterExternalTexture(TexHNeg);
	FRDGTextureRef Disp_x_RDG = GraphBuilder.RegisterExternalTexture(TexDisp_x);
	FRDGTextureRef Disp_y_RDG = GraphBuilder.RegisterExternalTexture(TexDisp_y);
	FRDGTextureRef DelH_x_RDG = GraphBuilder.RegisterExternalTexture(TexDelH_x);
	FRDGTextureRef DelH_y_RDG = GraphBuilder.RegisterExternalTexture(TexDelH_y);
	FRDGTextureRef Flow_x_RDG = GraphBuilder.RegisterExternalTexture(TexFlow_x);
	FRDGTextureRef Flow_y_RDG = GraphBuilder.RegisterExternalTexture(TexFlow_y);

	// Wind wave outputs
	FRDGTextureRef delH_x_RDG = GraphBuilder.RegisterExternalTexture(TexdelH_x);
	FRDGTextureRef delH_y_RDG = GraphBuilder.RegisterExternalTexture(TexdelH_y);
	FRDGTextureRef disp_x_RDG = GraphBuilder.RegisterExternalTexture(Texdisp_x);
	FRDGTextureRef disp_y_RDG = GraphBuilder.RegisterExternalTexture(Texdisp_y);

	// 3. Create Uniform Buffer
	TUniformBufferRef<FSimConstants> ConstantBuffer = CreateUniformBufferImmediate(Constants, EUniformBufferUsage::UniformBuffer_SingleFrame);

	// Get shader maps
	FGlobalShaderMap* ShaderMap = GetGlobalShaderMap(GMaxRHIFeatureLevel);

	FIntVector GridGroups(FMath::DivideAndRoundUp(GridSizeX, 16), FMath::DivideAndRoundUp(GridSizeY, 16), 1);
	FIntVector PaddedGroups(FMath::DivideAndRoundUp(PaddedSizeX, 16), FMath::DivideAndRoundUp(PaddedSizeY, 16), 1);
	FIntVector ComplexArrayGroups(FMath::DivideAndRoundUp(PaddedSizeX, 16), FMath::DivideAndRoundUp(PaddedSizeY, 16), Constants.depthNum);

	// ----------------------------------------------------
	// DECOMPOSITION STEP
	// ----------------------------------------------------

	// InitDecomp
	{
		TShaderMapRef<FInitDecompCS> Shader(ShaderMap);
		FInitDecompCS::FParameters* Params = GraphBuilder.AllocParameters<FInitDecompCS::FParameters>();
		Params->SimConstants = ConstantBuffer;
		Params->hIn = GraphBuilder.CreateSRV(h_RDG);
		Params->qIn_x = GraphBuilder.CreateSRV(qx_RDG);
		Params->qIn_y = GraphBuilder.CreateSRV(qy_RDG);
		Params->terrain = GraphBuilder.CreateSRV(Terrain_RDG);
		Params->H_Out = GraphBuilder.CreateUAV(H_RDG);
		Params->Q_Out_x = GraphBuilder.CreateUAV(Qx_RDG);
		Params->Q_Out_y = GraphBuilder.CreateUAV(Qy_RDG);

		FComputeShaderUtils::AddPass(GraphBuilder, RDG_EVENT_NAME("SWE_Decomp_Init"), ERDGPassFlags::Compute, Shader, Params, GridGroups);
	}

	// Copy initial fields to Orig fields before starting Jacobi solver (matching Sim2D.cpp CopyField)
	AddCopyTexturePass(GraphBuilder, H_RDG, HOrig_RDG);
	AddCopyTexturePass(GraphBuilder, Qx_RDG, QOrig_x_RDG);
	AddCopyTexturePass(GraphBuilder, Qy_RDG, QOrig_y_RDG);

	// CalcDiffusionCoeffs
	{
		TShaderMapRef<FCalcDiffusionCoeffsCS> Shader(ShaderMap);
		FCalcDiffusionCoeffsCS::FParameters* Params = GraphBuilder.AllocParameters<FCalcDiffusionCoeffsCS::FParameters>();
		Params->SimConstants = ConstantBuffer;
		Params->H_In = GraphBuilder.CreateSRV(H_RDG);
		Params->terrain = GraphBuilder.CreateSRV(Terrain_RDG);
		Params->alpha_HOut = GraphBuilder.CreateUAV(alpha_H_RDG);
		Params->alpha_QOut_x = GraphBuilder.CreateUAV(alpha_Q_x_RDG);
		Params->alpha_QOut_y = GraphBuilder.CreateUAV(alpha_Q_y_RDG);

		FComputeShaderUtils::AddPass(GraphBuilder, RDG_EVENT_NAME("SWE_Decomp_Coeffs"), ERDGPassFlags::Compute, Shader, Params, GridGroups);
	}

	// Diffusion loop - low pass filter H and Q using implicit Jacobi solver (ping-ponging H/HPast, Q_x/QPast_x, Q_y/QPast_y)
	{
		FRDGTextureRef H_Src = H_RDG;
		FRDGTextureRef H_Dst = HPast_RDG;
		FRDGTextureRef Qx_Src = Qx_RDG;
		FRDGTextureRef Qx_Dst = QPast_x_RDG;
		FRDGTextureRef Qy_Src = Qy_RDG;
		FRDGTextureRef Qy_Dst = QPast_y_RDG;

		TShaderMapRef<FDiffusionStepCS> Shader(ShaderMap);
		for (int32 j = 0; j < Constants.diffusionIterations; j++)
		{
			FDiffusionStepCS::FParameters* Params = GraphBuilder.AllocParameters<FDiffusionStepCS::FParameters>();
			Params->SimConstants = ConstantBuffer;
			Params->terrain = GraphBuilder.CreateSRV(Terrain_RDG);
			Params->H_Orig = GraphBuilder.CreateSRV(HOrig_RDG);
			Params->Q_Orig_x = GraphBuilder.CreateSRV(QOrig_x_RDG);
			Params->Q_Orig_y = GraphBuilder.CreateSRV(QOrig_y_RDG);
			Params->H_Past = GraphBuilder.CreateSRV(H_Src);
			Params->Q_Past_x = GraphBuilder.CreateSRV(Qx_Src);
			Params->Q_Past_y = GraphBuilder.CreateSRV(Qy_Src);
			Params->alpha_HIn = GraphBuilder.CreateSRV(alpha_H_RDG);
			Params->alpha_QIn_x = GraphBuilder.CreateSRV(alpha_Q_x_RDG);
			Params->alpha_QIn_y = GraphBuilder.CreateSRV(alpha_Q_y_RDG);
			Params->H_Out = GraphBuilder.CreateUAV(H_Dst);
			Params->Q_Out_x = GraphBuilder.CreateUAV(Qx_Dst);
			Params->Q_Out_y = GraphBuilder.CreateUAV(Qy_Dst);

			FComputeShaderUtils::AddPass(GraphBuilder, RDG_EVENT_NAME("SWE_Decomp_Diffusion"), ERDGPassFlags::Compute, Shader, Params, GridGroups);

			// Swap sources and destinations for next loop iteration
			Swap(H_Src, H_Dst);
			Swap(Qx_Src, Qx_Dst);
			Swap(Qy_Src, Qy_Dst);
		}

		// If Jacobi loop ended on Past buffer (odd iteration count), swap persistent pointers without copying
		if (H_Src != H_RDG)
		{
			Swap(TexH, TexHPast);
		}
		if (Qx_Src != Qx_RDG)
		{
			Swap(TexQ_x, TexQPast_x);
		}
		if (Qy_Src != Qy_RDG)
		{
			Swap(TexQ_y, TexQPast_y);
		}
		H_RDG = H_Src;
		Qx_RDG = Qx_Src;
		Qy_RDG = Qy_Src;
	}

	// DecomposeFields
	{
		TShaderMapRef<FDecomposeFieldsCS> Shader(ShaderMap);
		FDecomposeFieldsCS::FParameters* Params = GraphBuilder.AllocParameters<FDecomposeFieldsCS::FParameters>();
		Params->SimConstants = ConstantBuffer;
		Params->H_In = GraphBuilder.CreateSRV(H_RDG);
		Params->Q_In_x = GraphBuilder.CreateSRV(Qx_RDG);
		Params->Q_In_y = GraphBuilder.CreateSRV(Qy_RDG);
		Params->hIn = GraphBuilder.CreateSRV(h_RDG);
		Params->qIn_x = GraphBuilder.CreateSRV(qx_RDG);
		Params->qIn_y = GraphBuilder.CreateSRV(qy_RDG);
		Params->terrain = GraphBuilder.CreateSRV(Terrain_RDG);
		Params->hbarOut = GraphBuilder.CreateUAV(hbar_RDG);
		Params->qbarOut_x = GraphBuilder.CreateUAV(qbarx_RDG);
		Params->qbarOut_y = GraphBuilder.CreateUAV(qbary_RDG);
		Params->htildeOut = GraphBuilder.CreateUAV(htilde_RDG);
		Params->qtildeOut_x = GraphBuilder.CreateUAV(qtildex_RDG);
		Params->qtildeOut_y = GraphBuilder.CreateUAV(qtildey_RDG);

		FComputeShaderUtils::AddPass(GraphBuilder, RDG_EVENT_NAME("SWE_Decomp_Final"), ERDGPassFlags::Compute, Shader, Params, GridGroups);
	}

	// Recompute H
	{
		TShaderMapRef<FRecomputeHCS> Shader(ShaderMap);
		FRecomputeHCS::FParameters* Params = GraphBuilder.AllocParameters<FRecomputeHCS::FParameters>();
		Params->SimConstants = ConstantBuffer;
		Params->hIn = GraphBuilder.CreateSRV(h_RDG);
		Params->terrain = GraphBuilder.CreateSRV(Terrain_RDG);
		Params->H_Out = GraphBuilder.CreateUAV(H_RDG);

		FComputeShaderUtils::AddPass(GraphBuilder, RDG_EVENT_NAME("SWE_Decomp_ReH"), ERDGPassFlags::Compute, Shader, Params, GridGroups);
	}

	// ----------------------------------------------------
	// FFT WIND WAVE STEP
	// ----------------------------------------------------

	// Propagate wind waves
	{
		TShaderMapRef<FPropagateWavesCS> Shader(ShaderMap);
		FPropagateWavesCS::FParameters* Params = GraphBuilder.AllocParameters<FPropagateWavesCS::FParameters>();
		Params->SimConstants = ConstantBuffer;
		Params->HPosIn = GraphBuilder.CreateSRV(HPos_RDG);
		Params->HNegIn = GraphBuilder.CreateSRV(HNeg_RDG);
		Params->DispXOut = GraphBuilder.CreateUAV(Disp_x_RDG);
		Params->DispYOut = GraphBuilder.CreateUAV(Disp_y_RDG);
		Params->DelHXOut = GraphBuilder.CreateUAV(DelH_x_RDG);
		Params->DelHYOut = GraphBuilder.CreateUAV(DelH_y_RDG);
		Params->FlowXOut = GraphBuilder.CreateUAV(Flow_x_RDG);
		Params->FlowYOut = GraphBuilder.CreateUAV(Flow_y_RDG);

		FComputeShaderUtils::AddPass(GraphBuilder, RDG_EVENT_NAME("SWE_FFTWaves_Propagate"), ERDGPassFlags::Compute, Shader, Params, ComplexArrayGroups);
	}

	// ----------------------------------------------------
	// EWAVE DISPERSION STEP
	// ----------------------------------------------------

	// Transfer variables to fourier domain (reading htildeOld as SRV, writing htildeOldNext as UAV)
	{
		TShaderMapRef<FTransferToFFTCS> Shader(ShaderMap);
		FTransferToFFTCS::FParameters* Params = GraphBuilder.AllocParameters<FTransferToFFTCS::FParameters>();
		Params->SimConstants = ConstantBuffer;
		Params->htildeIn = GraphBuilder.CreateSRV(htilde_RDG);
		Params->htildeOldIn = GraphBuilder.CreateSRV(htildeOld_RDG);
		Params->qtildeIn_x = GraphBuilder.CreateSRV(qtildex_RDG);
		Params->qtildeIn_y = GraphBuilder.CreateSRV(qtildey_RDG);
		Params->htildeOldNext = GraphBuilder.CreateUAV(htildeOldNext_RDG);
		Params->hHat = GraphBuilder.CreateUAV(hHat_RDG);
		Params->qHat_x = GraphBuilder.CreateUAV(qHat_x_RDG);
		Params->qHat_y = GraphBuilder.CreateUAV(qHat_y_RDG);

		FComputeShaderUtils::AddPass(GraphBuilder, RDG_EVENT_NAME("SWE_eWave_Transfer"), ERDGPassFlags::Compute, Shader, Params, PaddedGroups);

		// Ping-pong swap persistent pointers for htildeOld
		Swap(TexhtildeOld, TexhtildeOldNext);
		Swap(htildeOld_RDG, htildeOldNext_RDG);
	}

	// Run Forward FFTs
	DispatchFFT_RenderThread(GraphBuilder, hHat_RDG, PaddedSizeX, PaddedSizeY, false, 1);
	DispatchFFT_RenderThread(GraphBuilder, qHat_x_RDG, PaddedSizeX, PaddedSizeY, false, 1);
	DispatchFFT_RenderThread(GraphBuilder, qHat_y_RDG, PaddedSizeX, PaddedSizeY, false, 1);

	// Compute eWave dispersion updates
	{
		TShaderMapRef<FCalcEWaveCS> Shader(ShaderMap);
		FCalcEWaveCS::FParameters* Params = GraphBuilder.AllocParameters<FCalcEWaveCS::FParameters>();
		Params->SimConstants = ConstantBuffer;
		Params->hhat = GraphBuilder.CreateSRV(hHat_RDG);
		Params->qhat_x = GraphBuilder.CreateSRV(qHat_x_RDG);
		Params->qhat_y = GraphBuilder.CreateSRV(qHat_y_RDG);
		Params->Flow_x = GraphBuilder.CreateSRV(Flow_x_RDG);
		Params->Flow_y = GraphBuilder.CreateSRV(Flow_y_RDG);
		Params->qhat_x_array = GraphBuilder.CreateUAV(qHat_x_array_RDG);
		Params->qhat_y_array = GraphBuilder.CreateUAV(qHat_y_array_RDG);

		FComputeShaderUtils::AddPass(GraphBuilder, RDG_EVENT_NAME("SWE_eWave_Calc"), ERDGPassFlags::Compute, Shader, Params, ComplexArrayGroups);
	}

	// Run Inverse FFTs
	DispatchFFT_RenderThread(GraphBuilder, qHat_x_array_RDG, PaddedSizeX, PaddedSizeY, true, Constants.depthNum);
	DispatchFFT_RenderThread(GraphBuilder, qHat_y_array_RDG, PaddedSizeX, PaddedSizeY, true, Constants.depthNum);

	// Interpolate between depths to get new qtilde
	{
		TShaderMapRef<FInterpQCS> Shader(ShaderMap);
		FInterpQCS::FParameters* Params = GraphBuilder.AllocParameters<FInterpQCS::FParameters>();
		Params->SimConstants = ConstantBuffer;
		Params->hbarIn = GraphBuilder.CreateSRV(hbar_RDG);
		Params->qHat_x_array = GraphBuilder.CreateSRV(qHat_x_array_RDG);
		Params->qHat_y_array = GraphBuilder.CreateSRV(qHat_y_array_RDG);
		Params->qtildeOut_x = GraphBuilder.CreateUAV(qtildex_RDG);
		Params->qtildeOut_y = GraphBuilder.CreateUAV(qtildey_RDG);

		FComputeShaderUtils::AddPass(GraphBuilder, RDG_EVENT_NAME("SWE_eWave_InterpQ"), ERDGPassFlags::Compute, Shader, Params, GridGroups);
	}

	// ----------------------------------------------------
	// SWE BULK STEP
	// ----------------------------------------------------

	// Run Inverse FFTs on Disp and DelH fields
	DispatchFFT_RenderThread(GraphBuilder, Disp_x_RDG, PaddedSizeX, PaddedSizeY, true, Constants.depthNum);
	DispatchFFT_RenderThread(GraphBuilder, Disp_y_RDG, PaddedSizeX, PaddedSizeY, true, Constants.depthNum);
	DispatchFFT_RenderThread(GraphBuilder, DelH_x_RDG, PaddedSizeX, PaddedSizeY, true, Constants.depthNum);
	DispatchFFT_RenderThread(GraphBuilder, DelH_y_RDG, PaddedSizeX, PaddedSizeY, true, Constants.depthNum);

	// Interpolate wind wave outputs (Disp -> disp, DelH -> delH) between depths
	{
		TShaderMapRef<FInterpCS> Shader(ShaderMap);
		FInterpCS::FParameters* Params = GraphBuilder.AllocParameters<FInterpCS::FParameters>();
		Params->SimConstants = ConstantBuffer;
		Params->HxIn = GraphBuilder.CreateSRV(DelH_x_RDG);
		Params->HyIn = GraphBuilder.CreateSRV(DelH_y_RDG);
		Params->DxIn = GraphBuilder.CreateSRV(Disp_x_RDG);
		Params->DyIn = GraphBuilder.CreateSRV(Disp_y_RDG);
		Params->hbarIn = GraphBuilder.CreateSRV(hbar_RDG);
		Params->HxOut = GraphBuilder.CreateUAV(delH_x_RDG);
		Params->HyOut = GraphBuilder.CreateUAV(delH_y_RDG);
		Params->DxOut = GraphBuilder.CreateUAV(disp_x_RDG);
		Params->DyOut = GraphBuilder.CreateUAV(disp_y_RDG);

		FComputeShaderUtils::AddPass(GraphBuilder, RDG_EVENT_NAME("SWE_FFTWaves_Interp"), ERDGPassFlags::Compute, Shader, Params, GridGroups);
	}

	// CalcUbar
	{
		TShaderMapRef<FCalcUbarCS> Shader(ShaderMap);
		FCalcUbarCS::FParameters* Params = GraphBuilder.AllocParameters<FCalcUbarCS::FParameters>();
		Params->SimConstants = ConstantBuffer;
		Params->qbarIn_x = GraphBuilder.CreateSRV(qbarx_RDG);
		Params->qbarIn_y = GraphBuilder.CreateSRV(qbary_RDG);
		Params->hbarIn = GraphBuilder.CreateSRV(hbarOld_RDG);
		Params->ubarOut_x = GraphBuilder.CreateUAV(ubarx_RDG);
		Params->ubarOut_y = GraphBuilder.CreateUAV(ubary_RDG);

		FComputeShaderUtils::AddPass(GraphBuilder, RDG_EVENT_NAME("SWE_CalcUbar"), ERDGPassFlags::Compute, Shader, Params, GridGroups);
	}

	// CalcSWE
	{
		TShaderMapRef<FCalcSWECS> Shader(ShaderMap);
		FCalcSWECS::FParameters* Params = GraphBuilder.AllocParameters<FCalcSWECS::FParameters>();
		Params->SimConstants = ConstantBuffer;
		Params->ubarIn_x = GraphBuilder.CreateSRV(ubarx_RDG);
		Params->ubarIn_y = GraphBuilder.CreateSRV(ubary_RDG);
		Params->hbarIn = GraphBuilder.CreateSRV(hbar_RDG);
		Params->H_In = GraphBuilder.CreateSRV(H_RDG);
		Params->delH_x = GraphBuilder.CreateSRV(delH_x_RDG);
		Params->delH_y = GraphBuilder.CreateSRV(delH_y_RDG);
		Params->terrain = GraphBuilder.CreateSRV(Terrain_RDG);
		Params->ubarNewOut_x = GraphBuilder.CreateUAV(ubarNew_x_RDG);
		Params->ubarNewOut_y = GraphBuilder.CreateUAV(ubarNew_y_RDG);
		Params->qbarOut_x = GraphBuilder.CreateUAV(qbarx_RDG);
		Params->qbarOut_y = GraphBuilder.CreateUAV(qbary_RDG);

		FComputeShaderUtils::AddPass(GraphBuilder, RDG_EVENT_NAME("SWE_CalcSWE"), ERDGPassFlags::Compute, Shader, Params, GridGroups);
	}

	// Ping-pong swap hbar and hbarOld pointers (matching Sim2D.cpp std::swap(hbar, hbarOld))
	Swap(Texhbar, TexhbarOld);
	Swap(hbar_RDG, hbarOld_RDG);

	// ----------------------------------------------------
	// TRANSPORT & COMBINE STEP
	// ----------------------------------------------------

	// Ping-pong swap pointers for qtilde and htilde before UpdateTilde (matching Sim2D.cpp std::swap(qtilde_x, qtildePast_x), etc.)
	Swap(Texqtilde_x, TexqtildePast_x);
	Swap(qtildex_RDG, qtildePast_x_RDG);
	Swap(Texqtilde_y, TexqtildePast_y);
	Swap(qtildey_RDG, qtildePast_y_RDG);
	Swap(Texhtilde, TexhtildePast);
	Swap(htilde_RDG, htildePast_RDG);

	// UpdateTilde (Advect wave height and flow rate)
	{
		TShaderMapRef<FUpdateTildeCS> Shader(ShaderMap);
		FUpdateTildeCS::FParameters* Params = GraphBuilder.AllocParameters<FUpdateTildeCS::FParameters>();
		Params->SimConstants = ConstantBuffer;
		Params->ubarNewIn_x = GraphBuilder.CreateSRV(ubarNew_x_RDG);
		Params->ubarIn_x = GraphBuilder.CreateSRV(ubarx_RDG);
		Params->ubarNewIn_y = GraphBuilder.CreateSRV(ubarNew_y_RDG);
		Params->ubarIn_y = GraphBuilder.CreateSRV(ubary_RDG);
		Params->qtildePast_x = GraphBuilder.CreateSRV(qtildePast_x_RDG);
		Params->qtildePast_y = GraphBuilder.CreateSRV(qtildePast_y_RDG);
		Params->hIn = GraphBuilder.CreateSRV(h_RDG);
		Params->htildePast = GraphBuilder.CreateSRV(htildePast_RDG);
		Params->htildeOut = GraphBuilder.CreateUAV(htilde_RDG);
		Params->qtildeOut_x = GraphBuilder.CreateUAV(qtildex_RDG);
		Params->qtildeOut_y = GraphBuilder.CreateUAV(qtildey_RDG);

		FComputeShaderUtils::AddPass(GraphBuilder, RDG_EVENT_NAME("SWE_UpdateTilde"), ERDGPassFlags::Compute, Shader, Params, GridGroups);
	}

	// CalcQAdvect
	{
		TShaderMapRef<FCalcQAdvectCS> Shader(ShaderMap);
		FCalcQAdvectCS::FParameters* Params = GraphBuilder.AllocParameters<FCalcQAdvectCS::FParameters>();
		Params->SimConstants = ConstantBuffer;
		Params->ubarNewIn_x = GraphBuilder.CreateSRV(ubarNew_x_RDG);
		Params->ubarNewIn_y = GraphBuilder.CreateSRV(ubarNew_y_RDG);
		Params->htildeIn = GraphBuilder.CreateSRV(htilde_RDG);
		Params->qAdvectOut_x = GraphBuilder.CreateUAV(qAdvect_x_RDG);
		Params->qAdvectOut_y = GraphBuilder.CreateUAV(qAdvect_y_RDG);

		FComputeShaderUtils::AddPass(GraphBuilder, RDG_EVENT_NAME("SWE_CalcQAdvect"), ERDGPassFlags::Compute, Shader, Params, GridGroups);
	}

	// Ping-pong swap h and hPast pointers before IntegrateH (matching Sim2D.cpp std::swap(h, hPast))
	Swap(Texh, TexhPast);
	Swap(h_RDG, hPast_RDG);

	// IntegrateH
	{
		TShaderMapRef<FIntegrateHCS> Shader(ShaderMap);
		FIntegrateHCS::FParameters* Params = GraphBuilder.AllocParameters<FIntegrateHCS::FParameters>();
		Params->SimConstants = ConstantBuffer;
		Params->qbarIn_x = GraphBuilder.CreateSRV(qbarx_RDG);
		Params->qtildeIn_x = GraphBuilder.CreateSRV(qtildex_RDG);
		Params->qAdvectIn_x = GraphBuilder.CreateSRV(qAdvect_x_RDG);
		Params->qbarIn_y = GraphBuilder.CreateSRV(qbary_RDG);
		Params->qtildeIn_y = GraphBuilder.CreateSRV(qtildey_RDG);
		Params->qAdvectIn_y = GraphBuilder.CreateSRV(qAdvect_y_RDG);
		Params->hPast = GraphBuilder.CreateSRV(hPast_RDG);
		Params->terrain = GraphBuilder.CreateSRV(Terrain_RDG);
		Params->hOut = GraphBuilder.CreateUAV(h_RDG);
		Params->qOut_x = GraphBuilder.CreateUAV(qx_RDG);
		Params->qOut_y = GraphBuilder.CreateUAV(qy_RDG);

		FComputeShaderUtils::AddPass(GraphBuilder, RDG_EVENT_NAME("SWE_IntegrateH"), ERDGPassFlags::Compute, Shader, Params, GridGroups);
	}

	// ----------------------------------------------------
	// EXPORT VISUAL OUTPUTS TO RENDER TARGETS
	// ----------------------------------------------------

	// Recompute total elevation H for rendering
	{
		TShaderMapRef<FRecomputeHCS> Shader(ShaderMap);
		FRecomputeHCS::FParameters* Params = GraphBuilder.AllocParameters<FRecomputeHCS::FParameters>();
		Params->SimConstants = ConstantBuffer;
		Params->hIn = GraphBuilder.CreateSRV(h_RDG);
		Params->terrain = GraphBuilder.CreateSRV(Terrain_RDG);
		Params->H_Out = GraphBuilder.CreateUAV(H_RDG);

		FComputeShaderUtils::AddPass(GraphBuilder, RDG_EVENT_NAME("SWE_Final_ReH"), ERDGPassFlags::Compute, Shader, Params, GridGroups);
	}

	// Copy current Displacement to DisplacementPast before updating
	if (DisplacementPastRT && DisplacementPastRT->GetRenderTargetResource() &&
		DisplacementRT && DisplacementRT->GetRenderTargetResource())
	{
		FRDGTextureRef SrcRDG = GraphBuilder.RegisterExternalTexture(CreateRenderTarget(DisplacementRT->GetRenderTargetResource()->GetTexture2DRHI(), TEXT("DispCurrent_CopySrc")));
		FRDGTextureRef DestRDG = GraphBuilder.RegisterExternalTexture(CreateRenderTarget(DisplacementPastRT->GetRenderTargetResource()->GetTexture2DRHI(), TEXT("DispPast_CopyDest")));
		AddCopyTexturePass(GraphBuilder, SrcRDG, DestRDG);
	}

	// Export Displacement: Combine dispX, dispY, and height H_RDG into a single PF_FloatRGBA Render Target
	if (DisplacementRT && DisplacementRT->GetRenderTargetResource() && H_RDG)
	{
		FRDGTextureRef ExportDispDest = GraphBuilder.RegisterExternalTexture(CreateRenderTarget(DisplacementRT->GetRenderTargetResource()->GetTexture2DRHI(), TEXT("DispExport")));

		TShaderMapRef<FScaleCopyDisplacementCS> ScaleCopyCS(GetGlobalShaderMap(GMaxRHIFeatureLevel));
		FScaleCopyDisplacementCS::FParameters* PassParams = GraphBuilder.AllocParameters<FScaleCopyDisplacementCS::FParameters>();
		PassParams->SimConstants = ConstantBuffer;
		PassParams->ScaleFactor = 100.0f; // m to cm
		PassParams->inDispX = GraphBuilder.CreateSRV(disp_x_RDG);
		PassParams->inDispY = GraphBuilder.CreateSRV(disp_y_RDG);
		PassParams->inHeight = GraphBuilder.CreateSRV(H_RDG);
		PassParams->outDisp4 = GraphBuilder.CreateUAV(ExportDispDest);

		FComputeShaderUtils::AddPass(
			GraphBuilder,
			RDG_EVENT_NAME("SWE_ExportDisp_Scale"),
			ERDGPassFlags::Compute,
			ScaleCopyCS,
			PassParams,
			FIntVector(FMath::DivideAndRoundUp(GridSizeX, 16), FMath::DivideAndRoundUp(GridSizeY, 16), 1)
		);
	}

	// Import persistent foam texture states
	FRDGTextureRef PreviousFoam_RDG = GraphBuilder.RegisterExternalTexture(TexFoam);
	FRDGTextureRef NewFoamRDG = GraphBuilder.RegisterExternalTexture(TexNewFoam);

	// Export Surface Normal, Foam & JacobianDet: Calculate combined fields using unified compute shader
	if (NormalRT && NormalRT->GetRenderTargetResource() &&
		FoamRT && FoamRT->GetRenderTargetResource() &&
		JacobianDetRT && JacobianDetRT->GetRenderTargetResource())
	{
		FRDGTextureRef ExportNormalDest = GraphBuilder.RegisterExternalTexture(CreateRenderTarget(NormalRT->GetRenderTargetResource()->GetTexture2DRHI(), TEXT("NormalExport")));
		FRDGTextureRef ExportFoamDest = GraphBuilder.RegisterExternalTexture(CreateRenderTarget(FoamRT->GetRenderTargetResource()->GetTexture2DRHI(), TEXT("FoamExport")));
		FRDGTextureRef ExportJacobianDest = GraphBuilder.RegisterExternalTexture(CreateRenderTarget(JacobianDetRT->GetRenderTargetResource()->GetTexture2DRHI(), TEXT("JacobianExport")));

		TShaderMapRef<FCalcSurfaceNormalAndFoamCS> NormalAndFoamCS(GetGlobalShaderMap(GMaxRHIFeatureLevel));
		FCalcSurfaceNormalAndFoamCS::FParameters* PassParams = GraphBuilder.AllocParameters<FCalcSurfaceNormalAndFoamCS::FParameters>();
		PassParams->SimConstants = ConstantBuffer;
		PassParams->inDispX = GraphBuilder.CreateSRV(disp_x_RDG);
		PassParams->inDispY = GraphBuilder.CreateSRV(disp_y_RDG);
		PassParams->inHeight = GraphBuilder.CreateSRV(H_RDG);
		PassParams->inPreviousFoam = GraphBuilder.CreateSRV(PreviousFoam_RDG);
		PassParams->outNormal = GraphBuilder.CreateUAV(ExportNormalDest);
		PassParams->outFoam = GraphBuilder.CreateUAV(NewFoamRDG);
		PassParams->outJacobianDet = GraphBuilder.CreateUAV(ExportJacobianDest);

		FComputeShaderUtils::AddPass(
			GraphBuilder,
			RDG_EVENT_NAME("SWE_CalcSurfaceNormalAndFoam"),
			ERDGPassFlags::Compute,
			NormalAndFoamCS,
			PassParams,
			FIntVector(FMath::DivideAndRoundUp(GridSizeX, 16), FMath::DivideAndRoundUp(GridSizeY, 16), 1)
		);

		// Copy the updated foam values to the user-facing FoamRT Render Target
		AddCopyTexturePass(GraphBuilder, NewFoamRDG, ExportFoamDest);

		// Ping-pong swap persistent foam textures for next frame
		Swap(TexFoam, TexNewFoam);
	}

	// ----------------------------------------------------
	// CALCULATE ROUGHNESS LOOK-UP TABLE (LUT)
	// ----------------------------------------------------
	if (TexRoughness.IsValid())
	{
		FRDGTextureRef PreviousRoughness_RDG = GraphBuilder.RegisterExternalTexture(TexRoughness);
		FRDGTextureRef NewRoughnessRDG = GraphBuilder.RegisterExternalTexture(TexNewRoughness);

		if (RoughnessRT && RoughnessRT->GetRenderTargetResource() && NormalRT && NormalRT->GetRenderTargetResource())
		{
			FRDGTextureRef ExportRoughnessDest = GraphBuilder.RegisterExternalTexture(CreateRenderTarget(RoughnessRT->GetRenderTargetResource()->GetTexture2DRHI(), TEXT("RoughnessExport")));
			FRDGTextureRef CurrentNormal_RDG = GraphBuilder.RegisterExternalTexture(CreateRenderTarget(NormalRT->GetRenderTargetResource()->GetTexture2DRHI(), TEXT("NormalForRoughness")));

			TShaderMapRef<FCalcRoughnessLUTCS> RoughnessCS(GetGlobalShaderMap(GMaxRHIFeatureLevel));
			FCalcRoughnessLUTCS::FParameters* PassParams = GraphBuilder.AllocParameters<FCalcRoughnessLUTCS::FParameters>();
			PassParams->SimConstants = ConstantBuffer;
			PassParams->IntegrationSamples = IntegrationSamples;
			PassParams->RoughnessPower = RoughnessPower;
			PassParams->inNormal = GraphBuilder.CreateSRV(CurrentNormal_RDG);
			PassParams->inPreviousRoughness = GraphBuilder.CreateSRV(PreviousRoughness_RDG);
			PassParams->BilinearWrapSampler = TStaticSamplerState<SF_Bilinear, AM_Wrap, AM_Wrap, AM_Wrap>::GetRHI();
			PassParams->outRoughness = GraphBuilder.CreateUAV(NewRoughnessRDG);

			FComputeShaderUtils::AddPass(
				GraphBuilder,
				RDG_EVENT_NAME("SWE_CalcRoughnessLUT"),
				ERDGPassFlags::Compute,
				RoughnessCS,
				PassParams,
				FIntVector(FMath::DivideAndRoundUp(GridSizeX, 16), 1, 1)
			);

			AddCopyTexturePass(GraphBuilder, NewRoughnessRDG, ExportRoughnessDest);

			// Ping-pong swap persistent roughness textures for next frame
			Swap(TexRoughness, TexNewRoughness);
		}
	}

	// ----------------------------------------------------
	// ASYNC READBACK OF WATER HEIGHT FIELD (TexH)
	// ----------------------------------------------------
	if (StagingTextures[StagingWriteIndex].IsValid() && H_RDG)
	{
		FTextureRHIRef CurrentStagingTexture = StagingTextures[StagingWriteIndex];
		
		TRefCountPtr<IPooledRenderTarget> StagingPooled = CreateRenderTarget(CurrentStagingTexture, TEXT("StagingTexture"));
		FRDGTextureRef StagingRDG = GraphBuilder.RegisterExternalTexture(StagingPooled);
		AddCopyTexturePass(GraphBuilder, H_RDG, StagingRDG);
	}

	if (StagingTextures[StagingReadIndex].IsValid())
	{
		FTextureRHIRef ReadStagingTexture = StagingTextures[StagingReadIndex];
		int32 TargetCPUIndex = 1 - ActiveCPUBufferIndex.load();

		GraphBuilder.AddPass(
			RDG_EVENT_NAME("SWE_ReadbackMap"),
			ERDGPassFlags::None,
			[ReadStagingTexture, TargetCPUIndex, this](FRHICommandListImmediate& RHICmdList)
			{
				void* LocalData = nullptr;
				int32 OutWidth = 0;
				int32 OutHeight = 0;
				RHICmdList.MapStagingSurface(ReadStagingTexture, LocalData, OutWidth, OutHeight);
				if (LocalData)
				{
					float* FloatData = (float*)LocalData;
					int32 Size = GridSizeX * GridSizeY;
					
					TArray<float>& DestArray = CPUHeightData[TargetCPUIndex];
					DestArray.SetNumUninitialized(Size);
					
					uint8* LocalByteData = (uint8*)LocalData;
					for (int32 y = 0; y < GridSizeY; ++y)
					{
						float* SrcRow = (float*)(LocalByteData + y * OutWidth);
						float* DestRow = DestArray.GetData() + y * GridSizeX;
						FMemory::Memcpy(DestRow, SrcRow, GridSizeX * sizeof(float));
					}
					
					RHICmdList.UnmapStagingSurface(ReadStagingTexture);
					
					// Swap indices
					ActiveCPUBufferIndex.store(TargetCPUIndex);
				}
			}
		);
	}

	// Swap indices for next frame
	StagingWriteIndex = (StagingWriteIndex + 1) % 2;
	StagingReadIndex = (StagingReadIndex + 1) % 2;

	GraphBuilder.Execute();
}

void UDispersiveSWESimulator::DispatchFFT_RenderThread(
	FRDGBuilder& GraphBuilder,
	FRDGTextureRef TargetTexture,
	int32 SizeX,
	int32 SizeY,
	bool bInverse,
	int32 NumLayers) 
{
	uint32 BitsX = 0, BitsY = 0;
	int32 TmpX = SizeX;
	while (TmpX >>= 1) ++BitsX;
	int32 TmpY = SizeY;
	while (TmpY >>= 1) ++BitsY;

	bool bIsArray = NumLayers > 1;

	// Row Pass: Transform each row
	{
		FFFTKernel1DCS::FPermutationDomain PermutationVector;
		PermutationVector.Set<FFFTKernel1DCS::FFFTSizeDim>(SizeX);
		PermutationVector.Set<FFFTKernel1DCS::FIsArrayDim>(bIsArray);

		TShaderMapRef<FFFTKernel1DCS> FFTShader(GetGlobalShaderMap(GMaxRHIFeatureLevel), PermutationVector);

		FFFTKernel1DCS::FParameters* PassParams = GraphBuilder.AllocParameters<FFFTKernel1DCS::FParameters>();
		PassParams->cb_Nx = SizeX;
		PassParams->cb_Ny = SizeY;
		PassParams->cb_BitsX = BitsX;
		PassParams->cb_BitsY = BitsY;
		PassParams->cb_Inverse = bInverse ? 1 : 0;
		PassParams->cb_IsRow = 1;
		PassParams->fft = GraphBuilder.CreateUAV(TargetTexture);

		FComputeShaderUtils::AddPass(
			GraphBuilder,
			RDG_EVENT_NAME("FFT_RowPass"),
			ERDGPassFlags::Compute,
			FFTShader,
			PassParams,
			FIntVector(1, SizeY, NumLayers)
		);
	}

	// Column Pass: Transform each column
	{
		FFFTKernel1DCS::FPermutationDomain PermutationVector;
		PermutationVector.Set<FFFTKernel1DCS::FFFTSizeDim>(SizeY);
		PermutationVector.Set<FFFTKernel1DCS::FIsArrayDim>(bIsArray);

		TShaderMapRef<FFFTKernel1DCS> FFTShader(GetGlobalShaderMap(GMaxRHIFeatureLevel), PermutationVector);

		FFFTKernel1DCS::FParameters* PassParams = GraphBuilder.AllocParameters<FFFTKernel1DCS::FParameters>();
		PassParams->cb_Nx = SizeX;
		PassParams->cb_Ny = SizeY;
		PassParams->cb_BitsX = BitsX;
		PassParams->cb_BitsY = BitsY;
		PassParams->cb_Inverse = bInverse ? 1 : 0;
		PassParams->cb_IsRow = 0;
		PassParams->fft = GraphBuilder.CreateUAV(TargetTexture);

		FComputeShaderUtils::AddPass(
			GraphBuilder,
			RDG_EVENT_NAME("FFT_ColPass"),
			ERDGPassFlags::Compute,
			FFTShader,
			PassParams,
			FIntVector(1, SizeX, NumLayers)
		);
	}
}

bool UDispersiveSWESimulator::LoadParametersFromJson(const FString& FilePath) 
{
	FString FinalPath = FilePath;
	if (FPaths::IsRelative(FinalPath)) {
		FinalPath = FPaths::Combine(FPaths::ProjectDir(), FilePath);
	}

	FString JsonString;
	if (!FFileHelper::LoadFileToString(JsonString, *FinalPath)) {
		// Try resolving relative to Content folder as a fallback
		FString FallbackPath = FPaths::Combine(FPaths::ProjectContentDir(), FilePath);
		if (!FFileHelper::LoadFileToString(JsonString, *FallbackPath)) {
			UE_LOG(LogTemp, Warning, TEXT("Failed to load JSON file from path: %s (or fallback %s)"), *FinalPath, *FallbackPath);
			return false;
		}
		FinalPath = FallbackPath;
	}

	TSharedPtr<FJsonObject> JsonObject;
	TSharedRef<TJsonReader<>> Reader = TJsonReaderFactory<>::Create(JsonString);

	if (!FJsonSerializer::Deserialize(Reader, JsonObject) || !JsonObject.IsValid()) {
		UE_LOG(LogTemp, Warning, TEXT("Failed to deserialize JSON string."));
		return false;
	}

	// Map JSON fields directly to component properties
	if (!FJsonObjectConverter::JsonObjectToUStruct(JsonObject.ToSharedRef(), GetClass(), this)) {
		UE_LOG(LogTemp, Warning, TEXT("Failed to map JSON object to component properties."));
		return false;
	}

	if (JsonObject->HasField(TEXT("CellSize"))) {
		double CustomCellSize = 0.0;
		if (JsonObject->TryGetNumberField(TEXT("CellSize"), CustomCellSize) && CustomCellSize > 0.0) {
			bAutoCalculateCellSize = false;
		}
	}

	UE_LOG(LogTemp, Log, TEXT("Successfully loaded simulation parameters from: %s"), *FinalPath);
	return true;
}

bool UDispersiveSWESimulator::SaveParametersToJson(const FString& FilePath) 
{
	FString FinalPath = FilePath;
	if (FPaths::IsRelative(FinalPath)) {
		FinalPath = FPaths::Combine(FPaths::ProjectDir(), FilePath);
	}

	TSharedRef<FJsonObject> JsonObject = MakeShared<FJsonObject>();
	if (!FJsonObjectConverter::UStructToJsonObject(GetClass(), this, JsonObject)) {
		UE_LOG(LogTemp, Warning, TEXT("Failed to convert component properties to JSON object."));
		return false;
	}

	// Remove runtime properties or input/output targets that shouldn't be serialized in a parameters config
	JsonObject->RemoveField(TEXT("terrainHeightInputRT"));
	JsonObject->RemoveField(TEXT("displacementRT"));
	JsonObject->RemoveField(TEXT("displacementPastRT"));
	JsonObject->RemoveField(TEXT("normalRT"));
	JsonObject->RemoveField(TEXT("foamRT"));
	JsonObject->RemoveField(TEXT("jacobianDetRT"));
	JsonObject->RemoveField(TEXT("roughnessRT"));
	JsonObject->RemoveField(TEXT("jsonConfigFilePath"));

	FString JsonString;
	TSharedRef<TJsonWriter<>> Writer = TJsonWriterFactory<>::Create(&JsonString);
	if (!FJsonSerializer::Serialize(JsonObject, Writer)) {
		UE_LOG(LogTemp, Warning, TEXT("Failed to serialize JSON object."));
		return false;
	}

	if (!FFileHelper::SaveStringToFile(JsonString, *FinalPath)) {
		UE_LOG(LogTemp, Warning, TEXT("Failed to save JSON string to path: %s"), *FinalPath);
		return false;
	}

	UE_LOG(LogTemp, Log, TEXT("Successfully saved simulation parameters to: %s"), *FinalPath);
	return true;
}

float UDispersiveSWESimulator::GetCachedHeight(int32 X, int32 Y) const 
{
	int32 Index = Y * GridSizeX + X;
	int32 ActiveIdx = ActiveCPUBufferIndex.load();
	if (CPUHeightData[ActiveIdx].IsValidIndex(Index)) {
		return CPUHeightData[ActiveIdx][Index];
	}
	return WaterLevel * 0.01f; // base water level in meters
}

float UDispersiveSWESimulator::GetWaterHeightAtLocation(const FVector& WorldLocation) const 
{
	if (GridSizeX <= 1 || GridSizeY <= 1 || CapturedWorldWidth <= 0.0f) {
		return WaterLevel;
	}

	AActor* Owner = GetOwner();
	if (!Owner) {
		return WaterLevel;
	}

	FVector ActorLoc = Owner->GetActorLocation();
	FVector LocalPos = WorldLocation - ActorLoc;

	// Map to UV space [-0.5, 0.5] -> [0, 1]
	float U = (LocalPos.X / CapturedWorldWidth) + 0.5f;
	float V = (LocalPos.Y / CapturedWorldWidth) + 0.5f;

	U = FMath::Clamp(U, 0.0f, 1.0f);
	V = FMath::Clamp(V, 0.0f, 1.0f);

	float GridX = U * (GridSizeX - 1);
	float GridY = V * (GridSizeY - 1);

	int32 X0 = FMath::FloorToInt(GridX);
	int32 Y0 = FMath::FloorToInt(GridY);
	int32 X1 = FMath::Min(X0 + 1, GridSizeX - 1);
	int32 Y1 = FMath::Min(Y0 + 1, GridSizeY - 1);

	float LerpX = GridX - X0;
	float LerpY = GridY - Y0;

	float H00 = GetCachedHeight(X0, Y0);
	float H10 = GetCachedHeight(X1, Y0);
	float H01 = GetCachedHeight(X0, Y1);
	float H11 = GetCachedHeight(X1, Y1);

	float H_Avg = FMath::BiLerp(H00, H10, H01, H11, LerpX, LerpY);

	// Convert meters to centimeters (Unreal units)
	return H_Avg * 100.0f;
}
