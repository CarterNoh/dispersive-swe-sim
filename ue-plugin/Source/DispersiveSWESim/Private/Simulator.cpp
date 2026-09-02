#include "Simulator.h"
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

USimulator::USimulator() 
{
	PrimaryComponentTick.bCanEverTick = true;
	PrimaryComponentTick.bStartWithTickEnabled = true;
}

void USimulator::BeginPlay() 
{
	UE_LOG(LogTemp, Warning, TEXT("USimulator::BeginPlay() called"));
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

void USimulator::InitializeSimulation() 
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

void USimulator::AllocatePersistentTargets(FRHICommandListImmediate& RHICmdList) 
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
	GRenderTargetPool.FindFreeElement(RHICmdList, Desc, Texhdot, TEXT("hdot"));
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
	FRHITextureCreateDesc DispDesc0 = FRHITextureCreateDesc::Create2D(TEXT("Staging0"), GridSizeX, GridSizeY, PF_R32_FLOAT)
		.SetFlags(TexCreate_CPUReadback)
		.SetNumMips(1);
	StagingTextures[0] = RHICreateTexture(DispDesc0);

	FRHITextureCreateDesc DispDesc1 = FRHITextureCreateDesc::Create2D(TEXT("Staging1"), GridSizeX, GridSizeY, PF_R32_FLOAT)
		.SetFlags(TexCreate_CPUReadback)
		.SetNumMips(1);
	StagingTextures[1] = RHICreateTexture(DispDesc1);

	FRHITextureCreateDesc VelDesc0 = FRHITextureCreateDesc::Create2D(TEXT("VelocityStaging0"), GridSizeX, GridSizeY, PF_FloatRGBA)
		.SetFlags(TexCreate_CPUReadback)
		.SetNumMips(1);
	VelocityStagingTextures[0] = RHICreateTexture(VelDesc0);

	FRHITextureCreateDesc VelDesc1 = FRHITextureCreateDesc::Create2D(TEXT("VelocityStaging1"), GridSizeX, GridSizeY, PF_FloatRGBA)
		.SetFlags(TexCreate_CPUReadback)
		.SetNumMips(1);
	VelocityStagingTextures[1] = RHICreateTexture(VelDesc1);

	FRHITextureCreateDesc AccelDesc0 = FRHITextureCreateDesc::Create2D(TEXT("AccelerationStaging0"), GridSizeX, GridSizeY, PF_FloatRGBA)
		.SetFlags(TexCreate_CPUReadback)
		.SetNumMips(1);
	AccelerationStagingTextures[0] = RHICreateTexture(AccelDesc0);

	FRHITextureCreateDesc AccelDesc1 = FRHITextureCreateDesc::Create2D(TEXT("AccelerationStaging1"), GridSizeX, GridSizeY, PF_FloatRGBA)
		.SetFlags(TexCreate_CPUReadback)
		.SetNumMips(1);
	AccelerationStagingTextures[1] = RHICreateTexture(AccelDesc1);

	CPUHeightData[0].Reset();
	CPUHeightData[0].SetNumZeroed(GridSizeX * GridSizeY);
	CPUHeightData[1].Reset();
	CPUHeightData[1].SetNumZeroed(GridSizeX * GridSizeY);

	CPUVelocityData[0].Reset();
	CPUVelocityData[0].SetNumZeroed(GridSizeX * GridSizeY);
	CPUVelocityData[1].Reset();
	CPUVelocityData[1].SetNumZeroed(GridSizeX * GridSizeY);

	CPUAccelerationData[0].Reset();
	CPUAccelerationData[0].SetNumZeroed(GridSizeX * GridSizeY);
	CPUAccelerationData[1].Reset();
	CPUAccelerationData[1].SetNumZeroed(GridSizeX * GridSizeY);

	ActiveCPUBufferIndex.store(0);
	StagingWriteIndex = 0;
	StagingReadIndex = 1;
}

void USimulator::AssignSimulationConstants(FSimConstants& Constants) const
{
	Constants.gridSizeX = GridSizeX;
	Constants.gridSizeY = GridSizeY;
	Constants.cellSize = CellSize * 0.01f; // Convert cm to meters
	Constants.timeStep = TimeStep;
	Constants.minWaterHeight = MinWaterHeight * 0.01f; // Convert cm to meters
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
	Constants.depthCutoff = DepthCutoff;
	Constants.paddedGridSizeX = PaddedSizeX;
	Constants.paddedGridSizeY = PaddedSizeY;
	Constants.maxSafeDepth = CalculatedMaxSafeDepth;
	for (int32 i = 0; i < 16; ++i) {
		float Val = i < DepthLevels.Num() ? DepthLevels[i] : 0.0f;
		Constants.depthLevels[i / 4][i % 4] = Val;
	}
}

void USimulator::AssignFFTWaveConstants(FFFTWaveConstants& Constants) const
{
	Constants.time = SimulationTime;
	Constants.gridSizeX = GridSizeX;
	Constants.gridSizeY = GridSizeY;
	Constants.cellSize = CellSize * 0.01f; // Convert cm to meters
	Constants.paddedGridSizeX = PaddedSizeX;
	Constants.paddedGridSizeY = PaddedSizeY;
	Constants.minWaterHeight = MinWaterHeight * 0.01f; // Convert cm to meters
	Constants.maxSafeDepth = CalculatedMaxSafeDepth;
	Constants.depthNum = DepthLevels.Num();
	Constants.surfaceTension = SurfaceTension;
	Constants.density = Density;
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
	for (int32 i = 0; i < 16; ++i) {
		float Val = i < DepthLevels.Num() ? DepthLevels[i] : 0.0f;
		Constants.depthLevels[i / 4][i % 4] = Val;
	}
}

void USimulator::AssignExportConstants(FExportConstants& Constants) const
{
	Constants.time = SimulationTime;
	Constants.gridSizeX = GridSizeX;
	Constants.gridSizeY = GridSizeY;
	Constants.cellSize = CellSize * 0.01f; // Convert cm to meters
	Constants.timeStep = TimeStep;
	Constants.foamThreshold = FoamThreshold;
	Constants.foamMultiplier = FoamMultiplier;
	Constants.foamFade = FoamFade;
	Constants.foamBlur = FoamBlur;
	Constants.minWaterHeight = MinWaterHeight * 0.01f; // Convert cm to meters
}

void USimulator::SetupInitialStates(FRHICommandListImmediate& RHICmdList) 
{
	// Setup typed uniform buffers
	FSimConstants CPUSimConstants = {};
	AssignSimulationConstants(CPUSimConstants);
	TUniformBufferRef<FSimConstants> SimConstantBuffer = CreateUniformBufferImmediate(CPUSimConstants, EUniformBufferUsage::UniformBuffer_SingleFrame);

	FFFTWaveConstants CPUFFTWaveConstants = {};
	AssignFFTWaveConstants(CPUFFTWaveConstants);
	TUniformBufferRef<FFFTWaveConstants> FFTWaveConstantBuffer = CreateUniformBufferImmediate(CPUFFTWaveConstants, EUniformBufferUsage::UniformBuffer_SingleFrame);

	// Clear all persistent simulation state textures to zero to avoid uninitialized GPU garbage (matching Sim2D.cpp)
	{
		FRDGBuilder GraphBuilder(RHICmdList);

		// 2D Scalar simulation fields
		if (TexQ_x.IsValid()) AddClearUAVPass(GraphBuilder, GraphBuilder.CreateUAV(GraphBuilder.RegisterExternalTexture(TexQ_x)), FLinearColor::Black);
		if (TexQ_y.IsValid()) AddClearUAVPass(GraphBuilder, GraphBuilder.CreateUAV(GraphBuilder.RegisterExternalTexture(TexQ_y)), FLinearColor::Black);
		if (Texq_x.IsValid()) AddClearUAVPass(GraphBuilder, GraphBuilder.CreateUAV(GraphBuilder.RegisterExternalTexture(Texq_x)), FLinearColor::Black);
		if (Texq_y.IsValid()) AddClearUAVPass(GraphBuilder, GraphBuilder.CreateUAV(GraphBuilder.RegisterExternalTexture(Texq_y)), FLinearColor::Black);
		if (Texhdot.IsValid()) AddClearUAVPass(GraphBuilder, GraphBuilder.CreateUAV(GraphBuilder.RegisterExternalTexture(Texhdot)), FLinearColor::Black);
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

		// 2D Complex float2 fields
		if (TexhHat.IsValid()) AddClearUAVPass(GraphBuilder, GraphBuilder.CreateUAV(GraphBuilder.RegisterExternalTexture(TexhHat)), FLinearColor::Black);
		if (TexqHat_x.IsValid()) AddClearUAVPass(GraphBuilder, GraphBuilder.CreateUAV(GraphBuilder.RegisterExternalTexture(TexqHat_x)), FLinearColor::Black);
		if (TexqHat_y.IsValid()) AddClearUAVPass(GraphBuilder, GraphBuilder.CreateUAV(GraphBuilder.RegisterExternalTexture(TexqHat_y)), FLinearColor::Black);

		// 2D Complex float2 array fields
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

		// Foam & Roughness
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

		FInitializeWaterInputs InitInputs;
		InitInputs.TerrainInput = TerrainInput_RDG;
		InitInputs.WaterLevel = WaterLevel * 0.01f; // Convert cm to meters
		InitInputs.TerrainCaptureCameraZ = TerrainCaptureCameraZ;
		InitInputs.GridSizeX = GridSizeX;
		InitInputs.GridSizeY = GridSizeY;

		FInitializeWaterOutputs InitOutputs;
		InitOutputs.hOut = GraphBuilder.RegisterExternalTexture(Texh);
		InitOutputs.H_Out = GraphBuilder.RegisterExternalTexture(TexH);
		InitOutputs.terrainOut = GraphBuilder.RegisterExternalTexture(TexTerrain);
		InitOutputs.terrainOutCM = TerrainExportRDG;
		InitOutputs.hbarOut = GraphBuilder.RegisterExternalTexture(Texhbar);
		InitOutputs.hbarOldOut = GraphBuilder.RegisterExternalTexture(TexhbarOld);
		if (TexFoam.IsValid()) InitOutputs.FoamOut = GraphBuilder.RegisterExternalTexture(TexFoam);
		if (TexRoughness.IsValid()) InitOutputs.RoughnessOut = GraphBuilder.RegisterExternalTexture(TexRoughness);

		AddInitializeWaterPass(GraphBuilder, SimConstantBuffer, InitInputs, InitOutputs);

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

		AddPopulateSpectrumPass(
			GraphBuilder,
			FFTWaveConstantBuffer,
			HPos_RDG,
			HNeg_RDG,
			PaddedSizeX,
			PaddedSizeY,
			DepthLevels.Num()
		);

		GraphBuilder.Execute();
	}

	// Export the captured terrain to a raw float32 file for standalone simulation testing if requested
	if (bExportCapturedTerrain && TexTerrain.IsValid())
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

void USimulator::TickComponent(float DeltaTime, ELevelTick TickType, FActorComponentTickFunction* ThisTickFunction) 
{
	Super::TickComponent(DeltaTime, TickType, ThisTickFunction);

	if (!bInitialized) return;

	SimulationTime += TimeStep;

	FSimConstants SimConstants = {};
	AssignSimulationConstants(SimConstants);

	FFFTWaveConstants FFTWaveConstants = {};
	AssignFFTWaveConstants(FFTWaveConstants);

	FExportConstants ExportConstants = {};
	AssignExportConstants(ExportConstants);

	ENQUEUE_RENDER_COMMAND(ExecuteSWESimulation)(
		[this, SimConstants, FFTWaveConstants, ExportConstants](FRHICommandListImmediate& RHICmdList) {
			ExecuteSimulation_RenderThread(RHICmdList, SimConstants, FFTWaveConstants, ExportConstants);
		}
	);
}

void USimulator::ExecuteSimulation_RenderThread(
	FRHICommandListImmediate& RHICmdList,
	const FSimConstants& SimConstants,
	const FFFTWaveConstants& FFTWaveConstants,
	const FExportConstants& ExportConstants) 
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
	FRDGTextureRef hdot_RDG = GraphBuilder.RegisterExternalTexture(Texhdot);
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

	// Create Uniform Buffers
	TUniformBufferRef<FSimConstants> SimConstantBuffer = CreateUniformBufferImmediate(SimConstants, EUniformBufferUsage::UniformBuffer_SingleFrame);
	TUniformBufferRef<FFFTWaveConstants> FFTWaveConstantBuffer = CreateUniformBufferImmediate(FFTWaveConstants, EUniformBufferUsage::UniformBuffer_SingleFrame);
	TUniformBufferRef<FExportConstants> ExportConstantBuffer = CreateUniformBufferImmediate(ExportConstants, EUniformBufferUsage::UniformBuffer_SingleFrame);

	// ----------------------------------------------------
	// DECOMPOSITION & JACOBI DIFFUSION STEP
	// ----------------------------------------------------
	FDecompositionInputs DecompInputs;
	DecompInputs.hIn = h_RDG;
	DecompInputs.qIn_x = qx_RDG;
	DecompInputs.qIn_y = qy_RDG;
	DecompInputs.terrain = Terrain_RDG;
	DecompInputs.HOrig = HOrig_RDG;
	DecompInputs.QOrig_x = QOrig_x_RDG;
	DecompInputs.QOrig_y = QOrig_y_RDG;
	DecompInputs.alpha_H = alpha_H_RDG;
	DecompInputs.alpha_Q_x = alpha_Q_x_RDG;
	DecompInputs.alpha_Q_y = alpha_Q_y_RDG;
	DecompInputs.GridSizeX = GridSizeX;
	DecompInputs.GridSizeY = GridSizeY;
	DecompInputs.DiffusionIterations = SimConstants.diffusionIterations;

	FDecompositionOutputs DecompOutputs;
	DecompOutputs.H_SrcDst = H_RDG;
	DecompOutputs.HPast_SrcDst = HPast_RDG;
	DecompOutputs.Qx_SrcDst = Qx_RDG;
	DecompOutputs.QPast_x_SrcDst = QPast_x_RDG;
	DecompOutputs.Qy_SrcDst = Qy_RDG;
	DecompOutputs.QPast_y_SrcDst = QPast_y_RDG;
	DecompOutputs.hbarOut = hbar_RDG;
	DecompOutputs.qbarOut_x = qbarx_RDG;
	DecompOutputs.qbarOut_y = qbary_RDG;
	DecompOutputs.htildeOut = htilde_RDG;
	DecompOutputs.qtildeOut_x = qtildex_RDG;
	DecompOutputs.qtildeOut_y = qtildey_RDG;

	AddDecompositionPasses(GraphBuilder, SimConstantBuffer, DecompInputs, DecompOutputs);

	// If Jacobi diffusion ended on Past buffer (odd iteration count), swap persistent pointers
	if (DecompOutputs.H_SrcDst != H_RDG) Swap(TexH, TexHPast);
	if (DecompOutputs.Qx_SrcDst != Qx_RDG) Swap(TexQ_x, TexQPast_x);
	if (DecompOutputs.Qy_SrcDst != Qy_RDG) Swap(TexQ_y, TexQPast_y);
	H_RDG = DecompOutputs.H_SrcDst;
	Qx_RDG = DecompOutputs.Qx_SrcDst;
	Qy_RDG = DecompOutputs.Qy_SrcDst;

	// ----------------------------------------------------
	// FFT WIND WAVE STEP
	// ----------------------------------------------------
	FPropagateFFTWavesInputs FFTWaveInputs;
	FFTWaveInputs.HPosIn = HPos_RDG;
	FFTWaveInputs.HNegIn = HNeg_RDG;
	FFTWaveInputs.Disp_x_Array = Disp_x_RDG;
	FFTWaveInputs.Disp_y_Array = Disp_y_RDG;
	FFTWaveInputs.DelH_x_Array = DelH_x_RDG;
	FFTWaveInputs.DelH_y_Array = DelH_y_RDG;
	FFTWaveInputs.Flow_x_Array = Flow_x_RDG;
	FFTWaveInputs.Flow_y_Array = Flow_y_RDG;
	FFTWaveInputs.hbarIn = hbar_RDG;
	FFTWaveInputs.PaddedSizeX = PaddedSizeX;
	FFTWaveInputs.PaddedSizeY = PaddedSizeY;
	FFTWaveInputs.GridSizeX = GridSizeX;
	FFTWaveInputs.GridSizeY = GridSizeY;
	FFTWaveInputs.DepthNum = SimConstants.depthNum;

	FPropagateFFTWavesOutputs FFTWaveOutputs;
	FFTWaveOutputs.disp_x_Out = disp_x_RDG;
	FFTWaveOutputs.disp_y_Out = disp_y_RDG;
	FFTWaveOutputs.delH_x_Out = delH_x_RDG;
	FFTWaveOutputs.delH_y_Out = delH_y_RDG;

	AddPropagateFFTWavesPasses(GraphBuilder, FFTWaveConstantBuffer, FFTWaveInputs, FFTWaveOutputs);

	// ----------------------------------------------------
	// EWAVE DISPERSION STEP
	// ----------------------------------------------------
	FEWaveInputs EWaveInputs;
	EWaveInputs.htildeIn = htilde_RDG;
	EWaveInputs.htildeOldIn = htildeOld_RDG;
	EWaveInputs.qtildeIn_x = qtildex_RDG;
	EWaveInputs.qtildeIn_y = qtildey_RDG;
	EWaveInputs.Flow_x = Flow_x_RDG;
	EWaveInputs.Flow_y = Flow_y_RDG;
	EWaveInputs.hbarIn = hbar_RDG;
	EWaveInputs.hHat = hHat_RDG;
	EWaveInputs.qHat_x = qHat_x_RDG;
	EWaveInputs.qHat_y = qHat_y_RDG;
	EWaveInputs.qHat_x_array = qHat_x_array_RDG;
	EWaveInputs.qHat_y_array = qHat_y_array_RDG;
	EWaveInputs.PaddedSizeX = PaddedSizeX;
	EWaveInputs.PaddedSizeY = PaddedSizeY;
	EWaveInputs.GridSizeX = GridSizeX;
	EWaveInputs.GridSizeY = GridSizeY;
	EWaveInputs.DepthNum = SimConstants.depthNum;

	FEWaveOutputs EWaveOutputs;
	EWaveOutputs.htildeOldNext = htildeOldNext_RDG;
	EWaveOutputs.qtildeOut_x = qtildex_RDG;
	EWaveOutputs.qtildeOut_y = qtildey_RDG;

	AddEWavePasses(GraphBuilder, SimConstantBuffer, EWaveInputs, EWaveOutputs);

	// Ping-pong swap persistent pointers for htildeOld
	Swap(TexhtildeOld, TexhtildeOldNext);
	Swap(htildeOld_RDG, htildeOldNext_RDG);

	// ----------------------------------------------------
	// BULK FLOW (SWE VELOCITY & MOMENTUM)
	// ----------------------------------------------------
	FSWEBulkInputs SWEInputs;
	SWEInputs.qbarIn_x = qbarx_RDG;
	SWEInputs.qbarIn_y = qbary_RDG;
	SWEInputs.hbarIn = hbar_RDG;
	SWEInputs.hbarOldIn = hbarOld_RDG;
	SWEInputs.H_In = H_RDG;
	SWEInputs.delH_x = delH_x_RDG;
	SWEInputs.delH_y = delH_y_RDG;
	SWEInputs.terrain = Terrain_RDG;
	SWEInputs.GridSizeX = GridSizeX;
	SWEInputs.GridSizeY = GridSizeY;

	FSWEBulkOutputs SWEOutputs;
	SWEOutputs.ubarOut_x = ubarx_RDG;
	SWEOutputs.ubarOut_y = ubary_RDG;
	SWEOutputs.ubarNewOut_x = ubarNew_x_RDG;
	SWEOutputs.ubarNewOut_y = ubarNew_y_RDG;
	SWEOutputs.qbarOut_x = qbarx_RDG;
	SWEOutputs.qbarOut_y = qbary_RDG;

	AddSWEBulkPasses(GraphBuilder, SimConstantBuffer, SWEInputs, SWEOutputs);

	// Ping-pong swap hbar and hbarOld pointers
	Swap(Texhbar, TexhbarOld);
	Swap(hbar_RDG, hbarOld_RDG);

	// ----------------------------------------------------
	// ADVECTIVE TRANSPORT & HEIGHT INTEGRATION
	// ----------------------------------------------------
	// Ping-pong swap pointers for qtilde, htilde, and h before transport
	Swap(Texqtilde_x, TexqtildePast_x);
	Swap(qtildex_RDG, qtildePast_x_RDG);
	Swap(Texqtilde_y, TexqtildePast_y);
	Swap(qtildey_RDG, qtildePast_y_RDG);
	Swap(Texhtilde, TexhtildePast);
	Swap(htilde_RDG, htildePast_RDG);
	Swap(Texh, TexhPast);
	Swap(h_RDG, hPast_RDG);

	FTransportAndIntegrateInputs TransportInputs;
	TransportInputs.ubarNewIn_x = ubarNew_x_RDG;
	TransportInputs.ubarIn_x = ubarx_RDG;
	TransportInputs.ubarNewIn_y = ubarNew_y_RDG;
	TransportInputs.ubarIn_y = ubary_RDG;
	TransportInputs.qtildePast_x = qtildePast_x_RDG;
	TransportInputs.qtildePast_y = qtildePast_y_RDG;
	TransportInputs.htildePast = htildePast_RDG;
	TransportInputs.hPast = hPast_RDG;
	TransportInputs.terrain = Terrain_RDG;
	TransportInputs.qbarIn_x = qbarx_RDG;
	TransportInputs.qbarIn_y = qbary_RDG;
	TransportInputs.GridSizeX = GridSizeX;
	TransportInputs.GridSizeY = GridSizeY;

	FRDGTextureDesc HAvgDesc = FRDGTextureDesc::Create2D(
		FIntPoint(GridSizeX, GridSizeY),
		PF_R32_FLOAT,
		FClearValueBinding::None,
		TexCreate_ShaderResource | TexCreate_UAV
	);
	FRDGTextureRef HAvg_RDG = GraphBuilder.CreateTexture(HAvgDesc, TEXT("HAvg_RDG"));

	FTransportAndIntegrateOutputs TransportOutputs;
	TransportOutputs.htildeOut = htilde_RDG;
	TransportOutputs.qtildeOut_x = qtildex_RDG;
	TransportOutputs.qtildeOut_y = qtildey_RDG;
	TransportOutputs.qAdvectOut_x = qAdvect_x_RDG;
	TransportOutputs.qAdvectOut_y = qAdvect_y_RDG;
	TransportOutputs.hOut = h_RDG;
	TransportOutputs.qOut_x = qx_RDG;
	TransportOutputs.qOut_y = qy_RDG;
	TransportOutputs.hdot_out = hdot_RDG;
	TransportOutputs.H_Out = H_RDG;
	TransportOutputs.HAvg_Out = HAvg_RDG;

	AddTransportAndIntegratePasses(GraphBuilder, SimConstantBuffer, TransportInputs, TransportOutputs);

	// ----------------------------------------------------
	// EXPORT VISUAL OUTPUTS TO RENDER TARGETS
	// ----------------------------------------------------
	// Copy current Displacement to DisplacementPast before updating
	if (DisplacementPastRT && DisplacementPastRT->GetRenderTargetResource() &&
		DisplacementRT && DisplacementRT->GetRenderTargetResource())
	{
		FRDGTextureRef SrcRDG = GraphBuilder.RegisterExternalTexture(CreateRenderTarget(DisplacementRT->GetRenderTargetResource()->GetTexture2DRHI(), TEXT("DispCurrent_CopySrc")));
		FRDGTextureRef DestRDG = GraphBuilder.RegisterExternalTexture(CreateRenderTarget(DisplacementPastRT->GetRenderTargetResource()->GetTexture2DRHI(), TEXT("DispPast_CopyDest")));
		AddCopyTexturePass(GraphBuilder, SrcRDG, DestRDG);
	}

	// Copy current Velocity to VelocityPast before updating
	FRDGTextureRef VelPast_RDG = nullptr;
	if (VelocityPastRT && VelocityPastRT->GetRenderTargetResource() &&
		VelocityRT && VelocityRT->GetRenderTargetResource())
	{
		FRDGTextureRef SrcRDG = GraphBuilder.RegisterExternalTexture(CreateRenderTarget(VelocityRT->GetRenderTargetResource()->GetTexture2DRHI(), TEXT("VelCurrent_CopySrc")));
		FRDGTextureRef DestRDG = GraphBuilder.RegisterExternalTexture(CreateRenderTarget(VelocityPastRT->GetRenderTargetResource()->GetTexture2DRHI(), TEXT("VelPast_CopyDest")));
		AddCopyTexturePass(GraphBuilder, SrcRDG, DestRDG);
		VelPast_RDG = DestRDG;
	}
	else if (VelocityPastRT && VelocityPastRT->GetRenderTargetResource())
	{
		VelPast_RDG = GraphBuilder.RegisterExternalTexture(CreateRenderTarget(VelocityPastRT->GetRenderTargetResource()->GetTexture2DRHI(), TEXT("VelPastExport")));
	}

	// Copy current Acceleration to AccelerationPast before updating
	if (AccelerationPastRT && AccelerationPastRT->GetRenderTargetResource() &&
		AccelerationRT && AccelerationRT->GetRenderTargetResource())
	{
		FRDGTextureRef SrcRDG = GraphBuilder.RegisterExternalTexture(CreateRenderTarget(AccelerationRT->GetRenderTargetResource()->GetTexture2DRHI(), TEXT("AccelCurrent_CopySrc")));
		FRDGTextureRef DestRDG = GraphBuilder.RegisterExternalTexture(CreateRenderTarget(AccelerationPastRT->GetRenderTargetResource()->GetTexture2DRHI(), TEXT("AccelPast_CopyDest")));
		AddCopyTexturePass(GraphBuilder, SrcRDG, DestRDG);
	}

	FVisualExportInputs ExportInputs;
	ExportInputs.inDispX = disp_x_RDG;
	ExportInputs.inDispY = disp_y_RDG;
	ExportInputs.inHeight = HAvg_RDG; // Use time-averaged surface height for displacement export & normals
	ExportInputs.inPreviousFoam = GraphBuilder.RegisterExternalTexture(TexFoam);
	if (TexRoughness.IsValid()) ExportInputs.inPreviousRoughness = GraphBuilder.RegisterExternalTexture(TexRoughness);
	ExportInputs.inQx = qx_RDG;
	ExportInputs.inQy = qy_RDG;
	ExportInputs.inHPast = hPast_RDG;
	ExportInputs.inHNew = h_RDG;
	ExportInputs.inHdot = hdot_RDG;

	if (DisplacementRT && DisplacementRT->GetRenderTargetResource())
	{
		ExportInputs.ExportDispDest = GraphBuilder.RegisterExternalTexture(CreateRenderTarget(DisplacementRT->GetRenderTargetResource()->GetTexture2DRHI(), TEXT("DispExport")));
	}

	FRDGTextureRef VelocityExportRDG = nullptr;
	if (VelocityRT && VelocityRT->GetRenderTargetResource())
	{
		VelocityExportRDG = GraphBuilder.RegisterExternalTexture(CreateRenderTarget(VelocityRT->GetRenderTargetResource()->GetTexture2DRHI(), TEXT("VelocityExport")));
	}
	else if (VelocityStagingTextures[0].IsValid())
	{
		FRDGTextureDesc VelDesc = FRDGTextureDesc::Create2D(
			FIntPoint(GridSizeX, GridSizeY),
			PF_FloatRGBA,
			FClearValueBinding::None,
			TexCreate_ShaderResource | TexCreate_UAV
		);
		VelocityExportRDG = GraphBuilder.CreateTexture(VelDesc, TEXT("VelocityTransient"));
	}
	ExportInputs.ExportVelocityDest = VelocityExportRDG;
	ExportInputs.inVel = VelocityExportRDG;
	ExportInputs.inVelPast = VelPast_RDG ? VelPast_RDG : VelocityExportRDG;

	FRDGTextureRef AccelExportRDG = nullptr;
	if (AccelerationRT && AccelerationRT->GetRenderTargetResource())
	{
		AccelExportRDG = GraphBuilder.RegisterExternalTexture(CreateRenderTarget(AccelerationRT->GetRenderTargetResource()->GetTexture2DRHI(), TEXT("AccelerationExport")));
	}
	else if (AccelerationStagingTextures[0].IsValid())
	{
		FRDGTextureDesc AccelDesc = FRDGTextureDesc::Create2D(
			FIntPoint(GridSizeX, GridSizeY),
			PF_FloatRGBA,
			FClearValueBinding::None,
			TexCreate_ShaderResource | TexCreate_UAV
		);
		AccelExportRDG = GraphBuilder.CreateTexture(AccelDesc, TEXT("AccelerationTransient"));
	}
	ExportInputs.ExportAccelDest = AccelExportRDG;

	if (NormalRT && NormalRT->GetRenderTargetResource())
	{
		ExportInputs.ExportNormalDest = GraphBuilder.RegisterExternalTexture(CreateRenderTarget(NormalRT->GetRenderTargetResource()->GetTexture2DRHI(), TEXT("NormalExport")));
	}
	if (FoamRT && FoamRT->GetRenderTargetResource())
	{
		ExportInputs.ExportFoamDest = GraphBuilder.RegisterExternalTexture(CreateRenderTarget(FoamRT->GetRenderTargetResource()->GetTexture2DRHI(), TEXT("FoamExport")));
	}
	if (JacobianDetRT && JacobianDetRT->GetRenderTargetResource())
	{
		ExportInputs.ExportJacobianDest = GraphBuilder.RegisterExternalTexture(CreateRenderTarget(JacobianDetRT->GetRenderTargetResource()->GetTexture2DRHI(), TEXT("JacobianExport")));
	}
	if (RoughnessRT && RoughnessRT->GetRenderTargetResource())
	{
		ExportInputs.ExportRoughnessDest = GraphBuilder.RegisterExternalTexture(CreateRenderTarget(RoughnessRT->GetRenderTargetResource()->GetTexture2DRHI(), TEXT("RoughnessExport")));
	}
	ExportInputs.ScaleFactor = 100.0f;
	ExportInputs.IntegrationSamples = IntegrationSamples;
	ExportInputs.RoughnessPower = RoughnessPower;
	ExportInputs.GridSizeX = GridSizeX;
	ExportInputs.GridSizeY = GridSizeY;

	FVisualExportOutputs ExportOutputs;
	ExportOutputs.outNewFoam = GraphBuilder.RegisterExternalTexture(TexNewFoam);
	ExportOutputs.outNewRoughness = TexNewRoughness.IsValid() ? GraphBuilder.RegisterExternalTexture(TexNewRoughness) : nullptr;

	AddVisualExportPasses(GraphBuilder, ExportConstantBuffer, ExportInputs, ExportOutputs);

	// Ping-pong swap persistent foam & roughness textures
	Swap(TexFoam, TexNewFoam);
	if (TexRoughness.IsValid() && TexNewRoughness.IsValid())
	{
		Swap(TexRoughness, TexNewRoughness);
	}

	// ----------------------------------------------------
	// ASYNC READBACK OF DISPLACEMENT & VELOCITY
	// ----------------------------------------------------
	if (StagingTextures[StagingWriteIndex].IsValid() && HAvg_RDG)
	{
		FTextureRHIRef CurrentStagingTexture = StagingTextures[StagingWriteIndex];
		TRefCountPtr<IPooledRenderTarget> StagingPooled = CreateRenderTarget(CurrentStagingTexture, TEXT("StagingTexture"));
		FRDGTextureRef StagingRDG = GraphBuilder.RegisterExternalTexture(StagingPooled);
		AddCopyTexturePass(GraphBuilder, HAvg_RDG, StagingRDG);
	}

	if (VelocityStagingTextures[StagingWriteIndex].IsValid() && VelocityExportRDG)
	{
		FTextureRHIRef CurrentVelStaging = VelocityStagingTextures[StagingWriteIndex];
		TRefCountPtr<IPooledRenderTarget> VelStagingPooled = CreateRenderTarget(CurrentVelStaging, TEXT("VelStagingTexture"));
		FRDGTextureRef VelStagingRDG = GraphBuilder.RegisterExternalTexture(VelStagingPooled);
		AddCopyTexturePass(GraphBuilder, VelocityExportRDG, VelStagingRDG);
	}

	if (AccelerationStagingTextures[StagingWriteIndex].IsValid() && AccelExportRDG)
	{
		FTextureRHIRef CurrentAccelStaging = AccelerationStagingTextures[StagingWriteIndex];
		TRefCountPtr<IPooledRenderTarget> AccelStagingPooled = CreateRenderTarget(CurrentAccelStaging, TEXT("AccelStagingTexture"));
		FRDGTextureRef AccelStagingRDG = GraphBuilder.RegisterExternalTexture(AccelStagingPooled);
		AddCopyTexturePass(GraphBuilder, AccelExportRDG, AccelStagingRDG);
	}

	if (StagingTextures[StagingReadIndex].IsValid() || VelocityStagingTextures[StagingReadIndex].IsValid() || AccelerationStagingTextures[StagingReadIndex].IsValid())
	{
		FTextureRHIRef ReadStagingTexture = StagingTextures[StagingReadIndex];
		FTextureRHIRef ReadVelStagingTexture = VelocityStagingTextures[StagingReadIndex];
		FTextureRHIRef ReadAccelStagingTexture = AccelerationStagingTextures[StagingReadIndex];
		int32 TargetCPUIndex = 1 - ActiveCPUBufferIndex.load();

		GraphBuilder.AddPass(
			RDG_EVENT_NAME("SWE_ReadbackMap"),
			ERDGPassFlags::None,
			[ReadStagingTexture, ReadVelStagingTexture, ReadAccelStagingTexture, TargetCPUIndex, this](FRHICommandListImmediate& RHICmdList)
			{
				if (ReadStagingTexture.IsValid())
				{
					void* LocalData = nullptr;
					int32 OutWidth = 0;
					int32 OutHeight = 0;
					RHICmdList.MapStagingSurface(ReadStagingTexture, LocalData, OutWidth, OutHeight);
					if (LocalData)
					{
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
					}
				}

				if (ReadVelStagingTexture.IsValid())
				{
					void* VelData = nullptr;
					int32 VelOutWidth = 0;
					int32 VelOutHeight = 0;
					RHICmdList.MapStagingSurface(ReadVelStagingTexture, VelData, VelOutWidth, VelOutHeight);
					if (VelData)
					{
						int32 Size = GridSizeX * GridSizeY;
						TArray<FVector>& DestVelArray = CPUVelocityData[TargetCPUIndex];
						DestVelArray.SetNumUninitialized(Size);

						uint8* VelByteData = (uint8*)VelData;
						for (int32 y = 0; y < GridSizeY; ++y)
						{
							FFloat16Color* SrcRow = (FFloat16Color*)(VelByteData + y * VelOutWidth);
							for (int32 x = 0; x < GridSizeX; ++x)
							{
								DestVelArray[y * GridSizeX + x] = FVector(
									SrcRow[x].R.GetFloat(),
									SrcRow[x].G.GetFloat(),
									SrcRow[x].B.GetFloat()
								);
							}
						}
						RHICmdList.UnmapStagingSurface(ReadVelStagingTexture);
					}
				}

				if (ReadAccelStagingTexture.IsValid())
				{
					void* AccelData = nullptr;
					int32 AccelOutWidth = 0;
					int32 AccelOutHeight = 0;
					RHICmdList.MapStagingSurface(ReadAccelStagingTexture, AccelData, AccelOutWidth, AccelOutHeight);
					if (AccelData)
					{
						int32 Size = GridSizeX * GridSizeY;
						TArray<FVector>& DestAccelArray = CPUAccelerationData[TargetCPUIndex];
						DestAccelArray.SetNumUninitialized(Size);

						uint8* AccelByteData = (uint8*)AccelData;
						for (int32 y = 0; y < GridSizeY; ++y)
						{
							FFloat16Color* SrcRow = (FFloat16Color*)(AccelByteData + y * AccelOutWidth);
							for (int32 x = 0; x < GridSizeX; ++x)
							{
								DestAccelArray[y * GridSizeX + x] = FVector(
									SrcRow[x].R.GetFloat(),
									SrcRow[x].G.GetFloat(),
									SrcRow[x].B.GetFloat()
								);
							}
						}
						RHICmdList.UnmapStagingSurface(ReadAccelStagingTexture);
					}
				}

				// Swap active buffer index
				ActiveCPUBufferIndex.store(TargetCPUIndex);
			}
		);
	}

	// Swap indices for next frame
	StagingWriteIndex = (StagingWriteIndex + 1) % 2;
	StagingReadIndex = (StagingReadIndex + 1) % 2;

	GraphBuilder.Execute();
}

bool USimulator::LoadParametersFromJson(const FString& FilePath) 
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

bool USimulator::SaveParametersToJson(const FString& FilePath) 
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
	JsonObject->RemoveField(TEXT("velocityRT"));
	JsonObject->RemoveField(TEXT("velocityPastRT"));
	JsonObject->RemoveField(TEXT("accelerationRT"));
	JsonObject->RemoveField(TEXT("accelerationPastRT"));
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

float USimulator::GetCachedHeight(int32 X, int32 Y) const 
{
	int32 Index = Y * GridSizeX + X;
	int32 ActiveIdx = ActiveCPUBufferIndex.load();
	if (CPUHeightData[ActiveIdx].IsValidIndex(Index)) {
		return CPUHeightData[ActiveIdx][Index];
	}
	return WaterLevel * 0.01f; // base water level in meters
}

FVector USimulator::GetCachedVelocity(int32 X, int32 Y) const 
{
	int32 Index = Y * GridSizeX + X;
	int32 ActiveIdx = ActiveCPUBufferIndex.load();
	if (CPUVelocityData[ActiveIdx].IsValidIndex(Index)) {
		return CPUVelocityData[ActiveIdx][Index];
	}
	return FVector::ZeroVector;
}

FVector USimulator::GetCachedAcceleration(int32 X, int32 Y) const 
{
	int32 Index = Y * GridSizeX + X;
	int32 ActiveIdx = ActiveCPUBufferIndex.load();
	if (CPUAccelerationData[ActiveIdx].IsValidIndex(Index)) {
		return CPUAccelerationData[ActiveIdx][Index];
	}
	return FVector::ZeroVector;
}

float USimulator::GetWaterHeightAtLocation(const FVector& WorldLocation) const 
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

FVector USimulator::GetWaterVelocityAtLocation(const FVector& WorldLocation) const 
{
	if (GridSizeX <= 1 || GridSizeY <= 1 || CapturedWorldWidth <= 0.0f) {
		return FVector::ZeroVector;
	}

	AActor* Owner = GetOwner();
	if (!Owner) {
		return FVector::ZeroVector;
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

	FVector V00 = GetCachedVelocity(X0, Y0);
	FVector V10 = GetCachedVelocity(X1, Y0);
	FVector V01 = GetCachedVelocity(X0, Y1);
	FVector V11 = GetCachedVelocity(X1, Y1);

	FVector V_Avg = FMath::BiLerp(V00, V10, V01, V11, LerpX, LerpY);

	// In meters/sec
	return V_Avg;
}

FVector USimulator::GetWaterAccelerationAtLocation(const FVector& WorldLocation) const 
{
	if (GridSizeX <= 1 || GridSizeY <= 1 || CapturedWorldWidth <= 0.0f) {
		return FVector::ZeroVector;
	}

	AActor* Owner = GetOwner();
	if (!Owner) {
		return FVector::ZeroVector;
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

	FVector A00 = GetCachedAcceleration(X0, Y0);
	FVector A10 = GetCachedAcceleration(X1, Y0);
	FVector A01 = GetCachedAcceleration(X0, Y1);
	FVector A11 = GetCachedAcceleration(X1, Y1);

	FVector A_Avg = FMath::BiLerp(A00, A10, A01, A11, LerpX, LerpY);

	// In m/s^2
	return A_Avg;
}
