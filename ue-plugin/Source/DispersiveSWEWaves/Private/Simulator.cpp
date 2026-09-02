#include "Simulator.h"
#include "SimConfigSerializer.h"
#include "Misc/Paths.h"
#include "HAL/FileManager.h"
#include "RenderGraphBuilder.h"
#include "RenderGraphUtils.h"
#include "RHI.h"
#include "RHIResources.h"
#include "RHIGPUReadback.h"
#include "RHIStaticStates.h"

USimulator::USimulator() {
	PrimaryComponentTick.bCanEverTick = true;
	PrimaryComponentTick.bStartWithTickEnabled = true;
}

void USimulator::BeginPlay() {
	UE_LOG(LogTemp, Warning, TEXT("USimulator::BeginPlay() called"));
	Super::BeginPlay();
	InitializeSimulation();
}

void USimulator::InitializeSimulation() {
	UE_LOG(LogTemp, Warning, TEXT("InitializeSimulation() running: GridSizeX=%d, GridSizeY=%d, CellSize=%f, CapturedWorldWidth=%f"), GridSizeX, GridSizeY, CellSize, CapturedWorldWidth);

	if (!JsonConfigFilePath.IsEmpty()) {
		LoadParametersFromJson(JsonConfigFilePath);
	}

	if (bAutoCalculateCellSize && GridSizeX > 0) {
		CellSize = CapturedWorldWidth / (float)GridSizeX;
	} 

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

	PaddedSizeX = FMath::RoundUpToPowerOfTwo(GridSizeX);
	PaddedSizeY = FMath::RoundUpToPowerOfTwo(GridSizeY);
	SimulationTime = 0.0f;

	ENQUEUE_RENDER_COMMAND(InitializeSWESimulation)(
		[this](FRHICommandListImmediate& RHICmdList) {
			StateBuffers.AllocateAll(RHICmdList, GridSizeX, GridSizeY, PaddedSizeX, PaddedSizeY, DepthLevels.Num());
			ReadbackHandler.Allocate(GridSizeX, GridSizeY);
			SetupInitialStates(RHICmdList);
			bInitialized = true;
		}
	);
}

void USimulator::AssignSimulationConstants(FSimConstants& Constants) const {
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

void USimulator::AssignFFTWaveConstants(FFFTWaveConstants& Constants) const {
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

void USimulator::AssignExportConstants(FExportConstants& Constants) const {
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

static FRDGTextureRef RegisterRenderTargetSafe(FRDGBuilder& GraphBuilder, UTextureRenderTarget2D* RT, const TCHAR* Name) {
	if (!RT) return nullptr;
	FTextureRenderTargetResource* Resource = RT->GetRenderTargetResource();
	if (!Resource) return nullptr;
	FTexture2DRHIRef TextureRHI = Resource->GetTexture2DRHI();
	if (!TextureRHI.IsValid()) return nullptr;
	TRefCountPtr<IPooledRenderTarget> Pooled = CreateRenderTarget(TextureRHI, Name);
	if (!Pooled.IsValid()) return nullptr;
	return GraphBuilder.RegisterExternalTexture(Pooled);
}

void USimulator::SetupInitialStates(FRHICommandListImmediate& RHICmdList) {
	FSimConstants CPUSimConstants = {};
	AssignSimulationConstants(CPUSimConstants);
	TUniformBufferRef<FSimConstants> SimConstantBuffer = CreateUniformBufferImmediate(CPUSimConstants, EUniformBufferUsage::UniformBuffer_SingleFrame);

	FFFTWaveConstants CPUFFTWaveConstants = {};
	AssignFFTWaveConstants(CPUFFTWaveConstants);
	TUniformBufferRef<FFFTWaveConstants> FFTWaveConstantBuffer = CreateUniformBufferImmediate(CPUFFTWaveConstants, EUniformBufferUsage::UniformBuffer_SingleFrame);

	// Clear all persistent simulation state textures to zero
	FRDGBuilder ClearGraphBuilder(RHICmdList);
	StateBuffers.ClearAll(ClearGraphBuilder);
	ClearGraphBuilder.Execute();

	// Initialize Terrain and Water Height fields
	FRDGBuilder GraphBuilder(RHICmdList);
	FRDGTextureRef TerrainInput_RDG = RegisterRenderTargetSafe(GraphBuilder, TerrainHeightInputRT, TEXT("TerrainInput"));

	if (TerrainInput_RDG) {
		FRDGTextureRef TerrainExportRDG = RegisterRenderTargetSafe(GraphBuilder, TerrainRT, TEXT("TerrainExport"));
		if (!TerrainExportRDG) {
			FRDGTextureDesc DummyDesc = FRDGTextureDesc::Create2D(
				FIntPoint(GridSizeX, GridSizeY),
				PF_R32_FLOAT,
				FClearValueBinding::None,
				TexCreate_ShaderResource | TexCreate_UAV
			);
			TerrainExportRDG = GraphBuilder.CreateTexture(DummyDesc, TEXT("TerrainExportDummy"));
		}

		FInitializeWaterInputs InitInputs;
		InitInputs.TerrainInput = TerrainInput_RDG;
		InitInputs.WaterLevel = WaterLevel * 0.01f; // Convert cm to meters
		InitInputs.TerrainCaptureCameraZ = TerrainCaptureCameraZ;
		InitInputs.GridSizeX = GridSizeX;
		InitInputs.GridSizeY = GridSizeY;

		FInitializeWaterOutputs InitOutputs;
		InitOutputs.hOut = StateBuffers.h.RegisterCurrent(GraphBuilder);
		InitOutputs.H_Out = StateBuffers.H.RegisterCurrent(GraphBuilder);
		InitOutputs.terrainOut = GraphBuilder.RegisterExternalTexture(StateBuffers.Terrain);
		InitOutputs.terrainOutCM = TerrainExportRDG;
		InitOutputs.hbarOut = StateBuffers.hbar.RegisterCurrent(GraphBuilder);
		InitOutputs.hbarOldOut = StateBuffers.hbar.RegisterTarget(GraphBuilder);
		if (StateBuffers.Foam.Current.IsValid()) InitOutputs.FoamOut = StateBuffers.Foam.RegisterCurrent(GraphBuilder);
		if (StateBuffers.Roughness.Current.IsValid()) InitOutputs.RoughnessOut = StateBuffers.Roughness.RegisterCurrent(GraphBuilder);

		AddInitializeWaterPass(GraphBuilder, SimConstantBuffer, InitInputs, InitOutputs);
		GraphBuilder.Execute();
	} else {
		TArray<float> TerrainData;
		TerrainData.SetNumZeroed(GridSizeX * GridSizeY);
		TArray<float> WaterData;
		WaterData.SetNumZeroed(GridSizeX * GridSizeY);
		TArray<float> HData;
		HData.SetNumZeroed(GridSizeX * GridSizeY);

		float TerrainHeight = -13.0f;
		float TerrainScale = 20.0f;
		float WaterLevelLocal = WaterLevel * 0.01f;

		for (int32 y = 0; y < GridSizeY; y++) {
			for (int32 x = 0; x < GridSizeX; x++) {
				float xf = (float)x / (GridSizeX - 1);
				float yf = (float)y / (GridSizeY - 1);
				int32 i = y * GridSizeX + x;

				float dunes = 0.05f * FMath::Sin(20.f * yf);
				TerrainData[i] = TerrainHeight + TerrainScale * (xf * (1.0f + dunes));
				WaterData[i] = FMath::Max(0.0f, WaterLevelLocal - TerrainData[i]);
				HData[i] = TerrainData[i] + WaterData[i];
			}
		}

		FUpdateTextureRegion2D Region(0, 0, 0, 0, GridSizeX, GridSizeY);
		RHIUpdateTexture2D(static_cast<FRHITexture2D*>(StateBuffers.Terrain->GetRHI()), 0, Region, GridSizeX * sizeof(float), (uint8*)TerrainData.GetData());
		RHIUpdateTexture2D(static_cast<FRHITexture2D*>(StateBuffers.H.Current->GetRHI()), 0, Region, GridSizeX * sizeof(float), (uint8*)HData.GetData());
		RHIUpdateTexture2D(static_cast<FRHITexture2D*>(StateBuffers.h.Current->GetRHI()), 0, Region, GridSizeX * sizeof(float), (uint8*)WaterData.GetData());
		RHIUpdateTexture2D(static_cast<FRHITexture2D*>(StateBuffers.hbar.Current->GetRHI()), 0, Region, GridSizeX * sizeof(float), (uint8*)WaterData.GetData());
		RHIUpdateTexture2D(static_cast<FRHITexture2D*>(StateBuffers.hbar.Target->GetRHI()), 0, Region, GridSizeX * sizeof(float), (uint8*)WaterData.GetData());
	}

	// Initialize Wave Spectrum
	FRDGBuilder SpectrumGraphBuilder(RHICmdList);
	FRDGTextureRef HPos_RDG = SpectrumGraphBuilder.RegisterExternalTexture(StateBuffers.HPos);
	FRDGTextureRef HNeg_RDG = SpectrumGraphBuilder.RegisterExternalTexture(StateBuffers.HNeg);

	AddPopulateSpectrumPass(SpectrumGraphBuilder, FFTWaveConstantBuffer, HPos_RDG, HNeg_RDG, PaddedSizeX, PaddedSizeY, DepthLevels.Num());
	SpectrumGraphBuilder.Execute();

	// Export captured terrain to raw file if requested
	if (bExportCapturedTerrain && StateBuffers.Terrain.IsValid()) {
		FRHITexture2D* TerrainRHI = static_cast<FRHITexture2D*>(StateBuffers.Terrain->GetRHI());
		FRHIGPUTextureReadback* GPUReadback = new FRHIGPUTextureReadback(TEXT("TerrainExportReadback"));
		GPUReadback->EnqueueCopy(RHICmdList, TerrainRHI);
		RHICmdList.BlockUntilGPUIdle();

		if (GPUReadback->IsReady()) {
			int32 OutWidth = 0, OutHeight = 0;
			float* FloatData = static_cast<float*>(GPUReadback->Lock(OutWidth, &OutHeight));
			if (FloatData) {
				TArray<float> RawFloatData;
				RawFloatData.SetNumUninitialized(GridSizeX * GridSizeY);
				for (int32 y = 0; y < GridSizeY; ++y) {
					FMemory::Memcpy(&RawFloatData[y * GridSizeX], &FloatData[y * OutWidth], GridSizeX * sizeof(float));
				}
				GPUReadback->Unlock();

				FString ExportFilePath = FPaths::ProjectDir() / TEXT("terrain_captured.raw");
				if (FArchive* Ar = IFileManager::Get().CreateFileWriter(*ExportFilePath)) {
					Ar->Serialize(RawFloatData.GetData(), RawFloatData.Num() * sizeof(float));
					delete Ar;
					UE_LOG(LogTemp, Warning, TEXT("DispersiveSWEWaves: Successfully exported captured terrain to %s"), *ExportFilePath);
				}
			}
		}
		delete GPUReadback;
	}
}

void USimulator::TickComponent(float DeltaTime, ELevelTick TickType, FActorComponentTickFunction* ThisTickFunction) {
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
	const FExportConstants& ExportConstants) {
	FRDGBuilder GraphBuilder(RHICmdList);

	FRDGTextureDesc FloatGridDesc = FRDGTextureDesc::Create2D(
		FIntPoint(GridSizeX, GridSizeY),
		PF_R32_FLOAT,
		FClearValueBinding::None,
		TexCreate_ShaderResource | TexCreate_UAV
	);

	FRDGTextureDesc ComplexArrayDesc = FRDGTextureDesc::Create2DArray(
		FIntPoint(PaddedSizeX, PaddedSizeY),
		PF_G32R32F,
		FClearValueBinding::None,
		TexCreate_ShaderResource | TexCreate_UAV,
		SimConstants.depthNum
	);

	// Import persistent state buffers
	FRDGTextureRef Terrain_RDG = GraphBuilder.RegisterExternalTexture(StateBuffers.Terrain);
	FRDGTextureRef H_RDG = StateBuffers.H.RegisterCurrent(GraphBuilder);
	FRDGTextureRef HPast_RDG = StateBuffers.H.RegisterTarget(GraphBuilder);
	FRDGTextureRef Qx_RDG = StateBuffers.Q_x.RegisterCurrent(GraphBuilder);
	FRDGTextureRef QPast_x_RDG = StateBuffers.Q_x.RegisterTarget(GraphBuilder);
	FRDGTextureRef Qy_RDG = StateBuffers.Q_y.RegisterCurrent(GraphBuilder);
	FRDGTextureRef QPast_y_RDG = StateBuffers.Q_y.RegisterTarget(GraphBuilder);
	FRDGTextureRef h_RDG = StateBuffers.h.RegisterCurrent(GraphBuilder);
	FRDGTextureRef hPast_RDG = StateBuffers.h.RegisterTarget(GraphBuilder);
	FRDGTextureRef qx_RDG = GraphBuilder.RegisterExternalTexture(StateBuffers.q_x);
	FRDGTextureRef qy_RDG = GraphBuilder.RegisterExternalTexture(StateBuffers.q_y);
	FRDGTextureRef hbar_RDG = StateBuffers.hbar.RegisterCurrent(GraphBuilder);
	FRDGTextureRef hbarOld_RDG = StateBuffers.hbar.RegisterTarget(GraphBuilder);
	FRDGTextureRef htildeOld_RDG = StateBuffers.htildeOld.RegisterCurrent(GraphBuilder);
	FRDGTextureRef htildeOldNext_RDG = StateBuffers.htildeOld.RegisterTarget(GraphBuilder);
	FRDGTextureRef HPos_RDG = GraphBuilder.RegisterExternalTexture(StateBuffers.HPos);
	FRDGTextureRef HNeg_RDG = GraphBuilder.RegisterExternalTexture(StateBuffers.HNeg);

	// Create uniform buffers
	TUniformBufferRef<FSimConstants> SimConstantBuffer = CreateUniformBufferImmediate(SimConstants, EUniformBufferUsage::UniformBuffer_SingleFrame);
	TUniformBufferRef<FFFTWaveConstants> FFTWaveConstantBuffer = CreateUniformBufferImmediate(FFTWaveConstants, EUniformBufferUsage::UniformBuffer_SingleFrame);
	TUniformBufferRef<FExportConstants> ExportConstantBuffer = CreateUniformBufferImmediate(ExportConstants, EUniformBufferUsage::UniformBuffer_SingleFrame);

	///// RUN SIMULATION STEPS /////

	// DECOMPOSITION

	FRDGTextureRef qbarx_RDG = GraphBuilder.CreateTexture(FloatGridDesc, TEXT("qbar_x"));
	FRDGTextureRef qbary_RDG = GraphBuilder.CreateTexture(FloatGridDesc, TEXT("qbar_y"));
	FRDGTextureRef htilde_RDG = GraphBuilder.CreateTexture(FloatGridDesc, TEXT("htilde"));
	FRDGTextureRef qtildex_RDG = GraphBuilder.CreateTexture(FloatGridDesc, TEXT("qtilde_x"));
	FRDGTextureRef qtildey_RDG = GraphBuilder.CreateTexture(FloatGridDesc, TEXT("qtilde_y"));

	FDecompositionInputs DecompInputs;
	DecompInputs.hIn = h_RDG;
	DecompInputs.qIn_x = qx_RDG;
	DecompInputs.qIn_y = qy_RDG;
	DecompInputs.terrain = Terrain_RDG;
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

	if (DecompOutputs.H_SrcDst != H_RDG) StateBuffers.H.Swap();
	if (DecompOutputs.Qx_SrcDst != Qx_RDG) StateBuffers.Q_x.Swap();
	if (DecompOutputs.Qy_SrcDst != Qy_RDG) StateBuffers.Q_y.Swap();
	H_RDG = DecompOutputs.H_SrcDst;
	Qx_RDG = DecompOutputs.Qx_SrcDst;
	Qy_RDG = DecompOutputs.Qy_SrcDst;


	// FFT WIND WAVES

	FRDGTextureRef Flow_x_RDG = GraphBuilder.CreateTexture(ComplexArrayDesc, TEXT("Flow_x"));
	FRDGTextureRef Flow_y_RDG = GraphBuilder.CreateTexture(ComplexArrayDesc, TEXT("Flow_y"));
	FRDGTextureRef disp_x_RDG = GraphBuilder.CreateTexture(FloatGridDesc, TEXT("disp_x"));
	FRDGTextureRef disp_y_RDG = GraphBuilder.CreateTexture(FloatGridDesc, TEXT("disp_y"));
	FRDGTextureRef delH_x_RDG = GraphBuilder.CreateTexture(FloatGridDesc, TEXT("delH_x"));
	FRDGTextureRef delH_y_RDG = GraphBuilder.CreateTexture(FloatGridDesc, TEXT("delH_y"));

	FPropagateFFTWavesInputs FFTWaveInputs;
	FFTWaveInputs.HPosIn = HPos_RDG;
	FFTWaveInputs.HNegIn = HNeg_RDG;
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
	FFTWaveOutputs.Flow_x_Out = Flow_x_RDG;
	FFTWaveOutputs.Flow_y_Out = Flow_y_RDG;

	AddPropagateFFTWavesPasses(GraphBuilder, FFTWaveConstantBuffer, FFTWaveInputs, FFTWaveOutputs);


	// EWAVE DISPERSION

	FEWaveInputs EWaveInputs;
	EWaveInputs.htildeIn = htilde_RDG;
	EWaveInputs.htildeOldIn = htildeOld_RDG;
	EWaveInputs.qtildeIn_x = qtildex_RDG;
	EWaveInputs.qtildeIn_y = qtildey_RDG;
	EWaveInputs.Flow_x = Flow_x_RDG;
	EWaveInputs.Flow_y = Flow_y_RDG;
	EWaveInputs.hbarIn = hbar_RDG;
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
	StateBuffers.htildeOld.Swap();
	Swap(htildeOld_RDG, htildeOldNext_RDG);


	// BULK FLOW (SWE VELOCITY & MOMENTUM)

	FRDGTextureRef ubarx_RDG = GraphBuilder.CreateTexture(FloatGridDesc, TEXT("ubar_x"));
	FRDGTextureRef ubary_RDG = GraphBuilder.CreateTexture(FloatGridDesc, TEXT("ubar_y"));
	FRDGTextureRef ubarNew_x_RDG = GraphBuilder.CreateTexture(FloatGridDesc, TEXT("ubarNew_x"));
	FRDGTextureRef ubarNew_y_RDG = GraphBuilder.CreateTexture(FloatGridDesc, TEXT("ubarNew_y"));

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
	StateBuffers.hbar.Swap();
	Swap(hbar_RDG, hbarOld_RDG);


	// ADVECTIVE TRANSPORT & HEIGHT INTEGRATION

	FRDGTextureRef htildeAdvect_RDG = GraphBuilder.CreateTexture(FloatGridDesc, TEXT("htildeAdvect"));
	FRDGTextureRef qtildeAdvect_x_RDG = GraphBuilder.CreateTexture(FloatGridDesc, TEXT("qtildeAdvect_x"));
	FRDGTextureRef qtildeAdvect_y_RDG = GraphBuilder.CreateTexture(FloatGridDesc, TEXT("qtildeAdvect_y"));
	FRDGTextureRef hdot_RDG = GraphBuilder.CreateTexture(FloatGridDesc, TEXT("hdot"));
	FRDGTextureRef HAvg_RDG = GraphBuilder.CreateTexture(FloatGridDesc, TEXT("HAvg"));

	FTransportAndIntegrateInputs TransportInputs;
	TransportInputs.ubarNewIn_x = ubarNew_x_RDG;
	TransportInputs.ubarIn_x = ubarx_RDG;
	TransportInputs.ubarNewIn_y = ubarNew_y_RDG;
	TransportInputs.ubarIn_y = ubary_RDG;
	TransportInputs.qtildePast_x = qtildex_RDG;
	TransportInputs.qtildePast_y = qtildey_RDG;
	TransportInputs.htildePast = htilde_RDG;
	TransportInputs.hPast = h_RDG;
	TransportInputs.terrain = Terrain_RDG;
	TransportInputs.qbarIn_x = qbarx_RDG;
	TransportInputs.qbarIn_y = qbary_RDG;
	TransportInputs.GridSizeX = GridSizeX;
	TransportInputs.GridSizeY = GridSizeY;

	FTransportAndIntegrateOutputs TransportOutputs;
	TransportOutputs.htildeOut = htildeAdvect_RDG;
	TransportOutputs.qtildeOut_x = qtildeAdvect_x_RDG;
	TransportOutputs.qtildeOut_y = qtildeAdvect_y_RDG;
	TransportOutputs.hOut = hPast_RDG;
	TransportOutputs.qOut_x = qx_RDG;
	TransportOutputs.qOut_y = qy_RDG;
	TransportOutputs.hdot_out = hdot_RDG;
	TransportOutputs.H_Out = H_RDG;
	TransportOutputs.HAvg_Out = HAvg_RDG;

	AddTransportAndIntegratePasses(GraphBuilder, SimConstantBuffer, TransportInputs, TransportOutputs);
	StateBuffers.h.Swap();
	Swap(h_RDG, hPast_RDG);


	// EXPORT VALUES

	FRDGTextureRef DispCurrentRDG = RegisterRenderTargetSafe(GraphBuilder, DisplacementRT, TEXT("DispCurrent_CopySrc"));
	FRDGTextureRef DispPastRDG = RegisterRenderTargetSafe(GraphBuilder, DisplacementPastRT, TEXT("DispPast_CopyDest"));
	if (DispCurrentRDG && DispPastRDG) {
		AddCopyTexturePass(GraphBuilder, DispCurrentRDG, DispPastRDG);
	}

	FRDGTextureRef VelCurrentRDG = RegisterRenderTargetSafe(GraphBuilder, VelocityRT, TEXT("VelCurrent_CopySrc"));
	FRDGTextureRef VelPast_RDG = RegisterRenderTargetSafe(GraphBuilder, VelocityPastRT, TEXT("VelPast_CopyDest"));
	if (VelCurrentRDG && VelPast_RDG) {
		AddCopyTexturePass(GraphBuilder, VelCurrentRDG, VelPast_RDG);
	}

	FRDGTextureRef AccelCurrentRDG = RegisterRenderTargetSafe(GraphBuilder, AccelerationRT, TEXT("AccelCurrent_CopySrc"));
	FRDGTextureRef AccelPastRDG = RegisterRenderTargetSafe(GraphBuilder, AccelerationPastRT, TEXT("AccelPast_CopyDest"));
	if (AccelCurrentRDG && AccelPastRDG) {
		AddCopyTexturePass(GraphBuilder, AccelCurrentRDG, AccelPastRDG);
	}

	FVisualExportInputs ExportInputs;
	ExportInputs.inDispX = disp_x_RDG;
	ExportInputs.inDispY = disp_y_RDG;
	ExportInputs.inHeight = HAvg_RDG;
	ExportInputs.inPreviousFoam = StateBuffers.Foam.RegisterCurrent(GraphBuilder);
	if (StateBuffers.Roughness.Current.IsValid()) ExportInputs.inPreviousRoughness = StateBuffers.Roughness.RegisterCurrent(GraphBuilder);
	ExportInputs.inQx = qx_RDG;
	ExportInputs.inQy = qy_RDG;
	ExportInputs.inHPast = hPast_RDG;
	ExportInputs.inHNew = h_RDG;
	ExportInputs.inHdot = hdot_RDG;
	ExportInputs.ScaleFactor = 100.0f;
	ExportInputs.IntegrationSamples = IntegrationSamples;
	ExportInputs.RoughnessPower = RoughnessPower;
	ExportInputs.GridSizeX = GridSizeX;
	ExportInputs.GridSizeY = GridSizeY;

	FVisualExportOutputs ExportOutputs;
	ExportOutputs.ExportDispDest = RegisterRenderTargetSafe(GraphBuilder, DisplacementRT, TEXT("DispExport"));

	FRDGTextureRef VelocityExportRDG = RegisterRenderTargetSafe(GraphBuilder, VelocityRT, TEXT("VelocityExport"));
	if (!VelocityExportRDG && ReadbackHandler.HasVelocityStaging()) {
		FRDGTextureDesc VelDesc = FRDGTextureDesc::Create2D(
			FIntPoint(GridSizeX, GridSizeY),
			PF_FloatRGBA,
			FClearValueBinding::None,
			TexCreate_ShaderResource | TexCreate_UAV
		);
		VelocityExportRDG = GraphBuilder.CreateTexture(VelDesc, TEXT("VelocityTransient"));
	}
	ExportOutputs.ExportVelocityDest = VelocityExportRDG;
	ExportInputs.inVel = VelocityExportRDG;
	ExportInputs.inVelPast = VelPast_RDG ? VelPast_RDG : VelocityExportRDG;

	FRDGTextureRef AccelExportRDG = RegisterRenderTargetSafe(GraphBuilder, AccelerationRT, TEXT("AccelerationExport"));
	if (!AccelExportRDG && ReadbackHandler.HasAccelerationStaging()) {
		FRDGTextureDesc AccelDesc = FRDGTextureDesc::Create2D(
			FIntPoint(GridSizeX, GridSizeY),
			PF_FloatRGBA,
			FClearValueBinding::None,
			TexCreate_ShaderResource | TexCreate_UAV
		);
		AccelExportRDG = GraphBuilder.CreateTexture(AccelDesc, TEXT("AccelerationTransient"));
	}
	ExportOutputs.ExportAccelDest = AccelExportRDG;

	ExportOutputs.ExportNormalDest = RegisterRenderTargetSafe(GraphBuilder, NormalRT, TEXT("NormalExport"));
	ExportOutputs.ExportFoamDest = RegisterRenderTargetSafe(GraphBuilder, FoamRT, TEXT("FoamExport"));
	ExportOutputs.ExportJacobianDest = RegisterRenderTargetSafe(GraphBuilder, JacobianDetRT, TEXT("JacobianExport"));
	ExportOutputs.ExportRoughnessDest = RegisterRenderTargetSafe(GraphBuilder, RoughnessRT, TEXT("RoughnessExport"));

	ExportOutputs.outNewFoam = StateBuffers.Foam.RegisterTarget(GraphBuilder);
	ExportOutputs.outNewRoughness = StateBuffers.Roughness.Target.IsValid() ? StateBuffers.Roughness.RegisterTarget(GraphBuilder) : nullptr;

	AddVisualExportPasses(GraphBuilder, ExportConstantBuffer, ExportInputs, ExportOutputs);


	// ASYNC GPU READBACK
	ReadbackHandler.EnqueueReadback(GraphBuilder, HAvg_RDG, VelocityExportRDG, AccelExportRDG, GridSizeX, GridSizeY);
	GraphBuilder.Execute();
}

bool USimulator::LoadParametersFromJson(const FString& FilePath) {
	return FSimConfigSerializer::LoadParametersFromJson(this, FilePath, bAutoCalculateCellSize);
}

bool USimulator::SaveParametersToJson(const FString& FilePath) {
	return FSimConfigSerializer::SaveParametersToJson(this, FilePath);
}

float USimulator::GetWaterHeightAtLocation(const FVector& WorldLocation) const {
	AActor* Owner = GetOwner();
	FVector ActorLoc = Owner ? Owner->GetActorLocation() : FVector::ZeroVector;
	return ReadbackHandler.GetWaterHeightAtLocation(WorldLocation, ActorLoc, CapturedWorldWidth, GridSizeX, GridSizeY, WaterLevel);
}

FVector USimulator::GetWaterVelocityAtLocation(const FVector& WorldLocation) const {
	AActor* Owner = GetOwner();
	FVector ActorLoc = Owner ? Owner->GetActorLocation() : FVector::ZeroVector;
	return ReadbackHandler.GetWaterVelocityAtLocation(WorldLocation, ActorLoc, CapturedWorldWidth, GridSizeX, GridSizeY);
}

FVector USimulator::GetWaterAccelerationAtLocation(const FVector& WorldLocation) const {
	AActor* Owner = GetOwner();
	FVector ActorLoc = Owner ? Owner->GetActorLocation() : FVector::ZeroVector;
	return ReadbackHandler.GetWaterAccelerationAtLocation(WorldLocation, ActorLoc, CapturedWorldWidth, GridSizeX, GridSizeY);
}
