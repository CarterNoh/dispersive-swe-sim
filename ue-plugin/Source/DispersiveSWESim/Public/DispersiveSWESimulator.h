#pragma once

#include "CoreMinimal.h"
#include "Components/ActorComponent.h"
#include "Engine/TextureRenderTarget2D.h"
#include "RenderGraphResources.h"
#include "DispersiveSWEShaders.h"
#include <atomic>
#include "DispersiveSWESimulator.generated.h"

UCLASS(ClassGroup=(Custom), meta=(BlueprintSpawnableComponent))
class DISPERSIVESWESIM_API UDispersiveSWESimulator : public UActorComponent
{
	GENERATED_BODY()

public:
	UDispersiveSWESimulator();
	virtual void BeginPlay() override;
	virtual void TickComponent(float DeltaTime, ELevelTick TickType, FActorComponentTickFunction* ThisTickFunction) override;


	//////// Parameters ////////

	// Grid Resolution in X dimension
	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "SWE Simulation")
	int32 GridSizeX = 512;

	// Grid Resolution in Y dimension
	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "SWE Simulation")
	int32 GridSizeY = 512;

	// Physical cell size in centimeters per cell
	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "SWE Simulation")
	float CellSize = 100.0f;

	// Time step of the simulation in seconds
	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "SWE Simulation")
	float TimeStep = 0.016666f;

	// CFL Condition boundary check
	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "SWE Simulation")
	float CFLCondition = 0.25f;

	// Sponge Layer Thickness (for wave absorption at edges)
	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "SWE Simulation")
	int32 SpongeThickness = 8;

	// Minimum water height in centimeters for terrain boundary
	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "SWE Simulation")
	float MinWaterHeight = 0.1f;

	// Maximum depth in centimeters allowed for safe, stable wave simulation. If <= 0, automatically calculated.
	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "SWE Simulation")
	float MaxSafeDepth = 0.0f;
	float CalculatedMaxSafeDepth;

	// Safety multiplier applied to the calculated max safe depth curve to prevent numerical explosion.
	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "SWE Simulation", meta = (ToolTip = "Safety multiplier applied to the calculated max safe depth curve to prevent numerical explosion near the stability boundary. Defaults to 0.8."))
	float StabilitySafetyFactor = 0.8f;

	// Initial water free surface level in centimeters
	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "SWE Simulation")
	float WaterLevel = 0.0f;

	// Surface tension coefficient
	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "SWE Simulation")
	float SurfaceTension = 0.001f;

	// Fluid density in kg/m3
	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "SWE Simulation")
	float Density = 999.0f;

	// Iterations of the diffusion step per frame
	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "SWE Simulation|Decomposition")
	int32 DiffusionIterations = 128;

	// Diffusion delta timestep
	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "SWE Simulation|Decomposition")
	float DiffusionDeltaT = 0.25f;

	// Penalty parameter damping gradients
	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "SWE Simulation|Decomposition")
	float DiffusionPenalty = 0.01f;

	// Slope limit of crashing waves
	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "SWE Simulation|SWE")
	float SlopeLimit = 1.00f;

	// Advection damping coefficient
	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "SWE Simulation|SWE")
	float GammaTransport = 0.25f;

	// Fetch in kilometers for JONSWAP wave spectrum
	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "SWE Simulation|Wind Wave")
	float Fetch = 200.0f;

	// Wind speed in m/s
	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "SWE Simulation|Wind Wave")
	float WindSpeed = 14.0f;

	// Wind angle in degrees
	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "SWE Simulation|Wind Wave")
	float WindAngle = 135.0f;

	// Swell factor
	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "SWE Simulation|Wind Wave")
	float Swell = 0.3f;

	// Swell direction angle
	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "SWE Simulation|Wind Wave")
	float SwellAngle = 135.0f;

	// Wave choppiness (horizontal displacement scaling)
	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "SWE Simulation|Wind Wave")
	float Choppiness = 1.0f;

	// Depth attenuation cutoff
	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "SWE Simulation|Wind Wave")
	float DepthCutoff = 8.0f;

	// Band pass filter parameters
	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "SWE Simulation|Wind Wave|Filter")
	float FilterSmall = 0.0f;

	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "SWE Simulation|Wind Wave|Filter")
	float FilterBig = 10000.0f;

	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "SWE Simulation|Wind Wave|Filter")
	float FilterWidth = 1.0f;

	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "SWE Simulation|Wind Wave|Filter")
	float FilterMin = 0.01f;

	// Discrete water depths for wave dispersion calculations
	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "SWE Simulation|eWave")
	TArray<float> DepthLevels = { 1.0f, 2.0f, 4.0f, 16.0f, 64.0f };

	// Foam threshold, multiplier, decay fade and blur rates
	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "SWE Simulation|eWave|Foam")
	float FoamThreshold = -0.25f;

	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "SWE Simulation|eWave|Foam")
	float FoamMultiplier = 1.0f;

	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "SWE Simulation|eWave|Foam")
	float FoamFade = 0.1f;

	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "SWE Simulation|eWave|Foam")
	float FoamBlur = 2.0f;

	// Roughness scale integration samples and power parameter
	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "SWE Simulation|eWave|Roughness")
	float IntegrationSamples = 100.0f;

	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "SWE Simulation|eWave|Roughness")
	float RoughnessPower = 1.0f;
	

	//////// Render Targets ////////

	// Input Render Target containing the level's terrain/height capture map
	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "SWE Inputs")
	UTextureRenderTarget2D* TerrainHeightInputRT = nullptr;

	// The absolute Z height of the terrain capture camera (used to reconstruct world Z from depth)
	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "SWE Inputs")
	float TerrainCaptureCameraZ = 5000.0f;

	// Auto calculate CellSize based on CapturedWorldWidth and GridSizeX
	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "SWE Inputs")
	bool bAutoCalculateCellSize = true;

	// The width of the captured terrain region in the world, in centimeters
	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "SWE Inputs", meta = (EditCondition = "bAutoCalculateCellSize"))
	float CapturedWorldWidth = 51200.0f;

	// File path to a JSON configuration file to load parameters from on startup
	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "SWE Configuration")
	FString JsonConfigFilePath = "";

	// Output target textures containing wave fields
	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "SWE Outputs")
	UTextureRenderTarget2D* TerrainRT = nullptr;

	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "SWE Outputs")
	UTextureRenderTarget2D* DisplacementRT = nullptr;

	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "SWE Outputs")
	UTextureRenderTarget2D* DisplacementPastRT = nullptr;

	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "SWE Outputs")
	UTextureRenderTarget2D* NormalRT = nullptr;

	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "SWE Outputs")
	UTextureRenderTarget2D* FoamRT = nullptr;

	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "SWE Outputs")
	UTextureRenderTarget2D* JacobianDetRT = nullptr;

	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "SWE Outputs")
	UTextureRenderTarget2D* RoughnessRT = nullptr;

	UFUNCTION(BlueprintCallable, Category = "SWE Configuration")
	bool LoadParametersFromJson(const FString& FilePath);

	UFUNCTION(BlueprintCallable, Category = "SWE Configuration")
	bool SaveParametersToJson(const FString& FilePath);

	UFUNCTION(BlueprintCallable, Category = "SWE | Buoyancy")
	float GetWaterHeightAtLocation(const FVector& WorldLocation) const;

private:
	float SimulationTime = 0.0f;
	int32 PaddedSizeX = 512;
	int32 PaddedSizeY = 512;
	bool bInitialized = false;

	// Staging textures for double-buffered async readback
	FTextureRHIRef StagingTextures[2];
	int32 StagingWriteIndex = 0;
	int32 StagingReadIndex = 1;

	// Thread-safe double-buffered CPU height cache (in meters)
	TArray<float> CPUHeightData[2];
	std::atomic<int32> ActiveCPUBufferIndex{0};

	float GetCachedHeight(int32 X, int32 Y) const;

	// Persistent graphics buffers for simulation states
	TRefCountPtr<IPooledRenderTarget> TexTerrain;
	TRefCountPtr<IPooledRenderTarget> TexTerrainBulk;
	TRefCountPtr<IPooledRenderTarget> TexH;
	TRefCountPtr<IPooledRenderTarget> TexQ_x;
	TRefCountPtr<IPooledRenderTarget> TexQ_y;
	TRefCountPtr<IPooledRenderTarget> Texh;
	TRefCountPtr<IPooledRenderTarget> Texq_x;
	TRefCountPtr<IPooledRenderTarget> Texq_y;
	TRefCountPtr<IPooledRenderTarget> Texhbar;
	TRefCountPtr<IPooledRenderTarget> TexhbarOld;
	TRefCountPtr<IPooledRenderTarget> Texqbar_x;
	TRefCountPtr<IPooledRenderTarget> Texqbar_y;
	TRefCountPtr<IPooledRenderTarget> Texhtilde;
	TRefCountPtr<IPooledRenderTarget> TexhtildeOld;
	TRefCountPtr<IPooledRenderTarget> Texqtilde_x;
	TRefCountPtr<IPooledRenderTarget> Texqtilde_y;
	TRefCountPtr<IPooledRenderTarget> Texubar_x;
	TRefCountPtr<IPooledRenderTarget> Texubar_y;
	TRefCountPtr<IPooledRenderTarget> TexFoam;
	TRefCountPtr<IPooledRenderTarget> TexRoughness;

	// Stateful complex textures array for wave FFT propagation
	TRefCountPtr<IPooledRenderTarget> TexHPos;
	TRefCountPtr<IPooledRenderTarget> TexHNeg;

	void InitializeSimulation();
	void AllocatePersistentTargets(FRHICommandListImmediate& RHICmdList);
	void SetupInitialStates(FRHICommandListImmediate& RHICmdList);
	void AssignConstants(FSimConstants& OutConstants) const;
	
	void ExecuteSimulation_RenderThread(FRHICommandListImmediate& RHICmdList, const FSimConstants& Constants);
	void DispatchFFT_RenderThread(FRDGBuilder& GraphBuilder, FRDGTextureRef TargetTexture, int32 SizeX, int32 SizeY, bool bInverse, int32 NumLayers);
};
