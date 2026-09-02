#pragma once

#include "CoreMinimal.h"
#include "Components/ActorComponent.h"
#include "Engine/TextureRenderTarget2D.h"
#include "MainSimShaders.h"
#include "FFTWaveShaders.h"
#include "ExportShaders.h"
#include "SimStateBuffers.h"
#include "SimReadback.h"
#include "Simulator.generated.h"

UCLASS(ClassGroup=(Custom), meta=(BlueprintSpawnableComponent))
class DISPERSIVESWESIM_API USimulator : public UActorComponent {
	GENERATED_BODY()

public:
	USimulator();
	virtual void BeginPlay() override;
	virtual void TickComponent(float DeltaTime, ELevelTick TickType, FActorComponentTickFunction* ThisTickFunction) override;

	//////// Simulation Parameters ////////

	// Number of cells in X dimension
	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "Sim|Grid")
	int32 GridSizeX = 512;

	// Number of cells in Y dimension
	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "Sim|Grid")
	int32 GridSizeY = 512;

	// Physical cell size in centimeters per cell
	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "Sim|Grid")
	float CellSize = 100.0f;

	// Simulation timestep (sec)
	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "Sim|Grid")
	float TimeStep = 0.016666f;

	// Maximum allowed CFL condition for stability
	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "Sim|Grid")
	float CFLCondition = 0.25f;

	// Minimum water height for stability (cm)
	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "Sim|Grid")
	float MinWaterHeight = 0.1f;

	// Maximum depth in centimeters allowed for safe, stable wave simulation. If <= 0, automatically calculated.
	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "Sim|Grid")
	float MaxSafeDepth = 0.0f;
	float CalculatedMaxSafeDepth;

	// Safety multiplier applied to the calculated max safe depth curve to prevent numerical explosion near stability boundary. Defaults to 0.8.
	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "Sim|Grid", meta = (ToolTip = "Safety multiplier applied to the calculated max safe depth curve to prevent numerical explosion near the stability boundary. Defaults to 0.8."))
	float StabilitySafetyFactor = 0.8f;

	// Level of water free surface at start (cm)
	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "Sim|Grid")
	float WaterLevel = 0.0f;

	// Surface tension coefficient of water (Newtons)
	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "Sim|Grid")
	float SurfaceTension = 0.001f;

	// Fluid density (kg/m^3)
	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "Sim|Grid")
	float Density = 999.0f;

	//////// Decomposition Parameters ////////

	// Number of iterations for diffusion step; more iterations means more stable but also more expensive
	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "Sim|Decomposition")
	int32 DiffusionIterations = 128;

	// Maximum height of diffusion stencil in cells; higher means more diffusion in deep water
	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "Sim|Decomposition")
	int32 MaxDiffusionCells = 8;

	// Penalty factor for diffusion; higher means more diffusion and more stability but also more damping of waves
	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "Sim|Decomposition")
	float DiffusionPenalty = 0.01f;

	//////// SWE & Transport Parameters ////////

	// Max slope of bulk flow surface before energy dissipation occurs
	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "Sim|SWE & Transport")
	float SlopeLimit = 1.0f;

	// Factor for damping converging flow, lower is more attenuation
	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "Sim|SWE & Transport")
	float GammaTransport = 0.25f;

	// (cells) Thickness of the sponge layer used to absorb waves at the boundaries
	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "Sim|SWE & Transport")
	int32 SpongeThickness = 8;

	// Damping factor for Laplacian smoothing to reduce spikes and unstable grid-scale ripples
	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "Sim|SWE & Transport")
	float LaplacianDamping = 0.01f;

	//////// FFT Wave Parameters ////////

	// Fetch in kilometers for JONSWAP wave spectrum
	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "Sim|FFT Waves")
	float Fetch = 200.0f;

	// Wind speed in m/s at 10 meters above surface
	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "Sim|FFT Waves")
	float WindSpeed = 14.0f;

	// Wind angle in degrees from x-axis
	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "Sim|FFT Waves")
	float WindAngle = 45.0f;

	// Swell factor in range [0, 1]
	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "Sim|FFT Waves")
	float Swell = 0.3f;

	// Swell direction in degrees from x-axis
	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "Sim|FFT Waves")
	float SwellAngle = 45.0f;

	// Amount of horizontal displacement in waves
	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "Sim|FFT Waves")
	float Choppiness = 1.0f;

	// Cutoff for small wavelengths in meters
	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "Sim|FFT Waves|Filter")
	float FilterSmall = 0.0f;

	// Cutoff for large wavelengths in meters
	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "Sim|FFT Waves|Filter")
	float FilterBig = 10000.0f;

	// Filter dropoff width in meters
	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "Sim|FFT Waves|Filter")
	float FilterWidth = 1.0f;

	// Minimum filter strength (dimensionless)
	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "Sim|FFT Waves|Filter")
	float FilterMin = 0.01f;

	// Depth to start attenuating FFT waves in meters
	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "Sim|FFT Waves")
	float DepthCutoff = 4.0f;

	//////// eWave Parameters ////////

	// Number of discrete water depth solutions to compute for eWave dispersion correction (in meters)
	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "Sim|eWave")
	TArray<float> DepthLevels = { 1.0f, 2.0f, 4.0f, 16.0f, 64.0f };

	// Foam threshold, multiplier, decay fade and blur rates
	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "Sim|eWave|Foam")
	float FoamThreshold = -0.25f;

	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "Sim|eWave|Foam")
	float FoamMultiplier = 1.0f;

	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "Sim|eWave|Foam")
	float FoamFade = 0.1f;

	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "Sim|eWave|Foam")
	float FoamBlur = 2.0f;

	// Roughness scale integration samples and power parameter
	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "Sim|eWave|Roughness")
	float IntegrationSamples = 100.0f;

	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "Sim|eWave|Roughness")
	float RoughnessPower = 1.0f;

	//////// Render Targets ////////

	// Input Render Target containing the level's terrain/height capture map
	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "UE|Inputs")
	UTextureRenderTarget2D* TerrainHeightInputRT = nullptr;

	// The absolute Z height of the terrain capture camera (used to reconstruct world Z from depth)
	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "UE|Inputs")
	float TerrainCaptureCameraZ = 5000.0f;

	// Auto calculate CellSize based on CapturedWorldWidth and GridSizeX
	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "UE|Inputs")
	bool bAutoCalculateCellSize = true;

	// The width of the captured terrain region in the world, in centimeters
	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "UE|Inputs", meta = (EditCondition = "bAutoCalculateCellSize"))
	float CapturedWorldWidth = 51200.0f;

	// File path to a JSON configuration file to load parameters from on startup
	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "UE|Configuration")
	FString JsonConfigFilePath = "";

	// If true, exports the captured terrain height map to terrain_captured.raw upon initialization. Defaults to false.
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "UE|Debug")
	bool bExportCapturedTerrain = false;

	// Output target textures containing wave fields
	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "UE|Outputs")
	UTextureRenderTarget2D* TerrainRT = nullptr;

	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "UE|Outputs")
	UTextureRenderTarget2D* DisplacementRT = nullptr;

	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "UE|Outputs")
	UTextureRenderTarget2D* DisplacementPastRT = nullptr;

	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "UE|Outputs")
	UTextureRenderTarget2D* VelocityRT = nullptr;

	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "UE|Outputs")
	UTextureRenderTarget2D* VelocityPastRT = nullptr;

	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "UE|Outputs")
	UTextureRenderTarget2D* AccelerationRT = nullptr;

	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "UE|Outputs")
	UTextureRenderTarget2D* AccelerationPastRT = nullptr;

	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "UE|Outputs")
	UTextureRenderTarget2D* NormalRT = nullptr;

	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "UE|Outputs")
	UTextureRenderTarget2D* FoamRT = nullptr;

	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "UE|Outputs")
	UTextureRenderTarget2D* JacobianDetRT = nullptr;

	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "UE|Outputs")
	UTextureRenderTarget2D* RoughnessRT = nullptr;

	UFUNCTION(BlueprintCallable, Category = "UE|Configuration")
	bool LoadParametersFromJson(const FString& FilePath);

	UFUNCTION(BlueprintCallable, Category = "UE|Configuration")
	bool SaveParametersToJson(const FString& FilePath);

	UFUNCTION(BlueprintCallable, Category = "UE|Buoyancy")
	float GetWaterHeightAtLocation(const FVector& WorldLocation) const;

	UFUNCTION(BlueprintCallable, Category = "UE|Buoyancy")
	FVector GetWaterVelocityAtLocation(const FVector& WorldLocation) const;

	UFUNCTION(BlueprintCallable, Category = "UE|Buoyancy")
	FVector GetWaterAccelerationAtLocation(const FVector& WorldLocation) const;

private:
	float SimulationTime = 0.0f;
	int32 PaddedSizeX = 512;
	int32 PaddedSizeY = 512;
	bool bInitialized = false;

	// Persistent GPU graphics state buffers
	FSimStateBuffers StateBuffers;

	// Double-buffered async GPU readback handler
	FSimReadbackHandler ReadbackHandler;

	void InitializeSimulation();
	void SetupInitialStates(FRHICommandListImmediate& RHICmdList);
	void AssignSimulationConstants(FSimConstants& OutConstants) const;
	void AssignFFTWaveConstants(FFFTWaveConstants& OutConstants) const;
	void AssignExportConstants(FExportConstants& OutConstants) const;
	
	void ExecuteSimulation_RenderThread(
		FRHICommandListImmediate& RHICmdList,
		const FSimConstants& SimConstants,
		const FFFTWaveConstants& FFTWaveConstants,
		const FExportConstants& ExportConstants);
};
