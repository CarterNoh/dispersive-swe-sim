#pragma once

#include "CoreMinimal.h"
#include "GameFramework/Actor.h"
#include "Components/SceneCaptureComponent2D.h"
#include "Engine/TextureRenderTarget2D.h"
#include "Simulator.h"
#include "DispersiveSWEActor.generated.h"

UCLASS()
class DISPERSIVESWEWAVES_API ADispersiveSWEActor : public AActor {
    GENERATED_BODY()

public:
    ADispersiveSWEActor();

    // Editor-time viewport preview and auto-fitting
    virtual void OnConstruction(const FTransform& Transform) override;

    UFUNCTION(BlueprintCallable, Category = "DispersiveSWE|Buoyancy")
    float GetWaterHeightAtLocation(const FVector& WorldLocation) const;

    UFUNCTION(BlueprintCallable, Category = "DispersiveSWE|Buoyancy")
    FVector GetWaterVelocityAtLocation(const FVector& WorldLocation) const;

    UFUNCTION(BlueprintCallable, Category = "DispersiveSWE|Buoyancy")
    FVector GetWaterAccelerationAtLocation(const FVector& WorldLocation) const;

protected:
    virtual void BeginPlay() override;

public:
    // --- Components ---
    UPROPERTY(VisibleAnywhere, BlueprintReadOnly, Category = "DispersiveSWE|Components")
    USceneCaptureComponent2D* TerrainCaptureComponent;

    UPROPERTY(VisibleAnywhere, BlueprintReadOnly, Category = "DispersiveSWE|Components")
    USimulator* SimComponent;

    UPROPERTY(VisibleAnywhere, BlueprintReadOnly, Category = "DispersiveSWE|Components")
    UStaticMeshComponent* WaterMeshComponent;

    // --- Configuration ---
    UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "DispersiveSWE|Configuration")
    UStaticMesh* WaterStaticMeshAsset;

    UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "DispersiveSWE|Configuration", meta=(ToolTip="The default size of the static mesh in centimeters. Unreal Engine's basic shape plane is 100.0."))
    float StaticMeshDefaultSize = 4200.0f;

    UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "DispersiveSWE|Configuration")
    int32 GridResolution = 512;

    UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "DispersiveSWE|Configuration", meta=(ToolTip="If set, the simulator will automatically position and scale itself to match this actor's bounds."))
    AActor* TerrainActor;

    UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "DispersiveSWE|Configuration", meta=(ToolTip="If true, the water mesh and ortho capture will automatically fit to the bounds of the Terrain Actor."))
    bool bAutoFitToTerrain = true;

    UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "DispersiveSWE|Configuration", meta=(ToolTip="If true, will automatically load the Engine basic plane mesh if none is set."))
    bool bAutoLoadDefaultAssets = true;

    UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "DispersiveSWE|Configuration")
    float CapturedWorldWidth = 51200.0f; // Centimeters

    UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "DispersiveSWE|Configuration")
    float WaterLevel = 0.0f; // Centimeters (Water plane height)

    UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "DispersiveSWE|Material")
    float FoamThreshold = -0.5f;

    UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "DispersiveSWE|Material")
    float FoamMultiplier = 1.0f;

    UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "DispersiveSWE|Material")
    float FoamFade = 0.1f;

    UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "DispersiveSWE|Material")
    float FoamBlur = 0.3f;

    UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "DispersiveSWE|Material")
    float IntegrationSamples = 100.0f;

    UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "DispersiveSWE|Material")
    float RoughnessPower = 1.0f;

    UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "DispersiveSWE|Material")
    UMaterialInterface* BaseWaterMaterial;

    // If true, exports the captured terrain height map to terrain_captured.raw upon initialization. Defaults to false.
    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "DispersiveSWE|Debug", meta=(ToolTip="If true, exports the captured terrain height map to terrain_captured.raw upon initialization."))
    bool bExportCapturedTerrain = false;

private:
    // --- Runtime Render Targets ---
    UPROPERTY(VisibleAnywhere, BlueprintReadOnly, Category = "DispersiveSWE|Debug", meta=(AllowPrivateAccess="true"))
    UTextureRenderTarget2D* TerrainCaptureRT;

    UPROPERTY(VisibleAnywhere, BlueprintReadOnly, Category = "DispersiveSWE|Debug", meta=(AllowPrivateAccess="true"))
    UTextureRenderTarget2D* TerrainRT;

    UPROPERTY(VisibleAnywhere, BlueprintReadOnly, Category = "DispersiveSWE|Debug", meta=(AllowPrivateAccess="true"))
    UTextureRenderTarget2D* DisplacementRT;

    UPROPERTY(VisibleAnywhere, BlueprintReadOnly, Category = "DispersiveSWE|Debug", meta=(AllowPrivateAccess="true"))
    UTextureRenderTarget2D* DisplacementPastRT;

    UPROPERTY(VisibleAnywhere, BlueprintReadOnly, Category = "DispersiveSWE|Debug", meta=(AllowPrivateAccess="true"))
    UTextureRenderTarget2D* VelocityRT;

    UPROPERTY(VisibleAnywhere, BlueprintReadOnly, Category = "DispersiveSWE|Debug", meta=(AllowPrivateAccess="true"))
    UTextureRenderTarget2D* VelocityPastRT;

    UPROPERTY(VisibleAnywhere, BlueprintReadOnly, Category = "DispersiveSWE|Debug", meta=(AllowPrivateAccess="true"))
    UTextureRenderTarget2D* AccelerationRT;

    UPROPERTY(VisibleAnywhere, BlueprintReadOnly, Category = "DispersiveSWE|Debug", meta=(AllowPrivateAccess="true"))
    UTextureRenderTarget2D* AccelerationPastRT;

    UPROPERTY(VisibleAnywhere, BlueprintReadOnly, Category = "DispersiveSWE|Debug", meta=(AllowPrivateAccess="true"))
    UTextureRenderTarget2D* NormalRT;

    UPROPERTY(VisibleAnywhere, BlueprintReadOnly, Category = "DispersiveSWE|Debug", meta=(AllowPrivateAccess="true"))
    UTextureRenderTarget2D* FoamRT;

    UPROPERTY(VisibleAnywhere, BlueprintReadOnly, Category = "DispersiveSWE|Debug", meta=(AllowPrivateAccess="true"))
    UTextureRenderTarget2D* JacobianDetRT;

    UPROPERTY(VisibleAnywhere, BlueprintReadOnly, Category = "DispersiveSWE|Debug", meta=(AllowPrivateAccess="true"))
    UTextureRenderTarget2D* RoughnessRT;

    // --- Material ---
    UPROPERTY()
    UMaterialInstanceDynamic* DynamicWaterMaterial;

    void FitToTerrain();
};
