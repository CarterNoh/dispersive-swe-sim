#pragma once

#include "CoreMinimal.h"
#include "GameFramework/Actor.h"
#include "Components/SceneCaptureComponent2D.h"
#include "Engine/TextureRenderTarget2D.h"
#include "Simulator.h"
#include "SimActor.generated.h"

UCLASS()
class DISPERSIVESWESIM_API ASimActor : public AActor
{
    GENERATED_BODY()

public:
    ASimActor();

    // Editor-time viewport preview and auto-fitting
    virtual void OnConstruction(const FTransform& Transform) override;

    UFUNCTION(BlueprintCallable, Category = "SWE | Buoyancy")
    float GetWaterHeightAtLocation(const FVector& WorldLocation) const;

protected:
    virtual void BeginPlay() override;

public:
    // --- Components ---
    UPROPERTY(VisibleAnywhere, BlueprintReadOnly, Category="SWE | Components")
    USceneCaptureComponent2D* TerrainCaptureComponent;

    UPROPERTY(VisibleAnywhere, BlueprintReadOnly, Category="SWE | Components")
    USimulator* SimComponent;

    UPROPERTY(VisibleAnywhere, BlueprintReadOnly, Category="SWE | Components")
    UStaticMeshComponent* WaterMeshComponent;

    // --- Configuration ---
    UPROPERTY(EditAnywhere, BlueprintReadOnly, Category="SWE | Configuration")
    UStaticMesh* WaterStaticMeshAsset;

    UPROPERTY(EditAnywhere, BlueprintReadOnly, Category="SWE | Configuration", meta=(ToolTip="The default size of the static mesh in centimeters. Unreal Engine's basic shape plane is 100.0."))
    float StaticMeshDefaultSize = 4200.0f;

    UPROPERTY(EditAnywhere, BlueprintReadOnly, Category="SWE | Configuration")
    int32 GridResolution = 512;

    UPROPERTY(EditAnywhere, BlueprintReadOnly, Category="SWE | Configuration", meta=(ToolTip="If set, the simulator will automatically position and scale itself to match this actor's bounds."))
    AActor* TerrainActor;

    UPROPERTY(EditAnywhere, BlueprintReadOnly, Category="SWE | Configuration", meta=(ToolTip="If true, the water mesh and ortho capture will automatically fit to the bounds of the Terrain Actor."))
    bool bAutoFitToTerrain = true;

    UPROPERTY(EditAnywhere, BlueprintReadOnly, Category="SWE | Configuration", meta=(ToolTip="If true, will automatically load the Engine basic plane mesh if none is set."))
    bool bAutoLoadDefaultAssets = true;

    UPROPERTY(EditAnywhere, BlueprintReadOnly, Category="SWE | Configuration")
    float CapturedWorldWidth = 51200.0f; // Centimeters

    UPROPERTY(EditAnywhere, BlueprintReadOnly, Category="SWE | Configuration")
    float WaterLevel = 0.0f; // Centimeters (Water plane height)

    UPROPERTY(EditAnywhere, BlueprintReadOnly, Category="SWE | Configuration")
    float FoamThreshold = -0.5f;

    UPROPERTY(EditAnywhere, BlueprintReadOnly, Category="SWE | Configuration")
    float FoamMultiplier = 1.0f;

    UPROPERTY(EditAnywhere, BlueprintReadOnly, Category="SWE | Configuration")
    float FoamFade = 0.1f;

    UPROPERTY(EditAnywhere, BlueprintReadOnly, Category="SWE | Configuration")
    float FoamBlur = 0.3f;

    UPROPERTY(EditAnywhere, BlueprintReadOnly, Category="SWE | Configuration")
    float IntegrationSamples = 100.0f;

    UPROPERTY(EditAnywhere, BlueprintReadOnly, Category="SWE | Configuration")
    float RoughnessPower = 1.0f;

    UPROPERTY(EditAnywhere, BlueprintReadOnly, Category="SWE | Configuration")
    UMaterialInterface* BaseWaterMaterial;

    // If true, exports the captured terrain height map to terrain_captured.raw upon initialization. Defaults to false.
    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category="SWE | Debug", meta=(ToolTip="If true, exports the captured terrain height map to terrain_captured.raw upon initialization."))
    bool bExportCapturedTerrain = false;

private:
    // --- Runtime Render Targets ---
    UPROPERTY(VisibleAnywhere, BlueprintReadOnly, Category="SWE | Debug", meta=(AllowPrivateAccess="true"))
    UTextureRenderTarget2D* TerrainCaptureRT;

    UPROPERTY(VisibleAnywhere, BlueprintReadOnly, Category="SWE | Debug", meta=(AllowPrivateAccess="true"))
    UTextureRenderTarget2D* TerrainRT;

    UPROPERTY(VisibleAnywhere, BlueprintReadOnly, Category="SWE | Debug", meta=(AllowPrivateAccess="true"))
    UTextureRenderTarget2D* DisplacementRT;

    UPROPERTY(VisibleAnywhere, BlueprintReadOnly, Category="SWE | Debug", meta=(AllowPrivateAccess="true"))
    UTextureRenderTarget2D* DisplacementPastRT;

    UPROPERTY(VisibleAnywhere, BlueprintReadOnly, Category="SWE | Debug", meta=(AllowPrivateAccess="true"))
    UTextureRenderTarget2D* NormalRT;

    UPROPERTY(VisibleAnywhere, BlueprintReadOnly, Category="SWE | Debug", meta=(AllowPrivateAccess="true"))
    UTextureRenderTarget2D* FoamRT;

    UPROPERTY(VisibleAnywhere, BlueprintReadOnly, Category="SWE | Debug", meta=(AllowPrivateAccess="true"))
    UTextureRenderTarget2D* JacobianDetRT;

    UPROPERTY(VisibleAnywhere, BlueprintReadOnly, Category="SWE | Debug", meta=(AllowPrivateAccess="true"))
    UTextureRenderTarget2D* RoughnessRT;

    // --- Material ---
    UPROPERTY()
    UMaterialInstanceDynamic* DynamicWaterMaterial;

    void FitToTerrain();
};
