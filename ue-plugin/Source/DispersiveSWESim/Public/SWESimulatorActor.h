#pragma once

#include "CoreMinimal.h"
#include "GameFramework/Actor.h"
#include "Components/SceneCaptureComponent2D.h"
#include "Engine/TextureRenderTarget2D.h"
#include "DispersiveSWESimulator.h"
#include "SWESimulatorActor.generated.h"

UCLASS()
class DISPERSIVESWESIM_API ASWESimulatorActor : public AActor
{
    GENERATED_BODY()

public:
    ASWESimulatorActor();

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
    UDispersiveSWESimulator* SimComponent;

    UPROPERTY(VisibleAnywhere, BlueprintReadOnly, Category="SWE | Components")
    class UProceduralMeshComponent* WaterMeshComponent;

    // --- Configuration ---
    UPROPERTY(EditAnywhere, BlueprintReadOnly, Category="SWE | Configuration")
    int32 GridResolution = 512;

    UPROPERTY(EditAnywhere, BlueprintReadOnly, Category="SWE | Configuration", meta=(ToolTip="If set, the simulator will automatically position and scale itself to match this actor's bounds."))
    AActor* TerrainActor;

    UPROPERTY(EditAnywhere, BlueprintReadOnly, Category="SWE | Configuration", meta=(ToolTip="If true, the water mesh and ortho capture will automatically fit to the bounds of the Terrain Actor."))
    bool bAutoFitToTerrain = true;

    UPROPERTY(EditAnywhere, BlueprintReadOnly, Category="SWE | Configuration", meta=(ToolTip="If true, will automatically load the Engine basic plane mesh if none is set."))
    bool bAutoLoadDefaultAssets = true;

    UPROPERTY(EditAnywhere, BlueprintReadOnly, Category="SWE | Configuration")
    float CapturedWorldWidth = 12800.0f; // Centimeters (128 meters)

    UPROPERTY(EditAnywhere, BlueprintReadOnly, Category="SWE | Configuration")
    float WaterLevel = 0.0f; // Centimeters (Water plane height)

    UPROPERTY(EditAnywhere, BlueprintReadOnly, Category="SWE | Configuration")
    UMaterialInterface* BaseWaterMaterial;

private:
    void GenerateWaterGrid();

    // --- Runtime Render Targets ---
    UPROPERTY(VisibleAnywhere, BlueprintReadOnly, Category="SWE | Debug", meta=(AllowPrivateAccess="true"))
    UTextureRenderTarget2D* TerrainCaptureRT;

    UPROPERTY(VisibleAnywhere, BlueprintReadOnly, Category="SWE | Debug", meta=(AllowPrivateAccess="true"))
    UTextureRenderTarget2D* HeightOutputRT;

    UPROPERTY(VisibleAnywhere, BlueprintReadOnly, Category="SWE | Debug", meta=(AllowPrivateAccess="true"))
    UTextureRenderTarget2D* DisplacementOutputRT;

    UPROPERTY(VisibleAnywhere, BlueprintReadOnly, Category="SWE | Debug", meta=(AllowPrivateAccess="true"))
    UTextureRenderTarget2D* FoamOutputRT;

    // --- Material ---
    UPROPERTY()
    UMaterialInstanceDynamic* DynamicWaterMaterial;
};
