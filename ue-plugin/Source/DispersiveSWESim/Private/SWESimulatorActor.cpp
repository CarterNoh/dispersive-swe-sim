#include "SWESimulatorActor.h"
#include "Materials/MaterialInstanceDynamic.h"
#include "Components/StaticMeshComponent.h"
#include "Components/SceneComponent.h"
#include "Engine/StaticMesh.h"
#include "UObject/ConstructorHelpers.h"
#include "Math/Float16Color.h"

#if WITH_EDITOR
#include "LandscapeProxy.h"
#include "LandscapeInfo.h"
#endif

ASWESimulatorActor::ASWESimulatorActor()
{
    PrimaryActorTick.bCanEverTick = false;

#if WITH_EDITORONLY_DATA
    bIsSpatiallyLoaded = false;
#endif

    // 1. Root Component
    USceneComponent* Root = CreateDefaultSubobject<USceneComponent>(TEXT("Root"));
    RootComponent = Root;

    // 2. Attach Scene Capture (positioned facing straight down)
    TerrainCaptureComponent = CreateDefaultSubobject<USceneCaptureComponent2D>(TEXT("TerrainCapture"));
    TerrainCaptureComponent->SetupAttachment(Root);
    TerrainCaptureComponent->SetRelativeLocation(FVector(0.0f, 0.0f, 5000.0f)); // default 50m above
    TerrainCaptureComponent->SetRelativeRotation(FRotator(-90.0f, 0.0f, 0.0f)); // Pointing down

    // 3. Attach water mesh representation (flat plane)
    WaterMeshComponent = CreateDefaultSubobject<UStaticMeshComponent>(TEXT("WaterMesh"));
    WaterMeshComponent->SetupAttachment(Root);
    WaterMeshComponent->SetCollisionEnabled(ECollisionEnabled::NoCollision);

    // 4. Attach SWE Simulation orchestrator
    SimComponent = CreateDefaultSubobject<UDispersiveSWESimulator>(TEXT("SWESimulator"));

    // Set defaults
    GridResolution = 512;
    WaterLevel = 0.0f;
    CapturedWorldWidth = 12800.0f;
    bAutoFitToTerrain = true;
    bAutoLoadDefaultAssets = true;

    // Attempt to resolve the default Engine plane mesh
    static ConstructorHelpers::FObjectFinder<UStaticMesh> PlaneMeshFinder(TEXT("/Engine/BasicShapes/Plane.Plane"));
    if (PlaneMeshFinder.Succeeded() && WaterMeshComponent)
    {
        WaterMeshComponent->SetStaticMesh(PlaneMeshFinder.Object);
    }

    // Attempt to resolve the default water material
    static ConstructorHelpers::FObjectFinder<UMaterialInterface> WaterMaterialFinder(TEXT("/DispersiveSWESim/M_SWEDefaultWater"));
    if (WaterMaterialFinder.Succeeded())
    {
        BaseWaterMaterial = WaterMaterialFinder.Object;
        if (WaterMeshComponent)
        {
            WaterMeshComponent->SetMaterial(0, BaseWaterMaterial);
        }
    }
}

void ASWESimulatorActor::OnConstruction(const FTransform& Transform)
{
    Super::OnConstruction(Transform);

    // Handle auto-fit logic in editor time so it's instantly visual to the developer
    if (bAutoFitToTerrain && TerrainActor)
    {
        FVector Origin = FVector::ZeroVector;
        FVector BoxExtent = FVector::ZeroVector;
        bool bBoundsFound = false;

#if WITH_EDITOR
        if (ALandscapeProxy* LandscapeProxy = Cast<ALandscapeProxy>(TerrainActor))
        {
            if (ULandscapeInfo* LandscapeInfo = LandscapeProxy->GetLandscapeInfo())
            {
                FBox CompleteBounds = LandscapeInfo->GetCompleteBounds();
                if (CompleteBounds.IsValid)
                {
                    Origin = CompleteBounds.GetCenter();
                    BoxExtent = CompleteBounds.GetExtent();
                    CapturedWorldWidth = FMath::Max(BoxExtent.X, BoxExtent.Y) * 2.0f;
                    bBoundsFound = true;
                }
            }
        }
#endif

        if (!bBoundsFound)
        {
            TerrainActor->GetActorBounds(false, Origin, BoxExtent);

            // Fall back to landscape pivot and CapturedWorldWidth if bounds are zero (World Partition proxy case)
            if (BoxExtent.X <= 10.0f || BoxExtent.Y <= 10.0f)
            {
                if (CapturedWorldWidth <= 10.0f)
                {
                    CapturedWorldWidth = 12800.0f; // Reset to default 128m
                }
                // Landscape pivot is at the bottom-left corner; offset by half-width to center it
                Origin = TerrainActor->GetActorLocation() + FVector(CapturedWorldWidth * 0.5f, CapturedWorldWidth * 0.5f, 0.0f);
                BoxExtent = FVector(CapturedWorldWidth * 0.5f, CapturedWorldWidth * 0.5f, 1000.0f);
            }
            else
            {
                CapturedWorldWidth = FMath::Max(BoxExtent.X, BoxExtent.Y) * 2.0f;
            }
        }

        // Center horizontally and set the vertical level to WaterLevel
        FVector NewLoc = Origin;
        NewLoc.Z = WaterLevel;
        SetActorLocation(NewLoc);

        // Scale the mesh component to match the bounding box extent.
        // The default Engine plane (/Engine/BasicShapes/Plane) is 100x100 units.
        // Extent is half-size, so scale factor is BoxExtent / 50.0f
        if (WaterMeshComponent)
        {
            WaterMeshComponent->SetWorldScale3D(FVector(BoxExtent.X / 50.0f, BoxExtent.Y / 50.0f, 1.0f));
        }
        if (TerrainCaptureComponent)
        {
            TerrainCaptureComponent->OrthoWidth = CapturedWorldWidth;
            
            // Put scene capture above the highest point of the terrain bounds
            float TargetRelativeZ = (Origin.Z + BoxExtent.Z + 1000.0f) - WaterLevel;
            TargetRelativeZ = FMath::Max(TargetRelativeZ, 5000.0f);
            TerrainCaptureComponent->SetRelativeLocation(FVector(0.0f, 0.0f, TargetRelativeZ));
        }

        if (SimComponent)
        {
            SimComponent->WaterLevel = WaterLevel;
            SimComponent->CapturedWorldWidth = CapturedWorldWidth;
        }
    }

    // Auto-resolve basic mesh fallback at construction time if none is assigned
    if (bAutoLoadDefaultAssets && WaterMeshComponent && !WaterMeshComponent->GetStaticMesh())
    {
        UStaticMesh* DefaultPlane = Cast<UStaticMesh>(StaticLoadObject(UStaticMesh::StaticClass(), nullptr, TEXT("/Engine/BasicShapes/Plane.Plane")));
        if (DefaultPlane)
        {
            WaterMeshComponent->SetStaticMesh(DefaultPlane);
        }
    }

    // Auto-assign the base water material so it is visible in the editor viewport
    if (BaseWaterMaterial && WaterMeshComponent)
    {
        WaterMeshComponent->SetMaterial(0, BaseWaterMaterial);
    }
}

void ASWESimulatorActor::BeginPlay()
{
    UE_LOG(LogTemp, Warning, TEXT("ASWESimulatorActor::BeginPlay() started"));

    // 1. Re-run bounds calculation to ensure runtime matches any runtime changes
    if (bAutoFitToTerrain && TerrainActor)
    {
        FVector Origin = FVector::ZeroVector;
        FVector BoxExtent = FVector::ZeroVector;
        bool bBoundsFound = false;

#if WITH_EDITOR
        if (ALandscapeProxy* LandscapeProxy = Cast<ALandscapeProxy>(TerrainActor))
        {
            if (ULandscapeInfo* LandscapeInfo = LandscapeProxy->GetLandscapeInfo())
            {
                FBox CompleteBounds = LandscapeInfo->GetCompleteBounds();
                if (CompleteBounds.IsValid)
                {
                    Origin = CompleteBounds.GetCenter();
                    BoxExtent = CompleteBounds.GetExtent();
                    CapturedWorldWidth = FMath::Max(BoxExtent.X, BoxExtent.Y) * 2.0f;
                    bBoundsFound = true;
                }
            }
        }
#endif

        if (!bBoundsFound)
        {
            TerrainActor->GetActorBounds(false, Origin, BoxExtent);

            // Fall back to landscape pivot and CapturedWorldWidth if bounds are zero (World Partition proxy case)
            if (BoxExtent.X <= 10.0f || BoxExtent.Y <= 10.0f)
            {
                if (CapturedWorldWidth <= 10.0f)
                {
                    CapturedWorldWidth = 12800.0f; // Reset to default 128m
                }
                // Landscape pivot is at the bottom-left corner; offset by half-width to center it
                Origin = TerrainActor->GetActorLocation() + FVector(CapturedWorldWidth * 0.5f, CapturedWorldWidth * 0.5f, 0.0f);
                BoxExtent = FVector(CapturedWorldWidth * 0.5f, CapturedWorldWidth * 0.5f, 1000.0f);
            }
            else
            {
                CapturedWorldWidth = FMath::Max(BoxExtent.X, BoxExtent.Y) * 2.0f;
            }
        }

        FVector NewLoc = Origin;
        NewLoc.Z = WaterLevel;
        SetActorLocation(NewLoc);

        if (WaterMeshComponent)
        {
            WaterMeshComponent->SetWorldScale3D(FVector(BoxExtent.X / 50.0f, BoxExtent.Y / 50.0f, 1.0f));
        }
    }

    // 2. Programmatically allocate 32-bit float Render Targets to avoid asset cluttering
    TerrainCaptureRT = NewObject<UTextureRenderTarget2D>(this);
    TerrainCaptureRT->bCanCreateUAV = true;
    TerrainCaptureRT->InitCustomFormat(GridResolution, GridResolution, PF_FloatRGBA, false);
    TerrainCaptureRT->AddressX = TA_Clamp;
    TerrainCaptureRT->AddressY = TA_Clamp;

    HeightOutputRT = NewObject<UTextureRenderTarget2D>(this);
    HeightOutputRT->bCanCreateUAV = true;
    HeightOutputRT->InitCustomFormat(GridResolution, GridResolution, PF_FloatRGBA, false);
    HeightOutputRT->AddressX = TA_Clamp;
    HeightOutputRT->AddressY = TA_Clamp;

    DisplacementOutputRT = NewObject<UTextureRenderTarget2D>(this);
    DisplacementOutputRT->bCanCreateUAV = true;
    DisplacementOutputRT->InitCustomFormat(GridResolution, GridResolution, PF_FloatRGBA, false);
    DisplacementOutputRT->AddressX = TA_Clamp;
    DisplacementOutputRT->AddressY = TA_Clamp;

    FoamOutputRT = NewObject<UTextureRenderTarget2D>(this);
    FoamOutputRT->bCanCreateUAV = true;
    FoamOutputRT->InitCustomFormat(GridResolution, GridResolution, PF_FloatRGBA, false);
    FoamOutputRT->AddressX = TA_Clamp;
    FoamOutputRT->AddressY = TA_Clamp;

    // 3. Setup Scene Capture Component properties
    float CameraZ = 5000.0f;
    if (TerrainCaptureComponent)
    {
        // Temporarily hide the water mesh component so it doesn't block the depth capture of the terrain below it
        if (WaterMeshComponent)
        {
            WaterMeshComponent->SetVisibility(false);
        }

        TerrainCaptureComponent->ProjectionType = ECameraProjectionMode::Orthographic;
        TerrainCaptureComponent->OrthoWidth = CapturedWorldWidth;
        TerrainCaptureComponent->TextureTarget = TerrainCaptureRT;
        TerrainCaptureComponent->CaptureSource = ESceneCaptureSource::SCS_SceneDepth; // Absolute scene depth in centimeters
        TerrainCaptureComponent->CaptureScene(); // Render the terrain depth once on start
        CameraZ = TerrainCaptureComponent->GetComponentLocation().Z;

        if (WaterMeshComponent)
        {
            WaterMeshComponent->SetVisibility(true);
        }

        // Diagnostic readback of captured depth pixels
        if (TerrainCaptureRT)
        {
            FTextureRenderTargetResource* Resource = TerrainCaptureRT->GameThread_GetRenderTargetResource();
            TArray<FFloat16Color> Pixels;
            if (Resource && Resource->ReadFloat16Pixels(Pixels))
            {
                float MinR = 1e20f;
                float MaxR = -1e20f;
                float SumR = 0.f;
                for (const FFloat16Color& Pixel : Pixels)
                {
                    float Val = Pixel.R.GetFloat();
                    if (Val < MinR) MinR = Val;
                    if (Val > MaxR) MaxR = Val;
                    SumR += Val;
                }
                float AvgR = Pixels.Num() > 0 ? SumR / Pixels.Num() : 0.f;
                UE_LOG(LogTemp, Warning, TEXT("TerrainCaptureRT Diagnostic: MinR=%f, MaxR=%f, AvgR=%f, NumPixels=%d, CameraZ=%f"), MinR, MaxR, AvgR, Pixels.Num(), CameraZ);
            }
            else
            {
                UE_LOG(LogTemp, Warning, TEXT("TerrainCaptureRT Diagnostic: Failed to read pixels."));
            }
        }
    }

    // 4. Configure and Bind Render Targets to Simulation Component
    if (SimComponent)
    {
        SimComponent->GridSizeX = GridResolution;
        SimComponent->GridSizeY = GridResolution;
        SimComponent->CapturedWorldWidth = CapturedWorldWidth;
        SimComponent->bAutoCalculateCellSize = true;
        SimComponent->WaterLevel = WaterLevel;
        SimComponent->TerrainCaptureCameraZ = CameraZ;

        SimComponent->TerrainHeightInputRT = TerrainCaptureRT;
        SimComponent->HeightOutputRT = HeightOutputRT;
        SimComponent->DisplacementOutputRT = DisplacementOutputRT;
        SimComponent->FoamOutputRT = FoamOutputRT;
    }

    // Call Super::BeginPlay() after components are fully configured to trigger correct InitializeSimulation grid sizes
    Super::BeginPlay();

    // 5. Try to load default material if not set
    if (!BaseWaterMaterial && bAutoLoadDefaultAssets)
    {
        BaseWaterMaterial = Cast<UMaterialInterface>(StaticLoadObject(UMaterialInterface::StaticClass(), nullptr, TEXT("/DispersiveSWESim/M_SWEDefaultWater")));
    }

    // 6. Create and bind dynamic material instance
    if (BaseWaterMaterial && WaterMeshComponent)
    {
        DynamicWaterMaterial = UMaterialInstanceDynamic::Create(BaseWaterMaterial, this);
        DynamicWaterMaterial->SetTextureParameterValue(FName("HeightMap"), HeightOutputRT);
        DynamicWaterMaterial->SetTextureParameterValue(FName("Height Map"), HeightOutputRT);
        DynamicWaterMaterial->SetTextureParameterValue(FName("Height"), HeightOutputRT);
        DynamicWaterMaterial->SetTextureParameterValue(FName("DisplacementMap"), DisplacementOutputRT);
        DynamicWaterMaterial->SetTextureParameterValue(FName("Displacement Map"), DisplacementOutputRT);
        DynamicWaterMaterial->SetTextureParameterValue(FName("Displacement"), DisplacementOutputRT);
        DynamicWaterMaterial->SetTextureParameterValue(FName("FoamMap"), FoamOutputRT);
        DynamicWaterMaterial->SetTextureParameterValue(FName("Foam Map"), FoamOutputRT);
        DynamicWaterMaterial->SetTextureParameterValue(FName("Foam"), FoamOutputRT);
        WaterMeshComponent->SetMaterial(0, DynamicWaterMaterial);
    }
}

float ASWESimulatorActor::GetWaterHeightAtLocation(const FVector& WorldLocation) const
{
    if (SimComponent)
    {
        return SimComponent->GetWaterHeightAtLocation(WorldLocation);
    }
    return WaterLevel;
}
