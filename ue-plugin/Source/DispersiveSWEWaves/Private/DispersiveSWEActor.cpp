#include "DispersiveSWEActor.h"
#include "Materials/MaterialInstanceDynamic.h"
#include "ProceduralMeshComponent.h"
#include "Components/SceneComponent.h"
#include "Engine/StaticMesh.h"
#include "UObject/ConstructorHelpers.h"
#include "Math/Float16Color.h"
#include "Components/StaticMeshComponent.h"
#include "TextureResource.h"


ADispersiveSWEActor::ADispersiveSWEActor() {
    PrimaryActorTick.bCanEverTick = false;

    // Root Component
    USceneComponent* Root = CreateDefaultSubobject<USceneComponent>(TEXT("Root"));
    RootComponent = Root;

    // Attach Scene Capture (positioned facing straight down)
    TerrainCaptureComponent = CreateDefaultSubobject<USceneCaptureComponent2D>(TEXT("TerrainCapture"));
    TerrainCaptureComponent->SetupAttachment(Root);
    TerrainCaptureComponent->bCaptureEveryFrame = false;
    TerrainCaptureComponent->bCaptureOnMovement = false;
    TerrainCaptureComponent->SetRelativeLocation(FVector(0.0f, 0.0f, 5000.0f)); // default 50m above
    TerrainCaptureComponent->SetRelativeRotation(FRotator(-90.0f, 0.0f, -90.0f)); // Pointing down

    // Attach water mesh representation (flat plane)
    WaterMeshComponent = CreateDefaultSubobject<UStaticMeshComponent>(TEXT("WaterMesh"));
    WaterMeshComponent->SetupAttachment(Root);
    WaterMeshComponent->SetCollisionEnabled(ECollisionEnabled::NoCollision);

    // Attach SWE Simulation orchestrator
    SimComponent = CreateDefaultSubobject<USimulator>(TEXT("SWESimulator"));

    // Set defaults
    GridResolution = 512;
    WaterLevel = 0.0f;
    CapturedWorldWidth = 51200.0f;
    bAutoFitToTerrain = true;
    bAutoLoadDefaultAssets = true;
    StaticMeshDefaultSize = 4200.0f;

    // Attempt to resolve default plane1024 mesh
    static ConstructorHelpers::FObjectFinder<UStaticMesh> PlaneMeshFinder(TEXT("/DispersiveSWEWaves/Meshes/plane1024"));
    if (PlaneMeshFinder.Succeeded()) {
        WaterStaticMeshAsset = PlaneMeshFinder.Object;
        WaterMeshComponent->SetStaticMesh(WaterStaticMeshAsset);
    }

    // Attempt to resolve the default water material
    static ConstructorHelpers::FObjectFinder<UMaterialInterface> WaterMaterialFinder(TEXT("/DispersiveSWEWaves/Materials/M_DispersiveSWEBaseWater"));
    if (WaterMaterialFinder.Succeeded()) {
        BaseWaterMaterial = WaterMaterialFinder.Object;
    }
}

void ADispersiveSWEActor::OnConstruction(const FTransform& Transform) {
    Super::OnConstruction(Transform);

    // Handle auto-fit logic in editor time so it's instantly visual to the developer
    FitToTerrain();

    // Setup static water mesh representation and scale it
    if (WaterMeshComponent) {
        WaterMeshComponent->SetVisibility(true);
        if (WaterStaticMeshAsset) {
            WaterMeshComponent->SetStaticMesh(WaterStaticMeshAsset);
            float ScaleFactor = CapturedWorldWidth / FMath::Max(1.0f, StaticMeshDefaultSize);
            WaterMeshComponent->SetWorldScale3D(FVector(ScaleFactor, ScaleFactor, 1.0f));
        } else if (bAutoLoadDefaultAssets) {
            UStaticMesh* DefaultPlane = Cast<UStaticMesh>(StaticLoadObject(UStaticMesh::StaticClass(), nullptr, TEXT("/Engine/BasicShapes/Plane")));
            if (DefaultPlane) {
                WaterMeshComponent->SetStaticMesh(DefaultPlane);
                float ScaleFactor = CapturedWorldWidth / 100.0f; // Basic shapes plane is 100x100
                WaterMeshComponent->SetWorldScale3D(FVector(ScaleFactor, ScaleFactor, 1.0f));
            }
        }

        if (BaseWaterMaterial) {
            WaterMeshComponent->SetMaterial(0, BaseWaterMaterial);
        }
    }
}

void ADispersiveSWEActor::BeginPlay() {
    UE_LOG(LogTemp, Warning, TEXT("ADispersiveSWEActor::BeginPlay() started"));

    // Re-run bounds calculation to ensure runtime matches any runtime changes
    FitToTerrain();

    // Set up static water mesh at runtime
    if (WaterMeshComponent) {
        if (WaterStaticMeshAsset) {
            WaterMeshComponent->SetStaticMesh(WaterStaticMeshAsset);
            float ScaleFactor = CapturedWorldWidth / FMath::Max(1.0f, StaticMeshDefaultSize);
            WaterMeshComponent->SetWorldScale3D(FVector(ScaleFactor, ScaleFactor, 1.0f));
        } else if (bAutoLoadDefaultAssets) {
            UStaticMesh* DefaultPlane = Cast<UStaticMesh>(StaticLoadObject(UStaticMesh::StaticClass(), nullptr, TEXT("/Engine/BasicShapes/Plane")));
            if (DefaultPlane) {
                WaterMeshComponent->SetStaticMesh(DefaultPlane);
                float ScaleFactor = CapturedWorldWidth / 100.0f;
                WaterMeshComponent->SetWorldScale3D(FVector(ScaleFactor, ScaleFactor, 1.0f));
            }
        }
    }

    // Programmatically allocate 32-bit float Render Targets to avoid asset cluttering
    TerrainCaptureRT = NewObject<UTextureRenderTarget2D>(this);
    TerrainCaptureRT->bCanCreateUAV = true;
    TerrainCaptureRT->InitCustomFormat(GridResolution, GridResolution, PF_FloatRGBA, false);
    TerrainCaptureRT->AddressX = TA_Clamp;
    TerrainCaptureRT->AddressY = TA_Clamp;

    TerrainRT = NewObject<UTextureRenderTarget2D>(this);
    TerrainRT->bCanCreateUAV = true;
    TerrainRT->InitCustomFormat(GridResolution, GridResolution, PF_FloatRGBA, false);
    TerrainRT->AddressX = TA_Clamp;
    TerrainRT->AddressY = TA_Clamp;

    DisplacementRT = NewObject<UTextureRenderTarget2D>(this);
    DisplacementRT->bCanCreateUAV = true;
    DisplacementRT->InitCustomFormat(GridResolution, GridResolution, PF_FloatRGBA, false);
    DisplacementRT->AddressX = TA_Clamp;
    DisplacementRT->AddressY = TA_Clamp;

    DisplacementPastRT = NewObject<UTextureRenderTarget2D>(this);
    DisplacementPastRT->bCanCreateUAV = true;
    DisplacementPastRT->InitCustomFormat(GridResolution, GridResolution, PF_FloatRGBA, false);
    DisplacementPastRT->AddressX = TA_Clamp;
    DisplacementPastRT->AddressY = TA_Clamp;

    VelocityRT = NewObject<UTextureRenderTarget2D>(this);
    VelocityRT->bCanCreateUAV = true;
    VelocityRT->InitCustomFormat(GridResolution, GridResolution, PF_FloatRGBA, false);
    VelocityRT->AddressX = TA_Clamp;
    VelocityRT->AddressY = TA_Clamp;

    VelocityPastRT = NewObject<UTextureRenderTarget2D>(this);
    VelocityPastRT->bCanCreateUAV = true;
    VelocityPastRT->InitCustomFormat(GridResolution, GridResolution, PF_FloatRGBA, false);
    VelocityPastRT->AddressX = TA_Clamp;
    VelocityPastRT->AddressY = TA_Clamp;

    AccelerationRT = NewObject<UTextureRenderTarget2D>(this);
    AccelerationRT->bCanCreateUAV = true;
    AccelerationRT->InitCustomFormat(GridResolution, GridResolution, PF_FloatRGBA, false);
    AccelerationRT->AddressX = TA_Clamp;
    AccelerationRT->AddressY = TA_Clamp;

    AccelerationPastRT = NewObject<UTextureRenderTarget2D>(this);
    AccelerationPastRT->bCanCreateUAV = true;
    AccelerationPastRT->InitCustomFormat(GridResolution, GridResolution, PF_FloatRGBA, false);
    AccelerationPastRT->AddressX = TA_Clamp;
    AccelerationPastRT->AddressY = TA_Clamp;

    NormalRT = NewObject<UTextureRenderTarget2D>(this);
    NormalRT->bCanCreateUAV = true;
    NormalRT->InitCustomFormat(GridResolution, GridResolution, PF_FloatRGBA, false);
    NormalRT->AddressX = TA_Clamp;
    NormalRT->AddressY = TA_Clamp;

    FoamRT = NewObject<UTextureRenderTarget2D>(this);
    FoamRT->bCanCreateUAV = true;
    FoamRT->InitCustomFormat(GridResolution, GridResolution, PF_FloatRGBA, false);
    FoamRT->AddressX = TA_Clamp;
    FoamRT->AddressY = TA_Clamp;

    JacobianDetRT = NewObject<UTextureRenderTarget2D>(this);
    JacobianDetRT->bCanCreateUAV = true;
    JacobianDetRT->InitCustomFormat(GridResolution, GridResolution, PF_FloatRGBA, false);
    JacobianDetRT->AddressX = TA_Clamp;
    JacobianDetRT->AddressY = TA_Clamp;

    RoughnessRT = NewObject<UTextureRenderTarget2D>(this);
    RoughnessRT->bCanCreateUAV = true;
    RoughnessRT->InitCustomFormat(GridResolution, 1, PF_FloatRGBA, false);
    RoughnessRT->AddressX = TA_Clamp;
    RoughnessRT->AddressY = TA_Clamp;

    // Setup Scene Capture Component properties
    float CameraZ = 5000.0f;
    if (TerrainCaptureComponent) {
        // Temporarily hide the water mesh component so it doesn't block the depth capture of the terrain below it
        if (WaterMeshComponent) {
            WaterMeshComponent->SetVisibility(false);
        }
        if (TerrainActor) {
            TerrainCaptureComponent->ShowOnlyActors.Add(TerrainActor);
            TerrainCaptureComponent->PrimitiveRenderMode = ESceneCapturePrimitiveRenderMode::PRM_UseShowOnlyList;
        }

        TerrainCaptureComponent->ProjectionType = ECameraProjectionMode::Orthographic;
        TerrainCaptureComponent->OrthoWidth = CapturedWorldWidth;
        TerrainCaptureComponent->TextureTarget = TerrainCaptureRT;
        TerrainCaptureComponent->CaptureSource = ESceneCaptureSource::SCS_SceneDepth; // Absolute scene depth in centimeters
        TerrainCaptureComponent->CaptureScene(); // Render the terrain depth once on start
        CameraZ = TerrainCaptureComponent->GetComponentLocation().Z;

        if (WaterMeshComponent) {
            WaterMeshComponent->SetVisibility(true);
        }

        // Diagnostic readback of captured depth pixels
        if (TerrainCaptureRT) {
            FTextureRenderTargetResource* Resource = TerrainCaptureRT->GameThread_GetRenderTargetResource();
            TArray<FFloat16Color> Pixels;
            if (Resource && Resource->ReadFloat16Pixels(Pixels)) {
                float MinR = 1e20f;
                float MaxR = -1e20f;
                float SumR = 0.f;
                for (const FFloat16Color& Pixel : Pixels) {
                    float Val = Pixel.R.GetFloat();
                    if (Val < MinR) MinR = Val;
                    if (Val > MaxR) MaxR = Val;
                    SumR += Val;
                }
                float AvgR = Pixels.Num() > 0 ? SumR / Pixels.Num() : 0.f;
                UE_LOG(LogTemp, Warning, TEXT("TerrainCaptureRT Diagnostic: MinR=%f, MaxR=%f, AvgR=%f, NumPixels=%d, CameraZ=%f"), MinR, MaxR, AvgR, Pixels.Num(), CameraZ);
            } else {
                UE_LOG(LogTemp, Warning, TEXT("TerrainCaptureRT Diagnostic: Failed to read pixels."));
            }
        }
    }

    // Configure and Bind Render Targets to Simulation Component
    if (SimComponent) {
        SimComponent->GridSizeX = GridResolution;
        SimComponent->GridSizeY = GridResolution;
        SimComponent->CapturedWorldWidth = CapturedWorldWidth;
        SimComponent->bAutoCalculateCellSize = true;
        SimComponent->WaterLevel = WaterLevel;
        SimComponent->TerrainCaptureCameraZ = CameraZ;

        SimComponent->FoamThreshold = FoamThreshold;
        SimComponent->FoamMultiplier = FoamMultiplier;
        SimComponent->FoamFade = FoamFade;
        SimComponent->FoamBlur = FoamBlur;

        SimComponent->IntegrationSamples = IntegrationSamples;
        SimComponent->RoughnessPower = RoughnessPower;

        SimComponent->TerrainHeightInputRT = TerrainCaptureRT;
        SimComponent->TerrainRT = TerrainRT;
        SimComponent->DisplacementRT = DisplacementRT;
        SimComponent->DisplacementPastRT = DisplacementPastRT;
        SimComponent->VelocityRT = VelocityRT;
        SimComponent->VelocityPastRT = VelocityPastRT;
        SimComponent->AccelerationRT = AccelerationRT;
        SimComponent->AccelerationPastRT = AccelerationPastRT;
        SimComponent->NormalRT = NormalRT;
        SimComponent->FoamRT = FoamRT;
        SimComponent->JacobianDetRT = JacobianDetRT;
        SimComponent->RoughnessRT = RoughnessRT;
        SimComponent->bExportCapturedTerrain = bExportCapturedTerrain;
    }

    // Call Super::BeginPlay() after components are fully configured to trigger correct InitializeSimulation grid sizes
    Super::BeginPlay();

    // Try to load default material if not set
    if (!BaseWaterMaterial && bAutoLoadDefaultAssets) {
        BaseWaterMaterial = Cast<UMaterialInterface>(StaticLoadObject(UMaterialInterface::StaticClass(), nullptr, TEXT("/DispersiveSWEWaves/Materials/M_DispersiveSWEBaseWater")));
    }

    // Create and bind dynamic material instance
    if (BaseWaterMaterial) {
        DynamicWaterMaterial = UMaterialInstanceDynamic::Create(BaseWaterMaterial, this);
        DynamicWaterMaterial->SetTextureParameterValue(FName("DisplacementMap"), DisplacementRT);
        DynamicWaterMaterial->SetTextureParameterValue(FName("Displacement Map"), DisplacementRT);
        DynamicWaterMaterial->SetTextureParameterValue(FName("Displacement"), DisplacementRT);

        DynamicWaterMaterial->SetTextureParameterValue(FName("DisplacementPastMap"), DisplacementPastRT);
        DynamicWaterMaterial->SetTextureParameterValue(FName("Displacement Past Map"), DisplacementPastRT);
        DynamicWaterMaterial->SetTextureParameterValue(FName("DisplacementPast"), DisplacementPastRT);

        DynamicWaterMaterial->SetTextureParameterValue(FName("VelocityMap"), VelocityRT);
        DynamicWaterMaterial->SetTextureParameterValue(FName("Velocity Map"), VelocityRT);
        DynamicWaterMaterial->SetTextureParameterValue(FName("Velocity"), VelocityRT);

        DynamicWaterMaterial->SetTextureParameterValue(FName("VelocityPastMap"), VelocityPastRT);
        DynamicWaterMaterial->SetTextureParameterValue(FName("Velocity Past Map"), VelocityPastRT);
        DynamicWaterMaterial->SetTextureParameterValue(FName("VelocityPast"), VelocityPastRT);

        DynamicWaterMaterial->SetTextureParameterValue(FName("AccelerationMap"), AccelerationRT);
        DynamicWaterMaterial->SetTextureParameterValue(FName("Acceleration Map"), AccelerationRT);
        DynamicWaterMaterial->SetTextureParameterValue(FName("Acceleration"), AccelerationRT);

        DynamicWaterMaterial->SetTextureParameterValue(FName("AccelerationPastMap"), AccelerationPastRT);
        DynamicWaterMaterial->SetTextureParameterValue(FName("Acceleration Past Map"), AccelerationPastRT);
        DynamicWaterMaterial->SetTextureParameterValue(FName("AccelerationPast"), AccelerationPastRT);

        DynamicWaterMaterial->SetTextureParameterValue(FName("NormalMap"), NormalRT);
        DynamicWaterMaterial->SetTextureParameterValue(FName("Normal Map"), NormalRT);
        DynamicWaterMaterial->SetTextureParameterValue(FName("SurfaceNormal"), NormalRT);
        DynamicWaterMaterial->SetTextureParameterValue(FName("Surface Normal"), NormalRT);

        DynamicWaterMaterial->SetTextureParameterValue(FName("FoamMap"), FoamRT);
        DynamicWaterMaterial->SetTextureParameterValue(FName("Foam Map"), FoamRT);
        DynamicWaterMaterial->SetTextureParameterValue(FName("Foam"), FoamRT);

        DynamicWaterMaterial->SetTextureParameterValue(FName("JacobianDetMap"), JacobianDetRT);
        DynamicWaterMaterial->SetTextureParameterValue(FName("Jacobian Det Map"), JacobianDetRT);
        DynamicWaterMaterial->SetTextureParameterValue(FName("JacobianDet"), JacobianDetRT);

        DynamicWaterMaterial->SetTextureParameterValue(FName("RoughnessMap"), RoughnessRT);
        DynamicWaterMaterial->SetTextureParameterValue(FName("Roughness Map"), RoughnessRT);
        DynamicWaterMaterial->SetTextureParameterValue(FName("Roughness"), RoughnessRT);
        DynamicWaterMaterial->SetTextureParameterValue(FName("RoughnessLUT"), RoughnessRT);
        DynamicWaterMaterial->SetTextureParameterValue(FName("Roughness LUT"), RoughnessRT);

        DynamicWaterMaterial->SetTextureParameterValue(FName("TerrainHeightMap"), TerrainRT);
        DynamicWaterMaterial->SetTextureParameterValue(FName("Terrain Height Map"), TerrainRT);
        DynamicWaterMaterial->SetTextureParameterValue(FName("TerrainHeight"), TerrainRT);
        DynamicWaterMaterial->SetTextureParameterValue(FName("Terrain Height"), TerrainRT);
        DynamicWaterMaterial->SetTextureParameterValue(FName("TerrainCapture"), TerrainRT);
        DynamicWaterMaterial->SetTextureParameterValue(FName("Terrain Capture"), TerrainRT);

        WaterMeshComponent->SetMaterial(0, DynamicWaterMaterial);
    }
}

float ADispersiveSWEActor::GetWaterHeightAtLocation(const FVector& WorldLocation) const {
    if (SimComponent) {
        return SimComponent->GetWaterHeightAtLocation(WorldLocation);
    }
    return WaterLevel;
}

FVector ADispersiveSWEActor::GetWaterVelocityAtLocation(const FVector& WorldLocation) const {
    if (SimComponent) {
        return SimComponent->GetWaterVelocityAtLocation(WorldLocation);
    }
    return FVector::ZeroVector;
}

FVector ADispersiveSWEActor::GetWaterAccelerationAtLocation(const FVector& WorldLocation) const {
    if (SimComponent) {
        return SimComponent->GetWaterAccelerationAtLocation(WorldLocation);
    }
    return FVector::ZeroVector;
}

void ADispersiveSWEActor::FitToTerrain() {
    if (bAutoFitToTerrain && TerrainActor) {
        FVector Origin = FVector::ZeroVector;
        FVector BoxExtent = FVector::ZeroVector;
        bool bBoundsFound = false;

        if (!bBoundsFound) {
            TerrainActor->GetActorBounds(false, Origin, BoxExtent);

            // Fall back to landscape pivot and CapturedWorldWidth if bounds are zero (World Partition proxy case)
            if (BoxExtent.X <= 10.0f || BoxExtent.Y <= 10.0f) {
                if (CapturedWorldWidth <= 10.0f) {
                    CapturedWorldWidth = 51200.0f; // Reset to default 512m
                }
                // Landscape pivot is at the bottom-left corner; offset by half-width to center it
                Origin = TerrainActor->GetActorLocation() + FVector(CapturedWorldWidth * 0.5f, CapturedWorldWidth * 0.5f, 0.0f);
                BoxExtent = FVector(CapturedWorldWidth * 0.5f, CapturedWorldWidth * 0.5f, 1000.0f);
            } else {
                CapturedWorldWidth = FMath::Max(BoxExtent.X, BoxExtent.Y) * 2.0f;
            }
        }

        // Center horizontally and place actor at Z = 0.0f
        FVector NewLoc = Origin;
        NewLoc.Z = 0.0f;
        SetActorLocation(NewLoc);

        if (WaterMeshComponent) {
            WaterMeshComponent->SetWorldScale3D(FVector(1.0f, 1.0f, 1.0f));
        }
        if (TerrainCaptureComponent) {
            TerrainCaptureComponent->OrthoWidth = CapturedWorldWidth;
            
            // Put scene capture above the highest point of the terrain bounds
            float TargetRelativeZ = (Origin.Z + BoxExtent.Z + 1000.0f) - WaterLevel;
            TargetRelativeZ = FMath::Max(TargetRelativeZ, 5000.0f);
            TerrainCaptureComponent->SetRelativeLocation(FVector(0.0f, 0.0f, TargetRelativeZ));
        }

        if (SimComponent) {
            SimComponent->WaterLevel = WaterLevel;
            SimComponent->CapturedWorldWidth = CapturedWorldWidth;
        }
    }
}
