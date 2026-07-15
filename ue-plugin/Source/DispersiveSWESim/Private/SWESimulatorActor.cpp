#include "SWESimulatorActor.h"
#include "Materials/MaterialInstanceDynamic.h"
#include "ProceduralMeshComponent.h"
#include "Components/SceneComponent.h"
#include "Engine/StaticMesh.h"
#include "UObject/ConstructorHelpers.h"
#include "Math/Float16Color.h"
#include "Components/StaticMeshComponent.h"
#include "TextureResource.h"


#if WITH_EDITOR
#include "LandscapeProxy.h"
#include "LandscapeInfo.h"
#endif

ASWESimulatorActor::ASWESimulatorActor()
{
    PrimaryActorTick.bCanEverTick = true;

#if WITH_EDITORONLY_DATA
    bIsSpatiallyLoaded = false;
#endif

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
    SimComponent = CreateDefaultSubobject<UDispersiveSWESimulator>(TEXT("SWESimulator"));

    // Set defaults
    GridResolution = 512;
    WaterLevel = 0.0f;
    CapturedWorldWidth = 51200.0f;
    bAutoFitToTerrain = true;
    bAutoLoadDefaultAssets = true;
    StaticMeshDefaultSize = 4200.0f;

    // Attempt to resolve default plane1024 mesh
    static ConstructorHelpers::FObjectFinder<UStaticMesh> PlaneMeshFinder(TEXT("/DispersiveSWESim/Meshes/plane1024"));
    if (PlaneMeshFinder.Succeeded())
    {
        WaterStaticMeshAsset = PlaneMeshFinder.Object;
        WaterMeshComponent->SetStaticMesh(WaterStaticMeshAsset);
    }

    // Attempt to resolve the default water material
    static ConstructorHelpers::FObjectFinder<UMaterialInterface> WaterMaterialFinder(TEXT("/DispersiveSWESim/Materials/M_PreviewOceanWater"));
    if (WaterMaterialFinder.Succeeded())
    {
        BaseWaterMaterial = WaterMaterialFinder.Object;
    }
}

void ASWESimulatorActor::OnConstruction(const FTransform& Transform)
{
    Super::OnConstruction(Transform);

    // Handle auto-fit logic in editor time so it's instantly visual to the developer
    FitToTerrain();

    // Setup static water mesh representation and scale it
    if (WaterMeshComponent)
    {
        WaterMeshComponent->SetVisibility(true);
        if (WaterStaticMeshAsset)
        {
            WaterMeshComponent->SetStaticMesh(WaterStaticMeshAsset);
            float ScaleFactor = CapturedWorldWidth / FMath::Max(1.0f, StaticMeshDefaultSize);
            WaterMeshComponent->SetWorldScale3D(FVector(ScaleFactor, ScaleFactor, 1.0f));
        }
        else if (bAutoLoadDefaultAssets)
        {
            UStaticMesh* DefaultPlane = Cast<UStaticMesh>(StaticLoadObject(UStaticMesh::StaticClass(), nullptr, TEXT("/Engine/BasicShapes/Plane")));
            if (DefaultPlane)
            {
                WaterMeshComponent->SetStaticMesh(DefaultPlane);
                float ScaleFactor = CapturedWorldWidth / 100.0f; // Basic shapes plane is 100x100
                WaterMeshComponent->SetWorldScale3D(FVector(ScaleFactor, ScaleFactor, 1.0f));
            }
        }

        if (BaseWaterMaterial)
        {
            WaterMeshComponent->SetMaterial(0, BaseWaterMaterial);
        }
    }
}

void ASWESimulatorActor::BeginPlay()
{
    UE_LOG(LogTemp, Warning, TEXT("ASWESimulatorActor::BeginPlay() started"));

    // Re-run bounds calculation to ensure runtime matches any runtime changes
    FitToTerrain();

    // Set up static water mesh at runtime
    if (WaterMeshComponent)
    {
        if (WaterStaticMeshAsset)
        {
            WaterMeshComponent->SetStaticMesh(WaterStaticMeshAsset);
            float ScaleFactor = CapturedWorldWidth / FMath::Max(1.0f, StaticMeshDefaultSize);
            WaterMeshComponent->SetWorldScale3D(FVector(ScaleFactor, ScaleFactor, 1.0f));
        }
        else if (bAutoLoadDefaultAssets)
        {
            UStaticMesh* DefaultPlane = Cast<UStaticMesh>(StaticLoadObject(UStaticMesh::StaticClass(), nullptr, TEXT("/Engine/BasicShapes/Plane")));
            if (DefaultPlane)
            {
                WaterMeshComponent->SetStaticMesh(DefaultPlane);
                float ScaleFactor = CapturedWorldWidth / 100.0f;
                WaterMeshComponent->SetWorldScale3D(FVector(ScaleFactor, ScaleFactor, 1.0f));
            }
        }
    }

    // Programmatically allocate 32-bit float Render Targets to avoid asset cluttering
    SimComponent->TerrainHeightInputRT = NewObject<UTextureRenderTarget2D>(this);
    SimComponent->TerrainHeightInputRT->bCanCreateUAV = true;
    SimComponent->TerrainHeightInputRT->InitCustomFormat(GridResolution, GridResolution, PF_FloatRGBA, false);
    SimComponent->TerrainHeightInputRT->AddressX = TA_Clamp;
    SimComponent->TerrainHeightInputRT->AddressY = TA_Clamp;

    SimComponent->TerrainRT = NewObject<UTextureRenderTarget2D>(this);
    SimComponent->TerrainRT->bCanCreateUAV = true;
    SimComponent->TerrainRT->InitCustomFormat(GridResolution, GridResolution, PF_FloatRGBA, false);
    SimComponent->TerrainRT->AddressX = TA_Clamp;
    SimComponent->TerrainRT->AddressY = TA_Clamp;

    SimComponent->DisplacementRT.SetNum(3);
    SimComponent->NormalRT.SetNum(2);
    SimComponent->FoamRT.SetNum(2);
    SimComponent->JacobianDetRT.SetNum(2);
    SimComponent->RoughnessRT.SetNum(2);

    for (int32 i = 0; i < 3; i++)
    {
        SimComponent->DisplacementRT[i] = NewObject<UTextureRenderTarget2D>(this);
        SimComponent->DisplacementRT[i]->bCanCreateUAV = true;
        SimComponent->DisplacementRT[i]->InitCustomFormat(GridResolution, GridResolution, PF_FloatRGBA, false);
        SimComponent->DisplacementRT[i]->AddressX = TA_Clamp;
        SimComponent->DisplacementRT[i]->AddressY = TA_Clamp;
    }

    for (int32 i = 0; i < 2; i++)
    {
        SimComponent->NormalRT[i] = NewObject<UTextureRenderTarget2D>(this);
        SimComponent->NormalRT[i]->bCanCreateUAV = true;
        SimComponent->NormalRT[i]->InitCustomFormat(GridResolution, GridResolution, PF_FloatRGBA, false);
        SimComponent->NormalRT[i]->AddressX = TA_Clamp;
        SimComponent->NormalRT[i]->AddressY = TA_Clamp;

        SimComponent->FoamRT[i] = NewObject<UTextureRenderTarget2D>(this);
        SimComponent->FoamRT[i]->bCanCreateUAV = true;
        SimComponent->FoamRT[i]->InitCustomFormat(GridResolution, GridResolution, PF_FloatRGBA, false);
        SimComponent->FoamRT[i]->AddressX = TA_Clamp;
        SimComponent->FoamRT[i]->AddressY = TA_Clamp;

        SimComponent->JacobianDetRT[i] = NewObject<UTextureRenderTarget2D>(this);
        SimComponent->JacobianDetRT[i]->bCanCreateUAV = true;
        SimComponent->JacobianDetRT[i]->InitCustomFormat(GridResolution, GridResolution, PF_FloatRGBA, false);
        SimComponent->JacobianDetRT[i]->AddressX = TA_Clamp;
        SimComponent->JacobianDetRT[i]->AddressY = TA_Clamp;

        SimComponent->RoughnessRT[i] = NewObject<UTextureRenderTarget2D>(this);
        SimComponent->RoughnessRT[i]->bCanCreateUAV = true;
        SimComponent->RoughnessRT[i]->InitCustomFormat(GridResolution, 1, PF_FloatRGBA, false);
        SimComponent->RoughnessRT[i]->AddressX = TA_Clamp;
        SimComponent->RoughnessRT[i]->AddressY = TA_Clamp;
    }

    // Setup Scene Capture Component properties
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
        TerrainCaptureComponent->TextureTarget = SimComponent->TerrainHeightInputRT;
        TerrainCaptureComponent->CaptureSource = ESceneCaptureSource::SCS_SceneDepth; // Absolute scene depth in centimeters
        TerrainCaptureComponent->CaptureScene(); // Render the terrain depth once on start
        CameraZ = TerrainCaptureComponent->GetComponentLocation().Z;

        if (WaterMeshComponent)
        {
            WaterMeshComponent->SetVisibility(true);
        }
    }

    if (SimComponent)
    {
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
    }

    // Call Super::BeginPlay() after components are fully configured to trigger correct InitializeSimulation grid sizes
    Super::BeginPlay();

    // Try to load default material if not set
    if (!BaseWaterMaterial && bAutoLoadDefaultAssets)
    {
        BaseWaterMaterial = Cast<UMaterialInterface>(StaticLoadObject(UMaterialInterface::StaticClass(), nullptr, TEXT("/DispersiveSWESim/Materials/M_PreviewOceanWater")));
    }

    // Create and bind dynamic material instance
    if (BaseWaterMaterial)
    {
        DynamicWaterMaterial = UMaterialInstanceDynamic::Create(BaseWaterMaterial, this);
        BindWaterMaterial();
        WaterMeshComponent->SetMaterial(0, DynamicWaterMaterial);
    }
}

void ASWESimulatorActor::BindWaterMaterial()
{
    if (!SimComponent || !DynamicWaterMaterial) return;

	int32 ReadIndex = SimComponent->GetReadIndex();
	int32 DispReadIndex = SimComponent->GetDispReadIndex();
	int32 DispPastIndex = SimComponent->GetDispPastIndex();

	DynamicWaterMaterial->SetTextureParameterValue(FName("DisplacementMap"), SimComponent->DisplacementRT[DispReadIndex]);
	DynamicWaterMaterial->SetTextureParameterValue(FName("Displacement Map"), SimComponent->DisplacementRT[DispReadIndex]);
	DynamicWaterMaterial->SetTextureParameterValue(FName("Displacement"), SimComponent->DisplacementRT[DispReadIndex]);

	DynamicWaterMaterial->SetTextureParameterValue(FName("DisplacementPastMap"), SimComponent->DisplacementRT[DispPastIndex]);
	DynamicWaterMaterial->SetTextureParameterValue(FName("Displacement Past Map"), SimComponent->DisplacementRT[DispPastIndex]);
	DynamicWaterMaterial->SetTextureParameterValue(FName("DisplacementPast"), SimComponent->DisplacementRT[DispPastIndex]);

	DynamicWaterMaterial->SetTextureParameterValue(FName("NormalMap"), SimComponent->NormalRT[ReadIndex]);
	DynamicWaterMaterial->SetTextureParameterValue(FName("Normal Map"), SimComponent->NormalRT[ReadIndex]);
	DynamicWaterMaterial->SetTextureParameterValue(FName("SurfaceNormal"), SimComponent->NormalRT[ReadIndex]);
	DynamicWaterMaterial->SetTextureParameterValue(FName("Surface Normal"), SimComponent->NormalRT[ReadIndex]);

	DynamicWaterMaterial->SetTextureParameterValue(FName("FoamMap"), SimComponent->FoamRT[ReadIndex]);
	DynamicWaterMaterial->SetTextureParameterValue(FName("Foam Map"), SimComponent->FoamRT[ReadIndex]);
	DynamicWaterMaterial->SetTextureParameterValue(FName("Foam"), SimComponent->FoamRT[ReadIndex]);

	DynamicWaterMaterial->SetTextureParameterValue(FName("JacobianDetMap"), SimComponent->JacobianDetRT[ReadIndex]);
	DynamicWaterMaterial->SetTextureParameterValue(FName("Jacobian Det Map"), SimComponent->JacobianDetRT[ReadIndex]);
	DynamicWaterMaterial->SetTextureParameterValue(FName("JacobianDet"), SimComponent->JacobianDetRT[ReadIndex]);

	DynamicWaterMaterial->SetTextureParameterValue(FName("RoughnessMap"), SimComponent->RoughnessRT[ReadIndex]);
	DynamicWaterMaterial->SetTextureParameterValue(FName("Roughness Map"), SimComponent->RoughnessRT[ReadIndex]);
	DynamicWaterMaterial->SetTextureParameterValue(FName("Roughness"), SimComponent->RoughnessRT[ReadIndex]);
	DynamicWaterMaterial->SetTextureParameterValue(FName("RoughnessLUT"), SimComponent->RoughnessRT[ReadIndex]);
	DynamicWaterMaterial->SetTextureParameterValue(FName("Roughness LUT"), SimComponent->RoughnessRT[ReadIndex]);

	// Bind TerrainRT parameter values (since terrain height does not change per-frame dynamically)
	DynamicWaterMaterial->SetTextureParameterValue(FName("TerrainHeightMap"), SimComponent->TerrainRT);
	DynamicWaterMaterial->SetTextureParameterValue(FName("Terrain Height Map"), SimComponent->TerrainRT);
	DynamicWaterMaterial->SetTextureParameterValue(FName("TerrainHeight"), SimComponent->TerrainRT);
	DynamicWaterMaterial->SetTextureParameterValue(FName("Terrain Height"), SimComponent->TerrainRT);
	DynamicWaterMaterial->SetTextureParameterValue(FName("TerrainCapture"), SimComponent->TerrainRT);
	DynamicWaterMaterial->SetTextureParameterValue(FName("Terrain Capture"), SimComponent->TerrainRT);
}

// ASWESimulatorActor.cpp
void ASWESimulatorActor::Tick(float DeltaTime)
{
	Super::Tick(DeltaTime);
    BindWaterMaterial();
}

float ASWESimulatorActor::GetWaterHeightAtLocation(const FVector& WorldLocation) const
{
    if (SimComponent)
    {
        return SimComponent->GetWaterHeightAtLocation(WorldLocation);
    }
    return WaterLevel;
}

void ASWESimulatorActor::FitToTerrain()
{
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
                    CapturedWorldWidth = 51200.0f; // Reset to default 512m
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

        // Center horizontally and place actor at Z = 0.0f
        FVector NewLoc = Origin;
        NewLoc.Z = 0.0f;
        SetActorLocation(NewLoc);

        if (WaterMeshComponent)
        {
            WaterMeshComponent->SetWorldScale3D(FVector(1.0f, 1.0f, 1.0f));
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
}


