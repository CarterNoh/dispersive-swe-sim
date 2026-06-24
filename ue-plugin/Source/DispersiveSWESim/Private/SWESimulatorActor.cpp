#include "SWESimulatorActor.h"
#include "Materials/MaterialInstanceDynamic.h"
#include "Components/StaticMeshComponent.h"
#include "Components/SceneComponent.h"
#include "Engine/StaticMesh.h"
#include "UObject/ConstructorHelpers.h"

ASWESimulatorActor::ASWESimulatorActor()
{
    PrimaryActorTick.bCanEverTick = false;

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

    // 4. Attach SWE Simulation orchestrator
    SimComponent = CreateDefaultSubobject<UDispersiveSWESimulator>(TEXT("SWESimulator"));

    // Set defaults
    GridResolution = 128;
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
}

void ASWESimulatorActor::OnConstruction(const FTransform& Transform)
{
    Super::OnConstruction(Transform);

    // Handle auto-fit logic in editor time so it's instantly visual to the developer
    if (bAutoFitToTerrain && TerrainActor)
    {
        FVector Origin;
        FVector BoxExtent;
        TerrainActor->GetActorBounds(false, Origin, BoxExtent);

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

        // Set width and position scene capture to cover the whole terrain
        CapturedWorldWidth = FMath::Max(BoxExtent.X, BoxExtent.Y) * 2.0f;
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
}

void ASWESimulatorActor::BeginPlay()
{
    Super::BeginPlay();

    // 1. Re-run bounds calculation to ensure runtime matches any runtime changes
    if (bAutoFitToTerrain && TerrainActor)
    {
        FVector Origin;
        FVector BoxExtent;
        TerrainActor->GetActorBounds(false, Origin, BoxExtent);

        FVector NewLoc = Origin;
        NewLoc.Z = WaterLevel;
        SetActorLocation(NewLoc);

        if (WaterMeshComponent)
        {
            WaterMeshComponent->SetWorldScale3D(FVector(BoxExtent.X / 50.0f, BoxExtent.Y / 50.0f, 1.0f));
        }

        CapturedWorldWidth = FMath::Max(BoxExtent.X, BoxExtent.Y) * 2.0f;
    }

    // 2. Programmatically allocate 32-bit float Render Targets to avoid asset cluttering
    TerrainCaptureRT = NewObject<UTextureRenderTarget2D>(this);
    TerrainCaptureRT->InitCustomFormat(GridResolution, GridResolution, PF_R32_FLOAT, false);
    TerrainCaptureRT->AddressX = TA_Clamp;
    TerrainCaptureRT->AddressY = TA_Clamp;

    HeightOutputRT = NewObject<UTextureRenderTarget2D>(this);
    HeightOutputRT->InitCustomFormat(GridResolution, GridResolution, PF_R32_FLOAT, false);
    HeightOutputRT->AddressX = TA_Clamp;
    HeightOutputRT->AddressY = TA_Clamp;

    DisplacementOutputRT = NewObject<UTextureRenderTarget2D>(this);
    DisplacementOutputRT->InitCustomFormat(GridResolution, GridResolution, PF_A2B10G10R10, false);
    DisplacementOutputRT->AddressX = TA_Clamp;
    DisplacementOutputRT->AddressY = TA_Clamp;

    FoamOutputRT = NewObject<UTextureRenderTarget2D>(this);
    FoamOutputRT->InitCustomFormat(GridResolution, GridResolution, PF_R32_FLOAT, false);
    FoamOutputRT->AddressX = TA_Clamp;
    FoamOutputRT->AddressY = TA_Clamp;

    // 3. Setup Scene Capture Component properties
    if (TerrainCaptureComponent)
    {
        TerrainCaptureComponent->ProjectionType = ECameraProjectionMode::Orthographic;
        TerrainCaptureComponent->OrthoWidth = CapturedWorldWidth;
        TerrainCaptureComponent->TextureTarget = TerrainCaptureRT;
        TerrainCaptureComponent->CaptureSource = ESceneCaptureSource::SCS_DeviceDepth; // Grayscale depth map
        TerrainCaptureComponent->CaptureScene(); // Render the terrain depth once on start
    }

    // 4. Configure and Bind Render Targets to Simulation Component
    if (SimComponent)
    {
        SimComponent->GridSizeX = GridResolution;
        SimComponent->GridSizeY = GridResolution;
        SimComponent->CapturedWorldWidth = CapturedWorldWidth;
        SimComponent->bAutoCalculateCellSize = true;
        SimComponent->WaterLevel = WaterLevel;

        SimComponent->TerrainHeightInputRT = TerrainCaptureRT;
        SimComponent->HeightOutputRT = HeightOutputRT;
        SimComponent->DisplacementOutputRT = DisplacementOutputRT;
        SimComponent->FoamOutputRT = FoamOutputRT;
    }

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
        DynamicWaterMaterial->SetTextureParameterValue(FName("DisplacementMap"), DisplacementOutputRT);
        DynamicWaterMaterial->SetTextureParameterValue(FName("FoamMap"), FoamOutputRT);
        WaterMeshComponent->SetMaterial(0, DynamicWaterMaterial);
    }
}
