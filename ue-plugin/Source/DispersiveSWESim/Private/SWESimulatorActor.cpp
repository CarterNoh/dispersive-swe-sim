#include "SWESimulatorActor.h"
#include "Materials/MaterialInstanceDynamic.h"
#include "ProceduralMeshComponent.h"
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
    TerrainCaptureComponent->SetRelativeRotation(FRotator(-90.0f, 0.0f, -90.0f)); // Pointing down

    // 3. Attach water mesh representation (flat plane)
    WaterMeshComponent = CreateDefaultSubobject<UProceduralMeshComponent>(TEXT("WaterMesh"));
    WaterMeshComponent->SetupAttachment(Root);
    WaterMeshComponent->SetCollisionEnabled(ECollisionEnabled::NoCollision);

    // 4. Attach SWE Simulation orchestrator
    SimComponent = CreateDefaultSubobject<UDispersiveSWESimulator>(TEXT("SWESimulator"));

    // Set defaults
    GridResolution = 512;
    WaterLevel = 0.0f;
    CapturedWorldWidth = 51200.0f;
    bAutoFitToTerrain = true;
    bAutoLoadDefaultAssets = true;

    // Attempt to resolve the default water material
    static ConstructorHelpers::FObjectFinder<UMaterialInterface> WaterMaterialFinder(TEXT("/DispersiveSWESim/M_SWEDefaultWater"));
    if (WaterMaterialFinder.Succeeded())
    {
        BaseWaterMaterial = WaterMaterialFinder.Object;
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

        // Center horizontally and set the vertical level to WaterLevel
        FVector NewLoc = Origin;
        NewLoc.Z = WaterLevel;
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

    // 5. Generate procedural water mesh
    if (WaterMeshComponent)
    {
        WaterMeshComponent->SetWorldScale3D(FVector(1.0f, 1.0f, 1.0f));
        GenerateWaterGrid();
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
            WaterMeshComponent->SetWorldScale3D(FVector(1.0f, 1.0f, 1.0f));
        }
    }

    // Generate procedural water mesh at runtime
    if (WaterMeshComponent)
    {
        WaterMeshComponent->SetWorldScale3D(FVector(1.0f, 1.0f, 1.0f));
        GenerateWaterGrid();
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

void ASWESimulatorActor::GenerateWaterGrid()
{
    if (!WaterMeshComponent) return;

    int32 N = GridResolution;
    if (N <= 0) return;

    TArray<FVector> Vertices;
    TArray<int32> Triangles;
    TArray<FVector> Normals;
    TArray<FVector2D> UV0;
    TArray<FProcMeshTangent> Tangents;
    TArray<FLinearColor> VertexColors;

    Vertices.Reserve((N + 1) * (N + 1));
    UV0.Reserve((N + 1) * (N + 1));
    Normals.Reserve((N + 1) * (N + 1));
    Tangents.Reserve((N + 1) * (N + 1));

    float HalfWidth = CapturedWorldWidth * 0.5f;

    for (int32 y = 0; y <= N; ++y)
    {
        float V = (float)y / N;
        float YPos = -HalfWidth + V * CapturedWorldWidth;

        for (int32 x = 0; x <= N; ++x)
        {
            float U = (float)x / N;
            float XPos = -HalfWidth + U * CapturedWorldWidth;

            Vertices.Add(FVector(XPos, YPos, 0.0f));
            UV0.Add(FVector2D(U, V));
            Normals.Add(FVector(0.0f, 0.0f, 1.0f));
            Tangents.Add(FProcMeshTangent(1.0f, 0.0f, 0.0f));
        }
    }

    Triangles.Reserve(N * N * 6);
    for (int32 y = 0; y < N; ++y)
    {
        for (int32 x = 0; x < N; ++x)
        {
            int32 IndexBL = y * (N + 1) * 1 + x;
            int32 IndexBR = IndexBL + 1;
            int32 IndexTL = (y + 1) * (N + 1) * 1 + x;
            int32 IndexTR = IndexTL + 1;

            // Winding order: Clockwise
            Triangles.Add(IndexBL);
            Triangles.Add(IndexTL);
            Triangles.Add(IndexBR);

            Triangles.Add(IndexBR);
            Triangles.Add(IndexTL);
            Triangles.Add(IndexTR);
        }
    }

    WaterMeshComponent->ClearAllMeshSections();
    WaterMeshComponent->CreateMeshSection_LinearColor(
        0,
        Vertices,
        Triangles,
        Normals,
        UV0,
        VertexColors,
        Tangents,
        true // Create collision
    );
}
