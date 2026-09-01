// Top5Shaders.cpp
#include "Top5Shaders.h"
#include "RenderGraphBuilder.h"
#include "RenderTargetPool.h"

IMPLEMENT_GLOBAL_SHADER(FGroupTop5CS,   "/Plugin/DispersiveSWESim/top5.usf", "CSGroupTop5", SF_Compute);
IMPLEMENT_GLOBAL_SHADER(FFinalTop5CS,   "/Plugin/DispersiveSWESim/top5.usf", "CSFinalTop5", SF_Compute);
IMPLEMENT_GLOBAL_SHADER(FIsolateTop5CS, "/Plugin/DispersiveSWESim/top5.usf", "CSIsolate",   SF_Compute);

FTop5Result AddIsolateTop5FrequenciesPass(
    FRDGBuilder& GraphBuilder,
    FRDGTextureRef InputSpectrum,
    int32 Width,
    int32 Height)
{
    const int32 GroupSizeX = FGroupTop5CS::ThreadGroupSizeX;
    const int32 GroupSizeY = FGroupTop5CS::ThreadGroupSizeY;
    const int32 GroupsX = FMath::DivideAndRoundUp(Width, GroupSizeX);
    const int32 GroupsY = FMath::DivideAndRoundUp(Height, GroupSizeY);
    const int32 NumGroups = GroupsX * GroupsY;
    const int32 NumCandidates = NumGroups * 5;

    check(NumCandidates > 0);

    // --- Intermediate buffers ---
    FRDGBufferRef CandidateVal = GraphBuilder.CreateBuffer(
        FRDGBufferDesc::CreateStructuredDesc(sizeof(float), NumCandidates),
        TEXT("Top5.CandidateVal"));
    FRDGBufferRef CandidateIdx = GraphBuilder.CreateBuffer(
        FRDGBufferDesc::CreateStructuredDesc(sizeof(FUintVector2), NumCandidates),
        TEXT("Top5.CandidateIdx"));
    FRDGBufferRef FinalVal = GraphBuilder.CreateBuffer(
        FRDGBufferDesc::CreateStructuredDesc(sizeof(float), 5),
        TEXT("Top5.FinalVal"));
    FRDGBufferRef FinalIdx = GraphBuilder.CreateBuffer(
        FRDGBufferDesc::CreateStructuredDesc(sizeof(FUintVector2), 5),
        TEXT("Top5.FinalIdx"));

    // --- Pass 1: per-group top-5 ---
    {
        FGroupTop5CS::FParameters* Params = GraphBuilder.AllocParameters<FGroupTop5CS::FParameters>();
        Params->InputTexture = InputSpectrum;
        Params->CandidateVal = GraphBuilder.CreateUAV(CandidateVal, PF_R32_FLOAT);
        Params->CandidateIdx = GraphBuilder.CreateUAV(CandidateIdx, PF_R32G32_UINT);
        Params->TextureSize = FIntPoint(Width, Height);
        Params->GroupsPerRow = GroupsX;

        TShaderMapRef<FGroupTop5CS> Shader(GetGlobalShaderMap(GMaxRHIFeatureLevel));
        FComputeShaderUtils::AddPass(
            GraphBuilder,
            RDG_EVENT_NAME("Top5.GroupReduce (%dx%d)", Width, Height),
            ERDGPassFlags::Compute,
            Shader,
            Params,
            FIntVector(GroupsX, GroupsY, 1));
    }

    // --- Pass 2: final reduction, single threadgroup of 256 ---
    {
        FFinalTop5CS::FParameters* Params = GraphBuilder.AllocParameters<FFinalTop5CS::FParameters>();
        Params->CandidateValIn = GraphBuilder.CreateSRV(CandidateVal, PF_R32_FLOAT);
        Params->CandidateIdxIn = GraphBuilder.CreateSRV(CandidateIdx, PF_R32G32_UINT);
        Params->NumCandidates = static_cast<uint32>(NumCandidates);
        Params->FinalVal = GraphBuilder.CreateUAV(FinalVal, PF_R32_FLOAT);
        Params->FinalIdx = GraphBuilder.CreateUAV(FinalIdx, PF_R32G32_UINT);

        TShaderMapRef<FFinalTop5CS> Shader(GetGlobalShaderMap(GMaxRHIFeatureLevel));
        FComputeShaderUtils::AddPass(
            GraphBuilder,
            RDG_EVENT_NAME("Top5.FinalReduce (%d candidates)", NumCandidates),
            ERDGPassFlags::Compute,
            Shader,
            Params,
            FIntVector(1, 1, 1)); // exactly one threadgroup of 256 threads
    }

    // --- Pass 3: clear + isolate ---
    FRDGTextureDesc ArrayDesc = FRDGTextureDesc::Create2DArray(
        FIntPoint(Width, Height),
        PF_G16R16F,
        FClearValueBinding::Black,
        TexCreate_ShaderResource | TexCreate_UAV,
        /*ArraySize=*/ 5);

    FRDGTextureRef IsolatedArray = GraphBuilder.CreateTexture(ArrayDesc, TEXT("Top5.IsolatedArray"));

    // Fast native clear — must happen before the scatter-write pass below.
    AddClearUAVPass(GraphBuilder, GraphBuilder.CreateUAV(IsolatedArray), FLinearColor::Black);

    {
        FIsolateTop5CS::FParameters* Params = GraphBuilder.AllocParameters<FIsolateTop5CS::FParameters>();
        Params->IsolateInputTexture = InputSpectrum;
        Params->IsolateFinalIdx = GraphBuilder.CreateSRV(FinalIdx, PF_R32G32_UINT);
        Params->OutputArray = GraphBuilder.CreateUAV(IsolatedArray);

        TShaderMapRef<FIsolateTop5CS> Shader(GetGlobalShaderMap(GMaxRHIFeatureLevel));
        FComputeShaderUtils::AddPass(
            GraphBuilder,
            RDG_EVENT_NAME("Top5.Isolate"),
            ERDGPassFlags::Compute,
            Shader,
            Params,
            FIntVector(1, 1, 1)); // single group of 5 threads is plenty
    }

    FTop5Result Result;
    Result.FinalVal = FinalVal;
    Result.FinalIdx = FinalIdx;
    Result.IsolatedArray = IsolatedArray;
    return Result;
}