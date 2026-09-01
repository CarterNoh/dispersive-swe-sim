// Top5Shaders.h
#pragma once

#include "CoreMinimal.h"
#include "GlobalShader.h"
#include "ShaderParameterStruct.h"
#include "RenderGraphResources.h"
#include "RenderGraphUtils.h"

// ------------------------------------------------------------------
// Pass 1: per-threadgroup top-5 reduction over the full texture
// ------------------------------------------------------------------
class FGroupTop5CS : public FGlobalShader
{
    DECLARE_GLOBAL_SHADER(FGroupTop5CS);
    SHADER_USE_PARAMETER_STRUCT(FGroupTop5CS, FGlobalShader);

    BEGIN_SHADER_PARAMETER_STRUCT(FParameters, )
        SHADER_PARAMETER_RDG_TEXTURE(Texture2D<float2>, InputTexture)
        SHADER_PARAMETER_RDG_BUFFER_UAV(RWStructuredBuffer<float>, CandidateVal)
        SHADER_PARAMETER_RDG_BUFFER_UAV(RWStructuredBuffer<uint2>, CandidateIdx)
        SHADER_PARAMETER(FIntPoint, TextureSize)
        SHADER_PARAMETER(uint32, GroupsPerRow)
    END_SHADER_PARAMETER_STRUCT()

    static const int32 ThreadGroupSizeX = 16;
    static const int32 ThreadGroupSizeY = 16;

    static bool ShouldCompilePermutation(const FGlobalShaderPermutationParameters& Parameters)
    {
        return true;
    }

    static void ModifyCompilationEnvironment(const FGlobalShaderPermutationParameters& Parameters, FShaderCompilerEnvironment& OutEnvironment)
    {
        FGlobalShader::ModifyCompilationEnvironment(Parameters, OutEnvironment);
        OutEnvironment.SetDefine(TEXT("GROUP_SIZE_X"), ThreadGroupSizeX);
        OutEnvironment.SetDefine(TEXT("GROUP_SIZE_Y"), ThreadGroupSizeY);
    }
};

// ------------------------------------------------------------------
// Pass 2: reduce all per-group candidates to the final top-5
// (single threadgroup of 256 threads, no CPU readback)
// ------------------------------------------------------------------
class FFinalTop5CS : public FGlobalShader
{
    DECLARE_GLOBAL_SHADER(FFinalTop5CS);
    SHADER_USE_PARAMETER_STRUCT(FFinalTop5CS, FGlobalShader);

    BEGIN_SHADER_PARAMETER_STRUCT(FParameters, )
        SHADER_PARAMETER_RDG_BUFFER_SRV(StructuredBuffer<float>, CandidateValIn)
        SHADER_PARAMETER_RDG_BUFFER_SRV(StructuredBuffer<uint2>, CandidateIdxIn)
        SHADER_PARAMETER(uint32, NumCandidates)
        SHADER_PARAMETER_RDG_BUFFER_UAV(RWStructuredBuffer<float>, FinalVal)
        SHADER_PARAMETER_RDG_BUFFER_UAV(RWStructuredBuffer<uint2>, FinalIdx)
    END_SHADER_PARAMETER_STRUCT()

    static bool ShouldCompilePermutation(const FGlobalShaderPermutationParameters& Parameters)
    {
        return true;
    }
};

// ------------------------------------------------------------------
// Pass 3: scatter-write the 5 winning texels into their own array slices
// ------------------------------------------------------------------
class FIsolateTop5CS : public FGlobalShader
{
    DECLARE_GLOBAL_SHADER(FIsolateTop5CS);
    SHADER_USE_PARAMETER_STRUCT(FIsolateTop5CS, FGlobalShader);

    BEGIN_SHADER_PARAMETER_STRUCT(FParameters, )
        SHADER_PARAMETER_RDG_TEXTURE(Texture2D<float2>, IsolateInputTexture)
        SHADER_PARAMETER_RDG_BUFFER_SRV(StructuredBuffer<uint2>, IsolateFinalIdx)
        SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2DArray<float2>, OutputArray)
    END_SHADER_PARAMETER_STRUCT()

    static bool ShouldCompilePermutation(const FGlobalShaderPermutationParameters& Parameters)
    {
        return true;
    }
};

// ------------------------------------------------------------------
// Entry point called from simulation code to add all 3 passes
// (plus the array clear) to the render graph for one frame.
// ------------------------------------------------------------------
struct FTop5Result
{
    FRDGBufferRef FinalVal;   // 5 floats  (squared magnitude, descending)
    FRDGBufferRef FinalIdx;   // 5 uint2s  (matching texel coordinates)
    FRDGTextureRef IsolatedArray; // 5-slice Texture2DArray, one isolated value per slice
};

FTop5Result AddIsolateTop5FrequenciesPass(
    FRDGBuilder& GraphBuilder,
    FRDGTextureRef InputSpectrum,
    int32 Width,
    int32 Height);