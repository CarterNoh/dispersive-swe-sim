#pragma once

#include "CoreMinimal.h"
#include "RHI.h"
#include "RHIResources.h"
#include "RenderGraphBuilder.h"
#include "RenderGraphUtils.h"
#include <atomic>

/**
 * Handles double-buffered async GPU readback of fluid height, velocity, and acceleration
 * and thread-safe CPU bilinear spatial queries.
 */
class DISPERSIVESWEWAVES_API FSimReadbackHandler {
public:
	FSimReadbackHandler();

	void Allocate(int32 GridSizeX, int32 GridSizeY);

	void EnqueueReadback(
		FRDGBuilder& GraphBuilder,
		FRDGTextureRef HeightRDG,
		FRDGTextureRef VelocityRDG,
		FRDGTextureRef AccelerationRDG,
		int32 GridSizeX,
		int32 GridSizeY
	);

	bool HasVelocityStaging() const { return VelocityStagingTextures[0].IsValid(); }
	bool HasAccelerationStaging() const { return AccelerationStagingTextures[0].IsValid(); }

	float GetCachedHeight(int32 X, int32 Y, int32 GridSizeX, float DefaultWaterLevel) const;
	FVector GetCachedVelocity(int32 X, int32 Y, int32 GridSizeX) const;
	FVector GetCachedAcceleration(int32 X, int32 Y, int32 GridSizeX) const;

	float GetWaterHeightAtLocation(
		const FVector& WorldLocation,
		const FVector& ActorLocation,
		float CapturedWorldWidth,
		int32 GridSizeX,
		int32 GridSizeY,
		float DefaultWaterLevel
	) const;

	FVector GetWaterVelocityAtLocation(
		const FVector& WorldLocation,
		const FVector& ActorLocation,
		float CapturedWorldWidth,
		int32 GridSizeX,
		int32 GridSizeY
	) const;

	FVector GetWaterAccelerationAtLocation(
		const FVector& WorldLocation,
		const FVector& ActorLocation,
		float CapturedWorldWidth,
		int32 GridSizeX,
		int32 GridSizeY
	) const;

private:
	FTextureRHIRef StagingTextures[2];
	FTextureRHIRef VelocityStagingTextures[2];
	FTextureRHIRef AccelerationStagingTextures[2];
	int32 StagingWriteIndex = 0;
	int32 StagingReadIndex = 1;
	int32 FramesSinceAlloc = 0;

	TArray<float> CPUHeightData[2];
	TArray<FVector> CPUVelocityData[2];
	TArray<FVector> CPUAccelerationData[2];
	mutable std::atomic<int32> ActiveCPUBufferIndex{0};
};
