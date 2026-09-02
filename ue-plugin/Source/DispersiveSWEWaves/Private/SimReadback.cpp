#include "SimReadback.h"
#include "RenderTargetPool.h"

FSimReadbackHandler::FSimReadbackHandler() {
	ActiveCPUBufferIndex.store(0);
	StagingWriteIndex = 0;
	StagingReadIndex = 1;
}

void FSimReadbackHandler::Allocate(int32 GridSizeX, int32 GridSizeY) {
	FRHITextureCreateDesc DispDesc0 = FRHITextureCreateDesc::Create2D(TEXT("Staging0"), GridSizeX, GridSizeY, PF_R32_FLOAT)
		.SetFlags(TexCreate_CPUReadback)
		.SetNumMips(1);
	StagingTextures[0] = RHICreateTexture(DispDesc0);

	FRHITextureCreateDesc DispDesc1 = FRHITextureCreateDesc::Create2D(TEXT("Staging1"), GridSizeX, GridSizeY, PF_R32_FLOAT)
		.SetFlags(TexCreate_CPUReadback)
		.SetNumMips(1);
	StagingTextures[1] = RHICreateTexture(DispDesc1);

	FRHITextureCreateDesc VelDesc0 = FRHITextureCreateDesc::Create2D(TEXT("VelocityStaging0"), GridSizeX, GridSizeY, PF_FloatRGBA)
		.SetFlags(TexCreate_CPUReadback)
		.SetNumMips(1);
	VelocityStagingTextures[0] = RHICreateTexture(VelDesc0);

	FRHITextureCreateDesc VelDesc1 = FRHITextureCreateDesc::Create2D(TEXT("VelocityStaging1"), GridSizeX, GridSizeY, PF_FloatRGBA)
		.SetFlags(TexCreate_CPUReadback)
		.SetNumMips(1);
	VelocityStagingTextures[1] = RHICreateTexture(VelDesc1);

	FRHITextureCreateDesc AccelDesc0 = FRHITextureCreateDesc::Create2D(TEXT("AccelerationStaging0"), GridSizeX, GridSizeY, PF_FloatRGBA)
		.SetFlags(TexCreate_CPUReadback)
		.SetNumMips(1);
	AccelerationStagingTextures[0] = RHICreateTexture(AccelDesc0);

	FRHITextureCreateDesc AccelDesc1 = FRHITextureCreateDesc::Create2D(TEXT("AccelerationStaging1"), GridSizeX, GridSizeY, PF_FloatRGBA)
		.SetFlags(TexCreate_CPUReadback)
		.SetNumMips(1);
	AccelerationStagingTextures[1] = RHICreateTexture(AccelDesc1);

	CPUHeightData[0].Reset();
	CPUHeightData[0].SetNumZeroed(GridSizeX * GridSizeY);
	CPUHeightData[1].Reset();
	CPUHeightData[1].SetNumZeroed(GridSizeX * GridSizeY);

	CPUVelocityData[0].Reset();
	CPUVelocityData[0].SetNumZeroed(GridSizeX * GridSizeY);
	CPUVelocityData[1].Reset();
	CPUVelocityData[1].SetNumZeroed(GridSizeX * GridSizeY);

	CPUAccelerationData[0].Reset();
	CPUAccelerationData[0].SetNumZeroed(GridSizeX * GridSizeY);
	CPUAccelerationData[1].Reset();
	CPUAccelerationData[1].SetNumZeroed(GridSizeX * GridSizeY);

	ActiveCPUBufferIndex.store(0);
	StagingWriteIndex = 0;
	StagingReadIndex = 1;
	FramesSinceAlloc = 0;
}

void FSimReadbackHandler::EnqueueReadback(
	FRDGBuilder& GraphBuilder,
	FRDGTextureRef HeightRDG,
	FRDGTextureRef VelocityRDG,
	FRDGTextureRef AccelerationRDG,
	int32 GridSizeX,
	int32 GridSizeY) {
	if (StagingTextures[StagingWriteIndex].IsValid() && HeightRDG) {
		FTextureRHIRef CurrentStagingTexture = StagingTextures[StagingWriteIndex];
		TRefCountPtr<IPooledRenderTarget> StagingPooled = CreateRenderTarget(CurrentStagingTexture, TEXT("StagingTexture"));
		FRDGTextureRef StagingRDG = GraphBuilder.RegisterExternalTexture(StagingPooled);
		AddCopyTexturePass(GraphBuilder, HeightRDG, StagingRDG);
	}

	if (VelocityStagingTextures[StagingWriteIndex].IsValid() && VelocityRDG) {
		FTextureRHIRef CurrentVelStaging = VelocityStagingTextures[StagingWriteIndex];
		TRefCountPtr<IPooledRenderTarget> VelStagingPooled = CreateRenderTarget(CurrentVelStaging, TEXT("VelStagingTexture"));
		FRDGTextureRef VelStagingRDG = GraphBuilder.RegisterExternalTexture(VelStagingPooled);
		AddCopyTexturePass(GraphBuilder, VelocityRDG, VelStagingRDG);
	}

	if (AccelerationStagingTextures[StagingWriteIndex].IsValid() && AccelerationRDG) {
		FTextureRHIRef CurrentAccelStaging = AccelerationStagingTextures[StagingWriteIndex];
		TRefCountPtr<IPooledRenderTarget> AccelStagingPooled = CreateRenderTarget(CurrentAccelStaging, TEXT("AccelStagingTexture"));
		FRDGTextureRef AccelStagingRDG = GraphBuilder.RegisterExternalTexture(AccelStagingPooled);
		AddCopyTexturePass(GraphBuilder, AccelerationRDG, AccelStagingRDG);
	}

	if (FramesSinceAlloc > 0 && (StagingTextures[StagingReadIndex].IsValid() || VelocityStagingTextures[StagingReadIndex].IsValid() || AccelerationStagingTextures[StagingReadIndex].IsValid())) {
		FTextureRHIRef ReadStagingTexture = StagingTextures[StagingReadIndex];
		FTextureRHIRef ReadVelStagingTexture = VelocityStagingTextures[StagingReadIndex];
		FTextureRHIRef ReadAccelStagingTexture = AccelerationStagingTextures[StagingReadIndex];
		int32 TargetCPUIndex = 1 - ActiveCPUBufferIndex.load();

		GraphBuilder.AddPass(
			RDG_EVENT_NAME("SWE_ReadbackMap"),
			ERDGPassFlags::None,
			[ReadStagingTexture, ReadVelStagingTexture, ReadAccelStagingTexture, TargetCPUIndex, GridSizeX, GridSizeY, this](FRHICommandListImmediate& RHICmdList) {
				if (ReadStagingTexture.IsValid()) {
					void* LocalData = nullptr;
					int32 OutWidth = 0;
					int32 OutHeight = 0;
					RHICmdList.MapStagingSurface(ReadStagingTexture, LocalData, OutWidth, OutHeight);
					if (LocalData) {
						int32 Size = GridSizeX * GridSizeY;
						TArray<float>& DestArray = CPUHeightData[TargetCPUIndex];
						DestArray.SetNumUninitialized(Size);
						
						uint8* LocalByteData = (uint8*)LocalData;
						for (int32 y = 0; y < GridSizeY; ++y) {
							float* SrcRow = (float*)(LocalByteData + y * OutWidth);
							float* DestRow = DestArray.GetData() + y * GridSizeX;
							FMemory::Memcpy(DestRow, SrcRow, GridSizeX * sizeof(float));
						}
						
						RHICmdList.UnmapStagingSurface(ReadStagingTexture);
					}
				}

				if (ReadVelStagingTexture.IsValid()) {
					void* VelData = nullptr;
					int32 VelOutWidth = 0;
					int32 VelOutHeight = 0;
					RHICmdList.MapStagingSurface(ReadVelStagingTexture, VelData, VelOutWidth, VelOutHeight);
					if (VelData) {
						int32 Size = GridSizeX * GridSizeY;
						TArray<FVector>& DestVelArray = CPUVelocityData[TargetCPUIndex];
						DestVelArray.SetNumUninitialized(Size);

						uint8* VelByteData = (uint8*)VelData;
						for (int32 y = 0; y < GridSizeY; ++y) {
							FFloat16Color* SrcRow = (FFloat16Color*)(VelByteData + y * VelOutWidth);
							for (int32 x = 0; x < GridSizeX; ++x) {
								DestVelArray[y * GridSizeX + x] = FVector(
									SrcRow[x].R.GetFloat(),
									SrcRow[x].G.GetFloat(),
									SrcRow[x].B.GetFloat()
								);
							}
						}
						RHICmdList.UnmapStagingSurface(ReadVelStagingTexture);
					}
				}

				if (ReadAccelStagingTexture.IsValid()) {
					void* AccelData = nullptr;
					int32 AccelOutWidth = 0;
					int32 AccelOutHeight = 0;
					RHICmdList.MapStagingSurface(ReadAccelStagingTexture, AccelData, AccelOutWidth, AccelOutHeight);
					if (AccelData) {
						int32 Size = GridSizeX * GridSizeY;
						TArray<FVector>& DestAccelArray = CPUAccelerationData[TargetCPUIndex];
						DestAccelArray.SetNumUninitialized(Size);

						uint8* AccelByteData = (uint8*)AccelData;
						for (int32 y = 0; y < GridSizeY; ++y) {
							FFloat16Color* SrcRow = (FFloat16Color*)(AccelByteData + y * AccelOutWidth);
							for (int32 x = 0; x < GridSizeX; ++x) {
								DestAccelArray[y * GridSizeX + x] = FVector(
									SrcRow[x].R.GetFloat(),
									SrcRow[x].G.GetFloat(),
									SrcRow[x].B.GetFloat()
								);
							}
						}
						RHICmdList.UnmapStagingSurface(ReadAccelStagingTexture);
					}
				}

				ActiveCPUBufferIndex.store(TargetCPUIndex);
			}
		);
	}

	StagingWriteIndex = (StagingWriteIndex + 1) % 2;
	StagingReadIndex = (StagingReadIndex + 1) % 2;
	FramesSinceAlloc++;
}

float FSimReadbackHandler::GetCachedHeight(int32 X, int32 Y, int32 GridSizeX, float DefaultWaterLevel) const {
	int32 Index = Y * GridSizeX + X;
	int32 ActiveIdx = ActiveCPUBufferIndex.load();
	if (CPUHeightData[ActiveIdx].IsValidIndex(Index)) {
		return CPUHeightData[ActiveIdx][Index];
	}
	return DefaultWaterLevel * 0.01f;
}

FVector FSimReadbackHandler::GetCachedVelocity(int32 X, int32 Y, int32 GridSizeX) const {
	int32 Index = Y * GridSizeX + X;
	int32 ActiveIdx = ActiveCPUBufferIndex.load();
	if (CPUVelocityData[ActiveIdx].IsValidIndex(Index)) {
		return CPUVelocityData[ActiveIdx][Index];
	}
	return FVector::ZeroVector;
}

FVector FSimReadbackHandler::GetCachedAcceleration(int32 X, int32 Y, int32 GridSizeX) const {
	int32 Index = Y * GridSizeX + X;
	int32 ActiveIdx = ActiveCPUBufferIndex.load();
	if (CPUAccelerationData[ActiveIdx].IsValidIndex(Index)) {
		return CPUAccelerationData[ActiveIdx][Index];
	}
	return FVector::ZeroVector;
}

float FSimReadbackHandler::GetWaterHeightAtLocation(
	const FVector& WorldLocation,
	const FVector& ActorLocation,
	float CapturedWorldWidth,
	int32 GridSizeX,
	int32 GridSizeY,
	float DefaultWaterLevel) const {
	if (GridSizeX <= 1 || GridSizeY <= 1 || CapturedWorldWidth <= 0.0f) {
		return DefaultWaterLevel;
	}

	FVector LocalPos = WorldLocation - ActorLocation;
	float U = FMath::Clamp((LocalPos.X / CapturedWorldWidth) + 0.5f, 0.0f, 1.0f);
	float V = FMath::Clamp((LocalPos.Y / CapturedWorldWidth) + 0.5f, 0.0f, 1.0f);

	float GridX = U * (GridSizeX - 1);
	float GridY = V * (GridSizeY - 1);

	int32 X0 = FMath::FloorToInt(GridX);
	int32 Y0 = FMath::FloorToInt(GridY);
	int32 X1 = FMath::Min(X0 + 1, GridSizeX - 1);
	int32 Y1 = FMath::Min(Y0 + 1, GridSizeY - 1);

	float LerpX = GridX - X0;
	float LerpY = GridY - Y0;

	float H00 = GetCachedHeight(X0, Y0, GridSizeX, DefaultWaterLevel);
	float H10 = GetCachedHeight(X1, Y0, GridSizeX, DefaultWaterLevel);
	float H01 = GetCachedHeight(X0, Y1, GridSizeX, DefaultWaterLevel);
	float H11 = GetCachedHeight(X1, Y1, GridSizeX, DefaultWaterLevel);

	float H_Avg = FMath::BiLerp(H00, H10, H01, H11, LerpX, LerpY);
	return H_Avg * 100.0f; // Return in cm
}

FVector FSimReadbackHandler::GetWaterVelocityAtLocation(
	const FVector& WorldLocation,
	const FVector& ActorLocation,
	float CapturedWorldWidth,
	int32 GridSizeX,
	int32 GridSizeY) const {
	if (GridSizeX <= 1 || GridSizeY <= 1 || CapturedWorldWidth <= 0.0f) {
		return FVector::ZeroVector;
	}

	FVector LocalPos = WorldLocation - ActorLocation;
	float U = FMath::Clamp((LocalPos.X / CapturedWorldWidth) + 0.5f, 0.0f, 1.0f);
	float V = FMath::Clamp((LocalPos.Y / CapturedWorldWidth) + 0.5f, 0.0f, 1.0f);

	float GridX = U * (GridSizeX - 1);
	float GridY = V * (GridSizeY - 1);

	int32 X0 = FMath::FloorToInt(GridX);
	int32 Y0 = FMath::FloorToInt(GridY);
	int32 X1 = FMath::Min(X0 + 1, GridSizeX - 1);
	int32 Y1 = FMath::Min(Y0 + 1, GridSizeY - 1);

	float LerpX = GridX - X0;
	float LerpY = GridY - Y0;

	FVector V00 = GetCachedVelocity(X0, Y0, GridSizeX);
	FVector V10 = GetCachedVelocity(X1, Y0, GridSizeX);
	FVector V01 = GetCachedVelocity(X0, Y1, GridSizeX);
	FVector V11 = GetCachedVelocity(X1, Y1, GridSizeX);

	return FMath::BiLerp(V00, V10, V01, V11, LerpX, LerpY); // In m/s
}

FVector FSimReadbackHandler::GetWaterAccelerationAtLocation(
	const FVector& WorldLocation,
	const FVector& ActorLocation,
	float CapturedWorldWidth,
	int32 GridSizeX,
	int32 GridSizeY) const {
	if (GridSizeX <= 1 || GridSizeY <= 1 || CapturedWorldWidth <= 0.0f) {
		return FVector::ZeroVector;
	}

	FVector LocalPos = WorldLocation - ActorLocation;
	float U = FMath::Clamp((LocalPos.X / CapturedWorldWidth) + 0.5f, 0.0f, 1.0f);
	float V = FMath::Clamp((LocalPos.Y / CapturedWorldWidth) + 0.5f, 0.0f, 1.0f);

	float GridX = U * (GridSizeX - 1);
	float GridY = V * (GridSizeY - 1);

	int32 X0 = FMath::FloorToInt(GridX);
	int32 Y0 = FMath::FloorToInt(GridY);
	int32 X1 = FMath::Min(X0 + 1, GridSizeX - 1);
	int32 Y1 = FMath::Min(Y0 + 1, GridSizeY - 1);

	float LerpX = GridX - X0;
	float LerpY = GridY - Y0;

	FVector A00 = GetCachedAcceleration(X0, Y0, GridSizeX);
	FVector A10 = GetCachedAcceleration(X1, Y0, GridSizeX);
	FVector A01 = GetCachedAcceleration(X0, Y1, GridSizeX);
	FVector A11 = GetCachedAcceleration(X1, Y1, GridSizeX);

	return FMath::BiLerp(A00, A10, A01, A11, LerpX, LerpY); // In m/s^2
}
