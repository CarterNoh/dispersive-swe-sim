#pragma once

#include "CoreMinimal.h"
#include "RenderTargetPool.h"
#include "RenderGraphBuilder.h"
#include "RHI.h"

/**
 * A pair of pooled render targets for ping-pong double buffering across simulation ticks.
 */
struct FPingPongBuffer {
	TRefCountPtr<IPooledRenderTarget> Current;
	TRefCountPtr<IPooledRenderTarget> Target;

	void Swap() {
		::Swap(Current, Target);
	}

	void Allocate(FRHICommandListImmediate& RHICmdList, const FPooledRenderTargetDesc& Desc, const TCHAR* NameCurrent, const TCHAR* NameTarget) {
		GRenderTargetPool.FindFreeElement(RHICmdList, Desc, Current, NameCurrent);
		if (NameTarget) {
			GRenderTargetPool.FindFreeElement(RHICmdList, Desc, Target, NameTarget);
		}
	}

	void Clear(FRDGBuilder& GraphBuilder, const FLinearColor& ClearColor = FLinearColor::Black) {
		if (Current.IsValid()) {
			AddClearUAVPass(GraphBuilder, GraphBuilder.CreateUAV(GraphBuilder.RegisterExternalTexture(Current)), ClearColor);
		}
		if (Target.IsValid()) {
			AddClearUAVPass(GraphBuilder, GraphBuilder.CreateUAV(GraphBuilder.RegisterExternalTexture(Target)), ClearColor);
		}
	}

	FRDGTextureRef RegisterCurrent(FRDGBuilder& GraphBuilder) const {
		return Current.IsValid() ? GraphBuilder.RegisterExternalTexture(Current) : nullptr;
	}

	FRDGTextureRef RegisterTarget(FRDGBuilder& GraphBuilder) const {
		return Target.IsValid() ? GraphBuilder.RegisterExternalTexture(Target) : nullptr;
	}
};

/**
 * Container holding all persistent simulation state render targets.
 */
struct FSimStateBuffers {
	// Diffusion ping-pong states
	FPingPongBuffer H;
	FPingPongBuffer Q_x;
	FPingPongBuffer Q_y;

	// Water depth & flux
	FPingPongBuffer h;
	TRefCountPtr<IPooledRenderTarget> q_x;
	TRefCountPtr<IPooledRenderTarget> q_y;

	// Low-pass filtered depth (SWE velocity computation)
	FPingPongBuffer hbar;

	// Wave dispersion second-order time integrator history
	FPingPongBuffer htildeOld;

	// Ocean directional spectrum (2D complex arrays)
	TRefCountPtr<IPooledRenderTarget> HPos;
	TRefCountPtr<IPooledRenderTarget> HNeg;

	// Terrain elevation
	TRefCountPtr<IPooledRenderTarget> Terrain;

	// Visual accumulation states
	FPingPongBuffer Foam;
	FPingPongBuffer Roughness;

	void AllocateAll(FRHICommandListImmediate& RHICmdList, int32 GridSizeX, int32 GridSizeY, int32 PaddedSizeX, int32 PaddedSizeY, int32 DepthNum) {
		FPooledRenderTargetDesc GridFloatDesc = FPooledRenderTargetDesc::Create2DDesc(
			FIntPoint(GridSizeX, GridSizeY),
			PF_R32_FLOAT,
			FClearValueBinding::None,
			TexCreate_None,
			TexCreate_ShaderResource | TexCreate_UAV,
			false
		);

		H.Allocate(RHICmdList, GridFloatDesc, TEXT("H"), TEXT("HPast"));
		Q_x.Allocate(RHICmdList, GridFloatDesc, TEXT("Q_x"), TEXT("QPast_x"));
		Q_y.Allocate(RHICmdList, GridFloatDesc, TEXT("Q_y"), TEXT("QPast_y"));
		h.Allocate(RHICmdList, GridFloatDesc, TEXT("h"), TEXT("hPast"));
		GRenderTargetPool.FindFreeElement(RHICmdList, GridFloatDesc, q_x, TEXT("q_x"));
		GRenderTargetPool.FindFreeElement(RHICmdList, GridFloatDesc, q_y, TEXT("q_y"));
		hbar.Allocate(RHICmdList, GridFloatDesc, TEXT("hbar"), TEXT("hbarOld"));
		htildeOld.Allocate(RHICmdList, GridFloatDesc, TEXT("htildeOld"), TEXT("htildeOldNext"));
		GRenderTargetPool.FindFreeElement(RHICmdList, GridFloatDesc, Terrain, TEXT("Terrain"));

		FPooledRenderTargetDesc FoamDesc = FPooledRenderTargetDesc::Create2DDesc(
			FIntPoint(GridSizeX, GridSizeY),
			PF_FloatRGBA,
			FClearValueBinding::None,
			TexCreate_None,
			TexCreate_ShaderResource | TexCreate_UAV,
			false
		);
		Foam.Allocate(RHICmdList, FoamDesc, TEXT("FoamState"), TEXT("NewFoamState"));

		FPooledRenderTargetDesc RoughnessDesc = FPooledRenderTargetDesc::Create2DDesc(
			FIntPoint(GridSizeX, 1),
			PF_FloatRGBA,
			FClearValueBinding::None,
			TexCreate_None,
			TexCreate_ShaderResource | TexCreate_UAV,
			false
		);
		Roughness.Allocate(RHICmdList, RoughnessDesc, TEXT("RoughnessState"), TEXT("NewRoughnessState"));

		FPooledRenderTargetDesc ComplexArrayDesc = FPooledRenderTargetDesc::Create2DArrayDesc(
			FIntPoint(PaddedSizeX, PaddedSizeY),
			PF_G32R32F,
			FClearValueBinding::None,
			TexCreate_None,
			TexCreate_ShaderResource | TexCreate_UAV,
			false,
			DepthNum
		);
		GRenderTargetPool.FindFreeElement(RHICmdList, ComplexArrayDesc, HPos, TEXT("HPos"));
		GRenderTargetPool.FindFreeElement(RHICmdList, ComplexArrayDesc, HNeg, TEXT("HNeg"));
	}

	void ClearAll(FRDGBuilder& GraphBuilder) {
		Q_x.Clear(GraphBuilder);
		Q_y.Clear(GraphBuilder);
		if (q_x.IsValid()) AddClearUAVPass(GraphBuilder, GraphBuilder.CreateUAV(GraphBuilder.RegisterExternalTexture(q_x)), FLinearColor::Black);
		if (q_y.IsValid()) AddClearUAVPass(GraphBuilder, GraphBuilder.CreateUAV(GraphBuilder.RegisterExternalTexture(q_y)), FLinearColor::Black);
		htildeOld.Clear(GraphBuilder);
		if (H.Target.IsValid()) AddClearUAVPass(GraphBuilder, GraphBuilder.CreateUAV(GraphBuilder.RegisterExternalTexture(H.Target)), FLinearColor::Black);
		if (h.Target.IsValid()) AddClearUAVPass(GraphBuilder, GraphBuilder.CreateUAV(GraphBuilder.RegisterExternalTexture(h.Target)), FLinearColor::Black);
		if (hbar.Target.IsValid()) AddClearUAVPass(GraphBuilder, GraphBuilder.CreateUAV(GraphBuilder.RegisterExternalTexture(hbar.Target)), FLinearColor::Black);
		Foam.Clear(GraphBuilder);
		Roughness.Clear(GraphBuilder);
	}
};
