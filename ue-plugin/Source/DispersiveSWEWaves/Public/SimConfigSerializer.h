#pragma once

#include "CoreMinimal.h"

/**
 * Utility for serializing and deserializing simulation component parameters to and from JSON.
 */
class DISPERSIVESWEWAVES_API FSimConfigSerializer {
public:
	static bool LoadParametersFromJson(UObject* TargetObject, const FString& FilePath, bool& bOutAutoCalculateCellSize);
	static bool SaveParametersToJson(const UObject* TargetObject, const FString& FilePath);
};
