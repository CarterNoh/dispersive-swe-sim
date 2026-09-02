#include "SimConfigSerializer.h"
#include "JsonObjectConverter.h"
#include "Misc/FileHelper.h"
#include "Misc/Paths.h"
#include "Serialization/JsonReader.h"
#include "Serialization/JsonSerializer.h"

bool FSimConfigSerializer::LoadParametersFromJson(UObject* TargetObject, const FString& FilePath, bool& bOutAutoCalculateCellSize) {
	if (!TargetObject) return false;

	FString FinalPath = FilePath;
	if (FPaths::IsRelative(FinalPath)) {
		FinalPath = FPaths::Combine(FPaths::ProjectDir(), FilePath);
	}

	FString JsonString;
	if (!FFileHelper::LoadFileToString(JsonString, *FinalPath)) {
		// Try resolving relative to Content folder as a fallback
		FString FallbackPath = FPaths::Combine(FPaths::ProjectContentDir(), FilePath);
		if (!FFileHelper::LoadFileToString(JsonString, *FallbackPath)) {
			UE_LOG(LogTemp, Warning, TEXT("Failed to load JSON file from path: %s (or fallback %s)"), *FinalPath, *FallbackPath);
			return false;
		}
		FinalPath = FallbackPath;
	}

	TSharedPtr<FJsonObject> JsonObject;
	TSharedRef<TJsonReader<>> Reader = TJsonReaderFactory<>::Create(JsonString);

	if (!FJsonSerializer::Deserialize(Reader, JsonObject) || !JsonObject.IsValid()) {
		UE_LOG(LogTemp, Warning, TEXT("Failed to deserialize JSON string."));
		return false;
	}

	// Map JSON fields directly to component properties
	if (!FJsonObjectConverter::JsonObjectToUStruct(JsonObject.ToSharedRef(), TargetObject->GetClass(), TargetObject)) {
		UE_LOG(LogTemp, Warning, TEXT("Failed to map JSON object to component properties."));
		return false;
	}

	if (JsonObject->HasField(TEXT("CellSize"))) {
		double CustomCellSize = 0.0;
		if (JsonObject->TryGetNumberField(TEXT("CellSize"), CustomCellSize) && CustomCellSize > 0.0) {
			bOutAutoCalculateCellSize = false;
		}
	}

	UE_LOG(LogTemp, Log, TEXT("Successfully loaded simulation parameters from: %s"), *FinalPath);
	return true;
}

bool FSimConfigSerializer::SaveParametersToJson(const UObject* TargetObject, const FString& FilePath) {
	if (!TargetObject) return false;

	FString FinalPath = FilePath;
	if (FPaths::IsRelative(FinalPath)) {
		FinalPath = FPaths::Combine(FPaths::ProjectDir(), FilePath);
	}

	TSharedRef<FJsonObject> JsonObject = MakeShared<FJsonObject>();
	if (!FJsonObjectConverter::UStructToJsonObject(TargetObject->GetClass(), TargetObject, JsonObject)) {
		UE_LOG(LogTemp, Warning, TEXT("Failed to convert component properties to JSON object."));
		return false;
	}

	// Remove runtime properties or input/output targets that shouldn't be serialized in a parameters config
	JsonObject->RemoveField(TEXT("terrainHeightInputRT"));
	JsonObject->RemoveField(TEXT("displacementRT"));
	JsonObject->RemoveField(TEXT("displacementPastRT"));
	JsonObject->RemoveField(TEXT("velocityRT"));
	JsonObject->RemoveField(TEXT("velocityPastRT"));
	JsonObject->RemoveField(TEXT("accelerationRT"));
	JsonObject->RemoveField(TEXT("accelerationPastRT"));
	JsonObject->RemoveField(TEXT("normalRT"));
	JsonObject->RemoveField(TEXT("foamRT"));
	JsonObject->RemoveField(TEXT("jacobianDetRT"));
	JsonObject->RemoveField(TEXT("roughnessRT"));
	JsonObject->RemoveField(TEXT("jsonConfigFilePath"));

	FString JsonString;
	TSharedRef<TJsonWriter<>> Writer = TJsonWriterFactory<>::Create(&JsonString);
	if (!FJsonSerializer::Serialize(JsonObject, Writer)) {
		UE_LOG(LogTemp, Warning, TEXT("Failed to serialize JSON object."));
		return false;
	}

	if (!FFileHelper::SaveStringToFile(JsonString, *FinalPath)) {
		UE_LOG(LogTemp, Warning, TEXT("Failed to save JSON string to path: %s"), *FinalPath);
		return false;
	}

	UE_LOG(LogTemp, Log, TEXT("Successfully saved simulation parameters to: %s"), *FinalPath);
	return true;
}
