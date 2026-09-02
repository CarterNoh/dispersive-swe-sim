#include "SimModule.h"
#include "Interfaces/IPluginManager.h"
#include "ShaderCore.h"
#include "Misc/Paths.h"

void FSimModule::StartupModule()
{
	FString PluginShaderDir = FPaths::Combine(IPluginManager::Get().FindPlugin(TEXT("DispersiveSWESim"))->GetBaseDir(), TEXT("Shaders"));
	AddShaderSourceDirectoryMapping(TEXT("/Plugin/DispersiveSWESim"), PluginShaderDir);
}

void FSimModule::ShutdownModule()
{
}

IMPLEMENT_MODULE(FSimModule, DispersiveSWESim)
