#include "DispersiveSWESimModule.h"
#include "Interfaces/IPluginManager.h"
#include "ShaderCore.h"
#include "Misc/Paths.h"

void FDispersiveSWESimModule::StartupModule()
{
	FString PluginShaderDir = FPaths::Combine(IPluginManager::Get().FindPlugin(TEXT("DispersiveSWESim"))->GetBaseDir(), TEXT("Shaders"));
	AddShaderSourceDirectoryMapping(TEXT("/Plugin/DispersiveSWESim"), PluginShaderDir);
}

void FDispersiveSWESimModule::ShutdownModule()
{
}

IMPLEMENT_MODULE(FDispersiveSWESimModule, DispersiveSWESim)
