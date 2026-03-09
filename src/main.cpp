#include <iostream>
#include <vector>
#include <chrono>
#include <string>
#include <fstream>
#include <filesystem>
#include "Sim2D.h"
namespace fs = std::filesystem;

// Parameters
int numTicks = 100; // Number of simulation steps to run
bool SWEonly = false; // Whether to run only the SWE step or include the full simulation steps
int terrainType = 0; // 0 = flat, 1 = ramp, 2 = bumps, 3 = basins, 4 = beach
int waterType = 1; // 0 = localized splash, 1 = step/dam break, 2 = basin flood
float waterLevel = 10.0f; // Initial water level

void SaveToCSV(const std::vector<float>& h, int size, int tick) {
    if (!fs::exists("data")) {
        fs::create_directory("data");
    }
    std::string filename = "data/frame_" + std::to_string(tick) + ".csv";
    if (tick < 0) {
        filename = "data/terrain.csv";
    }
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << " for writing!" << std::endl;
        return;
    }
    for (int y = 0; y < size; ++y) {
        for (int x = 0; x < size; ++x) {
            file << h[y * size + x] << (x == size - 1 ? "" : ",");
        }
        file << "\n";
    }
    file.close();
}

int main() {
    // Clear the data folder before starting the simulation
    if (fs::exists("data")) {  
        for (const auto& entry : fs::directory_iterator("data")) {
            fs::remove(entry.path());
        }
    }

    // Initialize the simulation
    Sim sim(terrainType, waterType, waterLevel);
    int size = sim.GRIDSIZE;
    SaveToCSV(sim.terrain, size, -1); // Save terrain data for reference
    SaveToCSV(sim.h, size, 0); // Save initial water height data

    sim.SimStep(SWEonly); // Run one step to initialize any necessary internal states before the main loop
    SaveToCSV(sim.hbar, size, 1); // Save bulk water height data
    SaveToCSV(sim.htilde, size, 2); // Save surface water height data

    // std::cout << "Starting simulation loop..." << std::endl;
    // auto totalStart = std::chrono::high_resolution_clock::now();
    // for (int tick = 1; tick < numTicks+1; ++tick) {
    //     auto start = std::chrono::high_resolution_clock::now();

    //     sim.SimStep(SWEonly);
    //     SaveToCSV(sim.h, size, tick);

    //     auto sim_time = std::chrono::high_resolution_clock::now();
    //     std::chrono::duration<double, std::milli> elapsed = sim_time - start;
    //     std::cout << "Tick: " << tick << " | Simulation Time: " << elapsed.count() << "ms" << std::endl;
    // }
    // auto totalEnd = std::chrono::high_resolution_clock::now();
    // std::chrono::duration<double, std::milli> totalElapsed = totalEnd - totalStart;
    // std::cout << "Total Simulation Time for " << numTicks << " ticks: " << totalElapsed.count() << "ms" << std::endl;
    // std::cout << "Average Speed: " << (numTicks / (totalElapsed.count() * 1000)) << " ticks/sec" << std::endl;

    return 0;
}