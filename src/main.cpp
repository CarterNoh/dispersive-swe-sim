#include <iostream>
#include <vector>
#include <chrono>
#include "Sim2D.h"
#include <matplot/matplot.h>

int main() {

    // Parameters
    bool visualize = true; // Set to false to disable visualization and just run the simulation loop
    int numTicks = 100; // Number of simulation steps to run
    int plotInterval = 1; // Update visualization every N ticks
    bool SWEonly = true; // Whether to run only the SWE step or include the full simulation steps
    int terrainType = 1; // 0 = flat, 1 = ramp, 2 = bumps, 3 =basins, 4 = beach
    int waterType = 2; // 0 = localized splash, 1 = step/dam break, 2 = basin flood
    float waterLevel = 5.0f; // Initial water level

    // Initialize the simulation
    Sim sim;

    // Matplotplusplus setup for 3D surface
    auto f = matplot::figure(true);
    auto ax = f->current_axes();

    // Create a 2D grid (vector of vectors) that the plotter understands
    std::vector<std::vector<double>> Z(GRIDSIZE, std::vector<double>(GRIDSIZE));
    std::vector<std::vector<double>> T(GRIDSIZE, std::vector<double>(GRIDSIZE));
    // Map flattened terrain to 2D grid T for visualization
        for (int y = 0; y < GRIDSIZE; ++y) {
            for (int x = 0; x < GRIDSIZE; ++x) {
                T[y][x] = static_cast<double>(sim.terrain[y * GRIDSIZE + x]);
            }
        }

    std::cout << "Starting simulation loop..." << std::endl;

    auto totalStart = std::chrono::high_resolution_clock::now();
    for (int tick = 0; tick < numTicks; ++tick) {
        auto start = std::chrono::high_resolution_clock::now();

        sim.SimStep(SWEonly);

        // Map flattened array h to 2D grid Z for visualization
        for (int y = 0; y < GRIDSIZE; ++y) {
            for (int x = 0; x < GRIDSIZE; ++x) {
                Z[y][x] = static_cast<double>(sim.h[y * GRIDSIZE + x]);
            }
        }

        auto sim_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> elapsed = sim_time - start;

        // Update visualization
        if (visualize && (tick % plotInterval == 0)) {

            ax->clear();
            ax->surf(T); // Plot terrain with some transparency
            ax->surf(Z, 0.5); // Overlay water surface with some transparency
            ax->zlim({-15, 25}); 
            ax->title("Tick: " + std::to_string(tick) + " | Frame Time: " + std::to_string(elapsed.count()) + "ms");
            f->draw();
        }
        auto end = std::chrono::high_resolution_clock::now();
        auto plot_time = std::chrono::duration<double, std::milli>(end - sim_time);
        auto total_time = std::chrono::duration<double, std::milli>(end - start);
        std::cout << "Tick: " << tick << " | Simulation Time: " << elapsed.count() << "ms | Plot Time: " << plot_time.count() << "ms | Total Time: " << total_time.count() << "ms" << std::endl;
    }
    auto totalEnd = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> totalElapsed = totalEnd - totalStart;
    std::cout << "Total Simulation Time for " << numTicks << " ticks: " << totalElapsed.count() << "ms" << std::endl;
    std::cout << "Average Speed: " << (numTicks / (totalElapsed.count() * 1000)) << " ticks/sec" << std::endl;

    if (visualize) {
        matplot::show();
    }
    return 0;
}