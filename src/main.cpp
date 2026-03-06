#include <iostream>
#include <vector>
#include <chrono>
#include <thread>
#include "Sim2D.h"
#include <matplot/matplot.h>

// Parameters
bool visualize = true; // Set to false to disable visualization and just run the simulation loop
int numTicks = 100; // Number of simulation steps to run
int plotInterval = 5; // Update visualization every N ticks
bool SWEonly = true; // Whether to run only the SWE step or include the full simulation steps
int terrainType = 1; // 0 = flat, 1 = ramp, 2 = bumps, 3 =basins, 4 = beach
int waterType = 2; // 0 = localized splash, 1 = step/dam break, 2 = basin flood
float waterLevel = 10.0f; // Initial water level

int main() {
    // Initialize the simulation
    Sim sim(terrainType, waterType, waterLevel);
    int size = sim.GRIDSIZE;
    int delay = static_cast<int>(sim.TIMESTEP * 1000.f); // Delay between visualization updates in milliseconds

    // Initialize 2D grid
    auto [X, Y] = matplot::meshgrid(matplot::linspace(0, size - 1, size), matplot::linspace(0, size - 1, size));
    std::vector<std::vector<double>> Z(size, std::vector<double>(size));
    std::vector<std::vector<double>> T(size, std::vector<double>(size));
    for (int y = 0; y < size; ++y) {
        for (int x = 0; x < size; ++x) {
            T[y][x] = static_cast<double>(sim.terrain[y * size + x]);
            Z[y][x] = static_cast<double>(sim.h[y * size + x]);
        }
    }
    // Matplotplusplus setup for 3D surface
    auto f = matplot::figure(true);
    auto ax = f->current_axes();
    ax->zlim({-15, 25}); // Set z-axis limits for better visualization
    auto w = ax->surf(X, Y, Z); // Plot initial water surface
    // auto t = ax->surf(X, Y, T); // Plot terrain
    // t->face_alpha(0.5); // Make terrain semi-transparent


    std::cout << "Starting simulation loop..." << std::endl;

    auto totalStart = std::chrono::high_resolution_clock::now();
    for (int tick = 0; tick < numTicks; ++tick) {
        auto start = std::chrono::high_resolution_clock::now();

        sim.SimStep(SWEonly);

        // Map flattened array h to 2D grid Z for visualization
        for (int y = 0; y < size; ++y) {
            for (int x = 0; x < size; ++x) {
                Z[y][x] = static_cast<double>(sim.h[y * size + x]);
            }
        }

        auto sim_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> elapsed = sim_time - start;

        // Update visualization
        if (visualize && (tick % plotInterval == 0)) {
            w->Z_data(Z); // Update water surface
            // t->Z_data(T); // Update terrain
            ax->title("Tick: " + std::to_string(tick) + " | Frame Time: " + std::to_string(elapsed.count()) + "ms");
            // f->draw();
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