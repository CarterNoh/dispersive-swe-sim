# dispersive-swe-sim
2D fluid simulation based off of Jeschke &amp; Wojtan's "Generalizing Shallow Water Simulations with Dispersive Surface Waves": https://dl.acm.org/doi/10.1145/3592098



testing the fft sim: 
g++ -std=c++17 -O2 fft_validate.cpp fft_reference.cpp fft_gpu_cpu_emu.cpp -o fft_validate
