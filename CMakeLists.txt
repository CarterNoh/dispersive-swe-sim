cmake_minimum_required(VERSION 3.15)
project(DispersiveSWESimulator)

# Set C++ standard (C++17 or 20 is standard for modern simulation/Unreal)
set(CMAKE_CXX_STANDARD 17)

# Tell CMake where your header files are
include_directories(include)

# Define the executable and the files needed to build it
add_executable(Simulator 
    src/main.cpp 
    src/Sim2D.cpp
    src/gpu.cpp
)

# Link DirectX 11 and the Shader Compiler
target_link_libraries(Simulator d3d11 d3dcompiler)

# Copy HLSL files to the output directory automatically
add_custom_command(TARGET Simulator POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
    "${CMAKE_CURRENT_SOURCE_DIR}/src/shaders"
    "$<TARGET_FILE_DIR:Simulator>/shaders")