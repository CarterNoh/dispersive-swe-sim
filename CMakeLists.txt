cmake_minimum_required(VERSION 3.15)
project(DispersiveSWESimulator)

# Set C++ standard (C++17 or 20 is standard for modern simulation/Unreal)
set(CMAKE_CXX_STANDARD 17)

# Automatically download Matplotplusplus
include(FetchContent)
FetchContent_Declare(
    matplotplusplus
    GIT_REPOSITORY https://github.com/alandefreitas/matplotplusplus
    GIT_TAG origin/master
)
FetchContent_MakeAvailable(matplotplusplus)

# Tell CMake where your header files are
include_directories(include)

# Define the executable and the files needed to build it
add_executable(Simulator 
    src/main.cpp 
    src/Sim2D.cpp
)

# "Link" the plotting library to your executable
target_link_libraries(Simulator PRIVATE matplot)