cmake_minimum_required(VERSION 3.15)
project(DispersiveSWESimulator)

# Set C++ standard (C++17 or 20 is standard for modern simulation/Unreal)
set(CMAKE_CXX_STANDARD 17)

# Tell CMake where your header files are
include_directories(include)

# Add DirectX SDK paths
include_directories("C:/Program Files (x86)/Windows Kits/10/Include/10.0.22621.0/um")
include_directories("C:/Program Files (x86)/Windows Kits/10/Include/10.0.22621.0/shared")

# Link DirectX 12 libraries
link_directories("C:/Program Files (x86)/Windows Kits/10/Lib/10.0.22621.0/um/x64")

# Function to compile HLSL shaders
function(add_hlsl_shader SHADER_FILE SHADER_ENTRY SHADER_TARGET)
    get_filename_component(SHADER_NAME ${SHADER_FILE} NAME_WE)
    set(COMPILED_SHADER "${CMAKE_CURRENT_BINARY_DIR}/${SHADER_NAME}.cso")
    
    add_custom_command(
        OUTPUT ${COMPILED_SHADER}
        COMMAND dxc.exe -E${SHADER_ENTRY} -Tcs_6_6 "${SHADER_FILE}" -Fo"${COMPILED_SHADER}"
        DEPENDS ${SHADER_FILE}
        COMMENT "Compiling HLSL shader: ${SHADER_FILE}"
    )
    
    set(${SHADER_TARGET} ${COMPILED_SHADER} PARENT_SCOPE)
endfunction()

# Compile HLSL shader
add_hlsl_shader(src/fft/fft_gpu.hlsl FFTKernel_1D FFT_SHADER_CSO)

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
