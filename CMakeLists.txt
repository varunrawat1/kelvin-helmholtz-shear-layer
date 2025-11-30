cmake_minimum_required(VERSION 3.16)

project(KelvinHelmholtz2D LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Library with the fluid solver core
add_library(fluidcore
    src/Field2D.cpp
    src/PoissonSolver.cpp
    src/Simulation.cpp
)

target_include_directories(fluidcore
    PUBLIC
        ${CMAKE_SOURCE_DIR}/include
)

# Main executable
add_executable(kh_sim src/main.cpp)

target_link_libraries(kh_sim PRIVATE fluidcore)
