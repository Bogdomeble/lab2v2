cmake_minimum_required(VERSION 3.28.3)
project(lab2)

set(CMAKE_CXX_STANDARD 17)
add_compile_options(-Wall -O3)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
set(CMAKE_BUILD_PARALLEL_LEVEL 2)
add_executable(lab2 src/main.cpp src/matrix.cpp src/ode_solver.cpp src/solution.cpp src/opt_alg.cpp src/user_funs.cpp)
