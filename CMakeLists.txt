cmake_minimum_required(VERSION 3.31)
project(lab2)

set(CMAKE_CXX_STANDARD 17)

add_executable(lab2 main.cpp matrix.cpp ode_solver.cpp solution.cpp opt_alg.cpp user_funs.cpp)
