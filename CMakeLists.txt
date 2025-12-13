cmake_minimum_required(VERSION 3.28.3)
project(lab2)

# Ustawienie standardu C++
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Flagi kompilatora
add_compile_options(-Wall -O3 -g)

# --- KLUCZOWA POPRAWKA ---
# Używamy \" aby przekazać ścieżkę jako string (np. "/home/user/projekt/data/")
# Bez tych ukośników ścieżka mogłaby być źle interpretowana.
add_compile_definitions(DATA_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}/data/\")
# -------------------------

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

set(SOURCES
    src/main.cpp
    src/matrix.cpp
    src/ode_solver.cpp
    src/solution.cpp
    src/opt_alg.cpp
    src/user_funs.cpp
)

add_executable(lab2 ${SOURCES})
