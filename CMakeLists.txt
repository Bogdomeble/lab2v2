cmake_minimum_required(VERSION 3.28.3)
project(lab2)

# Ustawienie standardu C++
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_compile_options(-Wall -O3)

# --- KLUCZOWA ZMIANA ---
# Definiujemy makro preprocesora DATA_DIR, które zawiera bezwzględną ścieżkę do folderu data
# Dzięki temu C++ będzie wiedział dokładnie gdzie zapisywać pliki
add_compile_definitions(DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/data/")
# -----------------------

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
