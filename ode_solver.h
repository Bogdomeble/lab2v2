//Ten plik nie powinien być edytowany

#pragma once

#include"matrix.h"
#include"user_funs.h"

/// Rozwiązuje równanie różniczkowe metodą Rungego-Kutty.
///
/// @param diff równanie różniczkowe
/// @param t0 poczatkowy czas
/// @param dt krok czasowy
/// @param tend koncowy czas
/// @param Y0 poczatkowe warunki brzegowe
/// @param ud1 dane użytkownika 1
/// @param ud2 dane użytkownika 2
matrix* solve_ode(matrix(*diff)(double, matrix, matrix, matrix), double t0, double dt, double tend, matrix Y0, matrix ud1, matrix ud2); // throw (string);
