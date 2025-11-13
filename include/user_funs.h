#pragma once

#include"ode_solver.h"

/// Funkcja celu dla przypadku testowego
///
/// @return matrix wektor wartości funkcji celu
matrix ff0T(matrix, matrix = NAN, matrix = NAN);
/// Funkcja celu dla problemu rzeczywistego
///
/// @return matrix wektor wartości funkcji celu
matrix ff0R(matrix, matrix = NAN, matrix = NAN);
/// Zwraca wektor pochodnych szukanych funkcji
matrix df0(double, matrix, matrix = NAN, matrix = NAN);
matrix ff1T(matrix x, matrix ud1, matrix ud2);
matrix ff1R(matrix x, matrix ud1, matrix ud2);
matrix dff1R(double t, matrix Y, matrix ud1, matrix ud2);
matrix ff2T(matrix x, matrix ud1 = NAN, matrix ud2 = NAN);
matrix ff2R(matrix x, matrix ud1 = NAN, matrix ud2 = NAN);
matrix df2R(double t, matrix Y, matrix ud1 = NAN, matrix ud2 = NAN);