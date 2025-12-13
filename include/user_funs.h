// include/user_funs.h
// user_funs.h
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
matrix ff3T(matrix x, matrix ud1=NAN, matrix ud2=NAN);
matrix ff3T_zew(matrix x, matrix ud1=NAN, matrix ud2=NAN);
matrix ff3T_wew(matrix x, matrix ud1=NAN, matrix ud2=NAN);
matrix dff3R(double t, matrix Y, matrix ud1=NAN, matrix ud2=NAN);
matrix ff3R_zew(matrix x, matrix ud1=NAN, matrix ud2=NAN);

// --- Funkcje Lab 4 ---
// Testowa funkcja celu
//

matrix read_matrix_semicolon(string path, int rows, int cols);

double calculate_accuracy(matrix theta, matrix X, matrix Y);

matrix ff4T(matrix x, matrix ud1 = NAN, matrix ud2 = NAN);
// Gradient testowej funkcji celu
matrix gf4T(matrix x, matrix ud1 = NAN, matrix ud2 = NAN);
// Hesjan testowej funkcji celu
matrix Hf4T(matrix x, matrix ud1 = NAN, matrix ud2 = NAN);

// Rzeczywista funkcja celu (Regresja logistyczna)
matrix ff4R(matrix theta, matrix ud1, matrix ud2);
// Gradient rzeczywistej funkcji celu
matrix gf4R(matrix theta, matrix ud1, matrix ud2);
// --- Funkcje Lab 5 ---

/// Testowa funkcja celu (Ważona suma)
/// ud1(0) = w (waga), ud2(0) = a (parametr)
matrix ff5T(matrix x, matrix ud1, matrix ud2);

/// Rzeczywista funkcja celu (Ważona suma z funkcją kary)
/// ud1(0) = w (waga)
matrix ff5R(matrix x, matrix ud1, matrix ud2);
