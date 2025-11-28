//Ten plik ~~nie~~ powinien być edytowany
//
// ^Akurat powinien, bo trzeba dodać dokumentację interfejsu,
// obecnie jest wygenerowana przez AI i za dokładność
// nie odpowiadam. Może się uzupełni na podstawie info z
// poprzednich lat czy coś.
//
// PS: Monte Carlo wzięte z pliku cpp, więc się powinno zgadzać,
// reszta XD że nam kazali implementować to na stałej liczbie
// argumentów bez wyjaśnienia, co one robią.
//
// --- Dawid P.

#pragma once

#include"solution.h"

/// Metoda Monte Carlo
///
/// @param ff Funkcja celu
/// @param N Liczba zmiennych funkcji celu
/// @param lb Dolne ograniczenie
/// @param ub Górne ograniczenie
/// @param epsilon Zakładana dokładność rozwiązania
/// @param Nmax Maksymalna liczba wywołań funkcji celu
/// @param ud1 Dane użytkownika 1
/// @param ud2 Dane użytkownika 2
solution MC(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, double epsilon, int Nmax, matrix ud1 = NAN, matrix ud2 = NAN);

/// Metoda rozszerzenia
///
/// @param ff Funkcja celu
/// @param x0 Początkowa wartość
/// @param d Kierunek
/// @param alpha Parametr
/// @param Nmax Maksymalna liczba wywołań funkcji celu
/// @param ud1 Dane użytkownika 1
/// @param ud2 Dane użytkownika 2
double* expansion(matrix(*ff)(matrix, matrix, matrix), double x0, double d, double alpha, int Nmax, matrix ud1 = NAN, matrix ud2 = NAN); // throw (string);

/// Metoda Fibonacciego
///
/// @param ff Funkcja celu
/// @param a Lewa granica przedziału
/// @param b Prawa granica przedziału
/// @param epsilon Zakładana dokładność rozwiązania
/// @param ud1 Dane użytkownika 1
/// @param ud2 Dane użytkownika 2
solution fib(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, matrix ud1 = NAN, matrix ud2 = NAN); // throw (string);

/// Metoda lagrangego
///
/// @param ff Funkcja celu
/// @param a Lewa granica przedziału
/// @param b Prawa granica przedziału
/// @param epsilon Zakładana dokładność rozwiązania
/// @param gamma Współczynnik kroku
/// @param Nmax Maksymalna liczba wywołań funkcji celu
/// @param ud1 Dane użytkownika 1
/// @param ud2 Dane użytkownika 2
solution lag(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, double gamma, int Nmax, matrix ud1 = NAN, matrix ud2 = NAN); // throw (string);

/// Metoda Hooke-Jeevesa
///
/// @param ff Funkcja celu
/// @param x0 Początkowy punkt
/// @param s Współczynnik kroku
/// @param alpha Współczynnik zmiany kroku
/// @param epsilon Zakładana dokładność rozwiązania
/// @param Nmax Maksymalna liczba wywołań funkcji celu
/// @param ud1 Dane użytkownika 1
/// @param ud2 Dane użytkownika 2
solution HJ(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double epsilon, int Nmax, matrix ud1 = NAN, matrix ud2 = NAN); // throw (string);

/// Metoda Hestenesa-Johnsona z próbą kroku
///
/// @param ff Funkcja celu
/// @param XB Rozwiązanie początkowe
/// @param s Współczynnik kroku
/// @param ud1 Dane użytkownika 1
/// @param ud2 Dane użytkownika 2
solution HJ_trial(matrix(*ff)(matrix, matrix, matrix), solution XB, double s, matrix ud1 = NAN, matrix ud2 = NAN); // throw (string);

/// Metoda Rosenbrocka
///
/// @param ff Funkcja celu
/// @param x0 Rozwiązanie początkowe
/// @param s0 Współczynnik kroku
/// @param alpha Współczynnik kroku
/// @param beta Współczynnik kroku
/// @param epsilon Współczynnik kroku
/// @param Nmax Współczynnik kroku
/// @param ud1 Dane użytkownika 1
/// @param ud2 Dane użytkownika 2
solution Rosen(matrix(*ff)(matrix, matrix, matrix), matrix x0, matrix s0, double alpha, double beta, double epsilon, int Nmax, matrix ud1 = NAN, matrix ud2 = NAN); // throw (string);

/// Metoda Newtona-Penaltyego
///
/// @param ff Funkcja celu
/// @param x0 Początkowy wektor
/// @param c Współczynnik kroku
/// @param dc Współczynnik kroku
/// @param epsilon Współczynnik kroku
/// @param Nmax Współczynnik kroku
/// @param ud1 Dane użytkownika 1
/// @param ud2 Dane użytkownika 2
solution pen(matrix(*ff)(matrix, matrix, matrix), matrix x0, double c, double dc, double epsilon, int Nmax, matrix ud1 = NAN, matrix ud2 = NAN); // throw (string);

/// Metoda symetryczna Newtona-Morego
///
/// @param ff Funkcja celu
/// @param x0 Początkowy wektor
/// @param s Współczynnik kroku
/// @param alpha Współczynnik kroku
/// @param beta Współczynnik kroku
/// @param gamma Współczynnik kroku
/// @param delta Współczynnik kroku
/// @param epsilon Zakładana dokładność rozwiazania
/// @param Nmax Maksymalna liczba wywołań funkcji celu
/// @param ud1 Dane użytkownika 1
/// @param ud2 Dane użytkownika 2
solution sym_NM(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double beta, double gamma, double delta, double epsilon, int Nmax, matrix ud1 = NAN, matrix ud2 = NAN); // throw (string);

/// Metoda Steepest Descent
///
/// @param ff Funkcja celu
/// @param gf Funkcja gradientu
/// @param x0 Początkowy wektor
/// @param h0 Współczynnik kroku
/// @param epsilon Zakładana dokładność rozwiazania
/// @param Nmax Maksymalna liczba wywołań funkcji celu
/// @param ud1 Dane użytkownika 1
/// @param ud2 Dane użytkownika 2
solution SD(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1 = NAN, matrix ud2 = NAN); // throw (string);

/// Metoda Conjugate Gradient
///
/// @param ff Funkcja celu
/// @param gf Funkcja gradientu
/// @param x0 Początkowy wektor
/// @param h0 Współczynnik kroku
/// @param epsilon Zakładana dokładność rozwiazania
/// @param Nmax Maksymalna liczba wywołań funkcji celu
/// @param ud1 Dane użytkownika 1
/// @param ud2 Dane użytkownika 2
solution CG(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1 = NAN, matrix ud2 = NAN); // throw (string);

/// Metoda Newtona
///
/// @param ff Funkcja celu
/// @param gf Funkcja gradientu
/// @param Hf Funkcja hessa
/// @param x0 Początkowy wektor
/// @param h0 Współczynnik kroku
/// @param epsilon Zakładana dokładność rozwiazania
/// @param Nmax Maksymalna liczba wywołań funkcji celu
/// @param ud1 Dane użytkownika 1
/// @param ud2 Dane użytkownika 2
solution Newton(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix(*Hf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1 = NAN, matrix ud2 = NAN); // throw (string);

/// Metoda Golden
///
/// @param ff Funkcja celu
/// @param a Lewa granica przedziału
/// @param b Prawa granica przedziału
/// @param epsilon Zakładana dokładność rozwiazania
/// @param Nmax Maksymalna liczba wywołań funkcji celu
/// @param ud1 Dane użytkownika 1
/// @param ud2 Dane użytkownika 2
solution golden(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, int Nmax, matrix ud1 = NAN, matrix ud2 = NAN); // throw (string);

/// Metoda Powell'a
///
/// @param ff Funkcja celu
/// @param x0 Początkowy wektor
/// @param epsilon Zakładana dokładność rozwiazania
/// @param Nmax Maksymalna liczba iteracji
/// @param ud1 Dane użytkownika 1
/// @param ud2 Dane użytkownika 2
solution Powell(matrix(*ff)(matrix, matrix, matrix), matrix x0, double epsilon, int Nmax, matrix ud1 = NAN, matrix ud2 = NAN); // throw (string);

/// Metoda EA
///
/// @param ff Funkcja celu
/// @param N Liczba populacji
/// @param lb Dolne ograniczenie
/// @param ub Górne ograniczenie
/// @param mi Liczba rodziców
/// @param lambda Liczba potomków
/// @param sigma0 Wielkość początkowa
/// @param epsilon Zakładana dokładność rozwiazania
/// @param Nmax Maksymalna liczba iteracji
/// @param ud1 Dane użytkownika 1
/// @param ud2 Dane użytkownika 2
solution EA(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, int mi, int lambda, matrix sigma0, double epsilon, int Nmax, matrix ud1 = NAN, matrix ud2 = NAN); // throw (string);


std::string resolvePath(std::string filename);
