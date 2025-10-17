//Ten plik nie powinien być edytowany

#pragma once

#include"ode_solver.h"

/// Rozwiązanie (obiekt klasy solution)
class solution
{
public:
    /// Macierz (najczęściej nx1) zawierająca
    /// współrzędne punktu stanowiące rozwiązanie
	matrix x;
	/// Macierz (najczęściej 1x1) zawierająca
	/// wartość funkcji celu w punkcie x
	matrix y;
	/// Macierz nx1 zawierająca gradient funkcji
	/// celu w punkcie x
	matrix g;
	/// Macierz nxn zawierająca hessjan funkcji
	/// celu w punkcie x
	matrix H;
	/// Macierz do dyspozycji użytkownika
	matrix ud;
	/// Liczba całkowita zawierająca informację o
	/// przyczynie zakończenia poszukiwania minimum
	/// funkcji celu
	int flag;
	/// Zmienna statyczna zawierająca liczbę wywołań
	/// funkcji celu
	static int f_calls;
	/// Zmienna statyczna zawierająca liczbę obliczeń
	/// gradientu funkcji celu
	static int g_calls;
	/// Zmienna statyczna zawierająca liczbę obliczeń
	/// hessjana funkcji celu
	static int H_calls;
	/// Funkcja zerująca zmienne `f_calls`, `g_calls`
	/// oraz `H_calls`
	static void clear_calls();
	/// Tworzy rozwiązanie o podanej współrzędnej
	solution(double = NAN);
	/// Tworzy rozwiązanie o przesłanych współrzędnych
	solution(const matrix&);
	/// Tworzy rozwiązanie o przesłanych współrzędnych
	///
	/// @param n liczba wymiarów rozwiązania
	/// @param M wskaźnik na tablicę wartości współrzędnych rozwiązania
	///
	/// @throw string wiadomość o błędzie
	solution(int n, double *M); // throw (string);
	/// Konstruktor kopiujący (macierz `ud` nie jest kopiowana
	/// jeżeli w obiekcie wzorcowym `ud(0, 0)` ma wartość `NAN`)
	solution(const solution&);
	/// Operator przypisania (macierz `ud` nie jest przypisywana
	/// jeżeli w obiekcie wzorcowym `ud(0, 0)` ma wartość `NAN`)
	solution& operator=(const solution&);
	/// Funkcja wyznacza wartość zmiennej `y`
	/// poprzez wywołanie funkcji, której adres
	/// jest pierwszym argumentem oraz zwiększa `f_calls`
	///
	/// @remarks
	/// Dwa kolejne argumenty, `ud1` i `ud2`, przekazywane są
	/// do `fit_fun` przez implementacje metod optymalizacji.
	///
	/// @throw string wiadomość o błędzie
	matrix fit_fun(matrix(*)(matrix, matrix, matrix), matrix ud1 = NAN, matrix ud2 = NAN); // throw (string);
	/// Funkcja wyznacza wartość zmiennej `g`
	/// poprzez wywołanie funkcji, której adres
	/// jest pierwszym argumentem oraz zwiększa `g_calls`
	///
	/// @remarks
	/// Dwa kolejne argumenty, `ud1` i `ud2`, przekazywane są
	/// do `fit_fun` przez implementacje metod optymalizacji.
	///
	/// @throw string wiadomość o błędzie
	matrix grad(matrix(*)(matrix, matrix, matrix), matrix ud1 = NAN, matrix ud2 = NAN); // throw (string);
	/// Funkcja wyznacza wartość zmiennej `H`
	/// poprzez wywołanie funkcji, której adres
	/// jest pierwszym argumentem oraz zwiększa `h_calls`
	///
	/// @remarks
	/// Dwa kolejne argumenty, `ud1` i `ud2`, przekazywane są
	/// do `fit_fun` przez implementacje metod optymalizacji.
	///
	/// @throw string wiadomość o błędzie
	matrix hess(matrix(*)(matrix, matrix, matrix), matrix ud1 = NAN, matrix ud2 = NAN); // throw (string);
};

/// Funkcja zwraca długość wektora (macierzy nx1) `x`
int get_dim(const solution&); // throw (string);
/// Wypisanie rozwiązania na ekran
/// (wartości g_calls i H_calls są wypisywane tylko wtedy,
/// gdy są większe od 0)
ostream& operator<<(ostream&, const solution&);
