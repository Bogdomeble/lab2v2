// include/matrix.h
//Ten plik nie powinien być edytowany

#pragma once

#include<iostream>
#include<string>
#include<sstream>
#include<fstream>
#include<random>
#include<chrono>
#include<memory>

using namespace std;


#ifndef SEP_SYMBOL
/// Domyślna wartość separatora dla wypisywania macierzy.
///
/// @remarks
/// Sprawiłem, by się dało zmieniać przez `#define` przed `#include`,
/// ale przy kompilacji wszystkiego innego jako bibliotekę może
/// hardcodować to za ostatnio ustawioną wartość przed kompilacją,
/// więc może lepiej tego nie robić. Kod w sumie słaby do w ogóle
/// kompilacji na bibliotekę i do rozdzielania na wiele plików programów
/// zależnych, nie ma zabezpieczeń przez wieloma dyrektywami `#include`
/// więc może wywalać błędy.
#define SEP_SYMBOL ','
#endif

class matrix
{
	/// Wymiary (NxM)
	int n, m;
	/// Tablica dwuwymiarowa (z komórkami macierzy)
	double** M;
	/// Zwraca wymiary macierzy N oraz M
	friend int* get_size(const matrix&);
	/// Zwraca długość wektora pionowego N (macierzy Nx1).
	friend int get_len(const matrix&); // throw (string);
public:
	/// Tworzy macierz 1x1 o podanej wartości. Jest również wykorzystywany jako
	/// konstruktor konwertujący (zamienia m. in. double i int na matrix).
	matrix(double = 0.0);
	/// Tworzy macierz nxm, której każdy element jest równy
	/// trzeciemu argumentowi.
	matrix(int, int, double = 0.0); // throw (string);
	/// Tworzy macierz nx1 (wektor pionowy) i wypełnia wartościami
	/// podanymi w tablicy 1D.
	///
	/// @throw string Wiadomość o błędzie
	matrix(int, double*); // throw (string);
	/// Tworzy macierz nxm i wypełnia wartościami podanymi w tablicy
	/// 2D.
	///
	/// @throw string Wiadomość o błędzie
	matrix(int, int, double**); // throw (string);
	/// Konstruktor kopiujący.
	matrix(const matrix&);
	/// Destruktor.
	~matrix();
	/// Operator przypisania.
	matrix& operator=(const matrix&);
	/// Zwraca wskazaną kolumnę macierzy w postaci macierzy nx1
    /// (rezultat zwracany jest przez wartość -> użycie operatora
	/// **nie jest l-wartością**
	///
	/// @throw string Wiadomość o błędzie
	matrix operator[](int) const; // throw (string);
	/// Operator zwraca wybrany element macierzy
	/// (mogą być wykorzystane do konwersji `matrix` na `double`)
	///
	/// @throw string Wiadomość o błędzie
	double& operator()(int = 0, int = 0); // throw (string);
	/// Operator zwraca wybrany element macierzy
	/// (mogą być wykorzystane do konwersji `matrix` na `double`)
	///
	/// @throw string Wiadomość o błędzie
	double operator()(int = 0, int = 0) const; // throw (string);
	/// Wstawia we wskazaną kolumnę wektor pionowy (macierz nx1)
	void set_col(const matrix&, int); // throw (string);
	/// Wstawia we wskazaną kolumnę wektor poziomy (macierz 1xm)
	///
	/// @throw string Wiadomość o błędzie
	void set_row(const matrix&, int); // throw (string);
	/// Dodaje do macierzy nową kolumnę i wypełnia ją podaną wartością
	///
	/// @throw string Wiadomość o błędzie
	void add_col(double = 0.0);
	/// Dodaje do macierzy nowy wiersz i wypełnia go podaną wartością
	///
	/// @throw string Wiadomość o błędzie
	void add_row(double = 0.0);
	/// Dodaje do macierzy nową kolumnę i wstawiaw nią podany wektor
	/// pionowy (macierz nx1)
	void add_col(const matrix&); // throw (string);
	/// Dodaje do macierzy nowy wiersz i wstawia w niego podany
	/// wektor poziomu (macierz 1xm)
	void add_row(const matrix&); // throw (string);
};

matrix operator+(const matrix&, const matrix&); // throw (string);
matrix operator-(const matrix&, const matrix&); // throw (string);
matrix operator*(const matrix&, const matrix&); // throw (string);
matrix operator/(const matrix&, const matrix&); // throw (string);
/// Macierz przeciwna.
matrix operator-(const matrix&);
/// Operatory relacji działają tylko dla macierzy 1x1!
bool operator<(const matrix&, const matrix&); // throw (string);
/// Operatory relacji działają tylko dla macierzy 1x1!
bool operator>(const matrix&, const matrix&); // throw (string);
/// Operatory relacji działają tylko dla macierzy 1x1!
bool operator<=(const matrix&, const matrix&); // throw (string);
/// Operatory relacji działają tylko dla macierzy 1x1!
///
/// @throw string
bool operator>=(const matrix&, const matrix&); // throw (string);
/// Operatory relacji działają tylko dla macierzy 1x1!
bool operator==(const matrix&, const matrix&); // throw (string);
/// Operatory relacji działają tylko dla macierzy 1x1!
bool operator!=(const matrix&, const matrix&); // throw (string);
/// Tworzy macierz jednostkową NxN.
matrix ident_mat(int = 1); // throw (string);
/// Tworzy macierz NxM i wypełnia wartościami losowymi
/// z przedziału [0,1] o rozkładzie jednostajnym.
matrix rand_mat(int = 1, int = 1); // throw (string);
/// Tworzy macierz NxM i wypełnia wartościami losowymi
/// o standardowym rozkładzie normalnym.
matrix randn_mat(int = 1, int = 1); // throw (string);
/// Zamiana macierzy 1x1 na `double`.
double m2d(const matrix&); // throw (string);
/// Wyznacznik macierzy.
double det(const matrix&); // throw (string);
/// Macierz odwrotna.
matrix inv(const matrix&); // throw (string);
/// Macierz transponowana (LGBT puns intended).
matrix trans(const matrix&);
/// Podnosi macierz do potęgi.
matrix pow(const matrix&, int = 2); // throw (string);
/// Oblicza normę z wektora pionowego (macierzy Nx1).
double norm(const matrix&); // throw (string);
/// Poziome połączenie dwóch macierzy.
matrix hcat(const matrix&, const matrix&); // throw (string);
/// Pionowe połączenie dwóch macierzy.
matrix vcat(const matrix&, const matrix&); // throw (string);
/// Zwraca wskazaną kolumnę macierzy w postaci wektora pionowego
/// (macierzy Nx1).
matrix get_col(const matrix&, int); // throw (string);
/// Zwraca wskazany wiersz macierzy w postaci wektora poziomego
/// (macierzy 1xM).
matrix get_row(const matrix&, int); // throw (string);
/// Wypisanie macierzy na ekran lub do pliku csv.
///
/// Podczas wypisywani elementów macierzy, część całkowita
/// oddzielona jest od części ułamkowej za pomocą zdefiniowanego
/// w linii 13 `SEP_SYMBOL`.
///
/// @param os strumień wyjściowy
/// @param m macierz do wypisania
/// @return referencja do strumienia wyjściowego
ostream& operator<<(ostream& os, const matrix&);
/// Odczyt macierzy z pliku csv, txt lub klawiatury.
///
/// Podając elementy macierzy należy pamiętać, że każda liczba
/// musi kończyć się znakiem `;`. Podczas wprowadzania elementów
/// macierzy, część całkowita może być oddzielona od części ułamkowej
/// za pomocą zdefiniowanego w linii 13 `SEP_SYMBOL` lub " ".
istream& operator>>(istream&, matrix&); // throw (string);
