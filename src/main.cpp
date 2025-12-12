#include "../include/opt_alg.h"
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <cstdlib>
#include <ctime>

using namespace std;

// Deklaracje funkcji
void lab0();
void lab1();
void lab1_real();
void lab2();
void lab3();
void lab4();
void lab5();
void test_real_problem_DA50();

int main() {
    try {
        // Odkomentuj odpowiednią linię, aby uruchomić dane laboratorium
        
        // lab0();
        // lab1();
        // lab1_real();
        // lab2();
        // lab3();
        // lab4();
        lab5(); 

    } catch (string EX_INFO) {
        cerr << "ERROR (CRITICAL):\n";
        cerr << EX_INFO << endl << endl;
    } catch (...) {
        cerr << "ERROR (UNKNOWN): Wystapil nieznany blad krytyczny." << endl;
    }
    return 0;
}

// ===========================================================================
// LABORATORIUM 0
// ===========================================================================
void lab0() {
    // Test function
    double epsilon = 1e-2;
    int Nmax = 10000;
    matrix lb(2, 1, -5), ub(2, 1, 5), a(2, 1);
    solution opt;
    a(0) = -1;
    a(1) = 2;
    opt = MC(ff0T, 2, lb, ub, epsilon, Nmax, a);
    cout << opt << endl << endl;
    solution::clear_calls();

    // Pendulum
    Nmax = 1000;
    epsilon = 1e-5;
    lb = 0; ub = 5;
    double teta_opt = 1;
    opt = MC(ff0R, 1, lb, ub, epsilon, Nmax, teta_opt);
    cout << opt << endl << endl;
    solution::clear_calls();

    matrix Y0 = matrix(2, 1);
    matrix MT = matrix(2, new double[2]{m2d(opt.x), 0.5});
    matrix *Y = solve_ode(df0, 0, 0.1, 10, Y0, NAN, MT);
    ofstream Sout("symulacja_lab0.csv");
    Sout << hcat(Y[0], Y[1]);
    Sout.close();
    delete[] Y;
}

// ===========================================================================
// LABORATORIUM 1
// ===========================================================================
void lab1() {
    srand(time(nullptr));
    ofstream file("lab1_results.csv");
    double d = 1.0, alpha = 1.7, epsilon = 1e-5, gamma = 1e-8;
    int Nmax = 1000;
    
    // Uproszczona wersja dla przykładu
    double x0 = -10.0 + (double)rand()/RAND_MAX * 20.0;
    
    // Expansion
    solution::clear_calls();
    double* interval = expansion(ff1T, x0, d, alpha, Nmax);
    cout << "Expansion interval: [" << interval[0] << ", " << interval[1] << "]" << endl;
    
    // Fibonacci
    solution::clear_calls();
    solution fib_sol = fib(ff1T, interval[0], interval[1], epsilon);
    cout << "Fibonacci: " << fib_sol << endl;
    
    // Lagrange
    solution::clear_calls();
    solution lag_sol = lag(ff1T, interval[0], interval[1], epsilon, gamma, Nmax);
    cout << "Lagrange: " << lag_sol << endl;
    
    delete[] interval;
    file.close();
}

void lab1_real() {
    // Implementacja problemu rzeczywistego lab1...
}

// ===========================================================================
// LABORATORIUM 2
// ===========================================================================
void lab2() {
    srand(time(nullptr));
    ofstream file("lab2_results.csv");
    double epsilon = 1e-5;
    int Nmax = 20000;
    double alpha_hj = 0.5;
    matrix s0_rosen(2, 1); s0_rosen(0) = 0.5; s0_rosen(1) = 0.5;
    double alpha_rosen = 2.0, beta_rosen = 0.5;
    
    matrix x0(2, 1);
    x0(0) = -1.0 + (double)rand() / RAND_MAX * 2.0;
    x0(1) = -1.0 + (double)rand() / RAND_MAX * 2.0;

    // Hooke-Jeeves
    solution::clear_calls();
    solution sol_hj = HJ(ff2T, x0, 0.5, alpha_hj, epsilon, Nmax);
    cout << "HJ:\n" << sol_hj << endl;

    // Rosenbrock
    solution::clear_calls();
    solution sol_rosen = Rosen(ff2T, x0, s0_rosen, alpha_rosen, beta_rosen, epsilon, Nmax);
    cout << "Rosen:\n" << sol_rosen << endl;
    
    file.close();
}

// ===========================================================================
// LABORATORIUM 3
// ===========================================================================
void lab3() {
    srand(time(nullptr));
    ofstream file("lab3_results.csv");
    
    double s = 0.5, alpha = 1.0, beta = 0.5, gamma = 2.0, delta = 0.5;
    double epsilon = 1e-3;
    int Nmax = 10000;
    
    // Parametry kary
    double a = 4.0;
    matrix ud(2, 1); ud(0) = a; ud(1) = 1.0; // a, c

    matrix x0(2, 1);
    x0(0) = 1.5; x0(1) = 1.5;

    // Simplex z funkcją kary zewnętrznej
    solution::clear_calls();
    solution sol = sym_NM(ff3T_zew, x0, s, alpha, beta, gamma, delta, epsilon, Nmax, ud);
    cout << "Simplex (Kara zew): " << sol << endl;
    
    file.close();
}

// ===========================================================================
// LABORATORIUM 4
// ===========================================================================
void lab4() {
    srand(time(nullptr));
    cout << "--- LABORATORIUM 4: Metody Gradientowe ---" << endl;

    // --- CZĘŚĆ A: Testowa funkcja celu ---
    cout << "\n=== CZESC A: Testowa funkcja celu ===" << endl;
    
    double epsilon = 1e-5;
    int Nmax = 10000;
    
    struct StepConfig {
        double val;
        string name;
    };
    vector<StepConfig> steps = {
        {0.05, "0.05"},
        {0.25, "0.25"},
        {NAN,  "Zmienny"}
    };
    
    ofstream file_SD("lab4_test_SD.csv");
    ofstream file_CG("lab4_test_CG.csv");
    ofstream file_Newton("lab4_test_Newton.csv");

    if (file_SD.is_open() && file_CG.is_open() && file_Newton.is_open()) {
        string header = "Krok,Iteracja,x1_start,x2_start,x1_end,x2_end,y_end,f_calls,Flag\n";
        file_SD << header; file_CG << header; file_Newton << header;
        
        for (const auto& step : steps) {
            cout << "Testowanie dla kroku: " << step.name << endl;
            for (int i = 0; i < 100; ++i) {
                matrix x0(2, 1);
                x0(0) = -2.0 + (double)rand()/RAND_MAX * 4.0;
                x0(1) = -2.0 + (double)rand()/RAND_MAX * 4.0;
                
                try {
                    solution::clear_calls();
                    solution sol = SD(ff4T, gf4T, x0, step.val, epsilon, Nmax);
                    file_SD << step.name << "," << i+1 << "," << x0(0) << "," << x0(1) << ","
                            << sol.x(0) << "," << sol.x(1) << "," << sol.y(0) << "," 
                            << solution::f_calls << "," << sol.flag << "\n";
                } catch(string ex) { file_SD << step.name << "," << i+1 << ",ERROR," << ex << "\n"; }

                try {
                    solution::clear_calls();
                    solution sol = CG(ff4T, gf4T, x0, step.val, epsilon, Nmax);
                    file_CG << step.name << "," << i+1 << "," << x0(0) << "," << x0(1) << ","
                            << sol.x(0) << "," << sol.x(1) << "," << sol.y(0) << "," 
                            << solution::f_calls << "," << sol.flag << "\n";
                } catch(string ex) { file_CG << step.name << "," << i+1 << ",ERROR," << ex << "\n"; }
                
                try {
                    solution::clear_calls();
                    solution sol = Newton(ff4T, gf4T, Hf4T, x0, step.val, epsilon, Nmax);
                    file_Newton << step.name << "," << i+1 << "," << x0(0) << "," << x0(1) << ","
                                << sol.x(0) << "," << sol.x(1) << "," << sol.y(0) << "," 
                                << solution::f_calls << "," << sol.flag << "\n";
                } catch(string ex) { file_Newton << step.name << "," << i+1 << ",ERROR," << ex << "\n"; }
            }
        }
        file_SD.close(); file_CG.close(); file_Newton.close();
    }

    // --- CZĘŚĆ B: Problem rzeczywisty ---
    cout << "\n=== CZESC B: Problem rzeczywisty ===" << endl;
    
    // Dane HARDCODED (zaszyte w kodzie)
    double y_raw[] = {
        1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 
        0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 0, 1, 1, 
        0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 
        0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 
        1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1
    };
    double x1_raw[] = {
        57, 86, 86, 71, 75, 84, 89, 68, 50, 90, 53, 40, 40, 28, 97, 35, 59, 29, 44, 54, 
        22, 53, 66, 26, 43, 70, 100, 64, 79, 48, 82, 44, 50, 95, 28, 37, 71, 26, 75, 35, 
        23, 87, 92, 37, 45, 42, 31, 47, 68, 76, 86, 60, 33, 56, 78, 51, 40, 36, 51, 46, 
        36, 38, 96, 97, 95, 54, 34, 94, 65, 97, 51, 37, 77, 31, 59, 35, 81, 83, 77, 21, 
        59, 52, 93, 65, 45, 78, 50, 91, 80, 93, 82, 98, 43, 38, 56, 46, 34, 52, 82, 74
    };
    double x2_raw[] = {
        68, 87, 46, 25, 77, 53, 93, 96, 34, 92, 87, 49, 68, 56, 88, 35, 80, 75, 30, 73, 
        77, 48, 42, 87, 67, 62, 96, 56, 82, 91, 79, 52, 45, 53, 45, 38, 81, 34, 54, 98, 
        49, 99, 76, 71, 74, 59, 45, 52, 66, 49, 51, 93, 70, 37, 25, 71, 99, 55, 38, 100, 
        49, 56, 46, 26, 21, 34, 87, 25, 27, 82, 98, 63, 93, 52, 60, 72, 42, 67, 75, 28, 
        82, 66, 90, 28, 85, 85, 62, 71, 36, 61, 37, 33, 57, 67, 26, 98, 86, 92, 100, 87
    };

    matrix X(3, 100);
    matrix Y(1, 100);
    for (int i = 0; i < 100; ++i) {
        Y(0, i) = y_raw[i];
        X(0, i) = 1.0;
        X(1, i) = x1_raw[i];
        X(2, i) = x2_raw[i];
    }

    matrix theta0(3, 1, 0.0);
    vector<double> real_steps = {0.01, 0.001, 0.0001};
    ofstream real_res("lab4_real_results.csv");
    if(real_res.is_open()) {
        real_res << "Krok,theta0,theta1,theta2,Koszt_Koncowy,Iteracje,Accuracy(%)\n";
        for (double h : real_steps) {
            cout << "Optymalizacja CG dla h=" << h << "..." << endl;
            solution::clear_calls();
            try {
                solution sol = CG(ff4R, gf4R, theta0, h, epsilon, Nmax, X, Y);
                int correct = 0, m = 100;
                for (int i = 0; i < m; ++i) {
                    matrix xi = get_col(X, i);
                    double z = m2d(trans(sol.x) * xi);
                    double h_val = 1.0 / (1.0 + exp(-z));
                    if ((h_val >= 0.5 ? 1 : 0) == (int)Y(0, i)) correct++;
                }
                real_res << h << "," << sol.x(0) << "," << sol.x(1) << "," << sol.x(2) << ","
                         << sol.y(0) << "," << solution::f_calls << "," << (double)correct/m*100.0 << "\n";
            } catch (string ex) { real_res << h << ",ERROR,ERROR,ERROR,ERROR,ERROR,ERROR\n"; }
        }
        real_res.close();
    }
}

// ===========================================================================
// LABORATORIUM 5
// ===========================================================================
void lab5() {
    srand(time(nullptr));
    cout << "--- LABORATORIUM 5: Optymalizacja Wielokryterialna ---" << endl;
    
    double epsilon = 1e-5;
    int Nmax = 50000; 

    // --- CZĘŚĆ A: Testowa funkcja celu ---
    cout << "\n=== CZESC A: Testowa funkcja celu ===" << endl;
    
    vector<double> a_values = {1.0, 10.0, 100.0};
    ofstream res_test("lab5_test_results.csv");
    
    if(res_test.is_open()) {
        res_test << "a,w,x1,x2,f1,f2\n";

        for (double a : a_values) {
            cout << "Obliczenia dla a = " << a << "..." << endl;
            
            // --- OPTYMALIZACJA: HOT START ---
            // Losujemy punkt startowy TYLKO RAZ dla danej wartości 'a'
            matrix x0(2, 1);
            x0(0) = -10.0 + (double)rand()/RAND_MAX * 20.0;
            x0(1) = -10.0 + (double)rand()/RAND_MAX * 20.0;

            for (double w = 0.0; w <= 1.01; w += 0.01) { // Krok co 0.01
                if(w > 1.0) w = 1.0; 

                matrix ud1(1, 1); ud1(0) = w;
                matrix ud2(1, 1); ud2(0) = a;

                solution::clear_calls();
                // Startujemy z x0, który jest wynikiem poprzedniej pętli (lub losowy dla pierwszej)
                solution sol = Powell(ff5T, x0, epsilon, Nmax, ud1, ud2);

                // Aktualizacja punktu startowego dla następnej iteracji (w + 0.01)
                // Dzięki temu algorytm ma bardzo blisko do nowego minimum
                x0 = sol.x; 

                double x1 = sol.x(0);
                double x2 = sol.x(1);
                double f1 = a * (pow(x1 - 3, 2) + pow(x2 - 3, 2));
                double f2 = (1.0 / a) * (pow(x1 + 3, 2) + pow(x2 + 3, 2));

                res_test << a << "," << w << "," << x1 << "," << x2 << "," << f1 << "," << f2 << "\n";
            }
        }
        res_test.close();
        cout << "Wyniki testowe zapisano do lab5_test_results.csv" << endl;
    } else {
        cerr << "Blad: Nie mozna otworzyc lab5_test_results.csv" << endl;
    }

    // --- CZĘŚĆ B: Problem rzeczywisty ---
    cout << "\n=== CZESC B: Problem rzeczywisty ===" << endl;
    
    ofstream res_real("lab5_real_results.csv");
    if(res_real.is_open()) {
        res_real << "w,l,d,mass,deflection,sigma,f_calls\n";

        double P = 2000.0, E = 1.2e11, ro = 8920.0;

        cout << "Optymalizacja belki..." << endl;
        
        // --- OPTYMALIZACJA: HOT START ---
        // Losujemy punkt startowy RAZ, w bezpiecznym obszarze
        matrix x0(2, 1);
        x0(0) = 0.5;  // l (środek przedziału)
        x0(1) = 0.03; // d (środek przedziału)

        for (double w = 0.0; w <= 1.01; w += 0.01) {
            if (w > 1.0) w = 1.0;

            matrix ud1(1, 1); ud1(0) = w;

            solution::clear_calls();
            solution sol = Powell(ff5R, x0, epsilon, Nmax, ud1);

            // Aktualizacja punktu startowego dla kolejnej wagi
            x0 = sol.x;

            double l = sol.x(0);
            double d = sol.x(1);

            double mass = ro * l * M_PI * pow(d / 2.0, 2);
            double deflection = (64.0 * P * pow(l, 3)) / (3.0 * E * M_PI * pow(d, 4));
            double sigma = (32.0 * P * l) / (M_PI * pow(d, 3));

            res_real << w << "," << l << "," << d << "," 
                     << mass << "," << deflection << "," << sigma << "," 
                     << solution::f_calls << "\n";
                     
            // Opcjonalnie wypisz postęp co 10 iteracji
            if (abs(fmod(w, 0.1)) < 1e-9) {
                cout << "  w = " << w << " ukonczono." << endl;
            }
        }
        res_real.close();
        cout << "Wyniki rzeczywiste zapisano do lab5_real_results.csv" << endl;
    }
}