//opt_alg.cpp
#include "../include/opt_alg.h"
#include <fstream>
#include <string>
#include <cmath>



std::string resolvePath(std::string filename) {
#ifdef DATA_DIR
    return std::string(DATA_DIR) + filename;
#else
    // Fallback, gdyby nie było kompilowane przez CMake
    return "data/" + filename;
#endif
}

solution MC(matrix (*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, double epsilon, int Nmax, matrix ud1,
            matrix ud2) {
    try {
        solution Xopt;
        while (true) {
            Xopt = rand_mat(N);
            for (int i = 0; i < N; ++i)
                Xopt.x(i) = (ub(i) - lb(i)) * Xopt.x(i) + lb(i);
            Xopt.fit_fun(ff, ud1, ud2);
            if (Xopt.y < epsilon) {
                Xopt.flag = 1;
                break;
            }
            if (solution::f_calls > Nmax) {
                Xopt.flag = 0;
                break;
            }
        }
        return Xopt;
    } catch (string ex_info) {
        throw ("solution MC(...):\n" + ex_info);
    }
}

double *expansion(matrix (*ff)(matrix, matrix, matrix), double x0, double d, double alpha, int Nmax, matrix ud1,
                  matrix ud2) {
    try {
        solution s0(x0);
        s0.fit_fun(ff, ud1, ud2);

        solution s1(x0 + d);
        s1.fit_fun(ff, ud1, ud2);


        if (s1.y == s0.y) {
            return new double[2]{m2d(s0.x), m2d(s1.x)};
        }


        if (s1.y > s0.y) {
            double original_s1_x = s1.x(0);
            d = -d;
            s1.x = x0 + d;
            s1.fit_fun(ff, ud1, ud2);


            if (s1.y >= s0.y) {
                if (d > 0) return new double[2]{original_s1_x, s1.x(0)};
                else return new double[2]{s1.x(0), original_s1_x};
            }
        }


        solution s_prev = s1;
        solution s_curr;
        int i = 1;

        while (true) {
            if (solution::f_calls >= Nmax) {
                throw string("Maximum number of calls to target function in method exceeded.");
            }

            s_curr.x = x0 + pow(alpha, i) * d;
            s_curr.fit_fun(ff, ud1, ud2);

            if (s_curr.y >= s_prev.y) {
                break;
            }

            s_prev = s_curr;
            i++;
        }

        double *p = new double[2];
        double x_prev_prev = (i == 1) ? x0 : (x0 + pow(alpha, i - 2) * d);

        if (d > 0) {
            p[0] = x_prev_prev;
            p[1] = m2d(s_curr.x);
        } else {
            p[0] = m2d(s_curr.x);
            p[1] = x_prev_prev;
        }

        return p;
    } catch (string ex_info) {
        throw ("double* expansion(...):\n" + ex_info);
    }
}

solution fib(matrix (*ff)(matrix, matrix, matrix), double a, double b, double epsilon, matrix ud1, matrix ud2) {
    try {
        std::vector<long long> F;
        F.push_back(1);
        F.push_back(1);
        int k = 1;
        while (F.back() <= (b - a) / epsilon) {
            F.push_back(F[k] + F[k - 1]);
            k++;
        }


        solution A(a), B(b);
        solution C(B.x(0) - static_cast<double>(F[k - 1]) / F[k] * (B.x(0) - A.x(0)));
        solution D(A.x(0) + B.x(0) - C.x(0));
        C.fit_fun(ff, ud1, ud2);
        D.fit_fun(ff, ud1, ud2);

        // --- ADDED FOR PLOTTING ---
        // These lines will print the data needed for the chart.
        //cout << "--- Fibonacci Interval Lengths ---" << endl;
        //cout << "Iteration,Interval_Length" << endl;
        //cout << 0 << "," << (B.x(0) - A.x(0)) << endl; // Print initial state
        // --------------------------

        for (int i = 0; i <= k - 3; ++i) {
            if (C.y < D.y) {
                B = D;
                D = C;

                C.x = B.x - static_cast<double>(F[k - i - 2]) / F[k - i - 1] * (B.x - A.x);
                C.fit_fun(ff, ud1, ud2);
            } else {
                A = C;
                C = D;

                D.x = A.x + B.x - C.x;
                D.fit_fun(ff, ud1, ud2);
            }
            // --- ADDED FOR PLOTTING ---
            // This prints the interval length at the end of each iteration.
            //cout << i + 1 << "," << (B.x(0) - A.x(0)) << endl;
            // --------------------------
        }


        solution final_solution = (C.y < D.y) ? C : D;
        final_solution.flag = 1;
        return final_solution;
    } catch (string ex_info) {
        throw ("solution fib(...):\n" + ex_info);
    }
}

solution lag(matrix (*ff)(matrix, matrix, matrix), double a, double b, double epsilon, double gamma, int Nmax,
             matrix ud1, matrix ud2) {
    try {
        solution A(a), B(b), C((a + b) / 2.0);
        A.fit_fun(ff, ud1, ud2);
        B.fit_fun(ff, ud1, ud2);
        C.fit_fun(ff, ud1, ud2);

        solution D;
        double d_prev = NAN;

        // --- ADDED FOR PLOTTING ---
        //int iteration = 0;
        //cout << "--- Lagrange Interval Lengths ---" << endl;
        //cout << "Iteration,Interval_Length" << endl;
        //cout << iteration << "," << (B.x(0) - A.x(0)) << endl; // Print initial state
        // --------------------------

        while (true) {
            if (solution::f_calls > Nmax) {
                D.flag = 0;
                return D;
            }

            double l = A.y(0) * (pow(B.x(0), 2) - pow(C.x(0), 2)) + B.y(0) * (pow(C.x(0), 2) - pow(A.x(0), 2)) + C.y(0)
                       * (pow(A.x(0), 2) - pow(B.x(0), 2));
            double m = A.y(0) * (B.x(0) - C.x(0)) + B.y(0) * (C.x(0) - A.x(0)) + C.y(0) * (A.x(0) - B.x(0));

            if (m <= 0) {
                throw string("Error in Lagrange's method: m <= 0. The parabola is directed in the wrong direction.");
            }

            double d_val = 0.5 * l / m;
            D.x = d_val;

            if (abs(B.x(0) - A.x(0)) < epsilon || (!isnan(d_prev) && abs(d_val - d_prev) < gamma)) {
                D.flag = 1;
                D.fit_fun(ff, ud1, ud2);
                return D;
            }

            d_prev = d_val;
            D.fit_fun(ff, ud1, ud2);


            if (d_val > A.x(0) && d_val < C.x(0)) {
                if (D.y < C.y) {
                    B = C;
                    C = D;
                } else {
                    A = D;
                }
            } else if (d_val > C.x(0) && d_val < B.x(0)) {
                if (D.y < C.y) {
                    A = C;
                    C = D;
                } else {
                    B = D;
                }
            } else {
                throw string("Error in Lagrange's method: the new point D fell outside the interval [A, B].");
            }

            // --- ADDED FOR PLOTTING ---
            //iteration++;
            //cout << iteration << "," << (B.x(0) - A.x(0)) << endl;
            // --------------------------
        }
    } catch (string ex_info) {
        throw ("solution lag(...):\n" + ex_info);
    }
}

solution HJ(matrix (*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double epsilon, int Nmax,
            matrix ud1, matrix ud2) {
    try {
		std::ofstream log_file("hj_iterations.csv");
		if (log_file.is_open()) {
			log_file << "Iteration,x1*,x2*\\n";
		}
		int iter = 0;

        solution XB(x0), X;
        XB.fit_fun(ff, ud1, ud2);

		if (log_file.is_open()) {
			log_file << iter << "," << XB.x(0) << "," << XB.x(1) << "\\n";
		}

        do {
            X = HJ_trial(ff, XB, s, ud1, ud2);
            if (X.y < XB.y) {
                solution X_prev;
                do {
                    X_prev = XB;
                    XB = X;
                    // Pattern move
                    X.x = 2.0 * XB.x - X_prev.x;
                    X = HJ_trial(ff, X, s, ud1, ud2);

                    if(solution::f_calls > Nmax) break;

                } while (X.y < XB.y);
                X = XB;
            } else {
                s *= alpha;
            }
			iter++;
			if (log_file.is_open()) {
				log_file << iter << "," << XB.x(0) << "," << XB.x(1) << "\\n";
			}
            if(solution::f_calls > Nmax) {
                X.flag = 0; // Timeout
                break;
            }
        } while (s > epsilon);

		if (log_file.is_open()) {
			log_file.close();
		}

        X.flag = 1; // Success
        return X;
    } catch (string ex_info) {
        throw ("solution HJ(...):\\n" + ex_info);
    }
}


solution HJ_trial(matrix (*ff)(matrix, matrix, matrix), solution XB, double s, matrix ud1, matrix ud2) {
    int n = get_len(XB.x);
    solution X = XB;

    for (int j = 0; j < n; ++j) {
        // Probe in positive direction
        X.x(j) += s;
        X.fit_fun(ff, ud1, ud2);
        if (X.y < XB.y) {
            XB = X;
        } else {
            // Probe in negative direction
            X.x(j) -= 2 * s;
            X.fit_fun(ff, ud1, ud2);
            if (X.y < XB.y) {
                XB = X;
            } else {
                // Return to original if no improvement
                X.x(j) += s;
            }
        }
    }
    return XB;
}

solution Rosen(matrix (*ff)(matrix, matrix, matrix), matrix x0, matrix s0, double alpha, double beta, double epsilon,
               int Nmax, matrix ud1, matrix ud2) {
    try {
		std::ofstream log_file("rosen_iterations.csv");
		if (log_file.is_open()) {
			log_file << "Iteration,x1*,x2*\\n";
		}
		int iter = 0;

        int n = get_len(x0);
        solution X(x0);
        X.fit_fun(ff, ud1, ud2);

		if (log_file.is_open()) {
			log_file << iter << "," << X.x(0) << "," << X.x(1) << "\\n";
		}

        matrix s = s0;
        matrix D = ident_mat(n);
        matrix lambda(n, 1);
        matrix p(n, 1);

        do {
            for (int j = 0; j < n; ++j) {
                solution X_temp = X;
                X_temp.x = X.x + s(j) * D[j];
                X_temp.fit_fun(ff, ud1, ud2);

                if (X_temp.y < X.y) {
                    X = X_temp;
                    lambda(j) += s(j);
                    s(j) *= alpha;
                } else {
                    p(j) += 1;
                    s(j) *= -beta;
                }
            }

            bool change_basis = true;
            for(int j=0; j<n; ++j) {
                if(lambda(j) == 0 || p(j) == 0) {
                    change_basis = false;
                    break;
                }
            }

            if (change_basis) {
                matrix Q(n, n);
                for(int j=0; j<n; ++j) {
                    matrix L(n, 1);
                    for(int k=j; k<n; ++k)
                        L(k) = lambda(k);
                    Q.set_col(D * L, j);
                }

                matrix v = Q[0];
                D.set_col(v / norm(v), 0);

                for (int j = 1; j < n; ++j) {
                    matrix temp = Q[j];
                    for (int k = 0; k < j; ++k) {
                        temp = temp - (trans(Q[j]) * D[k])() * D[k];
                    }
                    v = temp;
                    D.set_col(v / norm(v), j);
                }
                s = s0;
                lambda = matrix(n,1,0.0);
                p = matrix(n,1,0.0);
            }

			iter++;
			if (log_file.is_open()) {
				log_file << iter << "," << X.x(0) << "," << X.x(1) << "\\n";
			}

            if(solution::f_calls > Nmax) {
                X.flag = 0;
                break;
            }

        } while (norm(s) > epsilon);

		if (log_file.is_open()) {
			log_file.close();
		}

        X.flag = 1;
        return X;
    } catch (string ex_info) {
        throw ("solution Rosen(...):\\n" + ex_info);
    }
}

solution pen(matrix (*ff)(matrix, matrix, matrix), matrix x0, double c, double dc, double epsilon, int Nmax, matrix ud1,
             matrix ud2) {
    try {
        solution Xopt;


        return Xopt;
    } catch (string ex_info) {
        throw ("solution pen(...):\n" + ex_info);
    }
}

solution golden(matrix (*ff)(matrix, matrix, matrix), double a, double b, double epsilon, int Nmax, matrix ud1,
                matrix ud2) {
    try {
        solution Xopt;


        return Xopt;
    } catch (string ex_info) {
        throw ("solution golden(...):\n" + ex_info);
    }
}

solution EA(matrix (*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, int mi, int lambda, matrix sigma0,
            double epsilon, int Nmax, matrix ud1, matrix ud2) {
    try {
        solution Xopt;


        return Xopt;
    } catch (string ex_info) {
        throw ("solution EA(...):\n" + ex_info);
    }
}

solution sym_NM(matrix (*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double beta, double gamma,
                double delta, double epsilon, int Nmax, matrix ud1, matrix ud2) {
    try {
        int n = get_len(x0);
        solution S[n + 1];

        // 1. Inicjalizacja sympleksu
        S[0].x = x0;
        for (int i = 1; i <= n; ++i) {
            matrix e(n, 1);
            e(i - 1) = 1.0;
            S[i].x = x0 + s * e;
        }

        // Pętla główna
        do {
            // Obliczenie wartości funkcji celu dla wierzchołków
            for (int i = 0; i <= n; ++i) {
                S[i].fit_fun(ff, ud1, ud2);
            }

            // 2. Sortowanie wierzchołków
            int i_min = 0, i_max = 0;
            for (int i = 1; i <= n; ++i) {
                if (S[i].y < S[i_min].y) i_min = i;
                if (S[i].y > S[i_max].y) i_max = i;
            }

            // 3. Obliczenie środka ciężkości (bez najgorszego punktu)
            matrix p(n, 1);
            for (int i = 0; i <= n; ++i) {
                if (i != i_max) {
                    p = p + S[i].x;
                }
            }
            p = p / n;

            // 4. Odbicie (refleksja)
            solution p_odb = S[i_max];
            p_odb.x = p + alpha * (p - S[i_max].x);
            p_odb.fit_fun(ff, ud1, ud2);

            if (p_odb.y < S[i_min].y) {
                // 5. Ekpansja
                solution p_e = S[i_max];
                p_e.x = p + gamma * (p_odb.x - p);
                p_e.fit_fun(ff, ud1, ud2);
                if (p_e.y < p_odb.y) {
                    S[i_max] = p_e;
                } else {
                    S[i_max] = p_odb;
                }
            } else if (p_odb.y >= S[i_min].y) {
                bool better_than_second_worst = true;
                for (int i = 0; i <= n; i++) {
                    if (i != i_max && p_odb.y > S[i].y) {
                        better_than_second_worst = false;
                        break;
                    }
                }

                if (better_than_second_worst) {
                    S[i_max] = p_odb;
                } else {
                    // 6. Zawężenie
                    if (p_odb.y < S[i_max].y) {
                        S[i_max] = p_odb;
                    }
                    solution p_z = S[i_max];
                    p_z.x = p + beta * (S[i_max].x - p);
                    p_z.fit_fun(ff, ud1, ud2);

                    if (p_z.y >= S[i_max].y) {
                        // 7. Redukcja (ściągnięcie sympleksu)
                        for (int i = 0; i <= n; ++i) {
                            if (i != i_min) {
                                S[i].x = delta * (S[i].x + S[i_min].x);
                            }
                        }
                    } else {
                        S[i_max] = p_z;
                    }
                }
            }
            // Sprawdzenie warunku stopu
            double max_dist = 0;
            for (int i = 0; i <= n; ++i) {
                double dist = norm(S[i_min].x - S[i].x);
                if (dist > max_dist) {
                    max_dist = dist;
                }
            }
            if (max_dist < epsilon || solution::f_calls > Nmax) {
                S[0].flag = (solution::f_calls <= Nmax);
                break;
            }

        } while (true);

        // Znajdź ostateczne najlepsze rozwiązanie
        int final_best_idx = 0;
        for (int i = 1; i <= n; i++) {
            if (S[i].y < S[final_best_idx].y) {
                final_best_idx = i;
            }
        }
        return S[final_best_idx];

    } catch (string ex_info) {
        throw ("solution sym_NM(...):\n" + ex_info);
    }
}

// Funkcja pomocnicza: Złoty podział dla szukania kroku h
// Minimalizuje f(x0 + h * d) w przedziale [a, b]
double golden_search(matrix(*ff)(matrix, matrix, matrix), matrix x0, matrix d, matrix ud1, matrix ud2, double a, double b, double epsilon, int Nmax) {
    double alpha = (sqrt(5.0) - 1.0) / 2.0;
    double c = b - alpha * (b - a);
    double d_point = a + alpha * (b - a); // nazwa d_point, bo d to kierunek
    
    // Obliczenie wartości funkcji w punktach c i d_point
    matrix xc = x0 + c * d;
    matrix xd = x0 + d_point * d;
    double fc = m2d(ff(xc, ud1, ud2));
    double fd = m2d(ff(xd, ud1, ud2));
    
    int calls = 0;
    while ((b - a) > epsilon && calls < Nmax) {
        if (fc < fd) {
            b = d_point;
            d_point = c;
            fd = fc;
            c = b - alpha * (b - a);
            xc = x0 + c * d;
            fc = m2d(ff(xc, ud1, ud2));
        } else {
            a = c;
            c = d_point;
            fc = fd;
            d_point = a + alpha * (b - a);
            xd = x0 + d_point * d;
            fd = m2d(ff(xd, ud1, ud2));
        }
        calls++;
    }
    return (a + b) / 2.0;
}

solution SD(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2) {
    try {
        solution X(x0);
        solution X_prev;
        matrix d;
        double h;
        
        // Oblicz wartość początkową
        X.fit_fun(ff, ud1, ud2);
        
        while (solution::f_calls < Nmax) {
            // 1. Wyznacz kierunek (antygradient)
            d = -gf(X.x, ud1, ud2);
            
            // 2. Wyznacz długość kroku
            if (std::isnan(h0)) {
                // Zmiennokrokowy (Metoda Złotego Podziału)
                // Szukamy w przedziale np. [0, 1] - typowe dla normalizowanych problemów, lub większe
                h = golden_search(ff, X.x, d, ud1, ud2, 0.0, 1.0, epsilon, Nmax); 
            } else {
                // Stałokrokowy
                h = h0;
            }
            
            // 3. Nowy punkt
            X_prev = X;
            X.x = X.x + h * d;
            X.fit_fun(ff, ud1, ud2);
            
            // 4. Warunek stopu (zmiana argumentu)
            if (norm(X.x - X_prev.x) < epsilon) {
                X.flag = 1; // Sukces
                break;
            }
        }
        
        if (solution::f_calls >= Nmax) X.flag = 0;
        return X;
    } catch (string ex_info) {
        throw ("solution SD(...):\n" + ex_info);
    }
}

solution CG(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2) {
    try {
        solution X(x0);
        solution X_prev;
        matrix d, d_prev, g, g_prev;
        double h, beta;
        
        X.fit_fun(ff, ud1, ud2);
        g = gf(X.x, ud1, ud2);
        d = -g;
        
        while (solution::f_calls < Nmax) {
            // Wyznacz długość kroku
             if (std::isnan(h0)) {
                h = golden_search(ff, X.x, d, ud1, ud2, 0.0, 1.0, epsilon, Nmax);
            } else {
                h = h0;
            }
            
            X_prev = X;
            X.x = X.x + h * d;
            X.fit_fun(ff, ud1, ud2);
            
            if (norm(X.x - X_prev.x) < epsilon) {
                X.flag = 1;
                break;
            }
            
            g_prev = g;
            g = gf(X.x, ud1, ud2);
            
            // Obliczenie Beta (Fletcher-Reeves)
            double nom = pow(norm(g), 2);
            double den = pow(norm(g_prev), 2);
            beta = nom / den;
            
            d_prev = d;
            d = -g + beta * d_prev;
        }
        
        if (solution::f_calls >= Nmax) X.flag = 0;
        return X;
    } catch (string ex_info) {
        throw ("solution CG(...):\n" + ex_info);
    }
}

solution Newton(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix(*Hf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2) {
    try {
        solution X(x0);
        solution X_prev;
        matrix d, g, H;
        double h;
        
        X.fit_fun(ff, ud1, ud2);
        
        while (solution::f_calls < Nmax) {
            g = gf(X.x, ud1, ud2);
            H = Hf(X.x, ud1, ud2);
            
            // Kierunek Newtona: d = -H^-1 * g
            d = -inv(H) * g;
            
            if (std::isnan(h0)) {
                h = golden_search(ff, X.x, d, ud1, ud2, 0.0, 1.0, epsilon, Nmax);
            } else {
                h = h0;
            }
            
            X_prev = X;
            X.x = X.x + h * d;
            X.fit_fun(ff, ud1, ud2);
            
            if (norm(X.x - X_prev.x) < epsilon) {
                X.flag = 1;
                break;
            }
        }
        
        if (solution::f_calls >= Nmax) X.flag = 0;
        return X;
    } catch (string ex_info) {
        throw ("solution Newton(...):\n" + ex_info);
    }
}

// Zmienne statyczne
static matrix (*powell_ff)(matrix, matrix, matrix) = nullptr;
static matrix powell_x0;
static matrix powell_d;
static matrix powell_ud1;
static matrix powell_ud2;

// Wrapper
matrix ff_powell_wrapper(matrix h, matrix ud1, matrix ud2) {
    // h(0) to szukany skalar (krok)
    if (std::isnan(h(0))) return matrix(1e100); // Zabezpieczenie
    matrix x = powell_x0 + h(0) * powell_d;
    return powell_ff(x, powell_ud1, powell_ud2);
}

solution Powell(matrix(*ff)(matrix, matrix, matrix), matrix x0, double epsilon, int Nmax, matrix ud1, matrix ud2) {
    try {
        int n = get_len(x0);
        matrix D = ident_mat(n); // Macierz kierunków
        solution P(x0);
        P.fit_fun(ff, ud1, ud2); // Oblicz wartość startową
        
        powell_ff = ff;
        powell_ud1 = ud1;
        powell_ud2 = ud2;

        while (solution::f_calls < Nmax) {
            matrix P0 = P.x; // Zapamiętaj punkt startowy tej iteracji
            
            // Pętla po kierunkach
            for (int i = 0; i < n; ++i) {
                powell_x0 = P.x;
                powell_d = get_col(D, i);
                
                double* interval = nullptr;
                solution h_opt;
                bool success = false;

                try {
                    // 1. Ekspansja
                    interval = expansion(ff_powell_wrapper, 0.0, 0.1, 2.0, Nmax);
                    
                    // 2. Złoty podział
                    if (interval != nullptr && !std::isnan(interval[0]) && !std::isnan(interval[1])) {
                        h_opt = golden(ff_powell_wrapper, interval[0], interval[1], epsilon, Nmax);
                        success = true;
                    }
                } catch (...) {
                    // Ignorujemy błędy w pojedynczym kroku, próbujemy iść dalej lub kończymy
                    success = false;
                }

                // Sprzątanie po ekspansji
                if (interval != nullptr) {
                    delete[] interval;
                    interval = nullptr;
                }

                // Aktualizacja punktu TYLKO jeśli znaleziono poprawny krok
                if (success && !std::isnan(h_opt.x(0))) {
                    P.x = P.x + h_opt.x(0) * powell_d;
                }
                
                // Check limitu wywołań wewnątrz pętli
                if (solution::f_calls >= Nmax) {
                    P.fit_fun(ff, ud1, ud2);
                    P.flag = 0;
                    return P;
                }
            }

            // Sprawdzenie zbieżności: ||P_n - P_0|| < epsilon
            if (norm(P.x - P0) < epsilon) {
                P.fit_fun(ff, ud1, ud2);
                P.flag = 1;
                return P;
            }
            
            // Aktualizacja bazy kierunków
            matrix new_dir = P.x - P0;
            if (norm(new_dir) > 1e-9) { 
                // Przesuń kierunki: usuń pierwszy, dodaj nowy na koniec
                matrix D_new(n, n);
                for (int i = 0; i < n - 1; ++i) {
                    D_new.set_col(get_col(D, i + 1), i);
                }
                D_new.set_col(new_dir, n - 1);
                D = D_new;

                // Dodatkowy krok w nowym kierunku sprzężonym
                powell_x0 = P.x;
                powell_d = new_dir;
                
                double* interval = nullptr;
                try {
                    interval = expansion(ff_powell_wrapper, 0.0, 0.1, 2.0, Nmax);
                    if (interval != nullptr && !std::isnan(interval[0]) && !std::isnan(interval[1])) {
                         solution h_opt = golden(ff_powell_wrapper, interval[0], interval[1], epsilon, Nmax);
                         if (!std::isnan(h_opt.x(0))) {
                             P.x = P.x + h_opt.x(0) * new_dir;
                         }
                    }
                } catch (...) {
                    // Ignoruj błąd optymalizacji w nowym kierunku
                }
                if (interval != nullptr) delete[] interval;
            }
        }

        P.fit_fun(ff, ud1, ud2);
        P.flag = 0; // Przekroczono Nmax
        return P;
        
    } catch (string ex_info) {
        // W ostateczności zwróć to co mamy, zamiast zabijać program
        solution P(x0);
        P.fit_fun(ff, ud1, ud2);
        P.flag = -1;
        return P;
    }
}
