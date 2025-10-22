#include"opt_alg.h"

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
        }
    } catch (string ex_info) {
        throw ("solution lag(...):\n" + ex_info);
    }
}

solution HJ(matrix (*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double epsilon, int Nmax,
            matrix ud1, matrix ud2) {
    try {
        solution Xopt;


        return Xopt;
    } catch (string ex_info) {
        throw ("solution HJ(...):\n" + ex_info);
    }
}

solution HJ_trial(matrix (*ff)(matrix, matrix, matrix), solution XB, double s, matrix ud1, matrix ud2) {
    try {
        return XB;
    } catch (string ex_info) {
        throw ("solution HJ_trial(...):\n" + ex_info);
    }
}

solution Rosen(matrix (*ff)(matrix, matrix, matrix), matrix x0, matrix s0, double alpha, double beta, double epsilon,
               int Nmax, matrix ud1, matrix ud2) {
    try {
        solution Xopt;


        return Xopt;
    } catch (string ex_info) {
        throw ("solution Rosen(...):\n" + ex_info);
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

solution sym_NM(matrix (*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double beta, double gamma,
                double delta, double epsilon, int Nmax, matrix ud1, matrix ud2) {
    try {
        solution Xopt;


        return Xopt;
    } catch (string ex_info) {
        throw ("solution sym_NM(...):\n" + ex_info);
    }
}

solution SD(matrix (*ff)(matrix, matrix, matrix), matrix (*gf)(matrix, matrix, matrix), matrix x0, double h0,
            double epsilon, int Nmax, matrix ud1, matrix ud2) {
    try {
        solution Xopt;


        return Xopt;
    } catch (string ex_info) {
        throw ("solution SD(...):\n" + ex_info);
    }
}

solution CG(matrix (*ff)(matrix, matrix, matrix), matrix (*gf)(matrix, matrix, matrix), matrix x0, double h0,
            double epsilon, int Nmax, matrix ud1, matrix ud2) {
    try {
        solution Xopt;


        return Xopt;
    } catch (string ex_info) {
        throw ("solution CG(...):\n" + ex_info);
    }
}

solution Newton(matrix (*ff)(matrix, matrix, matrix), matrix (*gf)(matrix, matrix, matrix),
                matrix (*Hf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1,
                matrix ud2) {
    try {
        solution Xopt;


        return Xopt;
    } catch (string ex_info) {
        throw ("solution Newton(...):\n" + ex_info);
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

solution Powell(matrix (*ff)(matrix, matrix, matrix), matrix x0, double epsilon, int Nmax, matrix ud1, matrix ud2) {
    try {
        solution Xopt;


        return Xopt;
    } catch (string ex_info) {
        throw ("solution Powell(...):\n" + ex_info);
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
