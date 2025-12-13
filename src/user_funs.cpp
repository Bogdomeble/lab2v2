// src/user_funs.cpp
#include "../include/user_funs.h"
#include <cmath>
#include <vector>
#include <fstream>
#include <sstream>

#define _USE_MATH_DEFINES

// --- Lab 0 Functions ---
matrix ff0T(matrix x, matrix ud1, matrix ud2) {
    matrix y;
    y = pow(x(0) - ud1(0), 2) + pow(x(1) - ud1(1), 2);
    return y;
}

matrix ff0R(matrix x, matrix ud1, matrix ud2) {
    matrix y;
    matrix Y0 = matrix(2, 1),
            MT = matrix(2, new double[2]{m2d(x), 0.5});
    matrix *Y = solve_ode(df0, 0, 0.1, 10, Y0, ud1, MT);
    int n = get_len(Y[0]);
    double teta_max = Y[1](0, 0);
    for (int i = 1; i < n; ++i)
        if (teta_max < Y[1](i, 0))
            teta_max = Y[1](i, 0);
    y = abs(teta_max - m2d(ud1));
    delete[] Y;
    return y;
}

matrix df0(double t, matrix Y, matrix ud1, matrix ud2) {
    matrix dY(2, 1);
    double m = 1, l = 0.5, b = 0.5, g = 9.81;
    double I = m * pow(l, 2);
    dY(0) = Y(1);
    dY(1) = ((t <= ud2(1)) * ud2(0) - m * g * l * sin(Y(0)) - b * Y(1)) / I;
    return dY;
}

// --- Lab 1 Functions ---
matrix ff1T(matrix x, matrix ud1, matrix ud2) {
    double val = x(0);
    double term1 = -cos(0.1 * val) * exp(-pow(0.1 * val - 2 * M_PI, 2));
    double term2 = 0.002 * pow(0.1 * val, 2);
    matrix y(term1 + term2);
    return y;
}

matrix ff1R(matrix x, matrix ud1, matrix ud2) {
    double DA = x(0) * 1e-4;
    matrix Y0(3, 1);
    Y0(0) = 5.0; Y0(1) = 1.0; Y0(2) = 20.0;
    matrix params(1, 1);
    params(0) = DA;
    matrix *Y = solve_ode(dff1R, 0, 1, 2000, Y0, NAN, params);
    int n = get_len(Y[0]);
    double T_max = 0;
    for (int i = 0; i < n; ++i) {
        if (Y[1](i, 2) > T_max) {
            T_max = Y[1](i, 2);
        }
    }
    delete[] Y;
    return abs(T_max - 50.0);
}

matrix dff1R(double t, matrix Y, matrix ud1, matrix ud2) {
    double PA = 2.0; double PB = 1.0; double g = 9.81;
    double a = 0.98; double b = 0.63;
    double TA_in = 95.0; double TB_in_temp = 20.0;
    double Fin_B = 10.0 / 1000.0;
    double DB = 36.5665 * 1e-4;
    double DA = ud2(0);
    double VA = Y(0); double VB = Y(1); double TB = Y(2);
    double hA = VA / PA; double hB = VB / PB;
    double Fout_A = (hA > 0) ? a * b * DA * sqrt(2 * g * hA) : 0;
    double Fout_B = (hB > 0) ? a * b * DB * sqrt(2 * g * hB) : 0;
    matrix dY(3, 1);
    dY(0) = -Fout_A;
    dY(1) = Fout_A + Fin_B - Fout_B;
    if (abs(VB) < 1e-9) dY(2) = 0;
    else dY(2) = (Fout_A * (TA_in - TB) + Fin_B * (TB_in_temp - TB)) / VB;
    return dY;
}

matrix ff2R(matrix x, matrix ud1, matrix ud2) {
    double k1 = x(0); double k2 = x(1);
    // Usunieto nieuzywane zmienne b, I
    double aref = M_PI; double wref = 0;
    double t0 = 0.0; double tend = 100.0; double dt = 0.1;
    matrix Y0(2, 1); Y0(0) = 0.0; Y0(1) = 0.0;
    matrix ode_params(2, 1); ode_params(0) = k1; ode_params(1) = k2;
    matrix* Y = solve_ode(df2R, t0, dt, tend, Y0, NAN, ode_params);
    double Q = 0.0;
    int num_steps = get_len(Y[0]);
    for (int i = 0; i < num_steps; ++i) {
        double current_alpha = Y[1](i, 0);
        double current_omega = Y[1](i, 1);
        double M = k1 * (aref - current_alpha) + k2 * (wref - current_omega);
        double integrand = 10 * pow(aref - current_alpha, 2) + pow(wref - current_omega, 2) + pow(M, 2);
        Q += integrand * dt;
    }
    delete[] Y;
    return matrix(Q);
}

matrix ff2T(matrix x, matrix ud1, matrix ud2) {
    double x1 = x(0); double x2 = x(1);
    double val = pow(x1, 2) + pow(x2, 2) - cos(2.5 * M_PI * x1) - cos(2.5 * M_PI * x2) + 2;
    return matrix(val);
}

matrix df2R(double t, matrix Y, matrix ud1, matrix ud2) {
    double b = 0.25; double mr = 1.0; double mc = 5.0; double l = 2.0;
    double I = (1.0/3.0) * mr * pow(l, 2) + mc * pow(l, 2);
    double k1 = ud2(0); double k2 = ud2(1);
    double aref = M_PI; double wref = 0;
    double alpha = Y(0); double omega = Y(1);
    double M = k1 * (aref - alpha) + k2 * (wref - omega);
    matrix dY(2, 1);
    dY(0) = omega;
    dY(1) = (M - b * omega) / I;
    return dY;
}

// --- Lab 3 Functions ---
matrix ff3T(matrix x, matrix ud1, matrix ud2) {
    double x1 = x(0); double x2 = x(1);
    double term_in_sqrt = pow(x1 / M_PI, 2) + pow(x2 / M_PI, 2);
    double term = M_PI * sqrt(term_in_sqrt);
    if (abs(term) < 1e-9) return matrix(1.0);
    return matrix(sin(term) / term);
}

matrix ff3T_zew(matrix x, matrix ud1, matrix ud2) {
    double a = ud1(0); double c = ud1(1);
    matrix y_f = ff3T(x);
    double g1 = -x(0) + 1.0;
    double g2 = -x(1) + 1.0;
    double g3 = sqrt(pow(x(0), 2) + pow(x(1), 2)) - a;
    double kara = c * (pow(max(0.0, g1), 2) + pow(max(0.0, g2), 2) + pow(max(0.0, g3), 2));
    return y_f + kara;
}

matrix ff3T_wew(matrix x, matrix ud1, matrix ud2) {
    double a = ud1(0); double c = ud1(1);
    double g1 = -x(0) + 1.0;
    double g2 = -x(1) + 1.0;
    double g3 = sqrt(pow(x(0), 2) + pow(x(1), 2)) - a;
    if (g1 >= 0 || g2 >= 0 || g3 >= 0) return matrix(1e100);
    matrix y_f = ff3T(x);
    double kara = -c * (1.0 / g1 + 1.0 / g2 + 1.0 / g3);
    return y_f + kara;
}

matrix dff3R(double t, matrix Y, matrix ud1, matrix ud2) {
    // Usunieto nieuzywane zmienne x, y
    double vx = Y(1); double vy = Y(3);
    double m = 0.6; double r = 0.12; double C = 0.47; double rho = 1.2; double g = 9.81;
    double S = M_PI * r * r;
    double omega = ud1(0);
    double Dx = 0.5 * C * rho * S * vx * abs(vx);
    double Dy = 0.5 * C * rho * S * vy * abs(vy);
    double FMx = rho * vy * omega * M_PI * pow(r, 3);
    double FMy = rho * vx * omega * M_PI * pow(r, 3);
    matrix dY(4, 1);
    dY(0) = vx;
    dY(1) = -(Dx + FMx) / m;
    dY(2) = vy;
    dY(3) = -(Dy + FMy) / m - g;
    return dY;
}

matrix ff3R_zew(matrix x, matrix ud1, matrix ud2) {
    double v0x = x(0); double omega = x(1);
    double c = ud2(0);
    matrix Y0(4, 1); Y0(0) = 0.0; Y0(1) = v0x; Y0(2) = 100.0; Y0(3) = 0.0;
    matrix params(1, 1); params(0) = omega;
    matrix* Y = solve_ode(dff3R, 0, 0.01, 7, Y0, params, NAN);
    int n = get_len(Y[0]);
    double x_end = 0.0; double x_at_50 = 0.0;
    bool hit_ground = false; bool passed_50 = false;
    for (int i = 0; i < n; ++i) {
        double current_y = Y[1](i, 2);
        double current_x = Y[1](i, 0);
        if (!passed_50 && current_y <= 50.0) {
            x_at_50 = current_x;
            passed_50 = true;
        }
        if (!hit_ground && current_y <= 0.0) {
            x_end = current_x;
            hit_ground = true;
        }
    }
    if (!hit_ground) x_end = Y[1](n - 1, 0);
    if (!passed_50) x_at_50 = 0.0;
    delete[] Y;
    double f_val = -x_end;
    double g1 = abs(v0x) - 10.0;
    double g2 = abs(omega) - 10.0;
    double g3 = abs(x_at_50 - 5.0) - 2.0;
    double penalty = c * (pow(max(0.0, g1), 2) + pow(max(0.0, g2), 2) + pow(max(0.0, g3), 2));
    return matrix(f_val + penalty);
}

// --- LAB 4 ---

// Implementacja funkcji pomocniczych
matrix read_matrix_semicolon(string path, int rows, int cols) {
    ifstream file(path);
    if (!file.is_open()) throw string("Nie mozna otworzyc pliku: " + path);
    matrix M(rows, cols);
    string line, val_str;
    for (int i = 0; i < rows; ++i) {
        if (!getline(file, line)) break;
        stringstream ss(line);
        for (int j = 0; j < cols; ++j) {
            if (!getline(ss, val_str, ';')) break;
            try { M(i, j) = stod(val_str); } catch (...) { M(i, j) = 0.0; }
        }
    }
    return M;
}

double sigmoid(double z);

double calculate_accuracy(matrix theta, matrix X, matrix Y) {
    int m = get_size(Y)[1];
    int correct = 0;
    for (int i = 0; i < m; ++i) {
        matrix xi = get_col(X, i);
        double z = m2d(trans(theta) * xi);
        double h = sigmoid(z);
        int prediction = (h >= 0.5) ? 1 : 0;
        if (prediction == (int)Y(0, i)) correct++;
    }
    return (double)correct / m * 100.0;
}

matrix ff4T(matrix x, matrix ud1, matrix ud2) {
    double x1 = x(0); double x2 = x(1);
    double y = (1.0/6.0)*pow(x1, 6) - 1.05*pow(x1, 4) + 2.0*pow(x1, 2) + pow(x2, 2) + x1*x2;
    return matrix(y);
}

matrix gf4T(matrix x, matrix ud1, matrix ud2) {
    double x1 = x(0); double x2 = x(1);
    matrix g(2, 1);
    g(0) = pow(x1, 5) - 4.2*pow(x1, 3) + 4.0*x1 + x2;
    g(1) = 2.0*x2 + x1;
    return g;
}

matrix Hf4T(matrix x, matrix ud1, matrix ud2) {
    double x1 = x(0);
    // Usunieto nieuzywana zmienna x2
    matrix H(2, 2);
    H(0, 0) = 5.0*pow(x1, 4) - 12.6*pow(x1, 2) + 4.0;
    H(0, 1) = 1.0;
    H(1, 0) = 1.0;
    H(1, 1) = 2.0;
    return H;
}

double sigmoid(double z) {
    if (z > 20.0) return 1.0;
    if (z < -20.0) return 0.0;
    return 1.0 / (1.0 + exp(-z));
}

matrix ff4R(matrix theta, matrix ud1, matrix ud2) {
    matrix X = ud1; matrix Y = ud2;
    int* size = get_size(Y);
    int m = size[1];
    delete[] size;
    double J = 0.0;
    for (int i = 0; i < m; ++i) {
        matrix xi = get_col(X, i);
        double yi = Y(0, i);
        double z = m2d(trans(theta) * xi);
        double h = sigmoid(z);
        double eps = 1e-15;
        if (h < eps) h = eps;
        if (h > 1.0 - eps) h = 1.0 - eps;
        J += yi * log(h) + (1.0 - yi) * log(1.0 - h);
    }
    return matrix(-J / m);
}

matrix gf4R(matrix theta, matrix ud1, matrix ud2) {
    matrix X = ud1; matrix Y = ud2;
    int* size = get_size(Y);
    int m = size[1];
    delete[] size;
    matrix grad(3, 1);
    for (int i = 0; i < m; ++i) {
        matrix xi = get_col(X, i);
        double yi = Y(0, i);
        double z = m2d(trans(theta) * xi);
        double h = sigmoid(z);
        grad = grad + (h - yi) * xi;
    }
    return grad / m;
}

// --- Lab 5 Functions ---
matrix ff5T(matrix x, matrix ud1, matrix ud2) {
    double w = ud1(0); double a = ud2(0);
    if (std::isnan(x(0)) || std::isnan(x(1))) return matrix(1e100);
    double x1 = x(0); double x2 = x(1);
    double f1 = a * (pow(x1 - 3.0, 2) + pow(x2 - 3.0, 2));
    double f2 = (1.0 / a) * (pow(x1 + 3.0, 2) + pow(x2 + 3.0, 2));
    return matrix(w * f1 + (1.0 - w) * f2);
}

matrix ff5R(matrix x, matrix ud1, matrix ud2) {
    matrix y;
    double w = ud1(0);
    double l = x(0); double d = x(1);
    double P = 2000.0; double E = 1.2e11; double ro = 8920.0;
    double u_max = 0.0025; double sigma_max = 300e6;
    double mass = ro * l * M_PI * pow(d / 2.0, 2);
    double deflection = (64.0 * P * pow(l, 3)) / (3.0 * E * M_PI * pow(d, 4));
    double sigma = (32.0 * P * l) / (M_PI * pow(d, 3));
    double penalty = 0.0;
    double c = 1e12;
    double g1 = 0.2 - l;
    double g2 = l - 1.0;
    double g3 = 0.01 - d;
    double g4 = d - 0.05;
    // Usunieto nieuzywane zmienne g5, g6 (uzywamy wersji znormalizowanej)
    double g5_norm = (deflection / u_max) - 1.0;
    double g6_norm = (sigma / sigma_max) - 1.0;
    if (g1 > 0) penalty += c * pow(g1, 2);
    if (g2 > 0) penalty += c * pow(g2, 2);
    if (g3 > 0) penalty += c * pow(g3, 2);
    if (g4 > 0) penalty += c * pow(g4, 2);
    double c_phys = 1e6;
    if (g5_norm > 0) penalty += c_phys * pow(g5_norm, 2);
    if (g6_norm > 0) penalty += c_phys * pow(g6_norm, 2);
    y = w * mass + (1.0 - w) * deflection + penalty;
    return y;
}
