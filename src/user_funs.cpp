#include "../include/user_funs.h"
#include <cmath>
#include <vector>

#define _USE_MATH_DEFINES
// --- Lab 0 Functions ---

// Objective function for the test case
matrix ff0T(matrix x, matrix ud1, matrix ud2) {
    matrix y; // y contains the value of the objective function
    y = pow(x(0) - ud1(0), 2) + pow(x(1) - ud1(1), 2); // ud1 contains the coordinates of the sought optimum
    return y;
}

// Objective function for the real-world problem (pendulum)
matrix ff0R(matrix x, matrix ud1, matrix ud2) {
    matrix y; // y contains the value of the objective function
    matrix Y0 = matrix(2, 1), // Y0 contains initial conditions
            MT = matrix(2, new double[2]{m2d(x), 0.5}); // MT contains the torque on the pendulum and its duration
    matrix *Y = solve_ode(df0, 0, 0.1, 10, Y0, ud1, MT); // solve the differential equation
    int n = get_len(Y[0]); // length of the solution
    double teta_max = Y[1](0, 0); // find the maximum pendulum deflection
    for (int i = 1; i < n; ++i)
        if (teta_max < Y[1](i, 0))
            teta_max = Y[1](i, 0);
    y = abs(teta_max - m2d(ud1)); // value of the objective function (ud1 is the assumed maximum deflection)

    // Correctly free memory to prevent leaks
    delete[] Y;

    return y;
}

// Differential equations for the pendulum model
matrix df0(double t, matrix Y, matrix ud1, matrix ud2) {
    matrix dY(2, 1); // define the vector of derivatives of the sought functions
    double m = 1, l = 0.5, b = 0.5, g = 9.81; // define model parameters
    double I = m * pow(l, 2);
    dY(0) = Y(1); // derivative of position is velocity
    dY(1) = ((t <= ud2(1)) * ud2(0) - m * g * l * sin(Y(0)) - b * Y(1)) / I; // derivative of velocity is acceleration
    return dY;
}


// --- Lab 1 Functions ---

// Test objective function for Lab 1
matrix ff1T(matrix x, matrix ud1, matrix ud2) {
    double val = x(0); // Extract the value from the matrix

    // f(x) = -cos(0.1x) * exp(-(0.1x - 2*PI)^2) + 0.002 * (0.1x)^2
    double term1 = -cos(0.1 * val) * exp(-pow(0.1 * val - 2 * M_PI, 2));
    double term2 = 0.002 * pow(0.1 * val, 2);

    matrix y(term1 + term2); // The result must be a matrix
    return y;
}


// Objective function for the real-world problem
matrix ff1R(matrix x, matrix ud1, matrix ud2) {
    // Convert the cross-section area DA from cm^2 to m^2
    // x(0) <-> 50 cm^2
    double DA = x(0) * 1e-4;

    // Initial conditions: [Volume in A, Volume in B, Temperature in B]
    matrix Y0(3, 1);
    Y0(0) = 5.0; // VA_start = 5 m^3
    Y0(1) = 1.0; // VB_start = 1 m^3
    Y0(2) = 20.0; // TB_start = 20 C

    // Simulation and physical parameters
    matrix params(1, 1);
    params(0) = DA; // Pass DA to the differential function

    // Solve the ordinary differential equation
    matrix *Y = solve_ode(dff1R, 0, 1, 2000, Y0, NAN, params);

    // Find the maximum temperature in tank B
    int n = get_len(Y[0]);
    double T_max = 0;
    for (int i = 0; i < n; ++i) {
        if (Y[1](i, 2) > T_max) {
            T_max = Y[1](i, 2);
        }
    }

    // Free the dynamically allocated memory
    delete[] Y;

    // Objective function: minimize the difference from 50°C
    return abs(T_max - 50.0);
    // -50.0
}

// Differential equations for the real-world problem
// Differential equations for the real-world problem
matrix dff1R(double t, matrix Y, matrix ud1, matrix ud2) {
    // Physical parameters and constants
    double PA = 2.0; // Base area of tank A [m^2]
    double PB = 1.0; // Base area of tank B [m^2]
    double g = 9.81; // Gravitational acceleration [m/s^2]
    double a = 0.98; // Coefficient for fluid viscosity
    double b = 0.63; // Coefficient for stream constriction
    double TA_in = 95.0; // Temperature of water from tank A [C]
    double TB_in_temp = 20.0; // Temperature of external inflow to B [C]
    double Fin_B = 10.0 / 1000.0; // External inflow to B [10 l/s -> m^3/s]
    double DB = 36.5665 * 1e-4; // Outflow area from B [cm^2 -> m^2]

    double DA = ud2(0); // Cross-section area DA from parameters [m^2]

    // Current states
    double VA = Y(0);
    double VB = Y(1);
    double TB = Y(2);

    // Calculate water column heights
    double hA = VA / PA;
    double hB = VB / PB;

    // Flows based on Torricelli's law
    double Fout_A = (hA > 0) ? a * b * DA * sqrt(2 * g * hA) : 0;
    double Fout_B = (hB > 0) ? a * b * DB * sqrt(2 * g * hB) : 0;

    // Differential equations
    matrix dY(3, 1);


    dY(0) = -Fout_A; // Change in volume in A (water flows OUT)
    dY(1) = Fout_A + Fin_B - Fout_B; // Change in volume in B (water from A flows IN)
    // -----------------------

    if (abs(VB) < 1e-9) {
        // Avoid division by zero if tank B is empty
        dY(2) = 0;
    } else {
        dY(2) = (Fout_A * (TA_in - TB) + Fin_B * (TB_in_temp - TB)) / VB; // Change in temperature in B
    }

    return dY;
}


matrix ff2R(matrix x, matrix ud1, matrix ud2) {
    // Zmienne optymalizacji to współczynniki wzmocnienia regulatora
    double k1 = x(0);
    double k2 = x(1);

    // Stałe modelu i symulacji
    double b = 0.25;
    double mr = 1.0;
    double mc = 5.0;
    double l = 2.0;
    double I = (1.0/3.0) * mr * pow(l, 2) + mc * pow(l, 2);
    double aref = M_PI;
    double wref = 0;

    double t0 = 0.0;
    double tend = 100.0;
    double dt = 0.1;

    // Warunki początkowe dla równania różniczkowego: [alpha(0), omega(0)]
    matrix Y0(2, 1);
    Y0(0) = 0.0; // Początkowy kąt
    Y0(1) = 0.0; // Początkowa prędkość kątowa

    // Przekazanie współczynników k1 i k2 do funkcji różniczkowej przez ud2
    matrix ode_params(2, 1);
    ode_params(0) = k1;
    ode_params(1) = k2;

    // Rozwiązanie równania różniczkowego
    matrix* Y = solve_ode(df2R, t0, dt, tend, Y0, NAN, ode_params);

    // Obliczenie funkcjonału jakości Q(k1, k2) metodą prostokątów
    double Q = 0.0;
    int num_steps = get_len(Y[0]);

    for (int i = 0; i < num_steps; ++i) {
        double current_alpha = Y[1](i, 0);
        double current_omega = Y[1](i, 1);

        // Obliczenie momentu siły M(t) w danym kroku
        double M = k1 * (aref - current_alpha) + k2 * (wref - current_omega);

        // Obliczenie wartości funkcji podcałkowej
        double integrand = 10 * pow(aref - current_alpha, 2) + pow(wref - current_omega, 2) + pow(M, 2);

        // Dodanie do sumy całkowej
        Q += integrand * dt;
    }

    // Zwolnienie pamięci zaalokowanej przez solve_ode
    delete[] Y;

    // Funkcja celu zwraca obliczoną wartość Q
    return matrix(Q);
}


matrix ff2T(matrix x, matrix ud1, matrix ud2) {
    double x1 = x(0);

    double x2 = x(1);

    double val = pow(x1, 2) + pow(x2, 2) - cos(2.5 * M_PI * x1) - cos(2.5 * M_PI * x2) + 2;

    return matrix(val);
}

matrix df2R(double t, matrix Y, matrix ud1, matrix ud2)
{
    // Model constants
    double b = 0.25;  // Friction coefficient
    double mr = 1.0;  // Mass of the arm
    double mc = 5.0;  // Mass of the weight
    double l = 2.0;   // Length of the arm
    double I = (1.0/3.0) * mr * pow(l, 2) + mc * pow(l, 2); // Moment of inertia

    // Control parameters (k1, k2) are passed via ud2
    double k1 = ud2(0);
    double k2 = ud2(1);

    // Target values
    double aref = M_PI;
    double wref = 0;

    // Current state from Y vector
    double alpha = Y(0); // Current angle
    double omega = Y(1); // Current angular velocity

    // Calculate the torque M(t)
    double M = k1 * (aref - alpha) + k2 * (wref - omega);

    // Define the system of 2 first-order ODEs
    matrix dY(2, 1);
    dY(0) = omega; // d(alpha)/dt = omega
    dY(1) = (M - b * omega) / I; // d(omega)/dt = (M - b*omega)/I

    return dY;
}

// Testowa funkcja celu dla lab3
matrix ff3T(matrix x, matrix ud1, matrix ud2)
{
    double x1 = x(0);
    double x2 = x(1);

    double term_in_sqrt = pow(x1 / M_PI, 2) + pow(x2 / M_PI, 2);
    double term = M_PI * sqrt(term_in_sqrt);

    // Warunek zabezpieczający przed dzieleniem przez zero, gdy x1=0 i x2=0
    // Zgodnie z granicą sin(z)/z -> 1 dla z -> 0
    if (abs(term) < 1e-9) {
        return matrix(1.0);
    }

    return matrix(sin(term) / term);
}


// Funkcja celu z zewnętrzną funkcją kary
matrix ff3T_zew(matrix x, matrix ud1, matrix ud2)
{
    // ud1(0) przechowuje parametr 'a'
    // ud1(1) przechowuje współczynnik kary 'c'
    double a = ud1(0);
    double c = ud1(1);

    // Wartość oryginalnej funkcji celu
    matrix y_f = ff3T(x);

    // Ograniczenia g(x) <= 0
    double g1 = -x(0) + 1.0;
    double g2 = -x(1) + 1.0;
    double g3 = sqrt(pow(x(0), 2) + pow(x(1), 2)) - a;

    // Obliczenie kary
    double kara = c * (pow(max(0.0, g1), 2) + pow(max(0.0, g2), 2) + pow(max(0.0, g3), 2));

    return y_f + kara;
}

// Funkcja celu z wewnętrzną funkcją kary
matrix ff3T_wew(matrix x, matrix ud1, matrix ud2)
{
    // ud1(0) przechowuje parametr 'a'
    // ud1(1) przechowuje współczynnik kary 'c'
    double a = ud1(0);
    double c = ud1(1);

    // Ograniczenia g(x) <= 0
    double g1 = -x(0) + 1.0;
    double g2 = -x(1) + 1.0;
    double g3 = sqrt(pow(x(0), 2) + pow(x(1), 2)) - a;

    // Sprawdzenie, czy punkt leży w obszarze dopuszczalnym
    // Jeśli nie, zwróć bardzo dużą wartość
    if (g1 >= 0 || g2 >= 0 || g3 >= 0) {
        return matrix(1e100); // Wartość "nieskończona"
    }

    // Wartość oryginalnej funkcji celu
    matrix y_f = ff3T(x);

    // Obliczenie kary
    double kara = -c * (1.0 / g1 + 1.0 / g2 + 1.0 / g3);

    return y_f + kara;
}

matrix dff3R(double t, matrix Y, matrix ud1, matrix ud2) {
    // Y(0) = x, Y(1) = vx, Y(2) = y, Y(3) = vy
    double x = Y(0);
    double vx = Y(1);
    double y = Y(2);
    double vy = Y(3);

    // Parametry fizyczne
    double m = 0.6;   // kg
    double r = 0.12;  // m
    double C = 0.47;
    double rho = 1.2; // kg/m^3
    double g = 9.81;  // m/s^2
    double S = M_PI * r * r;

    // Omega jest przekazywana w ud1 (parametr sterujący, ale stały w czasie symulacji)
    double omega = ud1(0);

    // Siły oporu (zgodnie z instrukcją: 0.5 * C * rho * S * v * |v|)
    double Dx = 0.5 * C * rho * S * vx * abs(vx);
    double Dy = 0.5 * C * rho * S * vy * abs(vy);

    // Siły Magnusa
    double FMx = rho * vy * omega * M_PI * pow(r, 3);
    double FMy = rho * vx * omega * M_PI * pow(r, 3);

    // Równania ruchu
    matrix dY(4, 1);
    dY(0) = vx;
    dY(1) = -(Dx + FMx) / m;
    dY(2) = vy;
    dY(3) = -(Dy + FMy) / m - g; // Grawitacja działa w dół (przeciwnie do osi Y)

    return dY;
}


matrix ff3R_zew(matrix x, matrix ud1, matrix ud2) {
    // Zmienne decyzyjne: x(0) = v0x, x(1) = omega
    double v0x = x(0);
    double omega = x(1);

    // Parametry kary przekazywane w ud2 (ud2(0)=c)
    double c = ud2(0);

    // Warunki początkowe symulacji
    // [x=0, vx=v0x, y=100, vy=0]
    matrix Y0(4, 1);
    Y0(0) = 0.0;
    Y0(1) = v0x;
    Y0(2) = 100.0;
    Y0(3) = 0.0;

    // Przekazujemy omegę do równania różniczkowego
    matrix params(1, 1);
    params(0) = omega;

    // Symulacja: t_end = 7s, dt = 0.01s
    matrix* Y = solve_ode(dff3R, 0, 0.01, 7, Y0, params, NAN);

    // Analiza wyników symulacji
    int n = get_len(Y[0]);
    double x_end = 0.0; // Wartość x przy uderzeniu w ziemię (y=0)
    double x_at_50 = 0.0; // Wartość x przy y=50m
    bool hit_ground = false;
    bool passed_50 = false;

    // Szukamy momentu przejścia przez y=50 i uderzenia w ziemię
    for (int i = 0; i < n; ++i) {
        double current_y = Y[1](i, 2);
        double current_x = Y[1](i, 0);

        // Sprawdzenie przejścia przez 50m
        if (!passed_50 && current_y <= 50.0) {
            x_at_50 = current_x;
            passed_50 = true;
        }

        // Sprawdzenie uderzenia w ziemię
        if (!hit_ground && current_y <= 0.0) {
            x_end = current_x;
            hit_ground = true;
            // Możemy przerwać pętlę, jeśli interesuje nas tylko pierwsze uderzenie
             // break; // Opcjonalnie, ale solver i tak policzył całość
        }
    }

    // Jeśli nie uderzyła w ziemię w ciągu 7s, bierzemy ostatni x (jako fallback)
    if (!hit_ground) x_end = Y[1](n - 1, 0);
    // Jeśli nie minęła 50m (dziwne przy spadku ze 100m), bierzemy start
    if (!passed_50) x_at_50 = 0.0;

    delete[] Y;

    // --- Funkcja celu ---
    // Maksymalizacja x_end => Minimalizacja -x_end
    double f_val = -x_end;

    // --- Ograniczenia ---
    // 1. v0x w [-10, 10] => |v0x| - 10 <= 0
    double g1 = abs(v0x) - 10.0;

    // 2. omega w [-10, 10] => |omega| - 10 <= 0
    double g2 = abs(omega) - 10.0;

    // 3. Przy y=50m, x musi być w [3, 7] => odległość od 5 <= 2 => |x_at_50 - 5| - 2 <= 0
    double g3 = abs(x_at_50 - 5.0) - 2.0;

    // Kara zewnętrzna: c * suma(max(0, g_i)^2)
    double penalty = c * (pow(max(0.0, g1), 2) + pow(max(0.0, g2), 2) + pow(max(0.0, g3), 2));

    return matrix(f_val + penalty);
}

// --- IMPLEMENTACJA LAB 4 ---

// Testowa funkcja celu: f(x) = 1/6*x1^6 - 1.05*x1^4 + 2*x1^2 + x2^2 + x1*x2
matrix ff4T(matrix x, matrix ud1, matrix ud2) {
    double x1 = x(0);
    double x2 = x(1);
    double y = (1.0/6.0)*pow(x1, 6) - 1.05*pow(x1, 4) + 2.0*pow(x1, 2) + pow(x2, 2) + x1*x2;
    return matrix(y);
}

// Gradient testowej funkcji celu
matrix gf4T(matrix x, matrix ud1, matrix ud2) {
    double x1 = x(0);
    double x2 = x(1);
    matrix g(2, 1);
    // df/dx1 = x1^5 - 4.2*x1^3 + 4*x1 + x2
    g(0) = pow(x1, 5) - 4.2*pow(x1, 3) + 4.0*x1 + x2;
    // df/dx2 = 2*x2 + x1
    g(1) = 2.0*x2 + x1;
    return g;
}

// Hesjan testowej funkcji celu
matrix Hf4T(matrix x, matrix ud1, matrix ud2) {
    double x1 = x(0);
    double x2 = x(1);
    matrix H(2, 2);
    // d2f/dx1^2 = 5*x1^4 - 12.6*x1^2 + 4
    H(0, 0) = 5.0*pow(x1, 4) - 12.6*pow(x1, 2) + 4.0;
    // d2f/dx1dx2 = 1
    H(0, 1) = 1.0;
    // d2f/dx2dx1 = 1
    H(1, 0) = 1.0;
    // d2f/dx2^2 = 2
    H(1, 1) = 2.0;
    return H;
}

// Funkcja pomocnicza: sigmoid
double sigmoid(double z) {
    if (z > 20.0) return 1.0;
    if (z < -20.0) return 0.0;
    return 1.0 / (1.0 + exp(-z));
}

// Rzeczywista funkcja celu: Log-Loss
// theta: wektor parametrów (3x1)
// ud1: X (3x100) - macierz cech
// ud2: Y (1x100) - wektor etykiet
matrix ff4R(matrix theta, matrix ud1, matrix ud2) {
    matrix X = ud1;
    matrix Y = ud2;
    
    // ZMIANA: Zamiast get_len(trans(Y)), używamy get_size(Y)
    int* size = get_size(Y);
    int m = size[1]; // Y jest 1xm, więc liczba próbek to drugi wymiar
    delete[] size;
    
    double J = 0.0;
    for (int i = 0; i < m; ++i) {
        matrix xi = get_col(X, i);       // i-ta kolumna X (3x1)
        double yi = Y(0, i);             // i-ta etykieta
        
        double z = m2d(trans(theta) * xi);
        double h = sigmoid(z);
        
        // Zabezpieczenie logarytmu
        double eps = 1e-15;
        if (h < eps) h = eps;
        if (h > 1.0 - eps) h = 1.0 - eps;
        
        J += yi * log(h) + (1.0 - yi) * log(1.0 - h);
    }
    
    return matrix(-J / m);
}

// Gradient rzeczywistej funkcji celu
matrix gf4R(matrix theta, matrix ud1, matrix ud2) {
    matrix X = ud1;
    matrix Y = ud2;
    
    // ZMIANA: Zamiast get_len(trans(Y)), używamy get_size(Y)
    int* size = get_size(Y);
    int m = size[1]; // Y jest 1xm
    delete[] size;
    
    matrix grad(3, 1); // Wektor zerowy 3x1
    
    for (int i = 0; i < m; ++i) {
        matrix xi = get_col(X, i);
        double yi = Y(0, i);
        
        double z = m2d(trans(theta) * xi);
        double h = sigmoid(z);
        
        // Sumowanie: (h(xi) - yi) * xi
        grad = grad + (h - yi) * xi;
    }
    
    return grad / m;
}
