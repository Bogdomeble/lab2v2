

#include "../include/opt_alg.h"
#include <cstdlib>
#include <ctime>

void lab0();
void lab1();
void lab1_real();
void test_real_problem_DA50();
void lab2();
void lab3();
void lab4();
void lab5();
void lab6();

int main() {
    try {
        // Call the function for the test problem or the real problem
        //lab1();
        //test_real_problem_DA50();
        //lab1_real();
         //lab2();
         //lab3();
         lab4();
    } catch (string EX_INFO) {
        cerr << "ERROR:\n";
        cerr << EX_INFO << endl << endl;
    }
    return 0;
}

void lab0() {
    // Test function
    double epsilon = 1e-2; // precision
    int Nmax = 10000;      // maximum number of objective function calls
    matrix lb(2, 1, -5), ub(2, 1, 5), // lower and upper bounds
            a(2, 1);       // exact optimal solution
    solution opt;          // optimal solution found by the algorithm
    a(0) = -1;
    a(1) = 2;
    opt = MC(ff0T, 2, lb, ub, epsilon, Nmax, a); // call the optimization procedure
    cout << opt << endl << endl;                // print the result
    solution::clear_calls();                    // reset counters

    // Pendulum
    Nmax = 1000;
    epsilon = 1e-5;
    lb = 0, ub = 5;
    double teta_opt = 1; // maximum pendulum deflection
    opt = MC(ff0R, 1, lb, ub, epsilon, Nmax, teta_opt); // call the optimization procedure
    cout << opt << endl << endl;                       // print the result
    solution::clear_calls();                           // reset counters

    // Save simulation to a CSV file
    matrix Y0 = matrix(2, 1), // Y0 contains initial conditions
            MT = matrix(2, new double[2]{m2d(opt.x), 0.5}); // MT contains the torque on the pendulum and its duration
    matrix *Y = solve_ode(df0, 0, 0.1, 10, Y0, NAN, MT); // solve the differential equation
    ofstream Sout(resolvePath("symulacja_lab0.csv"));                // define a stream to the .csv file
    Sout << hcat(Y[0], Y[1]);                            // save results in the file
    Sout.close();                                       // close the stream
    delete[] Y;                                         // free the memory of the DE solution
}

void lab1() {
    // Initialize the random number generator to get different values each time
    srand(time(nullptr));

    // Open a CSV file to save the results.
    // The file will be created in the same folder as the executable.
    ofstream file(resolvePath("lab1_results.csv"));

    // Check if the file was opened correctly
    if (!file.is_open()) {
        cerr << "ERROR: Cannot open file lab1_results.csv for writing!" << endl;
        return;
    }

    // Write the header to the CSV file
    file << "Iteration,x0_start,"
         << "Exp_start,Exp_end,Exp_f_calls,"
         << "Fib_with_Exp_x,Fib_with_Exp_y,Fib_with_Exp_f_calls,"
         << "Lag_with_Exp_x,Lag_with_Exp_y,Lag_with_Exp_f_calls,"
         << "Fib_no_Exp_x,Fib_no_Exp_y,Fib_no_Exp_f_calls,"
         << "Lag_no_Exp_x,Lag_no_Exp_y,Lag_no_Exp_f_calls\n";

    // Define constant parameters
    double d = 1.0;
    double alpha = 1.7;
    int Nmax = 200;
    double epsilon = 1e-5;
    double gamma = 1e-8;
    double full_interval_a = -100.0;
    double full_interval_b = 100.0;

    // Main loop for 100 iterations
    for (int i = 0; i < 1; ++i) {
        cout << "Performing iteration no: " << i + 1 << "/100" << endl;

        // 1. Draw a random starting point x0 from the interval [-100, 100]
        double x0 = full_interval_a + (double) rand() / RAND_MAX * (full_interval_b - full_interval_a);

        // Variables to store results from this iteration
        double exp_start = NAN, exp_end = NAN;
        int exp_calls = 0;
        solution fib_exp_sol, lag_exp_sol, fib_full_sol, lag_full_sol;

        // 2. Optimization using initial expansion
        try {
            solution::clear_calls();
            double *interval = expansion(&ff1T, x0, d, alpha, Nmax);
            exp_start = interval[0];
            exp_end = interval[1];
            exp_calls = solution::f_calls;

            // Fibonacci on the narrowed interval
            solution::clear_calls();
            fib_exp_sol = fib(&ff1T, exp_start, exp_end, epsilon);
            fib_exp_sol.ud = solution::f_calls; // Store the number of calls in the 'ud' field

            // Lagrange on the narrowed interval
            solution::clear_calls();
            lag_exp_sol = lag(&ff1T, exp_start, exp_end, epsilon, gamma, Nmax);
            lag_exp_sol.ud = solution::f_calls;

            delete[] interval;
        } catch (string &ex) {
            // If expansion fails, save an error (NAN) and continue
            cerr << "   Iteration " << i + 1 << ": Error in expansion step - " << ex << endl;
        }

        // 3. Optimization without expansion (on the full interval [-100, 100])
        try {
            // Fibonacci on the full interval
            solution::clear_calls();
            fib_full_sol = fib(&ff1T, full_interval_a, full_interval_b, epsilon);
            fib_full_sol.ud = solution::f_calls;

            // Lagrange on the full interval
            solution::clear_calls();
            lag_full_sol = lag(&ff1T, full_interval_a, full_interval_b, epsilon, gamma, Nmax);
            lag_full_sol.ud = solution::f_calls;
        } catch (string &ex) {
            cerr << "   Iteration " << i + 1 << ": Error in full interval step - " << ex << endl;
        }

        // 4. Save all results from the iteration to the CSV file
        file << i + 1 << "," << x0 << ","
             << exp_start << "," << exp_end << "," << exp_calls << ","
             << fib_exp_sol.x(0) << "," << fib_exp_sol.y(0) << "," << fib_exp_sol.ud(0) << ","
             << lag_exp_sol.x(0) << "," << lag_exp_sol.y(0) << "," << lag_exp_sol.ud(0) << ","
             << fib_full_sol.x(0) << "," << fib_full_sol.y(0) << "," << fib_full_sol.ud(0) << ","
             << lag_full_sol.x(0) << "," << lag_full_sol.y(0) << "," << lag_full_sol.ud(0) << "\n";
    }

    // Close the file
    file.close();
    cout << "\nFinished. Results have been saved to lab1_results.csv" << endl;
}

void lab1_real()
{
    // Optimization parameters
    double epsilon = 1e-5;
    double gamma = 1e-9;
    int Nmax = 100;

    // Search interval for DA [1, 100] cm^2
    double a = 1.0;
    double b = 100.0;

    // --- Fibonacci Method ---
    solution::clear_calls();
    cout << "--- Fibonacci Method ---" << endl;
    solution sol_fib = fib(&ff1R, a, b, epsilon);
    cout << sol_fib << endl;

    // --- Lagrange Method ---
    solution::clear_calls();
    cout << "\n--- Lagrange Method ---" << endl;
    solution sol_lag = lag(&ff1R, a, b, epsilon, gamma, Nmax);
    cout << sol_lag << endl;

    // --- Simulation for the optimal DA from Fibonacci method ---
    cout << "\n--- Simulation for optimal DA (Fibonacci) = " << sol_fib.x(0) << " cm^2 ---" << endl;

    double DA_opt_fib = sol_fib.x(0) * 1e-4; // cm^2 -> m^2
    matrix Y0(3, 1);
    Y0(0) = 5.0;
    Y0(1) = 1.0;
    Y0(2) = 20.0;

    matrix params_fib(1, 1);
    params_fib(0) = DA_opt_fib;

    matrix* Y_fib = solve_ode(dff1R, 0, 1, 2000, Y0, NAN, params_fib);

    // Save Fibonacci simulation results to a file
    ofstream sim_file_fib(resolvePath("lab_1_real_fibonacci.csv"));
    sim_file_fib << "t,VA,VB,TB\n";
    for (int i = 0; i < get_len(Y_fib[0]); ++i)
    {
        sim_file_fib << Y_fib[0](i) << "," << Y_fib[1](i, 0) << "," << Y_fib[1](i, 1) << "," << Y_fib[1](i, 2) << "\n";
    }
    sim_file_fib.close();
    cout << "Fibonacci simulation results saved to lab_1_real_fibonacci.csv" << endl;
    delete[] Y_fib; // Free memory

    //  Simulation for the optimal DA from Lagrange method ---
    cout << "\n--- Simulation for optimal DA (Lagrange) = " << sol_lag.x(0) << " cm^2 ---" << endl;

    double DA_opt_lag = sol_lag.x(0) * 1e-4; // cm^2 -> m^2

    // Initial conditions are the same
    // Y0(0) = 5.0;
    // Y0(1) = 1.0;
    // Y0(2) = 20.0;

    matrix params_lag(1, 1);
    params_lag(0) = DA_opt_lag;

    matrix* Y_lag = solve_ode(dff1R, 0, 1, 2000, Y0, NAN, params_lag);

    // Save Lagrange simulation results to a file
    ofstream sim_file_lag(resolvePath("lab_1_real_lagrange.csv"));
    sim_file_lag << "t,VA,VB,TB\n";
    for (int i = 0; i < get_len(Y_lag[0]); ++i)
    {
        sim_file_lag << Y_lag[0](i) << "," << Y_lag[1](i, 0) << "," << Y_lag[1](i, 1) << "," << Y_lag[1](i, 2) << "\n";
    }
    sim_file_lag.close();
    cout << "Lagrange simulation results saved to lab_1_real_lagrange.csv" << endl;
    delete[] Y_lag; // Free memory
}

void lab2() {

    srand(time(nullptr));


    ofstream file(resolvePath("lab2_results.csv"));
    if (!file.is_open()) {
        cerr << "ERROR: Nie mozna otworzyc pliku lab2_results.csv do zapisu!" << endl;
        return;
    }


    file << "Iteracja,Krok_startowy,x0_start,x1_start,"
         << "HJ_x0_end,HJ_x1_end,HJ_y_end,HJ_f_calls,"
         << "Rosen_x0_end,Rosen_x1_end,Rosen_y_end,Rosen_f_calls\n";



    double epsilon = 1e-5;
    int Nmax = 20000;


    double alpha_hj = 0.5;


    matrix s0_rosen(2, 1);
    s0_rosen(0) = 0.5;
    s0_rosen(1) = 0.5;
    double alpha_rosen = 2.0;
    double beta_rosen = 0.5;


    vector<double> kroki_startowe = {0.5, 0.25, 0.1};

    int global_iteration_count = 1;


    for (double start_s : kroki_startowe) {
        cout << "--- Rozpoczynam testy dla kroku startowego s = " << start_s << " ---" << endl;


        for (int i = 0; i < 100; ++i) {
            cout << "   Wykonuje iteracje nr: " << i + 1 << "/100" << endl;


            matrix x0(2, 1);
            x0(0) = -1.0 + (double)rand() / RAND_MAX * 2.0;
            x0(1) = -1.0 + (double)rand() / RAND_MAX * 2.0;

            solution sol_hj, sol_rosen;
            int hj_calls = 0, rosen_calls = 0;


            try {
                solution::clear_calls();
                sol_hj = HJ(&ff2T, x0, start_s, alpha_hj, epsilon, Nmax);
                hj_calls = solution::f_calls;
            } catch (const string& ex) {
                cerr << "   Blad w metodzie Hooke-Jeevesa: " << ex << endl;

            }


            try {
                solution::clear_calls();
                sol_rosen = Rosen(&ff2T, x0, s0_rosen, alpha_rosen, beta_rosen, epsilon, Nmax);
                rosen_calls = solution::f_calls;
            } catch (const string& ex) {
                cerr << "   Blad w metodzie Rosenbrocka: " << ex << endl;
            }


            file << global_iteration_count++ << "," << start_s << ","
                 << x0(0) << "," << x0(1) << ","
                 << sol_hj.x(0) << "," << sol_hj.x(1) << "," << sol_hj.y(0) << "," << hj_calls << ","
                 << sol_rosen.x(0) << "," << sol_rosen.x(1) << "," << sol_rosen.y(0) << "," << rosen_calls << "\n";
        }
    }


    file.close();
    cout << "\nZakonczono. Wyniki zostaly zapisane do pliku lab2_results.csv" << endl;
    // --- Zadanie 5b: Problem rzeczywisty ---
    cout << "--- Laboratorium 2: Optymalizacja dla problemu rzeczywistego ---" << endl;

    ofstream real_problem_file(resolvePath("lab2_real_problem_results.csv"));
    if (!real_problem_file.is_open()) {
        cerr << "ERROR: Nie mozna otworzyc pliku lab2_real_problem_results.csv do zapisu!" << endl;
    } else {
        real_problem_file << "Dlugosc_kroku,HJ_k1,HJ_k2,HJ_Q,HJ_f_calls,Rosen_k1,Rosen_k2,Rosen_Q,Rosen_f_calls\n";
    }

    // Inicjalizacja generatora liczb losowych
    srand(time(nullptr));

    // Parametry optymalizacji

    // Przedział poszukiwań punktu startowego dla k1 i k2
    double k_min = 0.0;
    double k_max = 20.0;

    // Wygenerowanie pojedynczego losowego punktu startowego x0 = [k1, k2]
    matrix x0(2, 1);
    x0(0) = k_min + (double)rand() / RAND_MAX * (k_max - k_min);
    x0(1) = k_min + (double)rand() / RAND_MAX * (k_max - k_min);

    cout << "Losowy punkt startowy (k1, k2): (" << x0(0) << ", " << x0(1) << ")" << endl << endl;

    vector<double> step_sizes = { 5.0, 2.0, 1.0, 0.5, 0.1 };

    for (double current_step : step_sizes) {
        cout << "--- Testowanie dla kroku startowego: " << current_step << " ---" << endl;

        solution sol_hj, sol_rosen;
        int hj_calls = 0, rosen_calls = 0;

        // --- Metoda Hooke'a-Jeevesa ---
        try {
            solution::clear_calls();
            sol_hj = HJ(&ff2R, x0, current_step, alpha_hj, epsilon, Nmax);
            hj_calls = solution::f_calls;
            cout << "Optymalizacja metoda Hooke'a-Jeevesa zakonczona." << endl;
            cout << sol_hj << endl;
        } catch (const string& ex) {
            cerr << "Blad w metodzie Hooke'a-Jeevesa: " << ex << endl;
        }

        // --- Metoda Rosenbrocka ---
        matrix s0_rosen_loop(2, 1);
        s0_rosen_loop(0) = current_step;
        s0_rosen_loop(1) = current_step;
        try {
            solution::clear_calls();
            sol_rosen = Rosen(&ff2R, x0, s0_rosen_loop, alpha_rosen, beta_rosen, epsilon, Nmax);
            rosen_calls = solution::f_calls;
            cout << "Optymalizacja metoda Rosenbrocka zakonczona." << endl;
            cout << sol_rosen << endl;
        } catch (const string& ex) {
            cerr << "Blad w metodzie Rosenbrocka: " << ex << endl;
        }

        if (real_problem_file.is_open()) {
            real_problem_file << current_step << ","
                              << sol_hj.x(0) << "," << sol_hj.x(1) << "," << sol_hj.y(0) << "," << hj_calls << ","
                              << sol_rosen.x(0) << "," << sol_rosen.x(1) << "," << sol_rosen.y(0) << "," << rosen_calls << "\n";
        }

        // This part is for generating data for the second table
        if (current_step == 0.5) { // Run simulation only for one representative step size
            // --- Symulacja i zapis do pliku dla optymalnych parametrów z metody Hooke'a-Jeevesa ---
            cout << "--- Przeprowadzanie symulacji dla optymalnych parametrow z Hooke-Jeeves (krok=" << current_step << ") ---" << endl;
            if (sol_hj.flag > 0) { // Sprawdzenie, czy metoda HJ znalazła rozwiązanie
                double k1_opt = sol_hj.x(0);
                double k2_opt = sol_hj.x(1);
                cout << "Parametry optymalne: k1 = " << k1_opt << ", k2 = " << k2_opt << endl;

                matrix Y0(2, 1); Y0(0) = 0.0; Y0(1) = 0.0;
                matrix ode_params(2, 1); ode_params(0) = k1_opt; ode_params(1) = k2_opt;

                matrix* Y = solve_ode(df2R, 0.0, 0.1, 100.0, Y0, NAN, ode_params);

                ofstream sim_file(resolvePath("lab2_real_simulation_HJ.csv"));
                if (sim_file.is_open()) {
                    sim_file << "t,alpha,omega\n";
                    int num_steps = get_len(Y[0]);
                    for (int i = 0; i < num_steps; ++i) {
                        sim_file << Y[0](i) << "," << Y[1](i, 0) << "," << Y[1](i, 1) << "\n";
                    }
                    sim_file.close();
                    cout << "Wyniki symulacji zostaly zapisane do pliku lab2_real_simulation_HJ.csv" << endl << endl;
                } else {
                    cerr << "BLAD: Nie mozna otworzyc pliku lab2_real_simulation_HJ.csv do zapisu!" << endl;
                }
                delete[] Y;
            } else {
                cout << "Pominieto symulacje, poniewaz metoda Hooke-Jeevesa nie zwrocila poprawnego rozwiazania." << endl << endl;
            }

            // --- Symulacja i zapis do pliku dla optymalnych parametrów z metody Rosenbrocka ---
            cout << "--- Przeprowadzanie symulacji dla optymalnych parametrow z Rosenbrocka (krok=" << current_step << ") ---" << endl;
            if (sol_rosen.flag > 0) { // Sprawdzenie, czy metoda Rosenbrocka znalazła rozwiązanie
                double k1_opt = sol_rosen.x(0);
                double k2_opt = sol_rosen.x(1);
                cout << "Parametry optymalne: k1 = " << k1_opt << ", k2 = " << k2_opt << endl;

                matrix Y0(2, 1); Y0(0) = 0.0; Y0(1) = 0.0;
                matrix ode_params(2, 1); ode_params(0) = k1_opt; ode_params(1) = k2_opt;

                matrix* Y = solve_ode(df2R, 0.0, 0.1, 100.0, Y0, NAN, ode_params);

                ofstream sim_file(resolvePath("lab2_real_simulation_Rosenbrock.csv"));
                if (sim_file.is_open()) {
                    sim_file << "t,alpha,omega\n";
                    int num_steps = get_len(Y[0]);
                    for (int i = 0; i < num_steps; ++i) {
                        sim_file << Y[0](i) << "," << Y[1](i, 0) << "," << Y[1](i, 1) << "\n";
                    }
                    sim_file.close();
                    cout << "Wyniki symulacji zostaly zapisane do pliku lab2_real_simulation_Rosenbrock.csv" << endl << endl;
                } else {
                    cerr << "BLAD: Nie mozna otworzyc pliku lab2_real_simulation_Rosenbrock.csv do zapisu!" << endl;
                }
                delete[] Y;
            } else {
                cout << "Pominieto symulacje, poniewaz metoda Rosenbrocka nie zwrocila poprawnego rozwiazania." << endl << endl;
            }
        }
    }
    if (real_problem_file.is_open()) {
        real_problem_file.close();
        cout << "\nZakonczono testy dla problemu rzeczywistego. Wyniki zostaly zapisane do pliku lab2_real_problem_results.csv" << endl;
    }
}

void lab3() {
    srand(time(nullptr));

    // Otwarcie pliku do zapisu wyników
    
    ofstream file("../data/lab3_results.csv");
    if (!file.is_open()) {
        cerr << "Błąd: Nie można otworzyć pliku lab3_results.csv" << endl;
        return;
    }
    file << "Metoda,Parametr_a,Iteracja,x0_start,x1_start,x0_end,x1_end,y_end,f_calls,r\n";

    // Parametry dla metody Neldera-Meada
    double s = 0.5;       // Długość boku sympleksu
    double alpha = 1.0;   // Odbicie
    double beta = 0.5;    // Zawężenie
    double gamma = 2.0;   // Ekspansja
    double delta = 0.5;   // Redukcja
    double epsilon_nm = 1e-3; // Dokładność
    int Nmax = 10000;     // Maksymalna liczba wywołań funkcji celu

    // Parametry dla pętli metody kar
    double c_start = 1.0;
    double dc = 5.0;      // Współczynnik skalowania kary
    double epsilon_pen = 1e-4; // Warunek stopu dla pętli kar
    int max_iter_pen = 20;

    vector<double> a_params = {4.0, 4.4934, 5.0};;
    int N_opt = 100;

    for (double a : a_params) {
        cout << "--- Testowanie dla a = " << a << " ---" << endl;

        // --- ZEWNĘTRZNA FUNKCJA KARY ---
        cout << "  Metoda kary zewnętrznej..." << endl;
        for (int i = 0; i < N_opt; ++i) {
            solution::clear_calls();

            // Losowy punkt startowy z całego obszaru [-a, a] x [-a, a]
            matrix x0(2, 1);
            x0(0) = -a + (double)rand() / RAND_MAX * (2 * a);
            x0(1) = -a + (double)rand() / RAND_MAX * (2 * a);

            matrix x_prev = x0;
            double c = c_start;
            solution sol;

            for (int k = 0; k < max_iter_pen; ++k) {
                matrix ud(2, 1);
                ud(0) = a;
                ud(1) = c;

                sol = sym_NM(ff3T_zew, x_prev, s, alpha, beta, gamma, delta, epsilon_nm, Nmax - solution::f_calls, ud);

                if (norm(sol.x - x_prev) < epsilon_pen || solution::f_calls >= Nmax) {
                    break;
                }
                x_prev = sol.x;
                c *= dc;
            }

            double final_y = m2d(ff3T(sol.x));
            double r = m2d(norm(sol.x));

            file << "Zewnetrzna," << a << "," << i + 1 << ","
                 << x0(0) << "," << x0(1) << ","
                 << sol.x(0) << "," << sol.x(1) << "," << final_y << ","
                 << solution::f_calls << "," << r << "\n";
        }

        // --- WEWNĘTRZNA FUNKCJA KARY ---
        cout << "  Metoda kary wewnętrznej..." << endl;
        for (int i = 0; i < N_opt; ++i) {
            solution::clear_calls();

            // Losowy punkt startowy ze ściśle dopuszczalnego obszaru
            matrix x0(2, 1);
            do {
                x0(0) = 1.001 + (double)rand() / RAND_MAX * (a - 1.002);
                x0(1) = 1.001 + (double)rand() / RAND_MAX * (a - 1.002);
            } while (sqrt(pow(x0(0), 2) + pow(x0(1), 2)) >= a - 0.001);

            matrix x_prev = x0;
            double c = 10.0; // Wyższy c0 dla metody wewn.
            solution sol;

            for (int k = 0; k < max_iter_pen; ++k) {
                matrix ud(2, 1);
                ud(0) = a;
                ud(1) = c;

                sol = sym_NM(ff3T_wew, x_prev, s, alpha, beta, gamma, delta, epsilon_nm, Nmax - solution::f_calls, ud);

                if (norm(sol.x - x_prev) < epsilon_pen || solution::f_calls >= Nmax) {
                    break;
                }
                x_prev = sol.x;
                c /= dc; // Zmniejszanie 'c' dla metody wewnętrznej
            }

            double final_y = m2d(ff3T(sol.x));
            double r = m2d(norm(sol.x));

            file << "Wewnetrzna," << a << "," << i + 1 << ","
                 << x0(0) << "," << x0(1) << ","
                 << sol.x(0) << "," << sol.x(1) << "," << final_y << ","
                 << solution::f_calls << "," << r << "\n";
        }
    }

    file.close();
    cout << "\nZakończono. Wyniki zapisano do pliku lab3_results.csv" << endl;

    cout << "\n\n--- Lab 3: Problem Rzeczywisty (Pilka) ---" << endl;

        ofstream file_real(resolvePath("lab3_real_problem.csv"));
        file_real << "Iteracja,c,v0x_opt,omega_opt,x_end,f_calls\n";

        // Parametry Neldera-Meada
         s = 0.5;
         alpha = 1.0, beta = 0.5, gamma = 2.0, delta = 0.5;
         epsilon_nm = 1e-6;
         Nmax = 5000;

        // Parametry metody kar
        double c = 1.0;       // c startowe
         dc = 2.0;      // Krok zwiększania kary (można dobrać np. 2 lub 5)
         epsilon_pen = 1e-6;
        int max_pen_iter = 20;

        // Punkt startowy (losowy w dopuszczalnych granicach [-10, 10])
        matrix x0(2, 1);
        x0(0) = -10.0 + (double)rand() / RAND_MAX * 20.0; // v0x
        x0(1) = -10.0 + (double)rand() / RAND_MAX * 20.0; // omega

        cout << "Punkt startowy: v0x=" << x0(0) << ", omega=" << x0(1) << endl;

        matrix x_prev = x0;
        solution sol;

        // Pętla metody kar zewnętrznych
        for (int k = 0; k < max_pen_iter; ++k) {
            solution::clear_calls();
            matrix ud2(1, 1);
            ud2(0) = c;

            // Wywołanie NM dla funkcji z karą
            // Uwaga: ff3R_zew nie używa ud1, więc dajemy NAN, parametry kary idą w ud2
            sol = sym_NM(ff3R_zew, x_prev, s, alpha, beta, gamma, delta, epsilon_nm, Nmax, NAN, ud2);

            // Obliczamy rzeczywistą wartość x_end (bez kary) dla logów
            // Wywołujemy ff3R_zew z c=0, żeby dostać czystą wartość funkcji celu (-x_end)
            matrix ud2_zero(1, 1); ud2_zero(0) = 0.0;
            double real_obj = m2d(ff3R_zew(sol.x, NAN, ud2_zero));
            double x_end_val = -real_obj; // Bo minimalizowaliśmy -x_end

            cout << "Metoda kar it." << k+1 << " c=" << c
                 << " v0x=" << sol.x(0) << " w=" << sol.x(1)
                 << " x_end=" << x_end_val << endl;

            file_real << k + 1 << "," << c << ","
                      << sol.x(0) << "," << sol.x(1) << ","
                      << x_end_val << "," << solution::f_calls << "\n";

            // Warunek stopu (zmiana x jest mała)
            if (norm(sol.x - x_prev) < epsilon_pen) {
                cout << "Osiagnieto zbieznosc metody kar." << endl;
                break;
            }

            x_prev = sol.x;
            c *= dc;
        }
        file_real.close();

        // --- Symulacja dla znalezionych optymalnych parametrów ---
        cout << "Generowanie symulacji dla optymalnego rozwiazania..." << endl;

        double opt_v0x = sol.x(0);
        double opt_omega = sol.x(1);

        matrix Y0(4, 1);
        Y0(0) = 0.0; Y0(1) = opt_v0x; Y0(2) = 100.0; Y0(3) = 0.0;
        matrix params(1, 1);
        params(0) = opt_omega;

        matrix* Y = solve_ode(dff3R, 0, 0.01, 7, Y0, params, NAN);

        ofstream sim_file(resolvePath("lab3_real_simulation.csv"));
        sim_file << "t,x,vx,y,vy\n";
        int n = get_len(Y[0]);
        for(int i=0; i<n; ++i) {
            sim_file << Y[0](i) << ","
                     << Y[1](i, 0) << "," << Y[1](i, 1) << ","
                     << Y[1](i, 2) << "," << Y[1](i, 3) << "\n";
        }
        sim_file.close();
        delete[] Y;

        cout << "Zapisano wyniki symulacji do lab3_real_simulation.csv" << endl;
}

void lab4() {
    srand(time(nullptr));
    cout << "--- LABORATORIUM 4: Metody Gradientowe ---" << endl;

    // ==========================================================
    // CZĘŚĆ A: Testowa funkcja celu
    // ==========================================================
    cout << "\n=== CZESC A: Testowa funkcja celu ===" << endl;
    
    double epsilon = 1e-5;
    int Nmax = 10000;
    
    // Konfiguracja kroków
    struct StepConfig {
        double val;
        string name;
    };
    vector<StepConfig> steps = {
        {0.05, "0.05"},
        {0.25, "0.25"},
        {NAN,  "Zmienny"}
    };
    
    // Otwarcie trzech osobnych plików
    ofstream file_SD("lab4_test_SD.csv");
    ofstream file_CG("lab4_test_CG.csv");
    ofstream file_Newton("lab4_test_Newton.csv");

    if (file_SD.is_open() && file_CG.is_open() && file_Newton.is_open()) {
        // Zapis nagłówków do każdego pliku (usunąłem kolumnę 'Metoda', bo nazwa pliku o tym mówi)
        string header = "Krok,Iteracja,x1_start,x2_start,x1_end,x2_end,y_end,f_calls,Flag\n";
        file_SD << header;
        file_CG << header;
        file_Newton << header;
        
        for (const auto& step : steps) {
            cout << "Testowanie dla kroku: " << step.name << endl;
            
            for (int i = 0; i < 100; ++i) {
                // Losowy punkt startowy [-2, 2] x [-2, 2]
                matrix x0(2, 1);
                x0(0) = -2.0 + (double)rand()/RAND_MAX * 4.0;
                x0(1) = -2.0 + (double)rand()/RAND_MAX * 4.0;
                
                // --- 1. Najszybszy spadek (SD) ---
                try {
                    solution::clear_calls();
                    solution sol = SD(ff4T, gf4T, x0, step.val, epsilon, Nmax);
                    file_SD << step.name << "," << i+1 << "," 
                            << x0(0) << "," << x0(1) << ","
                            << sol.x(0) << "," << sol.x(1) << "," << sol.y(0) << "," 
                            << solution::f_calls << "," << sol.flag << "\n";
                } catch(string ex) { 
                    file_SD << step.name << "," << i+1 << ",ERROR," << ex << "\n"; 
                }

                // --- 2. Gradienty Sprzężone (CG) ---
                try {
                    solution::clear_calls();
                    solution sol = CG(ff4T, gf4T, x0, step.val, epsilon, Nmax);
                    file_CG << step.name << "," << i+1 << "," 
                            << x0(0) << "," << x0(1) << ","
                            << sol.x(0) << "," << sol.x(1) << "," << sol.y(0) << "," 
                            << solution::f_calls << "," << sol.flag << "\n";
                } catch(string ex) { 
                    file_CG << step.name << "," << i+1 << ",ERROR," << ex << "\n"; 
                }
                
                // --- 3. Metoda Newtona ---
                try {
                    solution::clear_calls();
                    solution sol = Newton(ff4T, gf4T, Hf4T, x0, step.val, epsilon, Nmax);
                    file_Newton << step.name << "," << i+1 << "," 
                                << x0(0) << "," << x0(1) << ","
                                << sol.x(0) << "," << sol.x(1) << "," << sol.y(0) << "," 
                                << solution::f_calls << "," << sol.flag << "\n";
                } catch(string ex) { 
                    file_Newton << step.name << "," << i+1 << ",ERROR," << ex << "\n"; 
                }
            }
        }
        
        file_SD.close();
        file_CG.close();
        file_Newton.close();
        cout << "Wyniki czesci A zapisano do osobnych plikow:\n"
             << " - lab4_test_SD.csv\n"
             << " - lab4_test_CG.csv\n"
             << " - lab4_test_Newton.csv" << endl;
    } else {
        cerr << "Blad: Nie udalo sie otworzyc plikow CSV do zapisu!" << endl;
    }

    // ==========================================================
    // CZĘŚĆ B: Problem rzeczywisty (Klasyfikacja)
    // ==========================================================
    cout << "\n=== CZESC B: Problem rzeczywisty (Klasyfikacja) ===" << endl;
    
    // Dane zaszyte w kodzie (Hardcoded Data)
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

    // Tworzenie macierzy
    matrix X(3, 100);
    matrix Y(1, 100);
    
    for (int i = 0; i < 100; ++i) {
        Y(0, i) = y_raw[i];
        X(0, i) = 1.0;       // Bias
        X(1, i) = x1_raw[i]; // Feature 1
        X(2, i) = x2_raw[i]; // Feature 2
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
                // Tutaj używamy metody CG (problem rzeczywisty wymagał jednej metody)
                solution sol = CG(ff4R, gf4R, theta0, h, epsilon, Nmax, X, Y);
                
                int correct = 0;
                int m = 100;
                for (int i = 0; i < m; ++i) {
                    matrix xi = get_col(X, i);
                    double yi = Y(0, i);
                    double z = m2d(trans(sol.x) * xi);
                    double h_val = 1.0 / (1.0 + exp(-z));
                    int prediction = (h_val >= 0.5) ? 1 : 0;
                    if (prediction == (int)yi) correct++;
                }
                double accuracy = (double)correct / m * 100.0;
                
                cout << "  Wynik: J=" << sol.y(0) << ", Iter=" << solution::f_calls 
                     << ", Acc=" << accuracy << "%" << endl;
                
                real_res << h << "," 
                         << sol.x(0) << "," << sol.x(1) << "," << sol.x(2) << ","
                         << sol.y(0) << "," << solution::f_calls << "," << accuracy << "\n";
            } catch (string ex) {
                cerr << "  Blad podczas optymalizacji: " << ex << endl;
                real_res << h << ",ERROR,ERROR,ERROR,ERROR,ERROR,ERROR\n";
            }
        }
        real_res.close();
        cout << "Wyniki czesci B zapisano do lab4_real_results.csv" << endl;
    } else {
        cerr << "Nie udalo sie otworzyc pliku do zapisu wynikow czesci B." << endl;
    }
}

void lab5() {
}

void lab6() {
}

void test_real_problem_DA50()
{
    cout << "\n--- Running Test for DA = 50 cm^2 ---" << endl;

    // 1. Define the fixed DA value (50 cm^2 -> m^2)
    double DA_test = 50.0 * 1e-4;

    // 2. Set initial conditions
    matrix Y0(3, 1);
    Y0(0) = 5.0;  // VA_start = 5 m^3
    Y0(1) = 1.0;  // VB_start = 1 m^3
    Y0(2) = 20.0; // TB_start = 20 C

    // 3. Set the parameters for the ODE solver
    matrix params(1, 1);
    params(0) = DA_test;

    // 4. Solve the ODE with a time step of 1 second
    matrix* Y = solve_ode(dff1R, 0, 1, 2000, Y0, NAN, params);

    // 5. Find the maximum temperature in tank B from the results
    int n = get_len(Y[0]);
    double T_max = 0;
    for (int i = 0; i < n; ++i)
    {
        if (Y[1](i, 2) > T_max)
        {
            T_max = Y[1](i, 2);
        }
    }

    // 6. Print the result
    cout << "Simulation finished." << endl;
    cout << "Maximum temperature in Tank B: " << T_max << " C" << endl;

    // 7. Free the dynamically allocated memory
    delete[] Y;
}