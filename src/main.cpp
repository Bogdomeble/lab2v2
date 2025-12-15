

#include "../include/opt_alg.h"
#include "../include/user_funs.h"
#include <cstdlib>
#include <ctime>

void lab0();
void lab1();
void lab1_real();
void test_real_problem_DA50();
void lab2();
void lab3();
void lab4();
void lab4_wykresy();
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
        // lab4();
        // lab4_wykresy();
        lab5();
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
    ofstream file(resolvePath("lab3_results.csv"));
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
    // lab4_wykresy();
    srand(time(nullptr));

    // ==========================================================
    // CZĘŚĆ 1: FUNKCJA TESTOWA - Generowanie CSV pod Excel
    // ==========================================================
    cout << "Generowanie wynikow dla funkcji testowej..." << endl;

    ofstream file(resolvePath("lab4_testowa_tabela.csv"));
    if (!file.is_open()) {
        cerr << "ERROR: Nie mozna otworzyc pliku lab4_testowa_tabela.csv" << endl;
        return;
    }

    // Naglowek pasujacy do Twojego Excela (US locale - kropka jako separator)
    // Uklad: Krok, Lp, x1(0), x2(0), [SD], [CG], [Newton]
    file << "Dlugosc_kroku,Lp,x1_0,x2_0,"
         // Metoda najszybszego spadku
         << "SD_x1,SD_x2,SD_y,SD_f_calls,SD_g_calls,SD_Globalne,"
         // Metoda gradientow sprzezonych
         << "CG_x1,CG_x2,CG_y,CG_f_calls,CG_g_calls,CG_Globalne,"
         // Metoda Newtona
         << "Newt_x1,Newt_x2,Newt_y,Newt_f_calls,Newt_g_calls,Newt_H_calls,Newt_Globalne\n";

    vector<double> steps = {0.05, 0.25, NAN}; // NAN = zmienny
    int Nmax = 10000;
    double epsilon = 1e-5;
    double global_threshold = 1e-3; // Prog uznania za minimum globalne (0)

    for (double h : steps) {
        string step_label = isnan(h) ? "Zmienny" : to_string(h);

        for (int i = 0; i < 100; ++i) {
            // Punkt startowy losowy [-2, 2]
            matrix x0(2, 1);
            x0(0) = -2.0 + (double)rand() / RAND_MAX * 4.0;
            x0(1) = -2.0 + (double)rand() / RAND_MAX * 4.0;

            // --- 1. SD (Najszybszy Spadek) ---
            solution::clear_calls();
            solution sol_sd = SD(ff4T, gf4T, x0, h, epsilon, Nmax);
            string sd_glob = (abs(m2d(sol_sd.y)) < global_threshold) ? "TAK" : "NIE";

            // --- 2. CG (Gradienty Sprzezone) ---
            solution::clear_calls();
            solution sol_cg = CG(ff4T, gf4T, x0, h, epsilon, Nmax);
            string cg_glob = (abs(m2d(sol_cg.y)) < global_threshold) ? "TAK" : "NIE";

            // --- 3. Newton ---
            solution::clear_calls();
            solution sol_newton = Newton(ff4T, gf4T, Hf4T, x0, h, epsilon, Nmax);
            string newt_glob = (abs(m2d(sol_newton.y)) < global_threshold) ? "TAK" : "NIE";

            // ZAPIS WIERSZA DO CSV
            file << step_label << "," << (i + 1) << ","
                 << x0(0) << "," << x0(1) << ","
                 // SD Cols
                 << sol_sd.x(0) << "," << sol_sd.x(1) << "," << sol_sd.y(0) << ","
                 << solution::f_calls << "," << solution::g_calls << "," << sd_glob << ","
                 // CG Cols
                 << sol_cg.x(0) << "," << sol_cg.x(1) << "," << sol_cg.y(0) << ","
                 << solution::f_calls << "," << solution::g_calls << "," << cg_glob << ","
                 // Newton Cols
                 << sol_newton.x(0) << "," << sol_newton.x(1) << "," << sol_newton.y(0) << ","
                 << solution::f_calls << "," << solution::g_calls << "," << solution::H_calls << "," << newt_glob
                 << "\n";
        }
    }
    file.close();
    cout << "Plik lab4_testowa_tabela.csv wygenerowany." << endl;

    // ==========================================================
    // CZĘŚĆ 2: PROBLEM RZECZYWISTY (XData, YData)
    // ==========================================================
    cout << "Przetwarzanie problemu rzeczywistego (LogReg)..." << endl;

    ofstream file_real(resolvePath("lab4_rzeczywisty_wyniki.csv"));
    file_real << "Krok,Theta0,Theta1,Theta2,f_min,f_calls,Dokladnosc_%\n";

    try {
        // Wczytujemy Twoje specyficzne pliki (format z srednikami)
        // XData ma 3 wiersze (Bias, Feature1, Feature2) i 100 kolumn
        // YData ma 1 wiersz i 100 kolumn
        matrix X = read_matrix_semicolon(resolvePath("XData.txt"), 3, 100);
        matrix Y = read_matrix_semicolon(resolvePath("YData.txt"), 1, 100);

        // Skalowanie danych (opcjonalne, ale pomaga przy duzych wartosciach jak 80-100)
        // Tutaj zostawiamy surowe, bo takie dostales polecenie, ale to przyczyna bledow.

        matrix theta0(3, 1, 0.0); // Start theta = [0, 0, 0]
        vector<double> real_steps = {0.01, 0.001, 0.0001};

        for (double h : real_steps) {
            solution::clear_calls();

            // Uzywamy CG (Gradienty Sprzezone) dla problemu rzeczywistego
            solution sol = CG(ff4R, gf4R, theta0, h, epsilon, Nmax, X, Y);

            double acc = calculate_accuracy(sol.x, X, Y);

            file_real << h << ","
                      << sol.x(0) << "," << sol.x(1) << "," << sol.x(2) << ","
                      << sol.y(0) << "," << solution::f_calls << "," << acc << "\n";

            cout << "Krok " << h << ": f_min=" << sol.y(0) << ", Acc=" << acc << "%" << endl;
        }

    } catch (string& ex) {
        cerr << "BLAD w czesci rzeczywistej: " << ex << endl;
        cerr << "Upewnij sie, ze pliki XData.txt i YData.txt sa w folderze z danymi." << endl;
    }

    file_real.close();
}

// do tabeli 2 i wykresow wersja co iteracje
//
/*
 * void lab4() {
     srand(time(nullptr));

     // ==========================================
     // TABELA 1: FUNKCJA TESTOWA
     // ==========================================
     cout << "Generowanie tabeli 1 (Funkcja testowa)..." << endl;

     ofstream t1(resolvePath("lab4_tabela1.csv"));
     // Nagłówek ułatwiający identyfikację kolumn (opcjonalnie do usunięcia przed wklejeniem)
     t1 << "Nr_iteracji,x1(0),x2(0),"
        << "SD_005_x1,SD_005_x2,SD_005_y,SD_005_fc,SD_005_gc,SD_005_glob,"
        << "SD_025_x1,SD_025_x2,SD_025_y,SD_025_fc,SD_025_gc,SD_025_glob,"
        << "SD_Var_x1,SD_Var_x2,SD_Var_y,SD_Var_fc,SD_Var_gc,SD_Var_glob,"
        << "CG_005_x1,CG_005_x2,CG_005_y,CG_005_fc,CG_005_gc,CG_005_glob,"
        << "CG_025_x1,CG_025_x2,CG_025_y,CG_025_fc,CG_025_gc,CG_025_glob,"
        << "CG_Var_x1,CG_Var_x2,CG_Var_y,CG_Var_fc,CG_Var_gc,CG_Var_glob,"
        << "Newt_005_x1,Newt_005_x2,Newt_005_y,Newt_005_fc,Newt_005_gc,Newt_005_glob,"
        << "Newt_025_x1,Newt_025_x2,Newt_025_y,Newt_025_fc,Newt_025_gc,Newt_025_glob,"
        << "Newt_Var_x1,Newt_Var_x2,Newt_Var_y,Newt_Var_fc,Newt_Var_gc,Newt_Var_glob\n";

     double steps[] = {0.05, 0.25, NAN};
     int Nmax = 10000;
     double epsilon = 1e-5;
     double glob_thresh = 1e-3;

     // Główna pętla po 100 iteracjach (wierszach w Excelu)
     for (int i = 0; i < 100; ++i) {
         // Losowy punkt startowy wspólny dla wiersza
         matrix x0(2, 1);
         x0(0) = -2.0 + (double)rand() / RAND_MAX * 4.0;
         x0(1) = -2.0 + (double)rand() / RAND_MAX * 4.0;

         // Wypisz początek wiersza: Nr iteracji, x1(0), x2(0)
         t1 << (i + 1) << "," << x0(0) << "," << x0(1);

         // Pętla po metodach: SD, CG, Newton
         // Kolejność w Excelu: Najszybszy Spadek -> Gradienty Sprzężone -> Newton
         for (int m = 0; m < 3; ++m) {
             // Pętla po krokach: 0.05 -> 0.25 -> M. zm.
             for (double h : steps) {
                 solution::clear_calls();
                 solution sol;
                 try {
                     if (m == 0) sol = SD(ff4T, gf4T, x0, h, epsilon, Nmax);
                     else if (m == 1) sol = CG(ff4T, gf4T, x0, h, epsilon, Nmax);
                     else sol = Newton(ff4T, gf4T, Hf4T, x0, h, epsilon, Nmax);
                 } catch (...) {
                     sol.y = matrix(NAN); // W razie błędu
                 }

                 string is_glob = (abs(m2d(sol.y)) < glob_thresh) ? "TAK" : "NIE";

                 // Zapisz wyniki dla danej komórki tabeli
                 // Format: x1*, x2*, y*, f_calls, g_calls, Globalne?
                 t1 << "," << sol.x(0) << "," << sol.x(1) << "," << m2d(sol.y) << ","
                    << solution::f_calls << "," << solution::g_calls << "," << is_glob;
             }
         }
         t1 << "\n"; // Koniec wiersza Excela
     }
     t1.close();
     cout << "Plik lab4_tabela1.csv gotowy." << endl;


     // ==========================================
     // TABELA 3: PROBLEM RZECZYWISTY
     // ==========================================
     cout << "Generowanie tabeli 3 (Problem rzeczywisty)..." << endl;

     ofstream t3(resolvePath("lab4_tabela3.csv"));
     t3 << "Dlugosc_kroku,theta0,theta1,theta2,J(theta),P(theta),g_calls\n";

     try {
         matrix X = read_matrix_semicolon(resolvePath("XData.txt"), 3, 100);
         matrix Y = read_matrix_semicolon(resolvePath("YData.txt"), 1, 100);

         // Punkt startowy [0, 0, 0]
         matrix theta0(3, 1, 0.0);
         vector<double> r_steps = {0.01, 0.001, 0.0001};

         for (double h : r_steps) {
             solution::clear_calls();

             // Metoda CG dla problemu rzeczywistego
             solution sol = CG(ff4R, gf4R, theta0, h, epsilon, Nmax, X, Y);

             double acc = calculate_accuracy(sol.x, X, Y); // Funkcja P(theta)

             t3 << h << ","
                << sol.x(0) << "," << sol.x(1) << "," << sol.x(2) << ","
                << m2d(sol.y) << "," << acc << "," << solution::g_calls << "\n";
         }
     } catch (string& ex) {
         cerr << "Blad w problemie rzeczywistym: " << ex << endl;
     }
     t3.close();
     cout << "Plik lab4_tabela3.csv gotowy." << endl;
 }
 */
 void lab5() {
     srand(time(nullptr));
     double epsilon = 1e-5;
     int Nmax = 2000;

     // ==========================================================
     // TABELA 1: FUNKCJA TESTOWA (wersja szeroka)
     // ==========================================================
     cout << "Generowanie danych do Tabeli 1 (funkcja testowa)..." << endl;

     ofstream t1(resolvePath("lab5_tabela1.csv"));
     if (!t1.is_open()) {
         cerr << "BLAD: Nie mozna otworzyc pliku lab5_tabela1.csv" << endl;
         return;
     }

     // w; x1(0); x2(0); | a=1 (x1,x2,f1,f2,calls) | a=10 (...) | a=100 (...)
     t1 << "w;x1(0);x2(0);"
        << "a1_x1*;a1_x2*;a1_f1*;a1_f2*;a1_calls;"
        << "a10_x1*;a10_x2*;a10_f1*;a10_f2*;a10_calls;"
        << "a100_x1*;a100_x2*;a100_f1*;a100_f2*;a100_calls\n";

     vector<double> A_values = {1.0, 10.0, 100.0};

     // Pętla po wadze w (od 0 do 1 co 0.01)
     // Uwaga: używamy integera do pętli, żeby uniknąć błędów zaokrągleń float
     for (int i = 0; i <= 100; ++i) {
         double w = i / 100.0;

         // 1. Losujemy punkt startowy (wspólny dla wszystkich 'a' w tym wierszu)
         matrix x0(2, 1);
         x0(0) = -10.0 + (double)rand() / RAND_MAX * 20.0;
         x0(1) = -10.0 + (double)rand() / RAND_MAX * 20.0;

         // Zapisujemy początek wiersza: w oraz punkt startowy
         t1 << w << ";" << x0(0) << ";" << x0(1);

         // 2. Obliczamy optymalizację dla a=1, a=10, a=100
         for (double a_val : A_values) {
             matrix ud1(1, 1); ud1(0) = w;
             matrix ud2(1, 1); ud2(0) = a_val;

             solution::clear_calls();
             // Startujemy zawsze z tego samego x0 w ramach wiersza
             solution sol = Powell(ff5T, x0, epsilon, Nmax, ud1, ud2);

             // Obliczenie składowych f1 i f2 w punkcie optymalnym
             double x1_opt = sol.x(0);
             double x2_opt = sol.x(1);
             double f1 = a_val * (pow(x1_opt - 3.0, 2) + pow(x2_opt - 3.0, 2));
             double f2 = (1.0 / a_val) * (pow(x1_opt + 3.0, 2) + pow(x2_opt + 3.0, 2));

             // Dopisujemy wyniki dla danego 'a' do wiersza
             t1 << ";" << x1_opt << ";" << x2_opt
                << ";" << f1 << ";" << f2
                << ";" << solution::f_calls;
         }
         // Koniec wiersza w pliku
         t1 << "\n";
     }
     t1.close();
     cout << "Zapisano lab5_tabela1.csv" << endl;

     // ==========================================================
     // TABELA 2: PROBLEM RZECZYWISTY (Belka)
     // ==========================================================
     cout << "Generowanie danych do Tabeli 2 (problem rzeczywisty)..." << endl;

     ofstream t2(resolvePath("lab5_tabela2.csv"));
     if (!t2.is_open()) {
         cerr << "BLAD: Nie mozna otworzyc pliku lab5_tabela2.csv" << endl;
         return;
     }

     // w; l(0)[mm]; d(0)[mm]; l*[mm]; d*[mm]; masa*[kg]; ugiecie*[mm]; naprezenie*[MPa]; calls
     t2 << "w;l(0)_mm;d(0)_mm;l*_mm;d*_mm;masa_kg;ugiecie_mm;naprezenie_MPa;calls\n";

     // Parametry fizyczne do obliczeń "wynikowych" (poza funkcją celu)
     double ro = 8920.0;
     double P = 2000.0;
     double E = 1.2e11;

     for (int i = 0; i <= 100; ++i) {
         double w = i / 100.0;

         // Losowy punkt startowy (w metrach, bo tak działa solver)
         // l: [0.2, 1.0], d: [0.01, 0.05]
         matrix x0(2, 1);
         x0(0) = 0.2 + (double)rand() / RAND_MAX * 0.8;
         x0(1) = 0.01 + (double)rand() / RAND_MAX * 0.04;

         // Dane do tabeli startowej (konwersja na mm)
         double l0_mm = x0(0) * 1000.0;
         double d0_mm = x0(1) * 1000.0;

         // Optymalizacja
         matrix ud1(1, 1); ud1(0) = w;
         matrix ud2 = NAN;

         solution::clear_calls();
         solution sol = Powell(ff5R, x0, epsilon, Nmax, ud1, ud2);

         // Wyniki optymalne (w metrach)
         double l_opt = sol.x(0);
         double d_opt = sol.x(1);

         // Obliczenie parametrów fizycznych dla tabeli
         double mass = ro * l_opt * M_PI * pow(d_opt / 2.0, 2);
         double deflection = (64.0 * P * pow(l_opt, 3)) / (3.0 * E * M_PI * pow(d_opt, 4));
         double stress = (32.0 * P * l_opt) / (M_PI * pow(d_opt, 3));

         // Konwersja jednostek do tabeli
         double l_star_mm = l_opt * 1000.0;
         double d_star_mm = d_opt * 1000.0;
         double u_star_mm = deflection * 1000.0;
         double stress_MPa = stress / 1e6;

         // Zapis wiersza
         t2 << w << ";"
            << l0_mm << ";" << d0_mm << ";"
            << l_star_mm << ";" << d_star_mm << ";"
            << mass << ";" << u_star_mm << ";" << stress_MPa << ";"
            << solution::f_calls << "\n";
     }

     t2.close();
     cout << "Zapisano lab5_tabela2.csv" << endl;
     cout << "Gotowe. Pliki mozna otworzyc w Excelu." << endl;
 }

void lab6() {
}
void lab4_wykresy() {
    cout << "Generowanie danych do wykresow (trajektorie)..." << endl;
    srand(time(nullptr));

    // 1. Wybór jednego punktu startowego dla wszystkich wykresów
    matrix x0(2, 1);
    // Możesz tu wpisać konkretne wartości, żeby ładnie wyglądało na wykresie
    // np. x0(0) = -1.5; x0(1) = 1.5;
    x0(0) = -0.5 + (double)rand()/RAND_MAX * 3.0; // Losowy z [-0.5, 2.5]
    x0(1) = -0.5 + (double)rand()/RAND_MAX * 3.0;

    cout << "Punkt startowy: [" << x0(0) << ", " << x0(1) << "]" << endl;

    ofstream file(resolvePath("lab4_trajektorie.csv"));
    // Nagłówek: Metoda, Krok, Iteracja, x1, x2
    file << "Metoda,Krok,Nr_Iteracji,x1,x2\n";

    double steps[] = {0.05, 0.25, NAN};
    string method_names[] = {"SD", "CG", "Newton"};
    int Nmax = 200; // Do wykresów nie potrzebujemy tysięcy punktów
    double epsilon = 1e-5;

    // Pętle po metodach i krokach
    for (int m = 0; m < 3; ++m) { // 0=SD, 1=CG, 2=Newton
        for (double h_val : steps) {

            string step_name = isnan(h_val) ? "Zmienny" : to_string(h_val);
            string method_name = method_names[m];

            // Inicjalizacja dla pojedynczego przebiegu
            matrix x = x0;
            matrix d, g, g_prev, H;
            double h;
            solution::clear_calls(); // Reset liczników (opcjonalne tutaj)

            // Zapis punktu startowego (iteracja 0)
            file << method_name << "," << step_name << ",0," << x(0) << "," << x(1) << "\n";

            // Symulacja algorytmów "ręcznie" z zapisem każdego kroku
            try {
                // Pierwszy gradient (potrzebny dla wszystkich)
                g = gf4T(x);
                d = -g; // Dla SD i pierwszego kroku CG/Newton kierunek to -gradient

                int k = 0;
                while (k < Nmax) {
                    // Wyznaczanie kierunku d w zależności od metody
                    if (m == 0) { // SD
                        d = -gf4T(x);
                    }
                    else if (m == 1) { // CG
                        if (k > 0) {
                            g_prev = g;
                            g = gf4T(x);
                            double beta = pow(norm(g), 2) / pow(norm(g_prev), 2);
                            d = -g + beta * d;
                        } else {
                            // Pierwszy krok CG to SD
                            g = gf4T(x);
                            d = -g;
                        }
                    }
                    else if (m == 2) { // Newton
                        g = gf4T(x);
                        H = Hf4T(x);
                        d = -inv(H) * g;
                    }

                    // Wyznaczanie długości kroku h
                    if (isnan(h_val)) {
                        h = golden_search_local(ff4T, x, d, NAN, NAN, epsilon, 100);
                    } else {
                        h = h_val;
                    }

                    // Wykonanie kroku
                    matrix x_prev = x;
                    x = x + h * d;

                    // Zapis do pliku
                    k++;
                    file << method_name << "," << step_name << "," << k << "," << x(0) << "," << x(1) << "\n";

                    // Warunek stopu
                    if (norm(x - x_prev) < epsilon) break;

                    // Zabezpieczenie dla stałego kroku (żeby nie uciekło w nieskończoność)
                    if (norm(x) > 100.0) break;
                }
            } catch (string& ex) {
                // Ignorujemy błędy (np. osobliwa macierz), przerywamy trasę
            } catch (...) {
                // Inne błędy
            }
        }
    }

    file.close();
    cout << "Plik lab4_trajektorie.csv zostal wygenerowany." << endl;
    cout << "Uzyj tego pliku do narysowania wykresow w Excelu (Wstaw -> Wykres Punktowy)." << endl;
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
