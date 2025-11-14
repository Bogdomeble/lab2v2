

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
         lab2();
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
    ofstream Sout("symulacja_lab0.csv");                // define a stream to the .csv file
    Sout << hcat(Y[0], Y[1]);                            // save results in the file
    Sout.close();                                       // close the stream
    delete[] Y;                                         // free the memory of the DE solution
}

void lab1() {
    // Initialize the random number generator to get different values each time
    srand(time(nullptr));

    // Open a CSV file to save the results.
    // The file will be created in the same folder as the executable.
    ofstream file("lab1_results.csv");

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
    ofstream sim_file_fib("lab_1_real_fibonacci.csv");
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
    ofstream sim_file_lag("lab_1_real_lagrange.csv");
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


    ofstream file("lab2_results.csv");
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
}

void lab3() {
}

void lab4() {
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
