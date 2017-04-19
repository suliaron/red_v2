#include <ctime>
#include <chrono>
#include <iostream>
#include <iomanip>      // std::setw
#include <fstream>

#include "tools.h"
#include "type.h"
#include "constants.h"
#include "redutil2.h"

using namespace std;
using namespace redutil2;

typedef enum benchmark_output_name
{
    BENCHMARK_OUTPUT_NAME_LOG,
    BENCHMARK_OUTPUT_NAME_RESULT,
    BENCHMARK_OUTPUT_NAME_SUMMARY,
    BENCHMARK_OUTPUT_NAME_N
} benchmark_output_name_t;

typedef enum job_name
{
    JOB_NAME_UNDEFINED,
    JOB_NAME_COMPARE,
    JOB_NAME_BENCMARK_CPU,
    JOB_NAME_BENCMARK_GPU,
    JOB_NAME_N
} job_name_t;

typedef struct option
{
    string     o_dir;
    string     base_fn;
    comp_dev_t comp_dev;
    bool       compare;
    var_t      tol;
    int        id_dev;
    int        n0;
    int        n1;
    int        dn;
    int        n_iter;
    job_name_t job_name;
} option_t;

static string method_name[] = { "naive", "naive_sym", "tile", "tile_advanced" };
static string param_name[] = { "n_body", "interaction_bound" };


void print(proc_unit_t proc_unit, string& method_name, string& param_name, interaction_bound int_bound, int n_body, int n_tpb, var_t Dt_CPU, var_t Dt_GPU, ofstream& sout, bool prn_to_scr)
{
    static const char* computing_device_name[] =
    {
        "CPU",
        "GPU"
    };
    static char sep = ',';

    if (prn_to_scr)
    {
        cout << tools::get_time_stamp() << sep
            << setw(4) << computing_device_name[proc_unit] << sep
            << setw(20) << method_name << sep
            << setw(20) << param_name << sep
            << setw(6) << int_bound.sink.y - int_bound.sink.x << sep
            << setw(6) << int_bound.source.y - int_bound.source.x << sep
            << setw(6) << n_body << sep
            << setw(5) << n_tpb << sep
            << scientific << setprecision(4) << setw(12) << Dt_CPU << sep
            << scientific << setprecision(4) << setw(12) << Dt_GPU << endl;
    }

    sout << tools::get_time_stamp() << SEP
        << setw(4) << computing_device_name[proc_unit] << SEP
        << setw(20) << method_name << SEP
        << setw(20) << param_name << SEP
        << setw(6) << int_bound.sink.y - int_bound.sink.x << SEP
        << setw(6) << int_bound.source.y - int_bound.source.x << SEP
        << setw(6) << n_body << SEP
        << setw(5) << n_tpb << SEP
        << scientific << setprecision(4) << setw(12) << Dt_CPU << SEP
        << scientific << setprecision(4) << setw(12) << Dt_GPU << endl;
}


void cpu_calc_grav_accel_naive(uint32_t n_body, const var3_t* r, const var_t* mass, var3_t* a)
{
    for (uint32_t i = 0; i < n_body; i++)
    {
        var3_t dVec = { 0.0, 0.0, 0.0 };
        for (uint32_t j = 0; j < n_body; j++)
        {
            if (i == j)
            {
                continue;
            }
            // r_ij 3 FLOP
            dVec.x = r[j].x - r[i].x;
            dVec.y = r[j].y - r[i].y;
            dVec.z = r[j].z - r[i].z;
            // 5 FLOP
            dVec.w = SQR(dVec.x) + SQR(dVec.y) + SQR(dVec.z);	// = r2

                                                                // 20 FLOP
            var_t d = sqrt(dVec.w);								// = r
                                                                // 2 FLOP
            var_t r_3 = 1.0 / (d*dVec.w);
            // 1 FLOP
            dVec.w = mass[j] * r_3;
            // 6 FLOP

            a[i].x += dVec.w * dVec.x;
            a[i].y += dVec.w * dVec.y;
            a[i].z += dVec.w * dVec.z;
        } // 36 FLOP
    }
}

void cpu_calc_grav_accel_naive_sym(uint32_t n_body, const var3_t* r, const var_t* mass, var3_t* a)
{
    for (uint32_t i = 0; i < n_body; i++)
    {
        var3_t dVec = { 0.0, 0.0, 0.0 };
        for (uint32_t j = i + 1; j < n_body; j++)
        {
            // r_ij 3 FLOP
            dVec.x = r[j].x - r[i].x;
            dVec.y = r[j].y - r[i].y;
            dVec.z = r[j].z - r[i].z;
            // 5 FLOP
            dVec.w = SQR(dVec.x) + SQR(dVec.y) + SQR(dVec.z);	// = r2

                                                                // 20 FLOP
            var_t d = sqrt(dVec.w);								// = r
                                                                // 2 FLOP
            var_t r_3 = 1.0 / (d*dVec.w);
            // 1 FLOP
            dVec.w = mass[j] * r_3;
            // 6 FLOP
            a[i].x += dVec.w * dVec.x;
            a[i].y += dVec.w * dVec.y;
            a[i].z += dVec.w * dVec.z;

            // 2 FLOP
            dVec.w = mass[i] * r_3;
            // 6 FLOP
            a[j].x -= dVec.w * dVec.x;
            a[j].y -= dVec.w * dVec.y;
            a[j].z -= dVec.w * dVec.z;
        } // 36 + 8 = 44 FLOP
    }
}


void benchmark_CPU(int n_body, const var3_t* h_x, const var_t* h_m, var3_t* h_a, ofstream& o_result, ofstream& o_summary)
{
    int i = 0;

    var_t Dt_GPU = 0.0;
    //Naive method
    {
        chrono::time_point<chrono::system_clock> t0 = chrono::system_clock::now();
        if (50 >= n_body)
        {
            for (i = 0; i < 1000; i++)
            {
                cpu_calc_grav_accel_naive(n_body, h_x, h_m, h_a);
            }
        }
        else if (50 < n_body && 200 >= n_body)
        {
            for (i = 0; i < 100; i++)
            {
                cpu_calc_grav_accel_naive(n_body, h_x, h_m, h_a);
            }
        }
        else if (200 < n_body && 2000 >= n_body)
        {
            for (i = 0; i < 10; i++)
            {
                cpu_calc_grav_accel_naive(n_body, h_x, h_m, h_a);
            }
        }
        else
        {
            cpu_calc_grav_accel_naive(n_body, h_x, h_m, h_a);
        }
        chrono::time_point<chrono::system_clock> t1 = chrono::system_clock::now();
        chrono::duration<var_t> total_time = t1 - t0;
        var_t Dt_CPU = total_time.count();

        interaction_bound int_bound;
        print(PROC_UNIT_CPU, method_name[0], param_name[0], int_bound, n_body, 1, Dt_CPU, Dt_GPU, o_result, false);
        //print(PROC_UNIT_CPU, method_name[0], param_name[0], int_bound, n_body, 1, Dt_CPU, Dt_GPU, o_summary, true);
    }

    //Naive symmetric method
    {
        chrono::time_point<chrono::system_clock> t0 = chrono::system_clock::now();
        if (50 >= n_body)
        {
            for (i = 0; i < 1000; i++)
            {
                cpu_calc_grav_accel_naive_sym(n_body, h_x, h_m, h_a);
            }
        }
        else if (50 < n_body && 200 >= n_body)
        {
            for (i = 0; i < 100; i++)
            {
                cpu_calc_grav_accel_naive_sym(n_body, h_x, h_m, h_a);
            }
        }
        else if (200 < n_body && 2000 >= n_body)
        {
            for (i = 0; i < 10; i++)
            {
                cpu_calc_grav_accel_naive_sym(n_body, h_x, h_m, h_a);
            }
        }
        else
        {
            cpu_calc_grav_accel_naive_sym(n_body, h_x, h_m, h_a);
        }
        chrono::time_point<chrono::system_clock> t1 = chrono::system_clock::now();
        chrono::duration<var_t> total_time = t1 - t0;
        var_t Dt_CPU = total_time.count();

        interaction_bound int_bound;
        print(PROC_UNIT_CPU, method_name[1], param_name[0], int_bound, n_body, 1, Dt_CPU, Dt_GPU, o_result, false);
        //print(COMPUTING_DEVICE_CPU, method_name[1], param_name[0], int_bound, n_body, 1, Dt_CPU, Dt_GPU, o_summary, true);
    }
}


int parse_options(int argc, const char **argv, option_t& opt, bool& verbose)
{
    int i = 1;

    while (i < argc)
    {
        string p = argv[i];
        string v;

        if (p == "-oDir")
        {
            i++;
            opt.o_dir = argv[i];
        }
        else if (p == "-bFile")
        {
            i++;
            opt.base_fn = argv[i];
        }
        else if (p == "-CPU")
        {
            opt.comp_dev.proc_unit = PROC_UNIT_CPU;
        }
        else if (p == "-GPU")
        {
            opt.comp_dev.proc_unit = PROC_UNIT_GPU;
        }
        else if (p == "-devId")
        {
            i++;
            if (!tools::is_number(argv[i]))
            {
                throw string("Invalid number at: " + p);
            }
            opt.id_dev = atoi(argv[i]);
        }
        else if (p == "-tol")
        {
            i++;
            if (!tools::is_number(argv[i]))
            {
                throw string("Invalid number at: " + p);
            }
            opt.tol = atof(argv[i]);
            opt.compare = true;
        }
        else if (p == "-n0")
        {
            i++;
            if (!tools::is_number(argv[i]))
            {
                throw string("Invalid number at: " + p);
            }
            opt.n0 = atoi(argv[i]);
        }
        else if (p == "-n1")
        {
            i++;
            if (!tools::is_number(argv[i]))
            {
                throw string("Invalid number at: " + p);
            }
            opt.n1 = atoi(argv[i]);
        }
        else if (p == "-dn")
        {
            i++;
            if (!tools::is_number(argv[i]))
            {
                throw string("Invalid number at: " + p);
            }
            opt.dn = atoi(argv[i]);
        }
        else if (p == "-n_iter")
        {
            i++;
            if (!tools::is_number(argv[i]))
            {
                throw string("Invalid number at: " + p);
            }
            opt.n_iter = atoi(argv[i]);
        }
        else if (p == "-v" || p == "--verbose")
        {
            i++;
            v = argv[i];
            if (v == "true")
            {
                verbose = true;
            }
            else if (v == "false")
            {
                verbose = false;
            }
            else
            {
                throw string("Invalid value at: " + p);
            }
        }
        else if (p == "-h")
        {
            printf("Usage:\n");
            printf("\n\t-CPU               : the benchmark will be carry on the CPU\n");
            printf("\n\t-GPU               : the benchmark will be carry on the GPU\n");
            printf("\n\t-devId <number>    : the id of the GPU to benchmark\n");
            printf("\n\t-n0 <number>       : the starting number of SI bodies\n");
            printf("\n\t-n1 <number>       : the end number of SI bodies\n");
            printf("\n\t-dn <number>       : at each iteration the number of bodies will be increased by dn\n");
            printf("\n\t-n_iter <number>   : after n_iter the value of dn will be multiplyed by a factor of 10\n");
            printf("\n\t-tol <number>      : the comparison will be done with the defined tolarance level (default value is 1.0e-16)\n");
            printf("\n\t-v <true|false>    : the detailed result of the comparison will be printed to the screen (default value is false)\n");
            printf("\n\t-oDir <filename>   : the output file will be stored in this directory\n");
            printf("\n\t-bFile <filename>  : the base filename of the output (without extension), it will be extended by CPU and GPU name\n");
            printf("\n\t-h                 : print this help\n");
            exit(EXIT_SUCCESS);
        }
        else
        {
            throw string("Invalid switch on command-line: " + p + ".");
        }
        i++;
    }

    if (opt.compare)
    {
        opt.job_name = JOB_NAME_COMPARE;
    }
    else
    {
        if (PROC_UNIT_CPU == opt.comp_dev.proc_unit)
        {
            opt.job_name = JOB_NAME_BENCMARK_CPU;
        }
        else if (PROC_UNIT_GPU == opt.comp_dev.proc_unit)
        {
            opt.job_name = JOB_NAME_BENCMARK_GPU;
        }
        else
        {
            throw string("Unknown computing device.");
        }
    }

    return i;
}

void create_default_option(option_t& opt)
{
    opt.base_fn = "";
    opt.compare = false;
    opt.comp_dev.id_dev = 0;
    opt.comp_dev.proc_unit = PROC_UNIT_CPU;
    opt.dn = 1;
    opt.id_dev = 0;
    opt.job_name = JOB_NAME_UNDEFINED;
    opt.n0 = 0;
    opt.n1 = 0;
    opt.n_iter = 10;
    opt.o_dir = "";
    opt.tol = 1.0e-16;
}

int main(int argc, const char** argv, const char** env)
{
    static const string header_str = "date       time     dev  method_name             param_name        n_snk  n_src  n_bdy  n_tpb Dt_CPU[ms]    Dt_GPU[ms]";

    chrono::time_point<chrono::system_clock> start = chrono::system_clock::now();
    ofstream* slog = NULL;

    try
    {
        option_t opt;
        create_default_option(opt);


    } /* try */
    catch (const string& msg)
    {
        if (NULL != slog)
        {
            file::log_message(*slog, "ERROR: " + msg, false);
        }
        cerr << "ERROR: " << msg << endl;
    }
    chrono::time_point<chrono::system_clock> end = chrono::system_clock::now();
    chrono::duration<var_t> total_time = end - start;

    if (NULL != slog)
    {
        file::log_message(*slog, "Total time: " + tools::convert_var_t(total_time.count()) + " s.", false);
    }
    cout << "Total time: " << total_time.count() << " s." << endl;

    return (EXIT_SUCCESS);
}
