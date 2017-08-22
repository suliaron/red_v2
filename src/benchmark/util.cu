#include <algorithm>
#include <iomanip>
#include <iostream>
#include <fstream>

#include "constants.h"
#include "redutil2.h"
#include "util.h"

#define NDIM   3        // Number of space dimension
#define NVPO   6        // Number of variables per object (3 space and 3 velocity coordinates)

using namespace std;
using namespace redutil2;

void create_filename(option_t& opt, string& filename)
{
    const char sep = '_';

    string cuda_dev_name;
    string cpu_name;

    if (PROC_UNIT_GPU == opt.comp_dev.proc_unit)
    {
        cuda_dev_name = get_device_name(opt.id_dev);
        std::replace(cuda_dev_name.begin(), cuda_dev_name.end(), ' ', '_');
    }
    else
    {
        cpu_name = "CPU";
    }
    filename += opt.base_fn;

    if (PROC_UNIT_GPU == opt.comp_dev.proc_unit)
    {
        filename += sep + cuda_dev_name;
    }
    else
    {
        filename += sep + cpu_name;
    }

    filename += ".txt";
}

void open_stream(string& o_dir, string& filename, ofstream** output, benchmark_output_name_t name)
{
    string path = file::combine_path(o_dir, filename);
    output[name] = new ofstream(path.c_str(), ios::out | ios::app);
    if (!*output[name])
    {
        throw string("Cannot open " + path + ".");
    }
}

void print(proc_unit_t proc_unit, string& method_name, string& param_name, uint2_t snk, uint2_t src, int n_body, int n_tpb, var_t Dt_CPU, var_t Dt_GPU, ofstream& sout, bool prn_to_scr)
{
    static const char*proc_unit_name[] =
    {
        "CPU",
        "GPU"
    };
    static char sep = ',';

    if (prn_to_scr)
    {
        cout << tools::get_time_stamp() << sep
            << setw(4) << proc_unit_name[proc_unit] << sep
            << setw(20) << method_name << sep
            << setw(20) << param_name << sep
            << setw(6) << snk.n2 - snk.n1 << sep
            << setw(6) << src.n2 - src.n1 << sep
            << setw(6) << n_body << sep
            << setw(5) << n_tpb << sep
            << scientific << setprecision(4) << setw(12) << Dt_CPU << sep
            << scientific << setprecision(4) << setw(12) << Dt_GPU << endl;
    }

    sout << tools::get_time_stamp() << SEP
        << setw(4) << proc_unit_name[proc_unit] << SEP
        << setw(20) << method_name << SEP
        << setw(20) << param_name << SEP
        << setw(6) << snk.n2 - snk.n1 << SEP
        << setw(6) << src.n2 - src.n1 << SEP
        << setw(6) << n_body << SEP
        << setw(5) << n_tpb << SEP
        << scientific << setprecision(4) << setw(12) << Dt_CPU << SEP
        << scientific << setprecision(4) << setw(12) << Dt_GPU << endl;
}

void allocate_host_storage(uint32_t n_obj, var_t** h_y, var_t** h_dy, var_t** h_p, nbp_t::metadata_t** h_md)
{
    const uint16_t n_ppo = sizeof(nbp_t::param_t) / sizeof(var_t);
    const uint32_t n_var = n_obj * NVPO;
    const uint32_t n_par = n_obj * n_ppo;

    ALLOCATE_HOST_ARRAY((void**)(h_y), n_var * sizeof(var_t));
    ALLOCATE_HOST_ARRAY((void**)(h_dy), n_var * sizeof(var_t));
    ALLOCATE_HOST_ARRAY((void**)(h_p), n_par * sizeof(var_t));
    ALLOCATE_HOST_ARRAY((void**)(h_md), n_obj * sizeof(nbp_t::metadata_t));

}

void allocate_device_storage(uint32_t n_obj, var_t** d_y, var_t** d_a, var_t** d_p, nbp_t::metadata_t** d_md)
{
    const uint16_t n_ppo = sizeof(nbp_t::param_t) / sizeof(var_t);
    const uint32_t n_var = n_obj * NVPO;
    const uint32_t n_par = n_obj * n_ppo;

    ALLOCATE_DEVICE_ARRAY((void**)(d_y), n_var * sizeof(var_t));
    ALLOCATE_DEVICE_ARRAY((void**)(d_a), n_var * sizeof(var_t));
    ALLOCATE_DEVICE_ARRAY((void**)(d_p), n_par * sizeof(var_t));
    ALLOCATE_DEVICE_ARRAY((void**)(d_md), n_obj * sizeof(nbp_t::metadata_t));
}

void deallocate_host_storage(var_t** h_y, var_t** h_dy, var_t** h_p, nbp_t::metadata_t** h_md)
{
    FREE_HOST_ARRAY((void **)(h_y));
    FREE_HOST_ARRAY((void **)(h_dy));
    FREE_HOST_ARRAY((void **)(h_p));
    FREE_HOST_ARRAY((void **)(h_md));
}

void deallocate_device_storage(var_t** d_y, var_t** d_a, var_t** d_p, nbp_t::metadata_t** d_md)
{
    FREE_DEVICE_ARRAY((void **)(d_y));
    FREE_DEVICE_ARRAY((void **)(d_a));
    FREE_DEVICE_ARRAY((void **)(d_p));
    FREE_DEVICE_ARRAY((void **)(d_md));
}

void populate(uint32_t n_obj, var_t* h_y, var_t* h_p, nbp_t::metadata_t* h_md)
{
    // Number of space and velocity coordinates
    const uint32_t nv = NDIM * n_obj;

    // Create aliases
    var3_t* r = (var3_t*)h_y;
    var3_t* v = (var3_t*)(h_y + nv);
    nbp_t::param_t* p = (nbp_t::param_t*)h_p;

    for (uint32_t i = 0; i < n_obj; i++)
    {
        r[i].x = i;     // [AU]
        r[i].y = i;
        r[i].z = i;

        v[i].y = i;  // [AU/day]
        v[i].z = i;
        v[i].x = i;

        p[i].mass = i + 1.0;  // [solar mass]
        p[i].density = (1.0 + 2.0 * ((var_t)rand() / RAND_MAX)) * constants::GramPerCm3ToSolarPerAu3;
        p[i].radius = tools::calc_radius(p[i].mass, p[i].density);

        h_md[i].active = true;
        h_md[i].body_type = BODY_TYPE_STAR;
        h_md[i].id = i + 1;
        h_md[i].mig_stop_at = 0.0;
        h_md[i].mig_type = MIGRATION_TYPE_NO;
        h_md[i].unused1 = h_md[i].unused2 = h_md[i].unused3 = false;
    }
}

void populate(uint32_t seed, uint32_t n_obj, var_t* h_y, var_t* h_p, nbp_t::metadata_t* h_md)
{
    // Number of space and velocity coordinates
    const uint32_t nv = NDIM * n_obj;

    // Create aliases
    var3_t* r = (var3_t*)h_y;
    var3_t* v = (var3_t*)(h_y + nv);
    nbp_t::param_t* p = (nbp_t::param_t*)h_p;

    srand(seed);
    for (uint32_t i = 0; i < n_obj; i++)
    {
        r[i].x = -50.0 + ((var_t)rand() / RAND_MAX) * 100.0;     // [AU]
        r[i].y = -50.0 + ((var_t)rand() / RAND_MAX) * 100.0;
        r[i].z = -50.0 + ((var_t)rand() / RAND_MAX) * 100.0;

        v[i].y = -1.0e-2 + ((var_t)rand() / RAND_MAX) * 2.0e-2;  // [AU/day]
        v[i].z = -1.0e-2 + ((var_t)rand() / RAND_MAX) * 2.0e-2;
        v[i].x = -1.0e-2 + ((var_t)rand() / RAND_MAX) * 2.0e-2;

        p[i].mass = ((var_t)rand() / RAND_MAX);  // [solar mass]
        p[i].density = (1.0 + 2.0 * ((var_t)rand() / RAND_MAX)) * constants::GramPerCm3ToSolarPerAu3;
        p[i].radius = tools::calc_radius(p[i].mass, p[i].density);

        h_md[i].active = true;
        h_md[i].body_type = BODY_TYPE_STAR;
        h_md[i].id = i + 1;
        h_md[i].mig_stop_at = 0.0;
        h_md[i].mig_type = MIGRATION_TYPE_NO;
        h_md[i].unused1 = h_md[i].unused2 = h_md[i].unused3 = false;
    }
}

void populate(uint32_t seed, uint32_t n_si, uint32_t n_nsi, uint32_t n_ni, var_t* h_y, var_t* h_p, nbp_t::metadata_t* h_md)
{
    const uint32_t n_obj = n_si + n_nsi + n_ni;
    // Number of space and velocity coordinates
    const uint32_t nv = NDIM * n_obj;

    // Create aliases
    var3_t* r = (var3_t*)h_y;
    var3_t* v = (var3_t*)(h_y + nv);
    nbp_t::param_t* p = (nbp_t::param_t*)h_p;

    srand(seed);
    for (uint32_t i = 0; i < n_si; i++)
    {
        r[i].x = -50.0 + ((var_t)rand() / RAND_MAX) * 100.0;     // [AU]
        r[i].y = -50.0 + ((var_t)rand() / RAND_MAX) * 100.0;
        r[i].z = -50.0 + ((var_t)rand() / RAND_MAX) * 100.0;

        v[i].y = -1.0e-2 + ((var_t)rand() / RAND_MAX) * 2.0e-2;  // [AU/day]
        v[i].z = -1.0e-2 + ((var_t)rand() / RAND_MAX) * 2.0e-2;
        v[i].x = -1.0e-2 + ((var_t)rand() / RAND_MAX) * 2.0e-2;

        p[i].mass = ((var_t)rand() / RAND_MAX);  // [solar mass]
        p[i].density = (1.0 + 2.0 * ((var_t)rand() / RAND_MAX)) * constants::GramPerCm3ToSolarPerAu3;
        p[i].radius = tools::calc_radius(p[i].mass, p[i].density);

        h_md[i].active = true;
        h_md[i].body_type = BODY_TYPE_STAR;
        h_md[i].id = i + 1;
        h_md[i].mig_stop_at = 0.0;
        h_md[i].mig_type = MIGRATION_TYPE_NO;
        h_md[i].unused1 = h_md[i].unused2 = h_md[i].unused3 = false;
    }

    for (uint32_t i = n_si; i < n_si + n_nsi; i++)
    {
        r[i].x = -50.0 + ((var_t)rand() / RAND_MAX) * 100.0;     // [AU]
        r[i].y = -50.0 + ((var_t)rand() / RAND_MAX) * 100.0;
        r[i].z = -50.0 + ((var_t)rand() / RAND_MAX) * 100.0;

        v[i].y = -1.0e-2 + ((var_t)rand() / RAND_MAX) * 2.0e-2;  // [AU/day]
        v[i].z = -1.0e-2 + ((var_t)rand() / RAND_MAX) * 2.0e-2;
        v[i].x = -1.0e-2 + ((var_t)rand() / RAND_MAX) * 2.0e-2;

        p[i].mass = ((var_t)rand() / RAND_MAX);  // [solar mass]
        p[i].density = (1.0 + 2.0 * ((var_t)rand() / RAND_MAX)) * constants::GramPerCm3ToSolarPerAu3;
        p[i].radius = tools::calc_radius(p[i].mass, p[i].density);

        h_md[i].active = true;
        h_md[i].body_type = BODY_TYPE_PLANETESIMAL;
        h_md[i].id = i + 1;
        h_md[i].mig_stop_at = 0.0;
        h_md[i].mig_type = MIGRATION_TYPE_NO;
        h_md[i].unused1 = h_md[i].unused2 = h_md[i].unused3 = false;
    }

    for (uint32_t i = n_si + n_nsi; i < n_obj; i++)
    {
        r[i].x = -50.0 + ((var_t)rand() / RAND_MAX) * 100.0;     // [AU]
        r[i].y = -50.0 + ((var_t)rand() / RAND_MAX) * 100.0;
        r[i].z = -50.0 + ((var_t)rand() / RAND_MAX) * 100.0;

        v[i].y = -1.0e-2 + ((var_t)rand() / RAND_MAX) * 2.0e-2;  // [AU/day]
        v[i].z = -1.0e-2 + ((var_t)rand() / RAND_MAX) * 2.0e-2;
        v[i].x = -1.0e-2 + ((var_t)rand() / RAND_MAX) * 2.0e-2;

        p[i].mass = 0.0;  // [solar mass]
        p[i].density = 0.0;
        p[i].radius = 0.0;

        h_md[i].active = true;
        h_md[i].body_type = BODY_TYPE_TESTPARTICLE;
        h_md[i].id = i + 1;
        h_md[i].mig_stop_at = 0.0;
        h_md[i].mig_type = MIGRATION_TYPE_NO;
        h_md[i].unused1 = h_md[i].unused2 = h_md[i].unused3 = false;
    }
}

bool compare(uint32_t n, var_t tol, const var3_t* y1, const var3_t* y2)
{
    bool result = true;

    for (uint32_t i = 0; i < n; i++)
    {
        var_t dx = y1[i].x - y2[i].x;
        var_t dy = y1[i].y - y2[i].y;
        var_t dz = y1[i].z - y2[i].z;

        if (tol < fabs(dx) || tol < fabs(dy) || tol < fabs(dz))
        {
            printf("Error: i = %6d (%24.16le, %24.16le, %24.16le)\n", i, dx, dy, dz);
            result = false;
        }
    }
    return result;
}

int parse_options(int argc, const char **argv, option_t& opt, bool& verbose)
{
    int i = 1;

    while (i < argc)
    {
        string p = argv[i];

        if (p == "-odir")
        {
            i++;
            opt.o_dir = argv[i];
        }
        else if (p == "-bFile")
        {
            i++;
            opt.base_fn = argv[i];
        }
        else if (p == "-cpu")
        {
            opt.comp_dev.proc_unit = PROC_UNIT_CPU;
        }
        else if (p == "-gpu")
        {
            opt.comp_dev.proc_unit = PROC_UNIT_GPU;
        }
        else if (p == "-id_dev")
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
        else if (p == "-n_si")
        {
            i++;
            if (!tools::is_number(argv[i]))
            {
                throw string("Invalid number at: " + p);
            }
            opt.n_si = atoi(argv[i]);
        }
        else if (p == "-n_nsi")
        {
            i++;
            if (!tools::is_number(argv[i]))
            {
                throw string("Invalid number at: " + p);
            }
            opt.n_nsi = atoi(argv[i]);
        }
        else if (p == "-n_ni")
        {
            i++;
            if (!tools::is_number(argv[i]))
            {
                throw string("Invalid number at: " + p);
            }
            opt.n_ni = atoi(argv[i]);
        }
        else if (p == "-v" || p == "--verbose")
        {
            verbose = true;
        }
        else if (p == "-h")
        {
            printf("Usage:\n");
            printf("\n\t-cpu               : the benchmark will be carry on the CPU\n");
            printf("\n\t-gpu               : the benchmark will be carry on the GPU\n");
            printf("\n\t-id_dev <number>    : the id of the GPU to benchmark\n");
            printf("\n\t-n0 <number>       : the starting number of SI bodies\n");
            printf("\n\t-n1 <number>       : the end number of SI bodies\n");
            printf("\n\t-dn <number>       : at each iteration the number of bodies will be increased by dn\n");
            printf("\n\t-n_iter <number>   : after n_iter the value of dn will be multiplyed by a factor of 10\n");
            printf("\n\t-tol <number>      : the comparison will be done with the defined tolarance level (default value is 1.0e-16)\n");
            printf("\n\t-v                 : the detailed result of the comparison will be printed to the screen (default value is false)\n");
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
            throw string("Unknown processing unit.");
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
    opt.n_ni = 0;
    opt.n_nsi = 0;
    opt.n_si = 0;
}
