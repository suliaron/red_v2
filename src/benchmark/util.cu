#include <algorithm>
#include <iomanip>
#include <iostream>
#include <fstream>

#ifdef _WIN32
#include <chrono>
#include <Windows.h>
#else
#include <sys/time.h>
#include <ctime>
#endif

#include "constants.h"
#include "redutil2.h"
#include "util.h"

#define NDIM   3        // Number of space dimension
#define NVPO   6        // Number of variables per object (3 space and 3 velocity coordinates)

using namespace std;
using namespace redutil2;

/*
 * Returns the amount of microseconds elapsed since the UNIX epoch.
 * Works on windows and linux.
 */
uint64 GetTimeMs64()
{
#ifdef _WIN32
    FILETIME ft;
    LARGE_INTEGER li;

    /* Get the amount of 100 nano seconds intervals elapsed since January 1, 1601 (UTC)
    * and copy it to a LARGE_INTEGER structure.
    */
    GetSystemTimeAsFileTime(&ft);
    li.LowPart = ft.dwLowDateTime;
    li.HighPart = ft.dwHighDateTime;

    uint64 ret = li.QuadPart;
    /* Convert from file time to UNIX epoch time. */
    ret -= 116444736000000000LL;
    ret /= 10;      /* From 100 nano seconds (10^-7) to 1 microsecond (10^-6) intervals */

    return ret;
#else
    // Linux
    struct timeval tv;

    gettimeofday(&tv, NULL);
    uint64 ret = tv.tv_usec;
    /* Convert from micro seconds (10^-6) to milliseconds (10^-3) */
    //ret /= 1000;

    /* Adds the seconds (10^0) after converting them to microseconds (10^-6) */
    ret += (tv.tv_sec * 1000000);

    return ret;
#endif
}

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

void print(proc_unit_t proc_unit, string& method_name, string& param_name, interaction_bound int_bound, int n_body, int n_tpb, var_t Dt_CPU, var_t Dt_GPU, ofstream& sout, bool prn_to_scr)
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
            << setw(6) << int_bound.sink.n2 - int_bound.sink.n1 << sep
            << setw(6) << int_bound.source.n2 - int_bound.source.n1 << sep
            << setw(6) << n_body << sep
            << setw(5) << n_tpb << sep
            << scientific << setprecision(4) << setw(12) << Dt_CPU << sep
            << scientific << setprecision(4) << setw(12) << Dt_GPU << endl;
    }

    sout << tools::get_time_stamp() << SEP
        << setw(4) << proc_unit_name[proc_unit] << SEP
        << setw(20) << method_name << SEP
        << setw(20) << param_name << SEP
        << setw(6) << int_bound.sink.n2 - int_bound.sink.n1 << SEP
        << setw(6) << int_bound.source.n2 - int_bound.source.n1 << SEP
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

    ALLOCATE_HOST_VECTOR((void**)(h_y), n_var * sizeof(var_t));
    ALLOCATE_HOST_VECTOR((void**)(h_dy), n_var * sizeof(var_t));
    ALLOCATE_HOST_VECTOR((void**)(h_p), n_par * sizeof(var_t));
    ALLOCATE_HOST_VECTOR((void**)(h_md), n_obj * sizeof(nbp_t::metadata_t));
}

void allocate_device_storage(uint32_t n_obj, var_t** d_y, var_t** d_a, var_t** d_p, nbp_t::metadata_t** d_md)
{
    const uint16_t n_ppo = sizeof(nbp_t::param_t) / sizeof(var_t);
    const uint32_t n_var = n_obj * NVPO;
    const uint32_t n_par = n_obj * n_ppo;

    ALLOCATE_DEVICE_VECTOR((void**)(d_y), n_var * sizeof(var_t));
    ALLOCATE_DEVICE_VECTOR((void**)(d_a), n_var * sizeof(var_t));
    ALLOCATE_DEVICE_VECTOR((void**)(d_p), n_par * sizeof(var_t));
    ALLOCATE_DEVICE_VECTOR((void**)(d_md), n_obj * sizeof(nbp_t::metadata_t));
}

void deallocate_host_storage(var_t** h_y, var_t** h_dy, var_t** h_p, nbp_t::metadata_t** h_md)
{
    FREE_HOST_VECTOR((void **)(h_y));
    FREE_HOST_VECTOR((void **)(h_dy));
    FREE_HOST_VECTOR((void **)(h_p));
    FREE_HOST_VECTOR((void **)(h_md));
}

void deallocate_device_storage(var_t** d_y, var_t** d_a, var_t** d_p, nbp_t::metadata_t** d_md)
{
    FREE_DEVICE_VECTOR((void **)(d_y));
    FREE_DEVICE_VECTOR((void **)(d_a));
    FREE_DEVICE_VECTOR((void **)(d_p));
    FREE_DEVICE_VECTOR((void **)(d_md));
}

void populate(uint32_t seed, uint32_t n_obj, var_t* h_y, var_t* h_p, nbp_t::metadata_t* h_md)
{
    srand(seed);

    // Number of space and velocity coordinates
    const uint32_t nv = NDIM * n_obj;

    // Create aliases
    var3_t* r = (var3_t*)h_y;
    var3_t* v = (var3_t*)(h_y + nv);
    nbp_t::param_t* p = (nbp_t::param_t*)h_p;

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
