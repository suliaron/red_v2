#include <algorithm>
#include <ctime>
#include <iostream>
#include <iomanip>
#include <fstream>

#ifdef _WIN32
#include <chrono>
#include <Windows.h>
#else
#include <sys/time.h>
#include <ctime>
#endif

#include "tools.h"
#include "type.h"
#include "constants.h"
#include "redutil2.h"

#define NDIM   3        // Number of space dimension
#define NVPO   6        // Number of variables per object (3 space and 3 velocity coordinates)

using namespace std;
using namespace redutil2;

// Remove if already defined
typedef long long int64;
typedef unsigned long long uint64;

typedef enum benchmark_output_name
{
    BENCHMARK_OUTPUT_NAME_RESULT,
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
    uint32_t   n0;
    uint32_t   n1;
    uint32_t   dn;
    uint32_t   n_iter;
    job_name_t job_name;
} option_t;

static string method_name[] = { "naive", "naive_sym", "tile", "tile_advanced" };
static string param_name[] = { "n_body", "interaction_bound" };

namespace kernel
{
inline __host__ __device__
void body_body_grav_accel(const var3_t& ri, const var3_t& rj, var_t mj, var3_t& ai)
{
    var3_t r_ij = { 0.0, 0.0, 0.0 };

    // compute r_ij = r_j - r_i [3 FLOPS] [6 read, 3 write]
    r_ij.x = rj.x - ri.x;
    r_ij.y = rj.y - ri.y;
    r_ij.z = rj.z - ri.z;

    // compute norm square of d vector [5 FLOPS] [3 read, 1 write]
    var_t d2 = SQR(r_ij.x) + SQR(r_ij.y) + SQR(r_ij.z);
    // compute m_j / d^3 []
    d2 = mj / (d2 * sqrt(d2));

    ai.x += d2 * r_ij.x;
    ai.y += d2 * r_ij.y;
    ai.z += d2 * r_ij.z;
}
} /* namespace kernel */

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

inline
void body_body_grav_accel(const var3_t& ri, const var3_t& rj, var_t mj, var3_t& ai)
{
    var3_t r_ij = { 0.0, 0.0, 0.0 };

    // compute r_ij = r_j - r_i [3 FLOPS] [6 read, 3 write]
    r_ij.x = rj.x - ri.x;
    r_ij.y = rj.y - ri.y;
    r_ij.z = rj.z - ri.z;

    // compute norm square of d vector [5 FLOPS] [3 read, 1 write]
    var_t d2 = SQR(r_ij.x) + SQR(r_ij.y) + SQR(r_ij.z);
    var_t d = sqrt(d2);
    // compute K2 * m_j / d^3 []
    var_t s = K2 * mj / (d2 * d);

    ai.x += s * r_ij.x;
    ai.y += s * r_ij.y;
    ai.z += s * r_ij.z;
}

inline
void body_body_grav_accel_sym(const var3_t& ri, const var3_t& rj, var_t mi, var_t mj, var3_t& ai, var3_t& aj)
{
    var3_t r_ij = { 0.0, 0.0, 0.0 };

    // compute r_ij = r_j - r_i [3 FLOPS] [6 read, 3 write]
    r_ij.x = rj.x - ri.x;
    r_ij.y = rj.y - ri.y;
    r_ij.z = rj.z - ri.z;

    // compute norm square of d vector [5 FLOPS] [3 read, 1 write]
    var_t d2 = SQR(r_ij.x) + SQR(r_ij.y) + SQR(r_ij.z);
    var_t d = sqrt(d2);
    // compute K2 / d^3 []
    var_t d_3 = K2 / (d * d2);

    var_t s = mj * d_3;
    ai.x += s * r_ij.x;
    ai.y += s * r_ij.y;
    ai.z += s * r_ij.z;

    s = mi * d_3;
    aj.x -= s * r_ij.x;
    aj.y -= s * r_ij.y;
    aj.z -= s * r_ij.z;
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

    ALLOCATE_HOST_VECTOR((void**)(d_y), n_var * sizeof(var_t));
    ALLOCATE_HOST_VECTOR((void**)(d_a), n_var * sizeof(var_t));
    ALLOCATE_HOST_VECTOR((void**)(d_p), n_par * sizeof(var_t));
    ALLOCATE_HOST_VECTOR((void**)(d_md), n_obj * sizeof(nbp_t::metadata_t));
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

void cpu_calc_grav_accel(uint32_t n_obj, const var_t* h_y, const var_t* h_p, var_t* h_dy, bool use_sym)
{
    // Number of space and velocity coordinates
    const uint32_t nv = NDIM * n_obj;

    // Create aliases
    const var3_t* r = (var3_t*)h_y;
    const nbp_t::param_t* p = (nbp_t::param_t*)h_p;
    var3_t* a = (var3_t*)(h_dy + nv);

    // Clear the acceleration array: the += op can be used
    memset(a, 0, nv * sizeof(var_t));

    if (use_sym)
    {
        for (uint32_t i = 0; i < n_obj; i++)
        {
            var3_t r_ij = { 0.0, 0.0, 0.0 };
            for (uint32_t j = i + 1; j < n_obj; j++)
            {
                body_body_grav_accel_sym(r[i], r[j], p[i].mass, p[j].mass, a[i], a[j]);
                //r_ij.x = r[j].x - r[i].x;
                //r_ij.y = r[j].y - r[i].y;
                //r_ij.z = r[j].z - r[i].z;

                //var_t d2 = SQR(r_ij.x) + SQR(r_ij.y) + SQR(r_ij.z);
                //var_t d = sqrt(d2);
                //var_t d_3 = K2 / (d * d2);

                //var_t s = p[j].mass * d_3;
                //a[i].x += s * r_ij.x;
                //a[i].y += s * r_ij.y;
                //a[i].z += s * r_ij.z;

                //s = p[i].mass * d_3;
                //a[j].x -= s * r_ij.x;
                //a[j].y -= s * r_ij.y;
                //a[j].z -= s * r_ij.z;
            }
        }
    }
    else
    {
        for (uint32_t i = 0; i < n_obj; i++)
        {
            var3_t r_ij = { 0.0, 0.0, 0.0 };
            for (uint32_t j = 0; j < n_obj; j++)
            {
                if (i == j) continue;
                body_body_grav_accel(r[i], r[j], p[j].mass, a[i]);
                //r_ij.x = r[j].x - r[i].x;
                //r_ij.y = r[j].y - r[i].y;
                //r_ij.z = r[j].z - r[i].z;

                //var_t d2 = SQR(r_ij.x) + SQR(r_ij.y) + SQR(r_ij.z);
                //var_t d = sqrt(d2);
                //var_t s = K2 * p[j].mass / (d * d2);

                //a[i].x += s * r_ij.x;
                //a[i].y += s * r_ij.y;
                //a[i].z += s * r_ij.z;
            }
        }
    }
}

void cpu_calc_grav_accel(uint32_t n_obj, uint2_t snk, uint2_t src, const var_t* h_y, const var_t* h_p, var_t* h_dy, bool use_sym)
{
    // Number of space and velocity coordinates
    const uint32_t nv = NDIM * n_obj;

    // Create aliases
    const var3_t* r = (var3_t*)h_y;
    const nbp_t::param_t* p = (nbp_t::param_t*)h_p;
    var3_t* a = (var3_t*)(h_dy + nv);

    // Clear the acceleration array: the += op can be used
    memset(a, 0, nv * sizeof(var_t));

    if (use_sym)
    {
        for (uint32_t i = snk.n1; i < snk.n2; i++)
        {
            var3_t r_ij = { 0.0, 0.0, 0.0 };
            for (uint32_t j = i + 1; j < src.n2; j++)
            {
                r_ij.x = r[j].x - r[i].x;
                r_ij.y = r[j].y - r[i].y;
                r_ij.z = r[j].z - r[i].z;

                var_t d2 = SQR(r_ij.x) + SQR(r_ij.y) + SQR(r_ij.z);
                var_t d = sqrt(d2);
                var_t d_3 = K2 / (d * d2);

                var_t s = p[j].mass * d_3;
                a[i].x += s * r_ij.x;
                a[i].y += s * r_ij.y;
                a[i].z += s * r_ij.z;

                s = p[i].mass * d_3;
                a[j].x -= s * r_ij.x;
                a[j].y -= s * r_ij.y;
                a[j].z -= s * r_ij.z;
            }
        }
    }
    else
    {
        for (uint32_t i = snk.n1; i < snk.n2; i++)
        {
            var3_t r_ij = { 0.0, 0.0, 0.0 };
            for (uint32_t j = src.n1; j < src.n2; j++)
            {
                if (i == j) continue;
                r_ij.x = r[j].x - r[i].x;
                r_ij.y = r[j].y - r[i].y;
                r_ij.z = r[j].z - r[i].z;

                var_t d2 = SQR(r_ij.x) + SQR(r_ij.y) + SQR(r_ij.z);
                var_t d = sqrt(d2);
                var_t s = K2 * p[j].mass / (d * d2);

                a[i].x += s * r_ij.x;
                a[i].y += s * r_ij.y;
                a[i].z += s * r_ij.z;
            }
        }
    }
}

void benchmark_CPU(uint32_t n_obj, const var_t* h_y, const var_t* h_p, var_t* h_dy, ofstream& o_result)
{
    interaction_bound int_bound;
    int i = 0;

    var_t Dt_GPU = 0.0;
#ifdef _WIN32
    chrono::time_point<chrono::system_clock> t0 = chrono::system_clock::now();
#else
    uint64 t0 = GetTimeMs64();
#endif
    //Naive method
    if (100 >= n_obj)
    {
        for (i = 0; i < 100; i++)
        {
            cpu_calc_grav_accel(n_obj, h_y, h_p, h_dy, false);
        }
    }
    else if (100 < n_obj && 1000 >= n_obj)
    {
        for (i = 0; i < 10; i++)
        {
            cpu_calc_grav_accel(n_obj, h_y, h_p, h_dy, false);
        }
    }
    else
    {
        cpu_calc_grav_accel(n_obj, h_y, h_p, h_dy, false);
    }
#ifdef _WIN32
    chrono::time_point<chrono::system_clock> t1 = chrono::system_clock::now();
    chrono::duration<var_t> total_time = t1 - t0;
    var_t Dt_CPU = total_time.count() / (var_t)(i == 0 ? 1 : i);
#else
    uint64 t1 = GetTimeMs64();
    var_t Dt_CPU = ((var_t)(t1 - t0)) / (var_t)(i == 0 ? 1 : i) / 1.0e6;  // [sec]
#endif

    print(PROC_UNIT_CPU, method_name[0], param_name[0], int_bound, n_obj, 1, Dt_CPU, Dt_GPU, o_result, true);

#ifdef _WIN32
    t0 = chrono::system_clock::now();
#else
    t0 = GetTimeMs64();
#endif
    //Naive symmetric method
    if (100 >= n_obj)
    {
        for (i = 0; i < 100; i++)
        {
            cpu_calc_grav_accel(n_obj, h_y, h_p, h_dy, true);
        }
        }
    else if (100 < n_obj && 1000 >= n_obj)
    {
        for (i = 0; i < 10; i++)
        {
            cpu_calc_grav_accel(n_obj, h_y, h_p, h_dy, true);
        }
    }
    else
    {
        cpu_calc_grav_accel(n_obj, h_y, h_p, h_dy, true);
    }
#ifdef _WIN32
    t1 = chrono::system_clock::now();
    total_time = t1 - t0;
    Dt_CPU = total_time.count() / (var_t)(i == 0 ? 1 : i);
#else
    t1 = GetTimeMs64();
    Dt_CPU = ((var_t)(t1 - t0)) / (var_t)(i == 0 ? 1 : i) / 1.0e6;  // [sec]
#endif

    print(PROC_UNIT_CPU, method_name[1], param_name[0], int_bound, n_obj, 1, Dt_CPU, Dt_GPU, o_result, true);
    }

void benchmark_CPU(uint32_t n_obj, uint2_t snk, uint2_t src, const var_t* h_y, const var_t* h_p, var_t* h_dy, ofstream& o_result)
{
    interaction_bound int_bound(snk, src);
    int i = 0;

    var_t Dt_GPU = 0.0;
#ifdef _WIN32
    chrono::time_point<chrono::system_clock> t0 = chrono::system_clock::now();
#else
    uint64 t0 = GetTimeMs64();
#endif
    //Naive method
    if (100 >= n_obj)
    {
        for (i = 0; i < 100; i++)
        {
            cpu_calc_grav_accel(n_obj, snk, src, h_y, h_p, h_dy, false);
        }
    }
    else if (100 < n_obj && 1000 >= n_obj)
    {
        for (i = 0; i < 10; i++)
        {
            cpu_calc_grav_accel(n_obj, snk, src, h_y, h_p, h_dy, false);
        }
    }
    else
    {
        cpu_calc_grav_accel(n_obj, snk, src, h_y, h_p, h_dy, false);
    }
#ifdef _WIN32
    chrono::time_point<chrono::system_clock> t1 = chrono::system_clock::now();
    chrono::duration<var_t> total_time = t1 - t0;
    var_t Dt_CPU = total_time.count() / (var_t)(i == 0 ? 1 : i);
#else
    uint64 t1 = GetTimeMs64();
    var_t Dt_CPU = ((var_t)(t1 - t0)) / (var_t)(i == 0 ? 1 : i) / 1.0e6;  // [sec]
#endif

    print(PROC_UNIT_CPU, method_name[0], param_name[1], int_bound, n_obj, 1, Dt_CPU, Dt_GPU, o_result, true);

#ifdef _WIN32
    t0 = chrono::system_clock::now();
#else
    t0 = GetTimeMs64();
#endif
    //Naive symmetric method
    if (100 >= n_obj)
    {
        for (i = 0; i < 100; i++)
        {
            cpu_calc_grav_accel(n_obj, snk, src, h_y, h_p, h_dy, true);
        }
    }
    else if (100 < n_obj && 1000 >= n_obj)
    {
        for (i = 0; i < 10; i++)
        {
            cpu_calc_grav_accel(n_obj, snk, src, h_y, h_p, h_dy, true);
        }
    }
    else
    {
        cpu_calc_grav_accel(n_obj, snk, src, h_y, h_p, h_dy, true);
    }
#ifdef _WIN32
    t1 = chrono::system_clock::now();
    total_time = t1 - t0;
    Dt_CPU = total_time.count() / (var_t)(i == 0 ? 1 : i);
#else
    t1 = GetTimeMs64();
    Dt_CPU = ((var_t)(t1 - t0)) / (var_t)(i == 0 ? 1 : i) / 1.0e6;  // [sec]
#endif

    print(PROC_UNIT_CPU, method_name[1], param_name[1], int_bound, n_obj, 1, Dt_CPU, Dt_GPU, o_result, true);
}

void benchmark(option& opt, ofstream& o_result)
{
    static const string header_str = "date       time     dev  method_name             param_name        n_snk  n_src  n_bdy  n_tpb Dt_CPU[s]     Dt_GPU[ms]";

    uint32_t seed = (uint32_t)time(NULL);
    cout << "The seed number is " << seed << endl;
    // The pseudo-random number generator is initialized using the argument passed as seed.
    // Used by the subsequent rand() function calls
    o_result << header_str << endl;

    srand(seed);
    if (PROC_UNIT_CPU == opt.comp_dev.proc_unit)
    {
        cout << "CPU Gravity acceleration:" << endl;

        var_t* h_y = 0x0;
        var_t* h_dy = 0x0;
        var_t* h_p = 0x0;
        nbp_t::metadata_t* h_md = 0x0;
        for (uint32_t i = opt.n0; i <= opt.n1; i *= opt.dn)
        {
            allocate_host_storage(i, &h_y, &h_dy, &h_p, &h_md);
            populate(seed, i, h_y, h_p, h_md);

            uint2_t snk;
            uint2_t src;
            snk.n1 = 0, snk.n2 = i;
            src.n1 = 0, src.n2 = i;
            //benchmark_CPU(i, h_y, h_p, h_dy, o_result);
            benchmark_CPU(i, snk, src, h_y, h_p, h_dy, o_result);

            deallocate_host_storage(&h_y, &h_dy, &h_p, &h_md);
        }
    }
    else
    {
        cout << "GPU Gravity acceleration:" << endl;
        throw string("GPU benchmark is not yet implemented.\n");
    }

    cout << "Done" << endl;
}

void compare(option& opt)
{
    const uint32_t n_obj = opt.n0;
    const uint32_t n_var = n_obj * NVPO;

    uint32_t seed = (uint32_t)time(NULL);
    cout << "The seed number is " << seed << endl;
    // The pseudo-random number generator is initialized using the argument passed as seed.
    // Used by the subsequent rand() function calls

    srand(seed);
    if (PROC_UNIT_CPU == opt.comp_dev.proc_unit)
    {
        cout << "Compare CPU gravity acceleration results:" << endl << endl;

        var_t* h_y = 0x0;
        var_t* h_dy1 = 0x0;
        var_t* h_dy2 = 0x0;
        var_t* h_p = 0x0;
        nbp_t::metadata_t* h_md = 0x0;

        allocate_host_storage(n_obj, &h_y, &h_dy1, &h_p, &h_md);
        ALLOCATE_HOST_VECTOR((void**)(&h_dy2), n_var * sizeof(var_t));

        populate(seed, n_obj, h_y, h_p, h_md);

        uint2_t snk;
        uint2_t src;
        snk.n1 = 0, snk.n2 = n_obj;
        src.n1 = 0, src.n2 = n_obj;

        //cpu_calc_grav_accel(n_obj, h_y, h_p, h_dy1, false);
        //cpu_calc_grav_accel(n_obj, h_y, h_p, h_dy2, true);
        cpu_calc_grav_accel(n_obj, snk, src, h_y, h_p, h_dy1, false);
        cpu_calc_grav_accel(n_obj, snk, src, h_y, h_p, h_dy2, true);

        // Number of space and velocity coordinates
        const uint32_t nv = NDIM * n_obj;
        const var3_t* a1 = (var3_t*)(h_dy1 + nv);
        const var3_t* a2 = (var3_t*)(h_dy2 + nv);

        bool success = true;
        for (uint32_t i = 0; i < n_obj; i++)
        {
            var_t dx = a1[i].x - a2[i].x;
            var_t dy = a1[i].y - a2[i].y;
            var_t dz = a1[i].z - a2[i].z;

            if (opt.tol < fabs(dx) || opt.tol < fabs(dy) || opt.tol < fabs(dz))
            {
                printf("Error: i = %6d (%24.16le, %24.16le, %24.16le)\n", i, dx, dy, dz);
                success = false;
            }
        }
        if (success)
        {
            printf("Success.\n");
        }

        deallocate_host_storage(&h_y, &h_dy1, &h_p, &h_md);
        FREE_HOST_VECTOR((void **)(h_dy2));
    }
    else
    {
        ;
    }

    cout << "Done" << endl;
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
            verbose = true;
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

/*
-n0 10 -n1 100000 -dn 10 -v -odir C:\Work\red.cuda.Results\v2.0\Benchmark\Test_01 -bFile benchmark
-n0 2000 -tol 1.0e-16 -v -odir C:\Work\red.cuda.Results\v2.0\Benchmark\Test_01 -bFile compare

*/
int main(int argc, const char** argv, const char** env)
{
    option_t opt;
    string result_filename;

    bool verbose = false;

    ofstream* output[BENCHMARK_OUTPUT_NAME_N];
    memset(output, 0x0, sizeof(output));

#ifdef _WIN32
    chrono::time_point<chrono::system_clock> start = chrono::system_clock::now();
#else
    uint64 start = GetTimeMs64();
#endif
    try
    {
        create_default_option(opt);
        parse_options(argc, argv, opt, verbose);
        create_filename(opt, result_filename);

        switch (opt.job_name)
        {
        case JOB_NAME_BENCMARK_CPU:
            open_stream(opt.o_dir, result_filename, output, BENCHMARK_OUTPUT_NAME_RESULT);
            benchmark(opt, *output[BENCHMARK_OUTPUT_NAME_RESULT]);
            break;
        case JOB_NAME_BENCMARK_GPU:
            throw string("GPU benchmark is not implemented.\n");
            break;
        case JOB_NAME_COMPARE:
            compare(opt);
            break;
        default:
            break;
        }
    } /* try */
    catch (const string& msg)
    {
        cerr << "ERROR: " << msg << endl;
    }
#ifdef _WIN32
    chrono::time_point<chrono::system_clock> end = chrono::system_clock::now();
    chrono::duration<var_t> total_time = end - start;
    cout << "Total time: " << total_time.count() << " s." << endl;
#else
    uint64 end = GetTimeMs64();
    var_t Dt_CPU = ((var_t)(end - start)) / 1.0e6;
    cout << "Total time: " << Dt_CPU << " s." << endl;
#endif

    return (EXIT_SUCCESS);
}
