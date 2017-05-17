#include <cstring>
#ifdef _WIN32
#include <chrono>
#include <Windows.h>
#endif
#include <math.h>

#include "macro.h"
#include "type.h"
#include "util.h"

#define NDIM   3        // Number of space dimension
#define NVPO   6        // Number of variables per object (3 space and 3 velocity coordinates)

using namespace std;

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
    var_t s = K2 * mj / (d * d2);

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

void cpu_calc_grav_accel(var_t t, uint32_t n_obj, const var3_t* r, const nbp_t::param_t* p, var3_t* a, bool use_sym)
{
    if (use_sym)
    {
        for (uint32_t i = 0; i < n_obj; i++)
        {
            var3_t r_ij = { 0.0, 0.0, 0.0 };
            for (uint32_t j = i + 1; j < n_obj; j++)
            {
                //body_body_grav_accel_sym(r[i], r[j], p[i].mass, p[j].mass, a[i], a[j]);
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
        for (uint32_t i = 0; i < n_obj; i++)
        {
            var3_t r_ij = { 0.0, 0.0, 0.0 };
            for (uint32_t j = 0; j < n_obj; j++)
            {
                if (i == j) continue;
                //body_body_grav_accel(r[i], r[j], p[j].mass, a[i]);
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

void cpu_calc_grav_accel(var_t t, uint32_t n_obj, uint2_t snk, uint2_t src, const var3_t* r, const nbp_t::param_t* p, var3_t* a, bool use_sym)
{
    if (use_sym)
    {
        for (uint32_t i = snk.n1; i < snk.n2; i++)
        {
            var3_t r_ij = { 0.0, 0.0, 0.0 };
            for (uint32_t j = i + 1; j < src.n2; j++)
            {
                //body_body_grav_accel_sym(r[i], r[j], p[i].mass, p[j].mass, a[i], a[j]);
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
                //body_body_grav_accel(r[i], r[j], p[j].mass, a[i]);
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
    static string method_name[] = { "base", "base_with_sym." };
    static string param_name[] = { "n_body", "snk_src" };

    // Number of space and velocity coordinates
    const uint32_t nv = NDIM * n_obj;

    // Create aliases
    const var3_t* r = (var3_t*)h_y;
    const nbp_t::param_t* p = (nbp_t::param_t*)h_p;
    var3_t* a = (var3_t*)(h_dy + nv);

    interaction_bound int_bound;
    var_t t = 0.0;
    int i = 0;

    var_t Dt_GPU = 0.0;
#ifdef _WIN32
    chrono::time_point<chrono::system_clock> t0 = chrono::system_clock::now();
#else
    uint64_t t0 = GetTimeMs64();
#endif
    // Base method
    memset(a, 0, n_obj * sizeof(var3_t));
    if (100 >= n_obj)
    {
        for (i = 0; i < 1000; i++)
        {
            cpu_calc_grav_accel(t, n_obj, r, p, a, false);
        }
    }
    else if (100 < n_obj && 1000 >= n_obj)
    {
        for (i = 0; i < 100; i++)
        {
            cpu_calc_grav_accel(t, n_obj, r, p, a, false);
        }
    }
    else if (1000 < n_obj && 10000 >= n_obj)
    {
        for (i = 0; i < 10; i++)
        {
            cpu_calc_grav_accel(t, n_obj, r, p, a, false);
        }
    }
    else
    {
        cpu_calc_grav_accel(t, n_obj, r, p, a, false);
    }
#ifdef _WIN32
    chrono::time_point<chrono::system_clock> t1 = chrono::system_clock::now();
    chrono::duration<var_t> total_time = t1 - t0;
    var_t Dt_CPU = total_time.count() / (var_t)(i == 0 ? 1 : i);
#else
    uint64_t t1 = GetTimeMs64();
    var_t Dt_CPU = ((var_t)(t1 - t0)) / (var_t)(i == 0 ? 1 : i) / 1.0e6;  // [sec]
#endif

    print(PROC_UNIT_CPU, method_name[0], param_name[0], int_bound, n_obj, 1, Dt_CPU, Dt_GPU, o_result, true);

#ifdef _WIN32
    t0 = chrono::system_clock::now();
#else
    t0 = GetTimeMs64();
#endif
    // Base symmetric method
    memset(a, 0, n_obj * sizeof(var3_t));
    if (100 >= n_obj)
    {
        for (i = 0; i < 1000; i++)
        {
            cpu_calc_grav_accel(t, n_obj, r, p, a, true);
        }
    }
    else if (100 < n_obj && 1000 >= n_obj)
    {
        for (i = 0; i < 100; i++)
        {
            cpu_calc_grav_accel(t, n_obj, r, p, a, true);
        }
    }
    else if (1000 < n_obj && 10000 >= n_obj)
    {
        for (i = 0; i < 10; i++)
        {
            cpu_calc_grav_accel(t, n_obj, r, p, a, true);
        }
    }
    else
    {
        cpu_calc_grav_accel(t, n_obj, r, p, a, true);
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
    static string method_name[] = { "base", "base_with_sym." };
    static string param_name[] = { "n_body", "snk_src" };

    // Number of space and velocity coordinates
    const uint32_t nv = NDIM * n_obj;

    // Create aliases
    const var3_t* r = (var3_t*)h_y;
    const nbp_t::param_t* p = (nbp_t::param_t*)h_p;
    var3_t* a = (var3_t*)(h_dy + nv);

    interaction_bound int_bound(snk, src);
    var_t t = 0.0;
    int i = 0;

    var_t Dt_GPU = 0.0;
#ifdef _WIN32
    chrono::time_point<chrono::system_clock> t0 = chrono::system_clock::now();
#else
    uint64_t t0 = GetTimeMs64();
#endif
    // Base method
    memset(a, 0, n_obj * sizeof(var3_t));
    if (100 >= n_obj)
    {
        for (i = 0; i < 1000; i++)
        {
            cpu_calc_grav_accel(t, n_obj, snk, src, r, p, a, false);
        }
    }
    else if (100 < n_obj && 1000 >= n_obj)
    {
        for (i = 0; i < 100; i++)
        {
            cpu_calc_grav_accel(t, n_obj, snk, src, r, p, a, false);
        }
    }
    else if (1000 < n_obj && 10000 >= n_obj)
    {
        for (i = 0; i < 10; i++)
        {
            cpu_calc_grav_accel(t, n_obj, snk, src, r, p, a, false);
        }
    }
    else
    {
        cpu_calc_grav_accel(t, n_obj, snk, src, r, p, a, false);
    }
#ifdef _WIN32
    chrono::time_point<chrono::system_clock> t1 = chrono::system_clock::now();
    chrono::duration<var_t> total_time = t1 - t0;
    var_t Dt_CPU = total_time.count() / (var_t)(i == 0 ? 1 : i);
#else
    uint64_t t1 = GetTimeMs64();
    var_t Dt_CPU = ((var_t)(t1 - t0)) / (var_t)(i == 0 ? 1 : i) / 1.0e6;  // [sec]
#endif

    print(PROC_UNIT_CPU, method_name[0], param_name[1], int_bound, n_obj, 1, Dt_CPU, Dt_GPU, o_result, true);

#ifdef _WIN32
    t0 = chrono::system_clock::now();
#else
    t0 = GetTimeMs64();
#endif
    // Base symmetric method
    memset(a, 0, n_obj * sizeof(var3_t));
    if (100 >= n_obj)
    {
        for (i = 0; i < 1000; i++)
        {
            cpu_calc_grav_accel(t, n_obj, snk, src, r, p, a, true);
        }
    }
    else if (100 < n_obj && 1000 >= n_obj)
    {
        for (i = 0; i < 100; i++)
        {
            cpu_calc_grav_accel(t, n_obj, snk, src, r, p, a, true);
        }
    }
    else if (1000 < n_obj && 10000 >= n_obj)
    {
        for (i = 0; i < 10; i++)
        {
            cpu_calc_grav_accel(t, n_obj, snk, src, r, p, a, true);
        }
    }
    else
    {
        cpu_calc_grav_accel(t, n_obj, snk, src, r, p, a, true);
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

#undef NDIM
#undef NVPO
