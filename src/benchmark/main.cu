#include <ctime>
#include <iostream>
#include <iomanip>
#include <cfloat>

#ifdef _WIN32
#include <chrono>
#include <Windows.h>
#else
#include <sys/time.h>
#include <ctime>
#endif

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include "cpu_base.h"
#include "util.h"

#include "tools.h"
#include "type.h"
#include "constants.h"
#include "redutil2.h"

#define NDIM   3        // Number of space dimension
#define NVPO   6        // Number of variables per object (3 space and 3 velocity coordinates)

using namespace std;
using namespace redutil2;

// Remove if already defined
typedef unsigned long long uint64;

static string method_name[] = { "base", "base with sym.", "tile", "tile advanced" };
static string param_name[] = { "n_body", "snk src" };

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

    //// compute norm square of d vector [5 FLOPS] [3 read, 1 write]
    //var_t d2 = SQR(r_ij.x) + SQR(r_ij.y) + SQR(r_ij.z);
    //var_t d = sqrt(d2);
    //var_t s = K2 * mj / (d2 * d);

    // compute norm square of d vector [5 FLOPS] [3 read, 1 write]
    var_t d2 = SQR(r_ij.x) + SQR(r_ij.y) + SQR(r_ij.z);
    d2 = K2 * mj / (sqrt(d2) * d2);

    ai.x += d2 * r_ij.x;
    ai.y += d2 * r_ij.y;
    ai.z += d2 * r_ij.z;
} /* body_body_grav_accel() */

__global__
void calc_grav_accel_naive(uint32_t n_obj, const var3_t* r, const nbp_t::param_t* p, var3_t* a)
{
    // i is the index of the SINK body
    const uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i < n_obj)
    {
        a[i].x = a[i].y = a[i].z = 0.0;

        var3_t r_ij = { 0.0, 0.0, 0.0 };
        // j is the index of the SOURCE body
        for (uint32_t j = 0; j < n_obj; j++)
        {
            if (i == j) continue;
            kernel::body_body_grav_accel(r[i], r[j], p[j].mass, a[i]);
            //r_ij.x = r[j].x - r[i].x;
            //r_ij.y = r[j].y - r[i].y;
            //r_ij.z = r[j].z - r[i].z;

            //// compute norm square of d vector [5 FLOPS] [3 read, 1 write]
            //var_t d2 = SQR(r_ij.x) + SQR(r_ij.y) + SQR(r_ij.z);
            //var_t d = sqrt(d2);
            //var_t s = K2 * p[j].mass / (d * d2);

            //// 6 FLOP
            //a[i].x += s * r_ij.x;
            //a[i].y += s * r_ij.y;
            //a[i].z += s * r_ij.z;
        } // 36 FLOP
    }
} /* calc_grav_accel_naive () */

__global__
void calc_grav_accel_tile_verbose(uint32_t n_obj, const var3_t* r, const nbp_t::param_t* p, var3_t* a)
{
    extern __shared__ var3_t sh_pos[];

    const uint32_t i = blockIdx.x*blockDim.x + threadIdx.x;

    var3_t acc = { 0.0, 0.0, 0.0 };
    var3_t my_pos;

    if (0 == i)
    {
        printf("[i = 0]   gridDim = [%3d, %3d %3d]\n", gridDim.x, gridDim.y, gridDim.z);
        printf("[i = 0]  blockDim = [%3d, %3d %3d]\n", blockDim.x, blockDim.y, blockDim.z);
        printf("[i = 0]: nThreads = gridDim.x * blockDim.x = %3d\n", gridDim.x * blockDim.x);
    }

    if (i < n_obj)
    {
        if (0 == threadIdx.x)
        {
            printf("[0 == threadIdx.x]: blockIdx.x = %3d\n", blockIdx.x);
        }
        if (0 == blockIdx.x)
        {
            printf("[0 == blockIdx.x]: threadIdx.x = %3d\n", threadIdx.x);
        }
    }
}
} /* namespace kernel */

void benchmark_CPU(uint32_t n_obj, const var_t* h_y, const var_t* h_p, var_t* h_dy, ofstream& o_result)
{
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
    uint64 t0 = GetTimeMs64();
#endif
    //Naive method
    if (100 >= n_obj)
    {
        for (i = 0; i < 100; i++)
        {
            cpu_calc_grav_accel(t, n_obj, r, p, a, false);
        }
    }
    else if (100 < n_obj && 1000 >= n_obj)
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
            cpu_calc_grav_accel(t, n_obj, r, p, a, true);
        }
    }
    else if (100 < n_obj && 1000 >= n_obj)
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
    uint64 t0 = GetTimeMs64();
#endif
    //Naive method
    if (100 >= n_obj)
    {
        for (i = 0; i < 100; i++)
        {
            cpu_calc_grav_accel(t, n_obj, snk, src, r, p, a, false);
        }
    }
    else if (100 < n_obj && 1000 >= n_obj)
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
            cpu_calc_grav_accel(t, n_obj, snk, src, r, p, a, true);
        }
    }
    else if (100 < n_obj && 1000 >= n_obj)
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

float2 gpu_calc_grav_accel_naive(uint32_t n_obj, unsigned int max_n_tpb, const var_t* d_y, const var_t* d_p, var_t* d_dy)
{
    static bool first_call = true;
    static uint32_t n_last;
    static unsigned int opt_n_tpb;

    float2 result = { 0.0f, FLT_MAX };

    cudaEvent_t start, stop;
    CUDA_SAFE_CALL(cudaEventCreate(&start));
    CUDA_SAFE_CALL(cudaEventCreate(&stop));

    if (first_call)
    {
        n_last = n_obj;
        opt_n_tpb = 16;
    }

    // Number of space and velocity coordinates
    const uint32_t nv = NDIM * n_obj;

    // Create aliases
    const var3_t* r = (var3_t*)d_y;
    const nbp_t::param_t* p = (nbp_t::param_t*)d_p;
    var3_t* a = (var3_t*)(d_dy + nv);

    if (first_call || n_last != n_obj)
    {
        for (unsigned int n_tpb = 16; n_tpb <= max_n_tpb; n_tpb += 16)
        {
            dim3 grid((n_obj + n_tpb - 1) / n_tpb);
            dim3 block(n_tpb);

            CUDA_SAFE_CALL(cudaEventRecord(start, 0));

            kernel::calc_grav_accel_naive<<< grid, block >>>(n_obj, r, p, a);
            CUDA_CHECK_ERROR();

            CUDA_SAFE_CALL(cudaEventRecord(stop, 0));
            CUDA_SAFE_CALL(cudaEventSynchronize(stop));

            float elapsed_time = 0.0f;
            // Computes the elapsed time between two events
            // (in milliseconds with a resolution of around 0.5 microseconds).
            CUDA_SAFE_CALL(cudaEventElapsedTime(&elapsed_time, start, stop));

            printf("    %4d %12.4e\n", n_tpb, elapsed_time);

            if (elapsed_time < result.y)
            {
                result.x = n_tpb;
                result.y = elapsed_time;
            }
        }
        opt_n_tpb = result.x;
    }
    else
    {
        dim3 grid((n_obj + opt_n_tpb - 1) / opt_n_tpb);
        dim3 block(opt_n_tpb);

        CUDA_SAFE_CALL(cudaEventRecord(start, 0));

        kernel::calc_grav_accel_naive<<< grid, block >>>(n_obj, r, p, a);
        CUDA_CHECK_ERROR();

        CUDA_SAFE_CALL(cudaEventRecord(stop, 0));
        CUDA_SAFE_CALL(cudaEventSynchronize(stop));

        float elapsed_time = 0.0f;
        // Computes the elapsed time between two events
        // (in milliseconds with a resolution of around 0.5 microseconds).
        CUDA_SAFE_CALL(cudaEventElapsedTime(&elapsed_time, start, stop));

        result.x = opt_n_tpb;
        result.y = elapsed_time;
    }
    first_call = false;

    return result;
}

void benchmark_GPU(int id_dev, uint32_t n_obj, const var_t* d_y, const var_t* d_p, var_t* d_dy, ofstream& o_result)
{
    interaction_bound int_bound;
    var_t Dt_CPU = 0.0;

    cudaDeviceProp deviceProp;
    CUDA_SAFE_CALL(cudaGetDeviceProperties(&deviceProp, id_dev));

    float2 result = gpu_calc_grav_accel_naive(n_obj, deviceProp.maxThreadsPerBlock, d_y, d_p, d_dy);
    int n_tpb = result.x;
    var_t Dt_GPU = result.y / 1.0e3;

    print(PROC_UNIT_GPU, method_name[0], param_name[0], int_bound, n_obj, n_tpb, Dt_CPU, Dt_GPU, o_result, true);
}

void benchmark(option& opt, ofstream& o_result)
{
    static const string header_str = "date       time     dev  method_name             param_name        n_snk  n_src  n_bdy  n_tpb Dt_CPU[s]     Dt_GPU[ms]";

    var_t* h_y = 0x0;
    var_t* h_dy = 0x0;
    var_t* h_p = 0x0;
    nbp_t::metadata_t* h_md = 0x0;

    // The pseudo-random number generator is initialized using the argument passed as seed.
    // Used by the subsequent rand() function calls
    o_result << header_str << endl;

    uint32_t seed = (uint32_t)time(NULL);
    cout << "The seed number is " << seed << endl;
    srand(seed);
    if (PROC_UNIT_CPU == opt.comp_dev.proc_unit)
    {
        cout << "CPU Gravity acceleration:" << endl;
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

        var_t* d_y = 0x0;
        var_t* d_dy = 0x0;
        var_t* d_p = 0x0;
        nbp_t::metadata_t* d_md = 0x0;

        for (uint32_t i = opt.n0; i <= opt.n1; i *= opt.dn)
        {
            const uint16_t n_ppo = sizeof(nbp_t::param_t) / sizeof(var_t);
            const uint32_t n_var = i * NVPO;
            const uint32_t n_par = i * n_ppo;

            allocate_host_storage(i, &h_y, &h_dy, &h_p, &h_md);
            populate(seed, i, h_y, h_p, h_md);
            allocate_device_storage(i, &d_y, &d_dy, &d_p, &d_md);

            redutil2::copy_vector_to_device(d_y, h_y, n_var * sizeof(var_t));
            redutil2::copy_vector_to_device(d_p, h_p, n_par * sizeof(var_t));
            redutil2::copy_vector_to_device(d_md, h_md, i * sizeof(nbp_t::metadata_t));

            //uint2_t snk;
            //uint2_t src;
            //snk.n1 = 0, snk.n2 = i;
            //src.n1 = 0, src.n2 = i;
            benchmark_GPU(opt.id_dev, i, d_y, d_p, d_dy, o_result);
            //benchmark_GPU(i, snk, src, h_y, h_p, h_dy, o_result);

            deallocate_host_storage(&h_y, &h_dy, &h_p, &h_md);
            deallocate_device_storage(&d_y, &d_dy, &d_p, &d_md);
        }
    }

    cout << "Done" << endl;
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

void compare(option& opt)
{
    const uint32_t n_obj = opt.n0;

    const uint16_t n_ppo = sizeof(nbp_t::param_t) / sizeof(var_t);
    const uint32_t n_var = n_obj * NVPO;
    const uint32_t n_par = n_obj * n_ppo;
    // Number of space and velocity coordinates
    const uint32_t nv = NDIM * n_obj;

    var_t* h_y = 0x0;
    var_t* h_dy1 = 0x0;
    var_t* h_dy2 = 0x0;
    var_t* h_p = 0x0;
    nbp_t::metadata_t* h_md = 0x0;

    uint32_t seed = (uint32_t)time(NULL);
    cout << "The seed number is " << seed << endl;
    // The pseudo-random number generator is initialized using the argument passed as seed.
    // Used by the subsequent rand() function calls
    srand(seed);
    if (PROC_UNIT_CPU == opt.comp_dev.proc_unit)
    {
        cout << "Compare CPU gravity acceleration results:" << endl;
        printf("n_obj = %6d\n\n", n_obj);

        allocate_host_storage(n_obj, &h_y, &h_dy1, &h_p, &h_md);
        ALLOCATE_HOST_VECTOR((void**)(&h_dy2), n_var * sizeof(var_t));

        populate(seed, n_obj, h_y, h_p, h_md);

        uint2_t snk;
        uint2_t src;
        snk.n1 = 0, snk.n2 = n_obj;
        src.n1 = 0, src.n2 = n_obj;

        // Create aliases
        const var3_t* r = (var3_t*)h_y;
        const nbp_t::param_t* p = (nbp_t::param_t*)h_p;
        var3_t* a1 = (var3_t*)(h_dy1 + nv);
        var3_t* a2 = (var3_t*)(h_dy2 + nv);

        //cpu_calc_grav_accel(n_obj, h_y, h_p, h_dy1, false);
        //cpu_calc_grav_accel(n_obj, h_y, h_p, h_dy2, true);
        cpu_calc_grav_accel(0.0, n_obj, snk, src, r, p, a1, false);
        cpu_calc_grav_accel(0.0, n_obj, snk, src, r, p, a2, true);

        bool success = compare(n_obj, opt.tol, a1, a2);
        if (success)
        {
            printf("Success.\n");
        }

        deallocate_host_storage(&h_y, &h_dy1, &h_p, &h_md);
        FREE_HOST_VECTOR((void **)(h_dy2));
    }
    else
    {
        cout << "Compare GPU gravity acceleration results:" << endl;
        printf("n_obj = %6d\n\n", n_obj);

        var_t* d_y = 0x0;
        var_t* d_dy = 0x0;
        var_t* d_p = 0x0;
        nbp_t::metadata_t* d_md = 0x0;

        allocate_host_storage(n_obj, &h_y, &h_dy1, &h_p, &h_md);
        ALLOCATE_HOST_VECTOR((void**)(&h_dy2), n_var * sizeof(var_t));
        populate(seed, n_obj, h_y, h_p, h_md);
        allocate_device_storage(n_obj, &d_y, &d_dy, &d_p, &d_md);

        redutil2::copy_vector_to_device(d_y, h_y, n_var * sizeof(var_t));
        redutil2::copy_vector_to_device(d_p, h_p, n_par * sizeof(var_t));
        redutil2::copy_vector_to_device(d_md, h_md, n_obj * sizeof(nbp_t::metadata_t));

        //uint2_t snk;
        //uint2_t src;
        //snk.n1 = 0, snk.n2 = n_obj;
        //src.n1 = 0, src.n2 = n_obj;

        {
            // Create aliases
            const var3_t* r = (var3_t*)h_y;
            const nbp_t::param_t* p = (nbp_t::param_t*)h_p;
            var3_t* a1 = (var3_t*)(h_dy1 + nv);

            cpu_calc_grav_accel(0.0, n_obj, r, p, a1, true);
            //cpu_calc_grav_accel(n_obj, snk, src, h_y, h_p, h_dy1, false);
            //cpu_calc_grav_accel(n_obj, snk, src, h_y, h_p, h_dy2, true);
        }
        {
            // Create aliases
            const var3_t* r = (var3_t*)d_y;
            const nbp_t::param_t* p = (nbp_t::param_t*)d_p;
            var3_t* a = (var3_t*)(d_dy + nv);

            dim3 grid((n_obj + 256 - 1) / 256);
            dim3 block(256);
            kernel::calc_grav_accel_naive << < grid, block >> > (n_obj, r, p, a);
            CUDA_CHECK_ERROR();
            redutil2::copy_vector_to_host(h_dy2, d_dy, n_var * sizeof(var_t));
        }
        const var3_t* a1 = (var3_t*)(h_dy1 + nv);
        const var3_t* a2 = (var3_t*)(h_dy2 + nv);

        bool success = compare(n_obj, opt.tol, a1, a2);
        if (success)
        {
            printf("Success.\n");
        }

        deallocate_host_storage(&h_y, &h_dy1, &h_p, &h_md);
        FREE_HOST_VECTOR((void **)(h_dy2));
        deallocate_device_storage(&d_y, &d_dy, &d_p, &d_md);
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
        else if (p == "-cpu")
        {
            opt.comp_dev.proc_unit = PROC_UNIT_CPU;
        }
        else if (p == "-gpu")
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
            printf("\n\t-cpu               : the benchmark will be carry on the CPU\n");
            printf("\n\t-gpu               : the benchmark will be carry on the GPU\n");
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
        case JOB_NAME_BENCMARK_GPU:
            open_stream(opt.o_dir, result_filename, output, BENCHMARK_OUTPUT_NAME_RESULT);
            benchmark(opt, *output[BENCHMARK_OUTPUT_NAME_RESULT]);
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
