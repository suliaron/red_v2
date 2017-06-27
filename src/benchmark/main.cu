#include <ctime>
#include <iostream>
#include <iomanip>

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
#include "gpu_base.h"
#include "util.h"

#include "tools.h"
#include "type.h"
#include "constants.h"
#include "redutil2.h"

#define NDIM   3        // Number of space dimension
#define NVPO   6        // Number of variables per object (3 space and 3 velocity coordinates)

using namespace std;
using namespace redutil2;

void benchmark(option& opt, ofstream& o_result)
{
    static const string header_str = "date       time     dev  method_name             param_name        n_snk  n_src  n_body n_tpb Dt_CPU[ms]    Dt_GPU[ms]";

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
            memset(h_dy, 0, i * NVPO * sizeof(var_t));
            populate(seed, i, h_y, h_p, h_md);
            benchmark_CPU(i, h_md, h_y, h_p, h_dy, o_result);
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

        const uint16_t n_ppo = sizeof(nbp_t::param_t) / sizeof(var_t);

        for (uint32_t i = opt.n0; i <= opt.n1; i *= opt.dn)
        {

            const uint32_t n_var = i * NVPO;
            const uint32_t n_par = i * n_ppo;

            allocate_host_storage(i, &h_y, &h_dy, &h_p, &h_md);
            populate(seed, i, h_y, h_p, h_md);
            allocate_device_storage(i, &d_y, &d_dy, &d_p, &d_md);
            memset(h_dy, 0, i * NVPO * sizeof(var_t));
            cudaMemset(d_dy, 0, i * NVPO * sizeof(var_t));

            redutil2::copy_vector_to_device(d_y, h_y, n_var * sizeof(var_t));
            redutil2::copy_vector_to_device(d_p, h_p, n_par * sizeof(var_t));
            redutil2::copy_vector_to_device(d_md, h_md, i * sizeof(nbp_t::metadata_t));

            benchmark_GPU(opt.id_dev, i, d_md, d_y, d_p, d_dy, o_result);

            deallocate_host_storage(&h_y, &h_dy, &h_p, &h_md);
            deallocate_device_storage(&d_y, &d_dy, &d_p, &d_md);
        }
    }

    cout << "Done" << endl;
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
    //cout << "The seed number is " << seed << endl;
    // The pseudo-random number generator is initialized using the argument passed as seed.
    // Used by the subsequent rand() function calls
    if (PROC_UNIT_CPU == opt.comp_dev.proc_unit)
    {
        printf("Compare gravity accelerations computed by CPU (n_obj = %6d):\n\n", n_obj);

        allocate_host_storage(n_obj, &h_y, &h_dy1, &h_p, &h_md);
        ALLOCATE_HOST_VECTOR((void**)(&h_dy2), n_var * sizeof(var_t));
        populate(seed, n_obj, h_y, h_p, h_md);

        // Create aliases
        const var3_t* r = (var3_t*)h_y;
        const nbp_t::param_t* p = (nbp_t::param_t*)h_p;
        var3_t* a1 = (var3_t*)(h_dy1 + nv);
        var3_t* a2 = (var3_t*)(h_dy2 + nv);
        memset(a1, 0, n_obj * sizeof(var3_t));
        memset(a2, 0, n_obj * sizeof(var3_t));

        printf("Comparing the base methods with n_obj parameter:\na) without the use of symmetry and\nb) with symmetry.\n");
#ifdef _WIN32
        chrono::time_point<chrono::system_clock> t0 = chrono::system_clock::now();
#else
        uint64_t t0 = redutil2::GetTimeMs64();
#endif
        cpu_calc_grav_accel(n_obj, h_md, r, p, a1, false);
#ifdef _WIN32
        chrono::time_point<chrono::system_clock> t1 = chrono::system_clock::now();
        chrono::duration<var_t> total_time = t1 - t0;
        var_t Dt_CPU = total_time.count() * 1.0e3;  // [ms]
#else
        uint64_t t1 = redutil2::GetTimeMs64();
        var_t Dt_CPU = (var_t)(t1 - t0) / 1.0e3;  // [ms]
#endif
        printf("The CPU computation without symmetry took: %10.6f [ms]\n", Dt_CPU);

#ifdef _WIN32
        t0 = chrono::system_clock::now();
#else
        t0 = redutil2::GetTimeMs64();
#endif
        cpu_calc_grav_accel(n_obj, h_md, r, p, a2, true);
#ifdef _WIN32
        t1 = chrono::system_clock::now();
        total_time = t1 - t0;
        Dt_CPU = total_time.count() * 1.0e3;  // [ms]
#else
        t1 = redutil2::GetTimeMs64();
        Dt_CPU = (var_t)(t1 - t0) / 1.0e3;  // [ms]
#endif
        printf("The CPU computation with symmetry took: %10.6f [ms]\n", Dt_CPU);

        bool success = compare(n_obj, opt.tol, a1, a2);
        if (success)
        {
            printf(" -- Success.\n\n");
        }
        memset(a1, 0, n_obj * sizeof(var3_t));
        memset(a2, 0, n_obj * sizeof(var3_t));

        printf("Comparing the base methods with snk and src parameters:\na) without the use of symmetry and\nb) with symmetry.\n");
        uint2_t snk = { 0, n_obj };
        uint2_t src = { 0, n_obj };
#ifdef _WIN32
        t0 = chrono::system_clock::now();
#else
        t0 = redutil2::GetTimeMs64();
#endif
        cpu_calc_grav_accel(snk, src, h_md, r, p, a1, false);
#ifdef _WIN32
        t1 = chrono::system_clock::now();
        total_time = t1 - t0;
        Dt_CPU = total_time.count() * 1.0e3;  // [ms]
#else
        t1 = redutil2::GetTimeMs64();
        Dt_CPU = (var_t)(t1 - t0) / 1.0e3;  // [ms]
#endif
        printf("The CPU computation without symmetry took: %10.6f [ms]\n", Dt_CPU);

#ifdef _WIN32
        t0 = chrono::system_clock::now();
#else
        t0 = redutil2::GetTimeMs64();
#endif
        cpu_calc_grav_accel(snk, src, h_md, r, p, a2, true);
#ifdef _WIN32
        t1 = chrono::system_clock::now();
        total_time = t1 - t0;
        Dt_CPU = total_time.count() * 1.0e3;  // [ms]
#else
        t1 = redutil2::GetTimeMs64();
        Dt_CPU = (var_t)(t1 - t0) / 1.0e3;  // [ms]
#endif
        printf("The CPU computation with symmetry took: %10.6f [ms]\n", Dt_CPU);

        success = compare(n_obj, opt.tol, a1, a2);
        if (success)
        {
            printf(" -- Success.\n");
        }

        deallocate_host_storage(&h_y, &h_dy1, &h_p, &h_md);
        FREE_HOST_VECTOR((void **)(h_dy2));
    }
    else
    {
        printf("Compare gravity accelerations computed by GPU (n_obj = %6d):\n\n", n_obj);

        cudaEvent_t start, stop;
        CUDA_SAFE_CALL(cudaEventCreate(&start));
        CUDA_SAFE_CALL(cudaEventCreate(&stop));

        var_t* d_y = 0x0;
        var_t* d_dy = 0x0;
        var_t* d_p = 0x0;
        nbp_t::metadata_t* d_md = 0x0;

        allocate_host_storage(n_obj, &h_y, &h_dy1, &h_p, &h_md);
        ALLOCATE_HOST_VECTOR((void**)(&h_dy2), n_var * sizeof(var_t));
        //populate(seed, n_obj, h_y, h_p, h_md);
        populate(n_obj, h_y, h_p, h_md);
        allocate_device_storage(n_obj, &d_y, &d_dy, &d_p, &d_md);

        redutil2::copy_vector_to_device(d_y, h_y, n_var * sizeof(var_t));
        redutil2::copy_vector_to_device(d_p, h_p, n_par * sizeof(var_t));
        redutil2::copy_vector_to_device(d_md, h_md, n_obj * sizeof(nbp_t::metadata_t));

        {
            // Create aliases
            const var3_t* r = (var3_t*)h_y;
            const nbp_t::param_t* p = (nbp_t::param_t*)h_p;
            var3_t* a1 = (var3_t*)(h_dy1 + nv);
            memset(a1, 0, n_obj * sizeof(var3_t));

            cpu_calc_grav_accel(n_obj, h_md, r, p, a1, true);
        }
        {
            // Create aliases
            const var3_t* r = (var3_t*)d_y;
            const nbp_t::param_t* p = (nbp_t::param_t*)d_p;
            var3_t* a = (var3_t*)(d_dy + nv);

            //float elapsed_time = gpu_calc_grav_accel_naive(n_obj, 256, start, stop, r, p, a);

            uint2_t snk = { 0, n_obj };
            uint2_t src = { 0, n_obj };
            float elapsed_time = gpu_calc_grav_accel_naive(snk, src, 256, start, stop, d_md, r, p, a);

            //float elapsed_time = gpu_calc_grav_accel_tile(n_obj, 256, start, stop, d_md, r, p, a);

            //uint2_t snk = { 0, n_obj };
            //uint2_t src = { 0, n_obj };
            //float elapsed_time = gpu_calc_grav_accel_tile(snk, src, 256, start, stop, d_md, r, p, a);

            redutil2::copy_vector_to_host(h_dy2, d_dy, n_var * sizeof(var_t));
            printf("The GPU computation took: %10.6f [ms]\n", elapsed_time);
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

void compare_part(option& opt, uint2_t snk, uint2_t src)
{
    const uint32_t n_obj = opt.n0;

    const uint16_t n_ppo = sizeof(nbp_t::param_t) / sizeof(var_t);
    const uint32_t n_var = n_obj * NVPO;
    const uint32_t n_par = n_obj * n_ppo;
    // Number of space and velocity coordinates
    const uint32_t nv = NDIM * n_obj;

    printf("Compare gravity accelerations computed by GPU on a block (snk = {%3d, %3d}, src = {%3d, %3d}):\n\n", snk.n1, snk.n2, src.n1, src.n2);

    var_t* h_y = 0x0;
    var_t* h_dy1 = 0x0;
    var_t* h_dy2 = 0x0;
    var_t* h_p = 0x0;
    nbp_t::metadata_t* h_md = 0x0;

    var_t* d_y = 0x0;
    var_t* d_dy = 0x0;
    var_t* d_p = 0x0;
    nbp_t::metadata_t* d_md = 0x0;

    allocate_host_storage(n_obj, &h_y, &h_dy1, &h_p, &h_md);
    ALLOCATE_HOST_VECTOR((void**)(&h_dy2), n_var * sizeof(var_t));
    populate(n_obj, h_y, h_p, h_md);
    allocate_device_storage(n_obj, &d_y, &d_dy, &d_p, &d_md);

    redutil2::copy_vector_to_device(d_y, h_y, n_var * sizeof(var_t));
    redutil2::copy_vector_to_device(d_p, h_p, n_par * sizeof(var_t));
    redutil2::copy_vector_to_device(d_md, h_md, n_obj * sizeof(nbp_t::metadata_t));

    {
        // Create aliases
        const var3_t* r = (var3_t*)h_y;
        const nbp_t::param_t* p = (nbp_t::param_t*)h_p;
        var3_t* a1 = (var3_t*)(h_dy1 + nv);
        memset(a1, 0, n_obj * sizeof(var3_t));
#ifdef _WIN32
        chrono::time_point<chrono::system_clock> t0 = chrono::system_clock::now();
#else
        uint64_t t0 = redutil2::GetTimeMs64();
#endif
        cpu_calc_grav_accel(snk, src, h_md, r, p, a1, false);
#ifdef _WIN32
        chrono::time_point<chrono::system_clock> t1 = chrono::system_clock::now();
        chrono::duration<var_t> total_time = t1 - t0;
        var_t Dt_CPU = total_time.count() * 1.0e3;  // [ms]
#else
        uint64_t t1 = redutil2::GetTimeMs64();
        var_t Dt_CPU = (var_t)(t1 - t0) / 1.0e3;  // [ms]
#endif
        printf("The CPU computation took: %10.6f [ms]\n", Dt_CPU);
    }
    {
        cudaEvent_t start, stop;
        CUDA_SAFE_CALL(cudaEventCreate(&start));
        CUDA_SAFE_CALL(cudaEventCreate(&stop));

        // Create aliases
        const var3_t* r = (var3_t*)d_y;
        const nbp_t::param_t* p = (nbp_t::param_t*)d_p;
        var3_t* a = (var3_t*)(d_dy + nv);
        cudaMemset(a, 0, n_obj * sizeof(var3_t));

        float elapsed_time = gpu_calc_grav_accel_naive(snk, src, 256, start, stop, d_md, r, p, a);

        //float elapsed_time = gpu_calc_grav_accel_tile(n_obj, 256, start, stop, r, p, a);

        //uint2_t snk = { 0, n_obj };
        //uint2_t src = { 0, n_obj };
        //float elapsed_time = gpu_calc_grav_accel_tile(snk, src, 256, start, stop, r, p, a);

        redutil2::copy_vector_to_host(h_dy2, d_dy, n_var * sizeof(var_t));
        printf("The GPU computation took: %10.6f [ms]\n", elapsed_time);
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

    cout << "Done" << endl;
}

void compare_test_case(option& opt, uint32_t n_si, uint32_t n_nsi, uint32_t n_ni)
{
    const uint32_t n_obj = n_si + n_nsi + n_ni;

    const uint16_t n_ppo = sizeof(nbp_t::param_t) / sizeof(var_t);
    const uint32_t n_var = n_obj * NVPO;
    const uint32_t n_par = n_obj * n_ppo;
    // Number of space and velocity coordinates
    const uint32_t nv = NDIM * n_obj;

    printf("Compare gravity accelerations computed by GPU using different body typpes (n_si = %4d, n_nsi = %4d, n_ni = %4d):\n\n", n_si, n_nsi, n_ni);

    var_t* h_y = 0x0;
    var_t* h_dy1 = 0x0;
    var_t* h_dy2 = 0x0;
    var_t* h_p = 0x0;
    nbp_t::metadata_t* h_md = 0x0;

    var_t* d_y = 0x0;
    var_t* d_dy = 0x0;
    var_t* d_p = 0x0;
    nbp_t::metadata_t* d_md = 0x0;

    uint32_t seed = (uint32_t)time(NULL);
    //cout << "The seed number is " << seed << endl;
    // The pseudo-random number generator is initialized using the argument passed as seed.
    // Used by the subsequent rand() function calls

    allocate_host_storage(n_obj, &h_y, &h_dy1, &h_p, &h_md);
    ALLOCATE_HOST_VECTOR((void**)(&h_dy2), n_var * sizeof(var_t));
    populate(seed, n_si, n_nsi, n_ni, h_y, h_p, h_md);
    allocate_device_storage(n_obj, &d_y, &d_dy, &d_p, &d_md);

    redutil2::copy_vector_to_device(d_y, h_y, n_var * sizeof(var_t));
    redutil2::copy_vector_to_device(d_p, h_p, n_par * sizeof(var_t));
    redutil2::copy_vector_to_device(d_md, h_md, n_obj * sizeof(nbp_t::metadata_t));

    uint2_t snk, src;
    // Compute the acceleration on the CPU with the naive method
    {
        // Create aliases
        const var3_t* r = (var3_t*)h_y;
        const nbp_t::param_t* p = (nbp_t::param_t*)h_p;
        var3_t* a = (var3_t*)(h_dy1 + nv);
        memset(a, 0, n_obj * sizeof(var3_t));

#ifdef _WIN32
        chrono::time_point<chrono::system_clock> t0 = chrono::system_clock::now();
#else
        uint64_t t0 = redutil2::GetTimeMs64();
#endif

        // 1. Compute the acceleration for the SI bodies:
        snk.n1 = 0, snk.n2 = n_si;
        src.n1 = 0, src.n2 = n_si + n_nsi;
        cpu_calc_grav_accel(snk, src, h_md, r, p, a, false);

        // 2. Compute the acceleration for the NSI bodies:
        snk.n1 = n_si, snk.n2 = n_si + n_nsi;
        src.n1 = 0, src.n2 = n_si;
        cpu_calc_grav_accel(snk, src, h_md, r, p, a, false);

        // 3. Compute the acceleration for the NI bodies:
        snk.n1 = n_si + n_nsi, snk.n2 = n_obj;
        src.n1 = 0, src.n2 = n_si + n_nsi;
        cpu_calc_grav_accel(snk, src, h_md, r, p, a, false);

#ifdef _WIN32
        chrono::time_point<chrono::system_clock> t1 = chrono::system_clock::now();
        chrono::duration<var_t> total_time = t1 - t0;
        var_t Dt_CPU = total_time.count() * 1.0e3;  // [ms]
#else
        uint64_t t1 = redutil2::GetTimeMs64();
        var_t Dt_CPU = ((var_t)(t1 - t0)) / 1.0e3;  // [ms]
#endif
        printf("The CPU computation took: %10.6f [ms]\n", Dt_CPU);
    }
    {
        cudaEvent_t start, stop;
        CUDA_SAFE_CALL(cudaEventCreate(&start));
        CUDA_SAFE_CALL(cudaEventCreate(&stop));
        float elapsed_time = 0.0f;

        // Create aliases
        const var3_t* r = (var3_t*)d_y;
        const nbp_t::param_t* p = (nbp_t::param_t*)d_p;
        var3_t* a = (var3_t*)(d_dy + nv);
        cudaMemset(a, 0, n_obj * sizeof(var3_t));

        // 1. Compute the acceleration for the SI bodies:
        snk.n1 = 0, snk.n2 = n_si;
        src.n1 = 0, src.n2 = n_si + n_nsi;
        if (0 != (snk.n2 - snk.n1) * (src.n2 - src.n1))
        {
            //elapsed_time = gpu_calc_grav_accel_naive(snk, src, 256, start, stop, r, p, a);
            elapsed_time = gpu_calc_grav_accel_tile(snk, src, 256, start, stop, d_md, r, p, a);
        }

        // 2. Compute the acceleration for the NSI bodies:
        snk.n1 = n_si, snk.n2 = n_si + n_nsi;
        src.n1 = 0, src.n2 = n_si;
        if (0 != (snk.n2 - snk.n1) * (src.n2 - src.n1))
        {
            //elapsed_time += gpu_calc_grav_accel_naive(snk, src, 256, start, stop, r, p, a);
            elapsed_time += gpu_calc_grav_accel_tile(snk, src, 256, start, stop, d_md, r, p, a);
        }

        // 3. Compute the acceleration for the NI bodies:
        snk.n1 = n_si + n_nsi, snk.n2 = n_obj;
        src.n1 = 0, src.n2 = n_si + n_nsi;
        if (0 != (snk.n2 - snk.n1) * (src.n2 - src.n1))
        {
            //elapsed_time += gpu_calc_grav_accel_naive(snk, src, 256, start, stop, r, p, a);
            elapsed_time += gpu_calc_grav_accel_tile(snk, src, 256, start, stop, d_md, r, p, a);
        }

        redutil2::copy_vector_to_host(h_dy2, d_dy, n_var * sizeof(var_t));
        printf("The GPU computation took: %10.6f [ms]\n", elapsed_time);
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

    cout << "Done" << endl;
}


/*
-gpu -id_dev -0 -n0 10 -n1 10000 -dn 10 -v -odir C:\Work\red.cuda.Results\v2.0\Benchmark\Test_01 -bFile inline_0618_1
-gpu -id_dev 0 -n0 5000 -tol 1.0e-16 -v -odir C:\Work\red.cuda.Results\v2.0\Benchmark\Test_01 -bFile compare
-n0 5000 -tol 1.0e-16 -v -odir C:\Work\red.cuda.Results\v2.0\Benchmark\Test_01 -bFile compare
-n_si 100 -n_nsi 1000 -n_ni 1000 -n0 100 -tol 1.0e-16 -v -odir C:\Work\red.cuda.Results\v2.0\Benchmark\Test_01 -bFile compare
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
    uint64_t start = GetTimeMs64();
#endif
    try
    {
        create_default_option(opt);
        parse_options(argc, argv, opt, verbose);
        if (0 > opt.id_dev || opt.id_dev >= get_n_cuda_device())
        {
            throw string("The requested GPU does not exist.");
        }
        create_filename(opt, result_filename);

        switch (opt.job_name)
        {
        case JOB_NAME_BENCMARK_CPU:
        case JOB_NAME_BENCMARK_GPU:
           
            device_query(cout, opt.id_dev, false);
            set_device(opt.id_dev, cout);
            open_stream(opt.o_dir, result_filename, output, BENCHMARK_OUTPUT_NAME_RESULT);
            benchmark(opt, *output[BENCHMARK_OUTPUT_NAME_RESULT]);
            break;
        case JOB_NAME_COMPARE:
            compare(opt);
            {
                opt.n0 = 1000;
                uint2_t snk = { 300, 600 };
                uint2_t src = { 500, 800 };
                compare_part(opt, snk, src);
            }
            compare_test_case(opt, opt.n_si, opt.n_nsi, opt.n_ni);
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
    uint64_t end = GetTimeMs64();
    var_t Dt_CPU = ((var_t)(end - start)) / 1.0e6;
    cout << "Total time: " << Dt_CPU << " s." << endl;
#endif

    return (EXIT_SUCCESS);
}

#undef NDIM
#undef NVPO
