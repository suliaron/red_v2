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
    static const string header_str = "date       time     dev  method_name             param_name        n_snk  n_src  n_body n_tpb Dt_CPU[s]     Dt_GPU[ms]";

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
            benchmark_CPU(i, h_y, h_p, h_dy, o_result);
            deallocate_host_storage(&h_y, &h_dy, &h_p, &h_md);
        }
        for (uint32_t i = opt.n0; i <= opt.n1; i *= opt.dn)
        {
            allocate_host_storage(i, &h_y, &h_dy, &h_p, &h_md);
            populate(seed, i, h_y, h_p, h_md);
            uint2_t snk = { 0, i };
            uint2_t src = { 0, i };
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

        const uint16_t n_ppo = sizeof(nbp_t::param_t) / sizeof(var_t);

        for (uint32_t i = opt.n0; i <= opt.n1; i *= opt.dn)
        {

            const uint32_t n_var = i * NVPO;
            const uint32_t n_par = i * n_ppo;

            allocate_host_storage(i, &h_y, &h_dy, &h_p, &h_md);
            populate(seed, i, h_y, h_p, h_md);
            allocate_device_storage(i, &d_y, &d_dy, &d_p, &d_md);

            redutil2::copy_vector_to_device(d_y, h_y, n_var * sizeof(var_t));
            redutil2::copy_vector_to_device(d_p, h_p, n_par * sizeof(var_t));
            redutil2::copy_vector_to_device(d_md, h_md, i * sizeof(nbp_t::metadata_t));

            benchmark_GPU(opt.id_dev, i, d_y, d_p, d_dy, o_result);

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
        memset(a1, 0, n_obj * sizeof(var3_t));
        memset(a2, 0, n_obj * sizeof(var3_t));

        //cpu_calc_grav_accel(n_obj, h_y, h_p, h_dy1, false);
        //cpu_calc_grav_accel(n_obj, h_y, h_p, h_dy2, true);
        cpu_calc_grav_accel(0.0, snk, src, r, p, a1, false);
        cpu_calc_grav_accel(0.0, snk, src, r, p, a2, true);

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

        cudaEvent_t start, stop;
        CUDA_SAFE_CALL(cudaEventCreate(&start));
        CUDA_SAFE_CALL(cudaEventCreate(&stop));

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
            memset(a1, 0, n_obj * sizeof(var3_t));

            cpu_calc_grav_accel(0.0, n_obj, r, p, a1, true);
            //cpu_calc_grav_accel(snk, src, h_y, h_p, h_dy1, false);
            //cpu_calc_grav_accel(snk, src, h_y, h_p, h_dy2, true);
        }
        {
            // Create aliases
            const var3_t* r = (var3_t*)d_y;
            const nbp_t::param_t* p = (nbp_t::param_t*)d_p;
            var3_t* a = (var3_t*)(d_dy + nv);

            //float elapsed_time = gpu_calc_grav_accel_naive(n_obj, 256, start, stop, r, p, a);
            float elapsed_time = gpu_calc_grav_accel_tile(n_obj, 256, start, stop, r, p, a);

            //dim3 grid((n_obj + 256 - 1) / 256);
            //dim3 block(256);
            //kernel::calc_grav_accel_naive<<< grid, block >>>(n_obj, r, p, a);
            //CUDA_CHECK_ERROR();
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


/*
-n0 10 -n1 100000 -dn 10 -v -odir C:\Work\red.cuda.Results\v2.0\Benchmark\Test_01 -bFile benchmark
-gpu -id_dev -0 -n0 10 -n1 10000 -dn 10 -v -odir C:\Work\red.cuda.Results\v2.0\Benchmark\Test_01 -bFile inline_0615_1
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
