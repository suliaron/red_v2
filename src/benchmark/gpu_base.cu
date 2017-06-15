#include <cfloat>
#include <string>

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include "macro.h"
#include "type.h"
#include "util.h"
#include "redutil2.h"


#define NDIM   3        // Number of space dimension
#define NVPO   6        // Number of variables per object (3 space and 3 velocity coordinates)

using namespace std;

namespace kernel
{
    inline __host__ __device__
        void body_body_grav_accel(const var3_t& ri, const var3_t& rj, var_t mj, var3_t& ai)
    {
        // compute r_ij = r_j - r_i [3 FLOPS] [6 read, 3 write]
        var3_t r_ij = { rj.x - ri.x, rj.y - ri.y, rj.z - ri.z };

        // compute square of r_ij vector [5 FLOPS] [3 read, 1 write]
        var_t d2 = SQR(r_ij.x) + SQR(r_ij.y) + SQR(r_ij.z);
        //var_t d = sqrt(d2);
        //var_t s = K2 * mj / (d2 * d);

        var_t s = K2 * mj / (d2 * sqrt(d2));

        ai.x += s * r_ij.x;
        ai.y += s * r_ij.y;
        ai.z += s * r_ij.z;
    } /* body_body_grav_accel() */

    __global__
        void calc_grav_accel_naive(uint32_t n_obj, const var3_t* r, const nbp_t::param_t* p, var3_t* a)
    {
        // i is the index of the SINK body
        const uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;

        if (i < n_obj)
        {
            a[i].x = a[i].y = a[i].z = 0.0;

            //var3_t r_ij = { 0.0, 0.0, 0.0 };
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
        void calc_grav_accel_naive(uint2_t snk, uint2_t src, const var3_t* r, const nbp_t::param_t* p, var3_t* a)
    {
        // i is the index of the SINK body
        const uint32_t i = snk.n1 + blockIdx.x * blockDim.x + threadIdx.x;
       
        if (i < snk.n2)
        {
            a[i].x = a[i].y = a[i].z = 0.0;

            //var3_t r_ij = { 0.0, 0.0, 0.0 };
            // j is the index of the SOURCE body
            for (uint32_t j = src.n1; j < src.n2; j++)
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



float gpu_calc_grav_accel_naive(uint32_t n_obj, unsigned int n_tpb, cudaEvent_t& start, cudaEvent_t& stop, const var3_t* r, const nbp_t::param_t* p, var3_t* a)
{
    float elapsed_time = 0.0f;

    dim3 grid((n_obj + n_tpb - 1) / n_tpb);
    dim3 block(n_tpb);

    CUDA_SAFE_CALL(cudaEventRecord(start, 0));

    kernel::calc_grav_accel_naive <<< grid, block >>>(n_obj, r, p, a);
    CUDA_CHECK_ERROR();

    CUDA_SAFE_CALL(cudaEventRecord(stop, 0));
    CUDA_SAFE_CALL(cudaEventSynchronize(stop));

    // Computes the elapsed time between two events in milliseconds with a resolution of around 0.5 microseconds.
    CUDA_SAFE_CALL(cudaEventElapsedTime(&elapsed_time, start, stop));

    return elapsed_time;
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
        for (unsigned int n_tpb = 16; n_tpb <= max_n_tpb/2; n_tpb += 16)
        {            
            float elapsed_time = gpu_calc_grav_accel_naive(n_obj, n_tpb, start, stop, r, p, a);
            printf("    %4d %12.4e [sec]\n", n_tpb, elapsed_time / 1.0e3);
            if (elapsed_time < result.y)
            {
                result.x = (float)n_tpb;
                result.y = elapsed_time;
            }
        }
        opt_n_tpb = (unsigned int)result.x;
        n_last = n_obj;
    }
    else
    {
        float elapsed_time = gpu_calc_grav_accel_naive(n_obj, opt_n_tpb, start, stop, r, p, a);
        result.x = (float)opt_n_tpb;
        result.y = elapsed_time;
    }
    first_call = false;

    CUDA_SAFE_CALL(cudaEventDestroy(stop));
    CUDA_SAFE_CALL(cudaEventDestroy(start));

    return result;
}

void benchmark_GPU(int id_dev, uint32_t n_obj, const var_t* d_y, const var_t* d_p, var_t* d_dy, ofstream& o_result)
{
    static string method_name[] = { "base", "base with sym.", "tile", "tile advanced" };
    static string param_name[] = { "n_body", "snk src" };

    uint2_t snk = { 0, 0 };
    uint2_t src = { 0, 0 };
    var_t Dt_CPU = 0.0;

    cudaDeviceProp deviceProp;
    CUDA_SAFE_CALL(cudaGetDeviceProperties(&deviceProp, id_dev));

    float2 result = gpu_calc_grav_accel_naive(n_obj, deviceProp.maxThreadsPerBlock, d_y, d_p, d_dy);
    unsigned int n_tpb = (unsigned int)result.x;
    var_t Dt_GPU = result.y / 1.0e3; // [sec]

    print(PROC_UNIT_GPU, method_name[0], param_name[0], snk, src, n_obj, n_tpb, Dt_CPU, Dt_GPU, o_result, true);
}

void benchmark_GPU(int id_dev, uint32_t n_obj, uint2_t snk, uint2_t src, const var_t* d_y, const var_t* d_p, var_t* d_dy, std::ofstream& o_result)
{
    static string method_name[] = { "base", "base with sym.", "tile", "tile advanced" };
    static string param_name[] = { "n_body", "snk src" };

    var_t Dt_CPU = 0.0;

    cudaDeviceProp deviceProp;
    CUDA_SAFE_CALL(cudaGetDeviceProperties(&deviceProp, id_dev));

    float2 result = gpu_calc_grav_accel_naive(n_obj, deviceProp.maxThreadsPerBlock, d_y, d_p, d_dy);
    unsigned int n_tpb = (unsigned int)result.x;
    var_t Dt_GPU = result.y / 1.0e3; // [sec]

    print(PROC_UNIT_GPU, method_name[0], param_name[0], snk, src, n_obj, n_tpb, Dt_CPU, Dt_GPU, o_result, true);
}

#undef NDIM
#undef NVPO
