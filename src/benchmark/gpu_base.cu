#include <cfloat>
#include <string>

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include "macro.h"
#include "type.h"
#include "util.h"
#include "redutil2.h"

#define VECTOR_FMT (%12.4le, %12.4le, %12.4le)

#define NDIM   3        // Number of space dimension
#define NVPO   6        // Number of variables per object (3 space and 3 velocity coordinates)

using namespace std;

namespace kernel
{
    // 36 FLOP
    inline __host__ __device__
        void body_body_grav_accel(const var3_t& ri, const var3_t& rj, var_t mj, var3_t& ai)
    {
        // compute r_ij = r_j - r_i [3 FLOPS] [6 read, 3 write]
        var3_t r_ij = { rj.x - ri.x, rj.y - ri.y, rj.z - ri.z };

        // compute square of r_ij vector [5 FLOPS + ] [3 read, 1 write]
        var_t d2 = SQR(r_ij.x) + SQR(r_ij.y) + SQR(r_ij.z);
        var_t s = K2 * mj / (d2 * sqrt(d2));

        // 6 FLOP
        ai.x += s * r_ij.x;
        ai.y += s * r_ij.y;
        ai.z += s * r_ij.z;
    } /* body_body_grav_accel() */

    __global__
        void calc_grav_accel_naive(uint32_t n_obj, const nbp_t::metadata_t* md, const var3_t* r, const nbp_t::param_t* p, var3_t* a)
    {
        // i is the index of the SINK body
        const uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;

        if (i < n_obj)
        {
            // j is the index of the SOURCE body
            for (uint32_t j = 0; j < n_obj; j++)
            {
                if (i == j) continue;
                kernel::body_body_grav_accel(r[i], r[j], p[j].mass, a[i]);
            }
        }
    } /* calc_grav_accel_naive () */

    __global__
        void calc_grav_accel_naive(uint2_t snk, uint2_t src, const nbp_t::metadata_t* md, const var3_t* r, const nbp_t::param_t* p, var3_t* a)
    {
        // i is the index of the SINK body
        const uint32_t i = snk.n1 + blockIdx.x * blockDim.x + threadIdx.x;
       
        if (snk.n2 > i)
        {
            // j is the index of the SOURCE body
            for (uint32_t j = src.n1; j < src.n2; j++)
            {
                if (i == j) continue;
                kernel::body_body_grav_accel(r[i], r[j], p[j].mass, a[i]);
            }
        }
    } /* calc_grav_accel_naive () */

    __global__
        void calc_grav_accel_tile(uint32_t n_obj, const nbp_t::metadata_t* md, const var3_t* r, const nbp_t::param_t* p, var3_t* a)
    {
        extern __shared__ var3_t sh_pos[];

        const uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
        //var3_t acc = { 0.0, 0.0, 0.0 };
        var3_t acc = a[i];
        var3_t my_pos;

        // To avoid overruning the r buffer
        if (n_obj > i)
        {
            my_pos = r[i];
        }
        // Note! : the for cycle must be outside the upper if clause, otherwise the sh_pos array will
        // not recive the input for the last tile! The reason is that some thread will be not considered
        // in the if (n_obj > idx) clause.
        for (uint32_t tile = 0; (tile * blockDim.x) < n_obj; tile++)
        {
            const uint32_t idx = tile * blockDim.x + threadIdx.x;
            // To avoid overruning the r and mass buffer
            if (n_obj > idx)
            {
                sh_pos[threadIdx.x] = r[idx];
            }
            __syncthreads();

            for (int j = 0; j < blockDim.x; j++)
            {
                // To avoid overrun the input arrays
                if (n_obj <= (tile * blockDim.x) + j)
                    break;
                // To avoid self-interaction
                if (i == (tile * blockDim.x) + j)
                    continue;
                body_body_grav_accel(my_pos, sh_pos[j], p[(tile * blockDim.x) + j].mass, acc);
            }
            __syncthreads();
        }
        if (n_obj > i)
        {
            a[i] = acc;
        }        
    }

    __global__
        void calc_grav_accel_tile(uint2_t snk, uint2_t src, const nbp_t::metadata_t* md, const var3_t* r, const nbp_t::param_t* p, var3_t* a)
    {
        extern __shared__ var3_t sh_pos[];

        const uint32_t i = snk.n1 + blockIdx.x * blockDim.x + threadIdx.x;
        //var3_t acc = { 0.0, 0.0, 0.0 };
        var3_t acc = a[i];
        var3_t my_pos;

        // To avoid overruning the r buffer
        if (snk.n2 > i)
        {
            my_pos = r[i];
        }
        // Note! : the for cycle must be outside the upper if clause, otherwise the sh_pos array will
        // not recive the input for the last tile! The reason is that some thread will be not considered
        // in the if (n_obj > idx) clause.
        for (uint32_t tile = 0; (tile * blockDim.x) < src.n2; tile++)
        {
            const uint32_t idx = src.n1 + tile * blockDim.x + threadIdx.x;
            // To avoid overruning the r and mass buffer
            if (src.n2 > idx)
            {
                sh_pos[threadIdx.x] = r[idx];
            }
            __syncthreads();

            for (int j = 0; j < blockDim.x; j++)
            {
                // To avoid overrun then input arrays
                if (src.n2 <= src.n1 + (tile * blockDim.x) + j)
                    break;
                // To avoid self-interaction
                if (i == src.n1 + (tile * blockDim.x) + j)
                    continue;
                body_body_grav_accel(my_pos, sh_pos[j], p[(tile * blockDim.x) + j].mass, acc);
            }
            __syncthreads();
        }
        if (snk.n2 > i)
        {
            a[i] = acc;
        }
    }
} /* namespace kernel */



float gpu_calc_grav_accel_naive(uint32_t n_obj, unsigned int n_tpb, cudaEvent_t& start, cudaEvent_t& stop, const nbp_t::metadata_t* md, const var3_t* r, const nbp_t::param_t* p, var3_t* a)
{
    float elapsed_time = 0.0f;

    dim3 grid((n_obj + n_tpb - 1) / n_tpb);
    dim3 block(n_tpb);

    CUDA_SAFE_CALL(cudaEventRecord(start, 0));

    kernel::calc_grav_accel_naive <<< grid, block >>>(n_obj, md, r, p, a);
    CUDA_CHECK_ERROR();

    CUDA_SAFE_CALL(cudaEventRecord(stop, 0));
    CUDA_SAFE_CALL(cudaEventSynchronize(stop));

    // Computes the elapsed time between two events in milliseconds with a resolution of around 0.5 microseconds.
    CUDA_SAFE_CALL(cudaEventElapsedTime(&elapsed_time, start, stop));

    return elapsed_time;
}

float gpu_calc_grav_accel_naive(uint2_t snk, uint2_t src, unsigned int n_tpb, cudaEvent_t& start, cudaEvent_t& stop, const nbp_t::metadata_t* md, const var3_t* r, const nbp_t::param_t* p, var3_t* a)
{
    float elapsed_time = 0.0f;

    dim3 grid(((snk.n2 - snk.n1) + n_tpb - 1) / n_tpb);
    dim3 block(n_tpb);

    CUDA_SAFE_CALL(cudaEventRecord(start, 0));

    kernel::calc_grav_accel_naive <<< grid, block >>>(snk, src, md, r, p, a);
    CUDA_CHECK_ERROR();

    CUDA_SAFE_CALL(cudaEventRecord(stop, 0));
    CUDA_SAFE_CALL(cudaEventSynchronize(stop));

    // Computes the elapsed time between two events in milliseconds with a resolution of around 0.5 microseconds.
    CUDA_SAFE_CALL(cudaEventElapsedTime(&elapsed_time, start, stop));

    return elapsed_time;
}

float gpu_calc_grav_accel_tile(uint32_t n_obj, unsigned int n_tpb, cudaEvent_t& start, cudaEvent_t& stop, const nbp_t::metadata_t* md, const var3_t* r, const nbp_t::param_t* p, var3_t* a)
{
    float elapsed_time = 0.0f;

    dim3 grid((n_obj + n_tpb - 1) / n_tpb);
    dim3 block(n_tpb);

    CUDA_SAFE_CALL(cudaEventRecord(start, 0));

    size_t sh_mem_size = n_tpb * sizeof(var3_t);
    kernel::calc_grav_accel_tile <<< grid, block, sh_mem_size >>>(n_obj, md, r, p, a);
    CUDA_CHECK_ERROR();

    CUDA_SAFE_CALL(cudaEventRecord(stop, 0));
    CUDA_SAFE_CALL(cudaEventSynchronize(stop));

    // Computes the elapsed time between two events in milliseconds with a resolution of around 0.5 microseconds.
    CUDA_SAFE_CALL(cudaEventElapsedTime(&elapsed_time, start, stop));

    return elapsed_time;
}

float gpu_calc_grav_accel_tile(uint2_t snk, uint2_t src, unsigned int n_tpb, cudaEvent_t& start, cudaEvent_t& stop, const nbp_t::metadata_t* md, const var3_t* r, const nbp_t::param_t* p, var3_t* a)
{
    float elapsed_time = 0.0f;

    dim3 grid(((snk.n2 - snk.n1) + n_tpb - 1) / n_tpb);
    dim3 block(n_tpb);

    CUDA_SAFE_CALL(cudaEventRecord(start, 0));

    size_t sh_mem_size = n_tpb * sizeof(var3_t);
    kernel::calc_grav_accel_tile << < grid, block, sh_mem_size >> >(snk, src, md, r, p, a);
    CUDA_CHECK_ERROR();

    CUDA_SAFE_CALL(cudaEventRecord(stop, 0));
    CUDA_SAFE_CALL(cudaEventSynchronize(stop));

    // Computes the elapsed time between two events in milliseconds with a resolution of around 0.5 microseconds.
    CUDA_SAFE_CALL(cudaEventElapsedTime(&elapsed_time, start, stop));

    return elapsed_time;
}


float2 gpu_calc_grav_accel_naive(uint32_t n_obj, unsigned int max_n_tpb, const nbp_t::metadata_t* d_md, const var_t* d_y, const var_t* d_p, var_t* d_dy)
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
            float elapsed_time = gpu_calc_grav_accel_naive(n_obj, n_tpb, start, stop, d_md, r, p, a);
            printf("    %4d %12.4e [sec]\n", n_tpb, elapsed_time / 1.0e3);
            if (elapsed_time < result.y)
            {
                result.x = (float)n_tpb;
                result.y = elapsed_time;  // [ms]
            }
        }
        opt_n_tpb = (unsigned int)result.x;
        n_last = n_obj;
    }
    else
    {
        float elapsed_time = gpu_calc_grav_accel_naive(n_obj, opt_n_tpb, start, stop, d_md, r, p, a);
        result.x = (float)opt_n_tpb;
        result.y = elapsed_time;  // [ms]
    }
    first_call = false;

    CUDA_SAFE_CALL(cudaEventDestroy(stop));
    CUDA_SAFE_CALL(cudaEventDestroy(start));

    return result;
}

float2 gpu_calc_grav_accel_naive(uint32_t n_obj, uint2_t snk, uint2_t src, unsigned int max_n_tpb, const nbp_t::metadata_t* d_md, const var_t* d_y, const var_t* d_p, var_t* d_dy)
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
        for (unsigned int n_tpb = 16; n_tpb <= max_n_tpb / 2; n_tpb += 16)
        {
            float elapsed_time = gpu_calc_grav_accel_naive(snk, src, n_tpb, start, stop, d_md, r, p, a);
            printf("    %4d %12.4e [sec]\n", n_tpb, elapsed_time / 1.0e3);
            if (elapsed_time < result.y)
            {
                result.x = (float)n_tpb;
                result.y = elapsed_time;  // [ms]
            }
        }
        opt_n_tpb = (unsigned int)result.x;
        n_last = n_obj;
    }
    else
    {
        float elapsed_time = gpu_calc_grav_accel_naive(snk, src, opt_n_tpb, start, stop, d_md, r, p, a);
        result.x = (float)opt_n_tpb;
        result.y = elapsed_time;  // [ms]
    }
    first_call = false;

    CUDA_SAFE_CALL(cudaEventDestroy(stop));
    CUDA_SAFE_CALL(cudaEventDestroy(start));

    return result;
}

float2 gpu_calc_grav_accel_tile(uint32_t n_obj, unsigned int max_n_tpb, const nbp_t::metadata_t* d_md, const var_t* d_y, const var_t* d_p, var_t* d_dy)
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
        for (unsigned int n_tpb = 16; n_tpb <= max_n_tpb / 2; n_tpb += 16)
        {
            float elapsed_time = gpu_calc_grav_accel_tile(n_obj, n_tpb, start, stop, d_md, r, p, a);
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
        float elapsed_time = gpu_calc_grav_accel_tile(n_obj, opt_n_tpb, start, stop, d_md, r, p, a);
        result.x = (float)opt_n_tpb;
        result.y = elapsed_time;
    }
    first_call = false;

    CUDA_SAFE_CALL(cudaEventDestroy(stop));
    CUDA_SAFE_CALL(cudaEventDestroy(start));

    return result;
}

float2 gpu_calc_grav_accel_tile(uint32_t n_obj, uint2_t snk, uint2_t src, unsigned int max_n_tpb, const nbp_t::metadata_t* d_md, const var_t* d_y, const var_t* d_p, var_t* d_dy)
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
        for (unsigned int n_tpb = 16; n_tpb <= max_n_tpb / 2; n_tpb += 16)
        {
            float elapsed_time = gpu_calc_grav_accel_tile(snk, src, n_tpb, start, stop, d_md, r, p, a);
            printf("    %4d %12.4e [sec]\n", n_tpb, elapsed_time / 1.0e3);
            if (elapsed_time < result.y)
            {
                result.x = (float)n_tpb;
                result.y = elapsed_time;  // [ms]
            }
        }
        opt_n_tpb = (unsigned int)result.x;
        n_last = n_obj;
    }
    else
    {
        float elapsed_time = gpu_calc_grav_accel_tile(snk, src, opt_n_tpb, start, stop, d_md, r, p, a);
        result.x = (float)opt_n_tpb;
        result.y = elapsed_time;  // [ms]
    }
    first_call = false;

    CUDA_SAFE_CALL(cudaEventDestroy(stop));
    CUDA_SAFE_CALL(cudaEventDestroy(start));

    return result;
}

void benchmark_GPU(int id_dev, uint32_t n_obj, const nbp_t::metadata_t* d_md, const var_t* d_y, const var_t* d_p, var_t* d_dy, ofstream& o_result)
{
    static string method_name[] = { "base", "base_with_sym", "tile", "tile_advanced" };
    static string param_name[] = { "n_body", "snk_src" };

    uint2_t snk = { 0, 0 };
    uint2_t src = { 0, 0 };
    var_t Dt_CPU = 0.0;

    cudaDeviceProp deviceProp;
    CUDA_SAFE_CALL(cudaGetDeviceProperties(&deviceProp, id_dev));

    // 1. Naive method on the GPU with n_obj parameter
    float2 result = gpu_calc_grav_accel_naive(n_obj, deviceProp.maxThreadsPerBlock, d_md, d_y, d_p, d_dy);
    unsigned int n_tpb = (unsigned int)result.x;
    var_t Dt_GPU = result.y;   // [ms]

    print(PROC_UNIT_GPU, method_name[0], param_name[0], snk, src, n_obj, n_tpb, Dt_CPU, Dt_GPU, o_result, true);

    // 2. Naive method on the GPU with snk and src parameters
    snk.n2 = n_obj;
    src.n2 = n_obj;
    result = gpu_calc_grav_accel_naive(n_obj, snk, src, deviceProp.maxThreadsPerBlock, d_md, d_y, d_p, d_dy);
    n_tpb = (unsigned int)result.x;
    Dt_GPU = result.y;   // [ms]

    print(PROC_UNIT_GPU, method_name[0], param_name[1], snk, src, n_obj, n_tpb, Dt_CPU, Dt_GPU, o_result, true);

    // 3. Tile method on the GPU with n_obj parameter
    result = gpu_calc_grav_accel_tile(n_obj, deviceProp.maxThreadsPerBlock, d_md, d_y, d_p, d_dy);
    n_tpb = (unsigned int)result.x;
    Dt_GPU = result.y;   // [ms]

    snk.n2 = 0;
    src.n2 = 0;
    print(PROC_UNIT_GPU, method_name[2], param_name[0], snk, src, n_obj, n_tpb, Dt_CPU, Dt_GPU, o_result, true);

    // 4. Tile method on the GPU with snk and src parameters
    snk.n2 = n_obj;
    src.n2 = n_obj;
    result = gpu_calc_grav_accel_tile(n_obj, snk, src, deviceProp.maxThreadsPerBlock, d_md, d_y, d_p, d_dy);
    n_tpb = (unsigned int)result.x;
    Dt_GPU = result.y;   // [ms]

    print(PROC_UNIT_GPU, method_name[2], param_name[1], snk, src, n_obj, n_tpb, Dt_CPU, Dt_GPU, o_result, true);
}

#undef NDIM
#undef NVPO
