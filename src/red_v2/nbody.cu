#include <cfloat>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include "nbody.h"

#include "redutil2.h"
#include "constants.h"

#define NDIM   3        // Number of space dimension
#define NVPO   6        // Number of variables per object (3 space and 3 velocity coordinates)


using namespace std;
using namespace redutil2;

namespace kernel_nbody
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
                body_body_grav_accel(r[i], r[j], p[j].mass, a[i]);
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
                body_body_grav_accel(r[i], r[j], p[j].mass, a[i]);
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
} /* kernel_nbody */

nbody::nbody(string& path_si, string& path_sd, uint32_t n_obj, uint16_t n_ppo, comp_dev_t comp_dev) :
	ode(NDIM, n_obj, NVPO, n_ppo, comp_dev)
{
	name = "Singular 3D n-body problem";
	
	initialize();
	allocate_storage();

    load_solution_info(path_si);
    load_solution_data(path_sd);

    calc_n_types();
	if (PROC_UNIT_GPU == comp_dev.proc_unit)
	{
		copy_vars(COPY_DIRECTION_TO_DEVICE);
		copy_params(COPY_DIRECTION_TO_DEVICE);
	}
	calc_integral();
	tout = t;
}

nbody::~nbody()
{
	deallocate_storage();
}

void nbody::initialize()
{
	h_md = NULL;
	d_md = NULL;
	md   = NULL;

    n_si  = 0;
    n_nsi = 0;
    n_ni  = 0;

    n_tpb_si  = 0;
    n_tpb_nsi = 0;
    n_tpb_ni  = 0;
}

void nbody::allocate_storage()
{
	allocate_host_storage();
	if (PROC_UNIT_GPU == comp_dev.proc_unit)
	{
		allocate_device_storage();
	}
}

void nbody::allocate_host_storage()
{
	ALLOCATE_HOST_VECTOR((void**)&(h_md), n_obj * sizeof(nbp_t::metadata_t));
}

void nbody::allocate_device_storage()
{
	ALLOCATE_DEVICE_VECTOR((void**)&(d_md), n_obj * sizeof(nbp_t::metadata_t));
}

void nbody::deallocate_storage()
{
	//NOTE : First always release the DEVICE memory
	if (PROC_UNIT_GPU == comp_dev.proc_unit)
	{
		deallocate_device_storage();
	}
	deallocate_host_storage();
}

void nbody::deallocate_host_storage()
{
	FREE_HOST_VECTOR((void **)&(h_md));
}

void nbody::deallocate_device_storage()
{
	FREE_DEVICE_VECTOR((void **)&(d_md));
}

void nbody::copy_metadata(copy_direction_t dir)
{
	switch (dir)
	{
	case COPY_DIRECTION_TO_DEVICE:
		copy_vector_to_device(d_md, h_md, n_obj*sizeof(nbp_t::metadata_t));
		break;
	case COPY_DIRECTION_TO_HOST:
		copy_vector_to_host(h_md, d_md, n_obj*sizeof(nbp_t::metadata_t));
		break;
	default:
		throw std::string("Parameter 'dir' is out of range.");
	}
}

void nbody::calc_n_types()
{
    for (uint32_t i = 0; i < n_obj; i++)
    {
        switch (h_md[i].body_type)
        {
        case BODY_TYPE_STAR:
        case BODY_TYPE_GIANTPLANET:
        case BODY_TYPE_ROCKYPLANET:
        case BODY_TYPE_PROTOPLANET:
            n_si++;
            break;
        case BODY_TYPE_SUPERPLANETESIMAL:
        case BODY_TYPE_PLANETESIMAL:
            n_nsi++;
            break;
        case BODY_TYPE_TESTPARTICLE:
            n_ni++;
            break;
        default:
            throw string("Unknown body type.");
        }
    }
}

// For Hermite4
void nbody::calc_dy(uint16_t stage, var_t curr_t, const var_t* y_temp, var_t* acc, var_t* jrk)
{
	if (PROC_UNIT_CPU == comp_dev.proc_unit)
	{
        // TODO: implement the symmetric version
		cpu_calc_dy(stage, curr_t, y_temp, acc, jrk, false);
	}
	else
	{
		throw string("The nbody::gpu_calc_dy for Hermite4 is not implemented.");
	}
}

// For Runge-Kutta type integrators
void nbody::calc_dy(uint16_t stage, var_t curr_t, const var_t* y_temp, var_t* dy)
{
	if (PROC_UNIT_CPU == comp_dev.proc_unit)
	{
		cpu_calc_dy(stage, curr_t, y_temp, dy);
	}
	else
	{
		gpu_calc_dy(stage, curr_t, y_temp, dy);
	}
}

// For Hermite4
void nbody::cpu_calc_dy(uint16_t stage, var_t curr_t, const var_t* y_temp, var_t* acc, var_t* jrk, bool use_symm_prop)
{
    // Number of space and velocity coordinates
    const uint32_t nv = NDIM * n_obj;

    // Create aliases
    const var3_t* r = (var3_t*)y_temp;
    const var3_t* v = (var3_t*)(y_temp + nv);
    const nbp_t::param_t* p = (nbp_t::param_t*)h_p;
    
    var3_t* _acc = (var3_t*)(acc);
	var3_t* _jrk = (var3_t*)(jrk);

    // Clear the acceleration and jerk arrays: the += op can be used
	memset(_acc, 0, n_obj * sizeof(var3_t));
	memset(_jrk, 0, n_obj * sizeof(var3_t));

	if (use_symm_prop)
	{
    	throw string("The symmetric version of nbody::cpu_calc_dy is not implemented.");
	}
	else
	{
		for (uint32_t i = 0; i < n_obj; i++)
		{
			var3_t r_ji = {0, 0, 0};
			var3_t v_ji = {0, 0, 0};
			for (uint32_t j = 0; j < n_obj; j++)
			{
				if (i == j)
				{
					continue;
				}
				r_ji.x = r[j].x - r[i].x;
				r_ji.y = r[j].y - r[i].y;
				r_ji.z = r[j].z - r[i].z;
				var_t d2 = SQR(r_ji.x) + SQR(r_ji.y) + SQR(r_ji.z);
				var_t d = sqrt(d2);
				var_t d_3 = K2 / (d*d2);
				var_t s = p[j].mass * d_3;

                _acc[i].x += s * r_ji.x;
				_acc[i].y += s * r_ji.y;
				_acc[i].z += s * r_ji.z;

                v_ji.x = v[j].x - v[i].x;
				v_ji.y = v[j].y - v[i].y;
				v_ji.z = v[j].z - v[i].z;
                var_t alpha = 3.0 * (r_ji.x * v_ji.x + r_ji.y * v_ji.y + r_ji.z * v_ji.z) / (d2 * d2 * d);

                _jrk[i].x += s * v_ji.x - alpha * _acc[i].x;
                _jrk[i].y += s * v_ji.y - alpha * _acc[i].y;
                _jrk[i].z += s * v_ji.z - alpha * _acc[i].z;
			}
		}
	}
}

// For Runge-Kutta type integrators
void nbody::cpu_calc_dy(uint16_t stage, var_t curr_t, const var_t* y_temp, var_t* dy)
{
    // Number of space and velocity coordinates
    const uint32_t nv = NDIM * n_obj;

    // Create aliases
    const var3_t* r = (var3_t*)y_temp;
    const var3_t* v = (var3_t*)(y_temp + nv);
    const nbp_t::param_t* p = (nbp_t::param_t*)h_p;

    var3_t* a = (var3_t*)(dy + nv);

    // Copy the velocities into dy
    memcpy(dy, v, nv * sizeof(var_t));

    // Clear the acceleration array: the += op can be used
	memset(a, 0, nv *sizeof(var_t));
	for (uint32_t i = 0; i < n_obj; i++)
	{
		//var3_t r_ij = {0, 0, 0};
		for (uint32_t j = i+1; j < n_obj; j++)
		{
            body_body_grav_accel(r[i], r[j], p[i].mass, p[j].mass, a[i], a[j]);
		}
	}
}

inline
void nbody::body_body_grav_accel(const var3_t& ri, const var3_t& rj, var_t mi, var_t mj, var3_t& ai, var3_t& aj)
{
    // compute r_ij = r_j - r_i [3 FLOPS] [6 read, 3 write]
    var3_t r_ij = { rj.x - ri.x, rj.y - ri.y, rj.z - ri.z };

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

void nbody::gpu_calc_dy(uint16_t stage, var_t curr_t, const var_t* y_temp, var_t* dy)
{
    static bool first_call = true;
    static uint32_t last_n_obj = n_obj;

    // Number of space and velocity coordinates
    const uint32_t nv = NDIM * n_obj;

    // Create aliases
    const var3_t* r = (var3_t*)y_temp;
    const var3_t* v = (var3_t*)(y_temp + nv);
    const nbp_t::param_t* p = (nbp_t::param_t*)d_p;

    var3_t* a = (var3_t*)(dy + nv);

    // 1. Calculate the acceleration from the gravitational forces
    if (first_call || last_n_obj != n_obj)
	{
        printf("Searching for the optimal thread number ");
        first_call = false;
        last_n_obj = n_obj;

        float min_dt = 1.0e10;
        for (unsigned int i = 16; i < 512; i += 16)
        {
            putc('.', stdout);
            cudaMemset(a, 0, n_obj * sizeof(var3_t));
            uint2_t snk = { 0, n_si };
            uint2_t src = { 0, n_si + n_nsi };
            //float dt_GPU = gpu_calc_grav_accel_naive(stage, snk, src, i, curr_t, r, p, a);
            float dt_GPU = gpu_calc_grav_accel_tile(stage, snk, src, i, curr_t, r, p, a);
            if (dt_GPU < min_dt)
            {
                min_dt = dt_GPU;
                n_tpb_si = i;
            }
        }
        printf(" Done.\nOptimal thread number = %3d.\n", n_tpb_si);
        // TODO: Move this line next to the kernel invocation since the kernel execution and data copy can be performed simultaneously
        // Copy the velocities into dy
        CUDA_SAFE_CALL(cudaMemcpy(dy, v, nv * sizeof(var_t), cudaMemcpyDeviceToDevice));
    }
    else
    {
        set_kernel_launch_param(n_obj, n_tpb_si, grid, block);

        cudaMemset(a, 0, n_obj * sizeof(var3_t));
        //kernel_nbody::calc_grav_accel_naive <<< grid, block >>>(n_obj, d_md, r, p, a);
        uint2_t snk = { 0, n_si };
        uint2_t src = { 0, n_si + n_nsi };
        float dt_GPU = gpu_calc_grav_accel_tile(stage, snk, src, n_tpb_si, curr_t, r, p, a);
        // TODO: Move this line next to the kernel invocation since the kernel execution and data copy can be performed simultaneously
        // Copy the velocities into dy
        CUDA_SAFE_CALL(cudaMemcpy(dy, v, nv * sizeof(var_t), cudaMemcpyDeviceToDevice));
        CUDA_CHECK_ERROR();
    }
    // ------- END -------

    // 2. Calculate the accelerations from other forces
    // ...
}

float nbody::gpu_calc_grav_accel_naive(uint16_t stage, uint2_t snk, uint2_t src, unsigned int n_tpb, var_t curr_t, const var3_t* r, const nbp_t::param_t* p, var3_t* a)
{
    cudaEvent_t start, stop;

    CUDA_SAFE_CALL(cudaEventCreate(&start));
    CUDA_SAFE_CALL(cudaEventCreate(&stop));

    set_kernel_launch_param(n_obj, n_tpb, grid, block);
    CUDA_SAFE_CALL(cudaEventRecord(start, 0));
    kernel_nbody::calc_grav_accel_naive <<< grid, block >>>(snk, src, md, r, p, a);

    CUDA_CHECK_ERROR();

    CUDA_SAFE_CALL(cudaEventRecord(stop, 0));
    CUDA_SAFE_CALL(cudaEventSynchronize(stop));

    // Computes the elapsed time between two events in milliseconds with a resolution of around 0.5 microseconds.
    float elapsed_time = 0.0f;
    CUDA_SAFE_CALL(cudaEventElapsedTime(&elapsed_time, start, stop));

    return elapsed_time;
}

float nbody::gpu_calc_grav_accel_tile(uint16_t stage, uint2_t snk, uint2_t src, unsigned int n_tpb, var_t curr_t, const var3_t* r, const nbp_t::param_t* p, var3_t* a)
{
    cudaEvent_t start, stop;

    CUDA_SAFE_CALL(cudaEventCreate(&start));
    CUDA_SAFE_CALL(cudaEventCreate(&stop));

    set_kernel_launch_param(n_obj, n_tpb, grid, block);
    CUDA_SAFE_CALL(cudaEventRecord(start, 0));
    size_t sh_mem_size = n_tpb * sizeof(var3_t);
    kernel_nbody::calc_grav_accel_tile <<< grid, block, sh_mem_size >>>(snk, src, md, r, p, a);

    CUDA_CHECK_ERROR();

    CUDA_SAFE_CALL(cudaEventRecord(stop, 0));
    CUDA_SAFE_CALL(cudaEventSynchronize(stop));

    // Computes the elapsed time between two events in milliseconds with a resolution of around 0.5 microseconds.
    float elapsed_time = 0.0f;
    CUDA_SAFE_CALL(cudaEventElapsedTime(&elapsed_time, start, stop));

    return elapsed_time;
}

void nbody::calc_integral()
{
	static bool first_call = true;

    // Number of space and velocity coordinates
    const uint32_t nv = NDIM * n_obj;
    // Create aliases
    const var3_t* r = (var3_t*)h_y;
    const var3_t* v = (var3_t*)(h_y + nv);
    const nbp_t::param_t* p = (nbp_t::param_t*)h_p;

	integral.R = tools::nbp::calc_position_of_bc(n_obj, p, r);
	integral.V = tools::nbp::calc_velocity_of_bc(n_obj, p, v);
	integral.c = tools::nbp::calc_angular_momentum(n_obj, p, r, v);
	integral.h = tools::nbp::calc_total_energy(n_obj, p, r, v);

	if (first_call)
	{
		first_call = false;

		integral.R0 = integral.R;
		integral.V0 = integral.V;
		integral.c0 = integral.c;
		integral.h0 = integral.h;
	}
}

void nbody::chk_coll(var_t a_ref)
{
    // Number of space and velocity coordinates
    const uint32_t nv = NDIM * n_obj;
    // Create aliases
    const var3_t* r = (var3_t*)h_y;
    const var3_t* v = (var3_t*)(h_y + nv);
    const nbp_t::param_t* p = (nbp_t::param_t*)h_p;

    for (uint32_t i = 0; i < n_obj; i++)
    {
        var3_t r_ij = { 0, 0, 0 };
        var_t min_d = 0.1;
        for (uint32_t j = i + 1; j < n_obj; j++)
        {
            r_ij.x = r[j].x - r[i].x;
            r_ij.y = r[j].y - r[i].y;
            r_ij.z = r[j].z - r[i].z;

            var_t d = sqrt(SQR(r_ij.x) + SQR(r_ij.y) + SQR(r_ij.z));
            if ((a_ref * (p[i].radius + p[j].radius)) > d)
            {
                printf("%*.16e Collision: %d %d (%*.16e)\n", VAR_T_W, t, h_md[i].id, h_md[j].id, VAR_T_W, d / (a_ref * (p[i].radius + p[j].radius)));
            }

            if (min_d > d)
            {
                // L suffix denotes the length of the corresponding vector
                // 2 suffix denotes the square of the length of the corresponding vector
                var3_t v_ij = {v[j].x - v[i].x, v[j].y - v[i].y, v[j].z - v[i].z};
                var_t r_ijL = d;
                var_t v_ijL = sqrt(SQR(v_ij.x) + SQR(v_ij.y) + SQR(v_ij.z));

                var_t rv = r_ij.x * v_ij.x + r_ij.y * v_ij.y + r_ij.z * v_ij.z;
                var_t alpha = PI - acos(rv / (r_ijL * v_ijL));
                var_t b = r_ijL * sin(alpha);
                var_t v_esc2 = 2.0 * K2 * p[i].mass / p[i].radius;
                var_t Gamma = 1.0 + v_esc2 / SQR(v_ijL);
                if (SQR(b) <= SQR(p[i].radius) * Gamma)
                {
                    printf("%*.16e Collision expected between %d %d (distance = %*.16e [R1+R2])\n", VAR_T_W, t, h_md[i].id, h_md[j].id, VAR_T_W, r_ijL / (p[i].radius + p[j].radius));

                    printf("ri_Vec.x  = %25.16e & ri_Vec.y  = %25.16e & ri_Vec.z  = %25.16e\n", r[i].x, r[i].y, r[i].z);
                    printf("rj_Vec.x  = %25.16e & rj_Vec.y  = %25.16e & rj_Vec.z  = %25.16e\n", r[j].x, r[j].y, r[j].z);
                    printf("rij_Vec.x = %25.16e & rij_Vec.y = %25.16e & rij_Vec.z = %25.16e\n", r_ij.x, r_ij.y, r_ij.z);
                    printf("vi_Vec.x  = %25.16e & vi_Vec.y  = %25.16e & vi_Vec.z  = %25.16e\n", v[i].x, v[i].y, v[i].z);
                    printf("vj_Vec.x  = %25.16e & vj_Vec.y  = %25.16e & vj_Vec.z  = %25.16e\n", v[j].x, v[j].y, v[j].z);
                    printf("vij_Vec.x = %25.16e & vij_Vec.y = %25.16e & vij_Vec.z = %25.16e\n", v_ij.x, v_ij.y, v_ij.z);
                    printf("m_i = %25.16e & radius_i = %25.16e\n", p[i].mass, p[i].radius);
                    printf("m_j = %25.16e & radius_j = %25.16e\n", p[j].mass, p[j].radius);
                    printf("rv = %25.16e\n", rv);
                    printf("alpha = %25.16e\n", alpha);
                    printf("b = %25.16e\n", b);
                    printf("v_esc2 = %25.16e\n", v_esc2);
                    printf("Gamma = %25.16e\n", Gamma);
                }
            }

        }
    }
}

void nbody::load_solution_info(string& path)
{
	ifstream input;

	cout << "Loading " << path << " ";

	data_rep_t repres = file::get_data_repres(path);
	switch (repres)
	{
	case DATA_REPRESENTATION_ASCII:
		input.open(path.c_str(), ios::in);
		if (input) 
		{
			input >> t >> dt >> t_wc >> n_obj;
		}
		else 
		{
			throw string("Cannot open " + path + ".");
		}
		break;
	case DATA_REPRESENTATION_BINARY:
		input.open(path.c_str(), ios::in | ios::binary);
		if (input) 
		{
    		input.read((char*)&t,     sizeof(var_t));
            input.read((char*)&dt,    sizeof(var_t));
            input.read((char*)&t_wc,  sizeof(var_t));
            input.read((char*)&n_obj, sizeof(uint32_t));
		}
		else 
		{
			throw string("Cannot open " + path + ".");
		}
		break;
	default:
		throw string("Parameter 'repres' is out of range.");
	}
	input.close();

	cout << " done" << endl;
}

void nbody::load_solution_data(string& path)
{
	ifstream input;

	cout << "Loading " << path << " ";

	data_rep_t repres = file::get_data_repres(path);
	switch (repres)
	{
	case DATA_REPRESENTATION_ASCII:
		input.open(path.c_str(), ios::in);
		if (input) 
		{
			load_ascii(input);
		}
		else 
		{
			throw string("Cannot open " + path + ".");
		}
		break;
	case DATA_REPRESENTATION_BINARY:
		input.open(path.c_str(), ios::in | ios::binary);
		if (input) 
		{
			load_binary(input);
		}
		else 
		{
			throw string("Cannot open " + path + ".");
		}
		break;
	default:
		throw string("Parameter 'repres' is out of range.");
	}
	input.close();

	cout << " done" << endl;
}

void nbody::load_ascii(ifstream& input)
{
    // Number of space and velocity coordinates
    const uint32_t nv = NDIM * n_obj;

    nbp_t::param_t* _p = (nbp_t::param_t*)h_p;

	for (uint32_t i = 0; i < n_obj; i++)
	{
        input >> h_md[i].id >> h_md[i].body_type >> h_md[i].active >> h_md[i].mig_type >> h_md[i].mig_stop_at >> h_md[i].unused1 >> h_md[i].unused2 >> h_md[i].unused3;
        input >> _p[i].density >> _p[i].mass >> _p[i].radius;
		uint32_t offset = NDIM * i;
		// position
		input >> h_y[offset+0] >> h_y[offset+1] >> h_y[offset+2];
		offset += nv;
		// velocity
		input >> h_y[offset+0] >> h_y[offset+1] >> h_y[offset+2];
	}
}

void nbody::load_binary(ifstream& input)
{
    input.read((char*)h_md, n_obj * sizeof(nbp_t::metadata_t));
    input.read((char*)h_p, n_obj * sizeof(nbp_t::param_t));
    input.read((char*)h_y, NVPO * n_obj * sizeof(var_t));
}

void nbody::print_solution(std::string& path_si, std::string& path_sd, data_rep_t repres)
{
	ofstream sout;

	switch (repres)
	{
	case DATA_REPRESENTATION_ASCII:
		sout.open(path_si.c_str(), ios::out | ios::app);
		break;
	case DATA_REPRESENTATION_BINARY:
		sout.open(path_si.c_str(), ios::out | ios::app | ios::binary);
		break;
	default:
		throw string("Parameter 'repres' is out of range.");
	}
	if (!sout)
	{
		throw string("Cannot open " + path_si + ".");
	}
	file::nbp::print_solution_info(sout, t, dt, t_wc, n_obj, repres);
	sout.close();

	switch (repres)
	{
	case DATA_REPRESENTATION_ASCII:
		sout.open(path_sd.c_str(), ios::out | ios::app);
		break;
	case DATA_REPRESENTATION_BINARY:
		sout.open(path_sd.c_str(), ios::out | ios::app | ios::binary);
		break;
	default:
		throw string("Parameter 'repres' is out of range.");
	}
	if (!sout)
	{
		throw string("Cannot open " + path_sd + ".");
	}
	file::nbp::print_solution_data(sout, n_obj, n_ppo, n_vpo, h_md, (nbp_t::param_t*)h_p, h_y, repres);
	sout.close();
}

void nbody::print_dump(std::string& path_si, std::string& path_sd)
{
    ofstream sout;

    sout.open(path_si.c_str(), ios::out | ios::binary);
    if (!sout)
    {
        throw string("Cannot open " + path_si + ".");
    }
    file::nbp::print_solution_info(sout, t, dt, t_wc, n_obj, DATA_REPRESENTATION_BINARY);
    sout.close();

    sout.open(path_sd.c_str(), ios::out | ios::binary);
    if (!sout)
    {
        throw string("Cannot open " + path_sd + ".");
    }
    file::nbp::print_solution_data(sout, n_obj, n_ppo, n_vpo, h_md, (nbp_t::param_t*)h_p, h_y, DATA_REPRESENTATION_BINARY);
    sout.close();
}

void nbody::print_integral(string& path)
{
	ofstream sout;

	sout.open(path.c_str(), ios::out | ios::app);
	if (sout)
	{
		sout.precision(16);
		sout.setf(ios::right);
		sout.setf(ios::scientific);

	    sout << setw(VAR_T_W) << t << SEP                    /* time of the record [day] (double)           */
			 << setw(VAR_T_W) << integral.R.x << SEP         /* x-position of the barycenter                */
			 << setw(VAR_T_W) << integral.R.y << SEP         /* y-position of the barycenter                */
			 << setw(VAR_T_W) << integral.R.z << SEP         /* z-position of the barycenter                */
			 << setw(VAR_T_W) << integral.V.x << SEP         /* x-velocity of the barycenter                */
			 << setw(VAR_T_W) << integral.V.y << SEP         /* y-velocity of the barycenter                */
			 << setw(VAR_T_W) << integral.V.z << SEP         /* z-velocity of the barycenter                */
			 << setw(VAR_T_W) << integral.c.x << SEP         /* x-angular momentum                          */
			 << setw(VAR_T_W) << integral.c.y << SEP         /* y-angular momentum                          */
			 << setw(VAR_T_W) << integral.c.z << SEP         /* z-angular momentum                          */
			 << setw(VAR_T_W) << integral.h << endl;         /* energy of the system                        */
	}
	else
	{
		throw string("Cannot open " + path + ".");
	}
	sout.close();
}

#undef NDIM
#undef NVPO
