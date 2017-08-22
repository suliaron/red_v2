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

__constant__ var_t dc_threshold[THRESHOLD_N];

namespace kernel_nbody
{
    __host__ __device__
        void store_event_data(event_name_t name, var_t t, var_t d, uint32_t surviv_idx, uint32_t merger_idx, const nbp_t::metadata_t* md, const nbp_t::param_t* p, const var3_t* r, const var3_t* v, nbp_t::event_data_t* evnt)
    {
        evnt->event_name = name;
        evnt->d = d;
        evnt->t = t;

        evnt->id0 = md[0].id;
        evnt->id1 = md[surviv_idx].id;
        evnt->id2 = md[merger_idx].id;

        evnt->idx0 = 0;
        evnt->idx1 = surviv_idx;
        evnt->idx2 = merger_idx;

        evnt->r0 = r[0];
        evnt->v0 = v[0];
        evnt->r1 = r[surviv_idx];
        evnt->v1 = v[surviv_idx];
        evnt->r2 = r[merger_idx];
        evnt->v2 = v[merger_idx];

        evnt->p0 = p[0];
        evnt->p1 = p[surviv_idx];
        evnt->p2 = p[merger_idx];

        if (EVENT_NAME_EJECTION == name)
        {
            evnt->rs = evnt->r1;
            evnt->vs = evnt->v1;
            evnt->ps = evnt->p1;
        }

        switch (name)
        {
        case EVENT_NAME_COLLISION:
            printf("  Collision t = %*.16e %5d %5d (%*.16e [R1 + R2])\n", VAR_T_W, t, md[surviv_idx].id, md[merger_idx].id, VAR_T_W, d / ((p[surviv_idx].radius + p[merger_idx].radius)));
            break;
        case EVENT_NAME_EJECTION:
            printf("   Ejection t = %*.16e %5d %5d\n", VAR_T_W, t, md[surviv_idx].id, md[merger_idx].id);
            break;
        case EVENT_NAME_HIT_CENTRUM:
            printf("Hit centrum t = %*.16e %5d %5d\n", VAR_T_W, t, md[surviv_idx].id, md[merger_idx].id);
            break;
        default:
            printf("Parameter 'name' is out of range.");
        }
    }

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

#if 0
    __global__
        void calc_grav_accel_naive(uint32_t n_obj, const nbp_t::metadata_t* md, const var3_t* r, const nbp_t::param_t* p, var3_t* a)
    {
        // i is the index of the SINK body
        const uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
        if (!md[i].active)
        {
            return;
        }

        if (i < n_obj)
        {
            // j is the index of the SOURCE body
            for (uint32_t j = 0; j < n_obj; j++)
            {
                if (i == j || !md[j].active) continue;
                body_body_grav_accel(r[i], r[j], p[j].mass, a[i]);
            }
        }
    } /* calc_grav_accel_naive () */

    __global__
        void calc_grav_accel_naive(uint2_t snk, uint2_t src, const nbp_t::metadata_t* md, const var3_t* r, const nbp_t::param_t* p, var3_t* a)
    {
        // i is the index of the SINK body
        const uint32_t i = snk.n1 + blockIdx.x * blockDim.x + threadIdx.x;
        if (!md[i].active)
        {
            return;
        }

        if (snk.n2 > i)
        {
            // j is the index of the SOURCE body
            for (uint32_t j = src.n1; j < src.n2; j++)
            {
                if (i == j || !md[j].active) continue;
                body_body_grav_accel(r[i], r[j], p[j].mass, a[i]);
            }
        }
    } /* calc_grav_accel_naive () */

    __global__
        void calc_grav_accel_tile(uint32_t n_obj, const nbp_t::metadata_t* md, const var3_t* r, const nbp_t::param_t* p, var3_t* a)
    {
        extern __shared__ var3_t sh_pos[];

        const uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
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
                // To avoid inactive bodies and overrun then input arrays
                if (!md[i].active || n_obj <= (tile * blockDim.x) + j)
                    break;
                // To avoid self-interaction and inactive bodies
                if (i == (tile * blockDim.x) + j || !md[(tile * blockDim.x) + j].active)
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
#endif

    __global__
        void calc_grav_accel_tile(uint2_t snk, uint2_t src, const nbp_t::metadata_t* md, const var3_t* r, const nbp_t::param_t* p, var3_t* a)
    {
        extern __shared__ var3_t sh_pos[];

        const uint32_t i = snk.n1 + blockIdx.x * blockDim.x + threadIdx.x;
        var3_t acc;
        var3_t my_pos;

        // To avoid overruning the r buffer
        if (i < snk.n2)
        {
            my_pos = r[i];
            acc = a[i];
        }
        // Note! : the for cycle must be outside the upper if clause, otherwise the sh_pos array will
        // not recive the input for the last tile! The reason is that some thread will be not considered
        // in the if (n_obj > idx) clause.
        for (uint32_t tile = 0; (tile * blockDim.x) < src.n2; tile++)
        {
            const uint32_t idx = src.n1 + tile * blockDim.x + threadIdx.x;
            // To avoid overruning the r and mass buffer
            if (idx < src.n2)
            {
                sh_pos[threadIdx.x] = r[idx];
            }
            __syncthreads();

            for (int j = 0; j < blockDim.x; j++)
            {
                // To avoid inactive bodies and overrun then input arrays
                if (i < snk.n2 && (!md[i].active || src.n2 <= src.n1 + (tile * blockDim.x) + j))
                    break;
                // To avoid self-interaction and inactive bodies
                if (i == (tile * blockDim.x) + j || !md[(tile * blockDim.x) + j].active)
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


    __global__
        void chk_coll(uint2_t snk, uint2_t src, var_t t, const nbp_t::metadata_t* md, const var3_t* r, const var3_t* v, const nbp_t::param_t* p, nbp_t::event_data_t* ed, uint32_t* n_event)
    {
        const uint32_t i = snk.n1 + blockIdx.x * blockDim.x + threadIdx.x;
        const var_t a_ref = dc_threshold[THRESHOLD_RADII_ENHANCE_FACTOR];

        if (i < snk.n2 && md[i].active)
        {
            for (uint32_t j = 1; j < src.n2; j++)
            {
                if (i == j || !md[j].active) continue;
                // compute r_ij = r_j - r_i [3 FLOPS] [6 read, 3 write]
                var3_t r_ij = { r[j].x - r[i].x, r[j].y - r[i].y, r[j].z - r[i].z };

                // compute square of r_ij vector [5 FLOPS + ] [3 read, 1 write]
                var_t d2 = SQR(r_ij.x) + SQR(r_ij.y) + SQR(r_ij.z);
                if (d2 < SQR(a_ref * (p[i].radius + p[j].radius)))
                {
                    uint32_t k = atomicAdd(n_event, 1);
                    uint32_t surviv_idx = i;
                    uint32_t merger_idx = j;
                    if (p[merger_idx].mass > p[surviv_idx].mass)
                    {
                        uint32_t idx = surviv_idx;
                        surviv_idx = merger_idx;
                        merger_idx = idx;
                    }
                    store_event_data(EVENT_NAME_COLLISION, t, sqrt(d2), surviv_idx, merger_idx, md, p, r, v, ed + k);
                }
            }
        }
    }

    __global__
        void chk_coll_tile(uint2_t snk, uint2_t src, var_t t, const nbp_t::metadata_t* md, const var3_t* r, const var3_t* v, const nbp_t::param_t* p, nbp_t::event_data_t* ed, uint32_t* n_event)
    {
        extern __shared__ var3_t sh_pos[];

        const uint32_t i = snk.n1 + blockIdx.x * blockDim.x + threadIdx.x;
        const var_t a_ref = dc_threshold[THRESHOLD_RADII_ENHANCE_FACTOR];
        var3_t my_pos;

        // To avoid overruning the r buffer
        if (i < snk.n2)
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
            if (idx < src.n2)
            {
                sh_pos[threadIdx.x] = r[idx];
            }
            __syncthreads();

            for (int j = 0; j < blockDim.x; j++)
            {
                // To avoid inactive bodies and overrun then input arrays
                if (i < snk.n2 && (!md[i].active || src.n2 <= src.n1 + (tile * blockDim.x) + j))
                    break;
                // To avoid self-interaction and inactive bodies
                if (i == (tile * blockDim.x) + j || !md[(tile * blockDim.x) + j].active)
                    continue;

                var3_t r_ij = { sh_pos[j].x - my_pos.x, sh_pos[j].y - my_pos.y, sh_pos[j].z - my_pos.z };
                // compute square of r_ij vector [5 FLOPS + ] [3 read, 1 write]
                var_t d2 = SQR(r_ij.x) + SQR(r_ij.y) + SQR(r_ij.z);
                if (d2 < SQR(a_ref * (p[i].radius + p[j].radius)))
                {
                    uint32_t k = atomicAdd(n_event, 1);
                    uint32_t surviv_idx = i;
                    uint32_t merger_idx = j;
                    if (p[merger_idx].mass > p[surviv_idx].mass)
                    {
                        uint32_t idx = surviv_idx;
                        surviv_idx = merger_idx;
                        merger_idx = idx;
                    }
                    store_event_data(EVENT_NAME_COLLISION, t, sqrt(d2), surviv_idx, merger_idx, md, p, r, v, ed + k);
                }
            }
            __syncthreads();
        }
    }

    __global__
        void chk_e_hc(uint2_t snk, var_t t, const nbp_t::metadata_t* md, const var3_t* r, const var3_t* v, const nbp_t::param_t* p, nbp_t::event_data_t* ed, uint32_t* n_event)
    {
        const uint32_t i = snk.n1 + blockIdx.x * blockDim.x + threadIdx.x;

        var_t ed2 = (0.0 < dc_threshold[THRESHOLD_EJECTION_DISTANCE] ? SQR(dc_threshold[THRESHOLD_EJECTION_DISTANCE]) : DBL_MAX);
        var_t hcd2 = SQR(dc_threshold[THRESHOLD_HIT_CENTRUM_DISTANCE]);

        // Ignore the inactive bodies (whose id < 0) and the star
        if (i < snk.n2 && md[i].active && md[i].body_type != BODY_TYPE_STAR)
        {
            uint32_t k = 0;

            // Calculate the distance from the barycenter
            var_t r2 = SQR(r[i].x) + SQR(r[i].y) + SQR(r[i].z);
            if (ed2 < r2)
            {
                k = atomicAdd(n_event, 1);
                store_event_data(EVENT_NAME_EJECTION, t, sqrt(r2), 0, i, md, p, r, v, ed + k);
            }
            if (hcd2 > r2)
            {
                k = atomicAdd(n_event, 1);
                store_event_data(EVENT_NAME_HIT_CENTRUM, t, sqrt(r2), 0, i, md, p, r, v, ed + k);
            }
        }
    }

} /* kernel_nbody */

nbody::nbody(string& path_si, string& path_sd, uint32_t n_obj, uint16_t n_ppo, size_t omd_size, comp_dev_t comp_dev) :
	ode(NDIM, n_obj, NVPO, n_ppo, omd_size, comp_dev)
{
	name = "Singular 3D n-body problem without collision";
	
	initialize();

    load_solution_info(path_si);
    load_solution_data(path_sd);

    calc_n_types();
	if (PROC_UNIT_GPU == comp_dev.proc_unit)
	{
		copy_vars(COPY_DIRECTION_TO_DEVICE);
		copy_params(COPY_DIRECTION_TO_DEVICE);
        copy_metadata(COPY_DIRECTION_TO_DEVICE);
	}
	calc_integral();
	tout = t;
}

nbody::nbody(string& path_si, string& path_sd, uint32_t n_obj, uint16_t n_ppo, size_t omd_size, size_t event_size, comp_dev_t comp_dev) :
    ode(NDIM, n_obj, NVPO, n_ppo, omd_size, event_size, comp_dev)
{
    name = "Singular 3D n-body problem including collision";

    initialize();

    load_solution_info(path_si);
    load_solution_data(path_sd);

    calc_n_types();
    if (PROC_UNIT_GPU == comp_dev.proc_unit)
    {
        copy_vars(COPY_DIRECTION_TO_DEVICE);
        copy_params(COPY_DIRECTION_TO_DEVICE);
        copy_metadata(COPY_DIRECTION_TO_DEVICE);

        ALLOCATE_DEVICE_ARRAY((void **)&d_n_event, sizeof(uint32_t));
    }
    calc_integral();
    tout = t;
}

nbody::~nbody()
{
    if (PROC_UNIT_GPU == comp_dev.proc_unit)
    {
        FREE_DEVICE_ARRAY((void **)&d_n_event);
    }
}

void nbody::initialize()
{
    print_oe = false;
    
    h_md = (nbp_t::metadata_t*)h_omd;
    d_md = (nbp_t::metadata_t*)d_omd;
      md = (nbp_t::metadata_t*)omd;

    h_ed = (nbp_t::event_data_t*)h_event;
    d_ed = (nbp_t::event_data_t*)d_event;
      ed = (nbp_t::event_data_t*)event;

    d_n_event = NULL;

    n_si  = 0;
    n_nsi = 0;
    n_ni  = 0;

    n_tpb_si  = 0;
    n_tpb_nsi = 0;
    n_tpb_ni  = 0;

    n_coll = 0;
    n_hitc = 0;
    n_ejec = 0;
}

void nbody::calc_n_types()
{
    n_si = n_nsi = n_ni = 0;

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

uint32_t nbody::calc_n_active()
{
    uint32_t n = 0;

    for (uint32_t i = 0; i < n_obj; i++)
    {
        if (h_md[i].active)
            n++;
    }
    return n;
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
    const nbp_t::param_t* p = (nbp_t::param_t*)h_par;
    
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
    const nbp_t::param_t* p = (nbp_t::param_t*)h_par;

    var3_t* a = (var3_t*)(dy + nv);

    // 1. Calculate the acceleration from the gravitational forces
    // Copy the velocities into dy
    memcpy(dy, v, nv * sizeof(var_t));

    // Clear the acceleration array: the += op can be used
	memset(a, 0, nv *sizeof(var_t));
    // -- Calculate the acceleration for the SI and NSI bodies
    uint2_t snk = { 0, n_si };
    uint2_t src = { 0, n_si + n_nsi };
    for (uint32_t i = snk.n1; i < snk.n2; i++)
	{
        if (!h_md[i].active) continue;
		for (uint32_t j = i+1; j < src.n2; j++)
		{
            if (!h_md[j].active) continue;
            body_body_grav_accel(r[i], r[j], p[i].mass, p[j].mass, a[i], a[j]);
        }
	}

    // -- Calculate the acceleration for the NI bodies
    snk.n1 = n_si + n_nsi, snk.n2 = n_obj;
    src.n1 = 0, src.n2 = n_si + n_nsi;
    for (uint32_t i = snk.n1; i < snk.n2; i++)
    {
        if (!h_md[i].active) continue;
        for (uint32_t j = src.n1; j < src.n2; j++)
        {
            if (!h_md[j].active) continue;
            body_body_grav_accel(r[i], r[j], p[j].mass, a[i]);
        }
    }
}

inline
void nbody::body_body_grav_accel(const var3_t& ri, const var3_t& rj, var_t mj, var3_t& ai)
{
    // compute r_ij = r_j - r_i [3 FLOPS] [6 read, 3 write]
    var3_t r_ij = { rj.x - ri.x, rj.y - ri.y, rj.z - ri.z };

    // compute norm square of d vector [5 FLOPS] [3 read, 1 write]
    var_t d2 = SQR(r_ij.x) + SQR(r_ij.y) + SQR(r_ij.z);
    var_t d = sqrt(d2);
    var_t s = K2 * mj / (d * d2);

    ai.x += s * r_ij.x;
    ai.y += s * r_ij.y;
    ai.z += s * r_ij.z;
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
    static uint32_t last_n_si  = 0;
    static uint32_t last_n_nsi = 0;
    static uint32_t last_n_ni  = 0;

    // Number of space and velocity coordinates
    const uint32_t nv = NDIM * n_obj;

    // Create aliases
    const var3_t* r = (var3_t*)y_temp;
    const var3_t* v = (var3_t*)(y_temp + nv);
    const nbp_t::param_t* p = (nbp_t::param_t*)d_par;

    var3_t* a = (var3_t*)(dy + nv);
    // Clear the acceleration array: the += op can be used
    cudaMemset(a, 0, n_obj * sizeof(var3_t));

    // 1. Calculate the acceleration from the gravitational forces
    // Copy the velocities into dy
    CUDA_SAFE_CALL(cudaMemcpy(dy, v, nv * sizeof(var_t), cudaMemcpyDeviceToDevice));

    uint2_t snk = { 0, n_si };
    uint2_t src = { 0, n_si + n_nsi };
    if ((snk.n2 - snk.n1) * (src.n2 - src.n1) > 1)
    {
        if (last_n_si != n_si)
        {
            last_n_si = n_si;
            printf("Searching for the optimal thread number [%s] (SI type bodies) ", __FUNCTION__);
            n_tpb_si = calc_opt_tpb(stage, snk, src, curr_t, r, p);
            printf(" [%3d] done.\n", n_tpb_si);
        }
        // Compute the acceleration and store it in a
        float dt_GPU = gpu_calc_grav_accel_tile(stage, snk, src, n_tpb_si, curr_t, r, p, a);
    }
//    redutil2::print_array("", " SI a", nv, (var_t*)a, (PROC_UNIT_CPU == comp_dev.proc_unit ? MEM_LOC_HOST : MEM_LOC_DEVICE));

    snk.n1 = n_si, snk.n2 = n_si + n_nsi;
    src.n1 = 0, src.n2 = n_si;
    if ((snk.n2 - snk.n1) * (src.n2 - src.n1) > 1)
    {
        if (last_n_nsi != n_nsi)
        {
            last_n_nsi = n_nsi;
            printf("Searching for the optimal thread number [%s] (NSI type bodies) ", __FUNCTION__);
            n_tpb_nsi = calc_opt_tpb(stage, snk, src, curr_t, r, p);
            printf(" [%3d] done.\n", n_tpb_nsi);
        }
        float dt_GPU = gpu_calc_grav_accel_tile(stage, snk, src, n_tpb_nsi, curr_t, r, p, a);
    }
//    redutil2::print_array("", "NSI a", nv, (var_t*)a, (PROC_UNIT_CPU == comp_dev.proc_unit ? MEM_LOC_HOST : MEM_LOC_DEVICE));

    snk.n1 = n_si + n_nsi, snk.n2 = n_obj;
    src.n1 = 0, src.n2 = n_si + n_nsi;
    if ((snk.n2 - snk.n1) * (src.n2 - src.n1) > 1)
    {
        if (last_n_ni != n_ni)
        {
            last_n_ni = n_ni;
            printf("Searching for the optimal thread number [%s] (NI type bodies) ", __FUNCTION__);
            n_tpb_ni = calc_opt_tpb(stage, snk, src, curr_t, r, p);
            printf(" [%3d] done.\n", n_tpb_ni);
        }
        float dt_GPU = gpu_calc_grav_accel_tile(stage, snk, src, n_tpb_ni, curr_t, r, p, a);
    }
//    redutil2::print_array("", " NI a", nv, (var_t*)a, (PROC_UNIT_CPU == comp_dev.proc_unit ? MEM_LOC_HOST : MEM_LOC_DEVICE));
    // ------- END -------

    // 2. Calculate the accelerations from other forces
    // ...
}

uint16_t nbody::calc_opt_tpb(uint16_t stage, uint2_t snk, uint2_t src, var_t curr_t, const var3_t* r, const nbp_t::param_t* p)
{
    uint16_t opt_n_tpb = 0;
    float min_dt = 1.0e10;

    // Temporary array to store the result of the gpu_calc_grav_accel_tile() function
    var3_t *d_a = NULL;
    ALLOCATE_DEVICE_ARRAY((void**)&(d_a), n_obj * sizeof(var3_t));
    for (unsigned int i = 16; i <= 512; i += 16)
    {
        putc('.', stdout);
        // Clear the acceleration array: the += op can be used
        cudaMemset(d_a, 0, n_obj * sizeof(var3_t));
        float dt_GPU = gpu_calc_grav_accel_tile(stage, snk, src, i, curr_t, r, p, d_a);
        if (dt_GPU < min_dt)
        {
            min_dt = dt_GPU;
            opt_n_tpb = i;
        }
    }
    FREE_DEVICE_ARRAY((void **)&d_a);

    return opt_n_tpb;
}

//float nbody::gpu_calc_grav_accel_naive(uint16_t stage, uint2_t snk, uint2_t src, unsigned int n_tpb, var_t curr_t, const var3_t* r, const nbp_t::param_t* p, var3_t* a)
//{
//    cudaEvent_t start, stop;
//
//    CUDA_SAFE_CALL(cudaEventCreate(&start));
//    CUDA_SAFE_CALL(cudaEventCreate(&stop));
//
//    set_kernel_launch_param(snk.n2 - snk.n1, n_tpb, grid, block);
//
//    CUDA_SAFE_CALL(cudaEventRecord(start, 0));
//    kernel_nbody::calc_grav_accel_naive <<< grid, block >>>(snk, src, d_md, r, p, a);
//    CUDA_CHECK_ERROR();
//
//    CUDA_SAFE_CALL(cudaEventRecord(stop, 0));
//    CUDA_SAFE_CALL(cudaEventSynchronize(stop));
//
//    // Computes the elapsed time between two events in milliseconds with a resolution of around 0.5 microseconds.
//    float elapsed_time = 0.0f;
//    CUDA_SAFE_CALL(cudaEventElapsedTime(&elapsed_time, start, stop));
//
//    return elapsed_time;
//}

float nbody::gpu_calc_grav_accel_tile(uint16_t stage, uint2_t snk, uint2_t src, unsigned int n_tpb, var_t curr_t, const var3_t* r, const nbp_t::param_t* p, var3_t* a)
{
    cudaEvent_t start, stop;

    CUDA_SAFE_CALL(cudaEventCreate(&start));
    CUDA_SAFE_CALL(cudaEventCreate(&stop));

    set_kernel_launch_param(snk.n2 - snk.n1, n_tpb, grid, block);
    size_t sh_mem_size = n_tpb * sizeof(var3_t);

    CUDA_SAFE_CALL(cudaEventRecord(start, 0));
    kernel_nbody::calc_grav_accel_tile <<< grid, block, sh_mem_size >>>(snk, src, d_md, r, p, a);
    CUDA_CHECK_ERROR();

    CUDA_SAFE_CALL(cudaEventRecord(stop, 0));
    CUDA_SAFE_CALL(cudaEventSynchronize(stop));

    // Computes the elapsed time between two events in milliseconds with a resolution of around 0.5 microseconds.
    float elapsed_time = 0.0f;
    CUDA_SAFE_CALL(cudaEventElapsedTime(&elapsed_time, start, stop));

    return elapsed_time;
}

uint32_t nbody::chk_event(string path, var_t* threshold)
{
    uint32_t n_event = 0;

    if (PROC_UNIT_CPU == comp_dev.proc_unit)
    {
        n_event = cpu_chk_event(threshold);
    }
    else
    {
        n_event = gpu_chk_event(threshold);
    }
    if (n_event)
    {
        handle_event(n_event);
        print_event_data(path, n_event);
    }

    return n_event;
}

uint32_t nbody::cpu_chk_event(var_t* threshold)
{
    // Number of space and velocity coordinates
    const uint32_t nv = NDIM * n_obj;
    // Create aliases
    const var3_t* r = (var3_t*)h_y;
    const var3_t* v = (var3_t*)(h_y + nv);
    const nbp_t::param_t* p = (nbp_t::param_t*)h_par;

    uint32_t n_event = 0;
    if (0.0 < threshold[THRESHOLD_HIT_CENTRUM_DISTANCE] || 0.0 < threshold[THRESHOLD_EJECTION_DISTANCE])
    {
        for (uint32_t j = 1; j < n_obj; j++)
        {
            if (!h_md[j].active) continue;
            // compute r_ij = r_j - r_i [3 FLOPS] [6 read, 3 write]
            var3_t r_ij = { r[j].x - r[0].x, r[j].y - r[0].y, r[0].z - r[0].z };
            // compute norm square of d vector [5 FLOPS] [3 read, 1 write]
            var_t d2 = SQR(r_ij.x) + SQR(r_ij.y) + SQR(r_ij.z);

            if (d2 > SQR(threshold[THRESHOLD_EJECTION_DISTANCE]))
            {
                kernel_nbody::store_event_data(EVENT_NAME_EJECTION, t, sqrt(d2), 0, j, h_md, p, r, v, h_ed + n_event);
                n_event++; n_ejec++;
            }
            if (d2 < SQR(threshold[THRESHOLD_HIT_CENTRUM_DISTANCE]))
            {
                kernel_nbody::store_event_data(EVENT_NAME_HIT_CENTRUM, t, sqrt(d2), 0, j, h_md, p, r, v, h_ed + n_event);
                n_event++; n_hitc++;
            }
        }
    }

    if (0.0 < threshold[THRESHOLD_RADII_ENHANCE_FACTOR])
    {
        var_t a_ref = threshold[THRESHOLD_RADII_ENHANCE_FACTOR];
        for (uint32_t i = 1; i < n_obj; i++)
        {
            if (!h_md[i].active) continue;
            //var_t min_d = 0.1;
            for (uint32_t j = i + 1; j < n_obj; j++)
            {
                if (!h_md[j].active) continue;
                // compute r_ij = r_j - r_i [3 FLOPS] [6 read, 3 write]
                var3_t r_ij = { r[j].x - r[i].x, r[j].y - r[i].y, r[j].z - r[i].z };
                // compute norm square of d vector [5 FLOPS] [3 read, 1 write]
                var_t d2 = SQR(r_ij.x) + SQR(r_ij.y) + SQR(r_ij.z);

                if (d2 < SQR(a_ref * (p[i].radius + p[j].radius)))
                {
                    uint32_t surviv_idx = i;
                    uint32_t merger_idx = j;
                    if (p[merger_idx].mass > p[surviv_idx].mass)
                    {
                        uint32_t idx = surviv_idx;
                        surviv_idx = merger_idx;
                        merger_idx = idx;
                    }
                    kernel_nbody::store_event_data(EVENT_NAME_COLLISION, t, sqrt(d2), surviv_idx, merger_idx, h_md, p, r, v, h_ed + n_event);
                    n_event++; n_coll++;
                }
                //if (min_d > d)
                //{
                //    // L suffix denotes the length of the corresponding vector
                //    // 2 suffix denotes the square of the length of the corresponding vector
                //    var3_t v_ij = {v[j].x - v[i].x, v[j].y - v[i].y, v[j].z - v[i].z};
                //    var_t r_ijL = d;
                //    var_t v_ijL = sqrt(SQR(v_ij.x) + SQR(v_ij.y) + SQR(v_ij.z));

                //    var_t rv = r_ij.x * v_ij.x + r_ij.y * v_ij.y + r_ij.z * v_ij.z;
                //    var_t alpha = PI - acos(rv / (r_ijL * v_ijL));
                //    var_t b = r_ijL * sin(alpha);
                //    var_t v_esc2 = 2.0 * K2 * p[i].mass / p[i].radius;
                //    var_t Gamma = 1.0 + v_esc2 / SQR(v_ijL);
                //    if (SQR(b) <= SQR(p[i].radius) * Gamma)
                //    {
                //        printf("%*.16e Collision expected between %d %d (distance = %*.16e [R1+R2])\n", VAR_T_W, t, h_md[i].id, h_md[j].id, VAR_T_W, r_ijL / (p[i].radius + p[j].radius));

                //        printf("ri_Vec.x  = %25.16e & ri_Vec.y  = %25.16e & ri_Vec.z  = %25.16e\n", r[i].x, r[i].y, r[i].z);
                //        printf("rj_Vec.x  = %25.16e & rj_Vec.y  = %25.16e & rj_Vec.z  = %25.16e\n", r[j].x, r[j].y, r[j].z);
                //        printf("rij_Vec.x = %25.16e & rij_Vec.y = %25.16e & rij_Vec.z = %25.16e\n", r_ij.x, r_ij.y, r_ij.z);
                //        printf("vi_Vec.x  = %25.16e & vi_Vec.y  = %25.16e & vi_Vec.z  = %25.16e\n", v[i].x, v[i].y, v[i].z);
                //        printf("vj_Vec.x  = %25.16e & vj_Vec.y  = %25.16e & vj_Vec.z  = %25.16e\n", v[j].x, v[j].y, v[j].z);
                //        printf("vij_Vec.x = %25.16e & vij_Vec.y = %25.16e & vij_Vec.z = %25.16e\n", v_ij.x, v_ij.y, v_ij.z);
                //        printf("m_i = %25.16e & radius_i = %25.16e\n", p[i].mass, p[i].radius);
                //        printf("m_j = %25.16e & radius_j = %25.16e\n", p[j].mass, p[j].radius);
                //        printf("rv = %25.16e\n", rv);
                //        printf("alpha = %25.16e\n", alpha);
                //        printf("b = %25.16e\n", b);
                //        printf("v_esc2 = %25.16e\n", v_esc2);
                //        printf("Gamma = %25.16e\n", Gamma);
                //    }
                //}
            }
        }
    }

    return n_event;
}

uint32_t nbody::gpu_chk_event(var_t* threshold)
{
    static bool first_call = true;

    // Number of space and velocity coordinates
    const uint32_t nv = NDIM * n_obj;
    // Create aliases
    const var3_t* r = (var3_t*)d_y;
    const var3_t* v = (var3_t*)(d_y + nv);
    const nbp_t::param_t* p = (nbp_t::param_t*)d_par;

    if (first_call)
    {
        first_call = false;
        copy_constant_to_device(dc_threshold, threshold, THRESHOLD_N * sizeof(var_t));
    }

    uint32_t n_event = 0;
    // Set the *d_n_event to zero
    CUDA_SAFE_CALL(cudaMemset(d_n_event, 0, sizeof(uint32_t)));

    if (0.0 < threshold[THRESHOLD_HIT_CENTRUM_DISTANCE] || 0.0 < threshold[THRESHOLD_EJECTION_DISTANCE])
    {
        uint2_t snk = { 1, n_obj };
        set_kernel_launch_param(snk.n2 - snk.n1, 256, grid, block);
        kernel_nbody::chk_e_hc<< < grid, block >> >(snk, t, d_md, r, v, p, d_ed, d_n_event);
        CUDA_CHECK_ERROR();
    }

    if (0.0 < threshold[THRESHOLD_RADII_ENHANCE_FACTOR])
    {
        uint2_t snk = { 1, n_si };
        uint2_t src = { 0, n_obj };

        unsigned int _n_tpb = (n_tpb_si == 0 ? 256 : n_tpb_si);
        set_kernel_launch_param(snk.n2 - snk.n1, _n_tpb, grid, block);
        //size_t sh_mem_size = _n_tpb * sizeof(var3_t);
        //kernel_nbody::chk_coll << < grid, block, sh_mem_size >> >(snk, src, t, d_md, r, v, p, d_ed, d_n_event);
        kernel_nbody::chk_coll << < grid, block >> >(snk, src, t, d_md, r, v, p, d_ed, d_n_event);
        CUDA_CHECK_ERROR();
    }

    copy_array_to_host((void *)&n_event, (void *)d_n_event, sizeof(uint32_t));
    if (n_event)
    {
        copy_array_to_host((void *)h_ed, (void *)d_ed, n_event * sizeof(nbp_t::event_data_t));
        n_event = redutil2::tools::nbp::remove_duplicate(h_ed, n_event);
    }

    return n_event;
}

void nbody::handle_event(uint32_t n_event)
{
    static uint32_t id = 1;

    // Number of space and velocity coordinates
    const uint32_t nv = NDIM * n_obj;
    // Create aliases
    const nbp_t::param_t* p = (nbp_t::param_t*)h_par;
    const var3_t* r = (var3_t*)h_y;
    const var3_t* v = (var3_t*)(h_y + nv);

    uint32_t _n_hc = 0;
    var3_t r0, v0;
    nbp_t::param_t p0;
    for (uint32_t i = 0; i < n_event; i++)
    {
        h_ed[i].id = id++;

        // If execution is on device you have no guarantee that the h_y contains the right data, therefore 
        // in the first iteration we get the data from the event data, than updateted data are written into
        // p, r, and v array so we take those from there.
        if (0 == _n_hc)
        {
            p0 = h_ed[i].p0;
            r0 = h_ed[i].r0;
            v0 = h_ed[i].v0;
        }
        else
        {
            p0 = p[0];
            r0 = r[0];
            v0 = v[0];
        }

        // If more than one hit centrum were detetcted the result of the previous events has to take into account to correctly (?) compute the next hit centrum
        if (EVENT_NAME_COLLISION == h_ed[i].event_name)
        {
            h_ed[i].p0 = p0;
            h_ed[i].r0 = r0;
            h_ed[i].v0 = v0;
            handle_collision_pair(h_ed + i);
        }
        else if (EVENT_NAME_EJECTION == h_ed[i].event_name)
        {
            h_ed[i].p0 = h_ed[i].p1 = p[0];
            h_ed[i].r0 = h_ed[i].r1 = r[0];
            h_ed[i].v0 = h_ed[i].v1 = v[0];
            handle_ejection(h_ed + i);
        }
        else
        {
            h_ed[i].p0 = h_ed[i].p1 = p[0];
            h_ed[i].r0 = h_ed[i].r1 = r[0];
            h_ed[i].v0 = h_ed[i].v1 = v[0];
            handle_collision_pair(h_ed + i);
            _n_hc++;
        }
    }
}

void nbody::handle_collision_pair(nbp_t::event_data_t *collision)
{
    // Number of space and velocity coordinates
    const uint32_t nv = NDIM * n_obj;
    // Create aliases
    var3_t* r = (var3_t*)h_y;
    var3_t* v = (var3_t*)(h_y + nv);
    nbp_t::param_t* p = (nbp_t::param_t*)h_par;

    uint32_t surviv_idx = collision->idx1;
    uint32_t merger_idx = collision->idx2;

    // Compute the kinetic energy of the two bodies before the collision
    //var_t T0 = 0.5 * (collision->p1.mass * (SQR(collision->v1.x) + SQR(collision->v1.y) + SQR(collision->v1.z)) + 
    //	              collision->p2.mass * (SQR(collision->v2.x) + SQR(collision->v2.y) + SQR(collision->v2.z)));

    // Calculate position and velocitiy of the new object
    tools::calc_position_after_collision(collision->p1.mass, collision->p2.mass, &(collision->r1), &(collision->r2), collision->rs);
    tools::calc_velocity_after_collision(collision->p1.mass, collision->p2.mass, &(collision->v1), &(collision->v2), collision->vs);
    // Calculate physical properties of the new object
    tools::calc_physical_properties(collision->p1, collision->p2, collision->ps);
    // Update position, velocity and physical properties of survivor
    r[surviv_idx] = collision->rs;
    v[surviv_idx] = collision->vs;
    p[surviv_idx] = collision->ps;

    h_md[merger_idx].active = false;
    v[merger_idx].x = v[merger_idx].y = v[merger_idx].z = 0.0;

    if (PROC_UNIT_GPU == comp_dev.proc_unit)
    {
        // Create aliases
        var3_t* d_r = (var3_t*)d_y;
        var3_t* d_v = (var3_t*)(d_y + nv);
        nbp_t::param_t* d_p = (nbp_t::param_t*)d_par;

        copy_array_to_device((void **)&d_r[surviv_idx], (void **)&r[surviv_idx], sizeof(var3_t));
        //printf("d_r[surviv_idx]: "); redutil2::prn(&d_r[surviv_idx]);
        copy_array_to_device((void **)&d_v[surviv_idx], (void **)&v[surviv_idx], sizeof(var3_t));
        //printf("d_v[surviv_idx]: "); redutil2::prn(&d_v[surviv_idx]);
        copy_array_to_device((void **)&d_p[surviv_idx], (void **)&p[surviv_idx], sizeof(nbp_t::param_t));
        //printf("d_p[surviv_idx]: "); redutil2::prn(&d_p[surviv_idx]);
        copy_array_to_device((void **)&d_v[merger_idx], (void **)&v[merger_idx], sizeof(var3_t));
        //printf("d_v[merger_idx]: "); redutil2::prn(&d_v[merger_idx]);
        copy_array_to_device((void **)&d_md[merger_idx], (void **)&h_md[merger_idx], sizeof(nbp_t::metadata_t));
        //printf("d_md[merger_idx]: "); redutil2::prn(&d_md[merger_idx]);
    }

    // Compute the kinetic energy of the surviver
    //var_t T1 = 0.5 * (collision->ps.mass * (SQR(collision->vs.x) + SQR(collision->vs.y) + SQR(collision->vs.z)));
    //var_t dT = T1 - T0;
}

void nbody::handle_ejection(nbp_t::event_data_t *ejection)
{
    // Number of space and velocity coordinates
    const uint32_t nv = NDIM * n_obj;
    // Create aliases
    var3_t* v = (var3_t*)(h_y + nv);

    uint32_t idx = ejection->idx2;

    h_md[idx].active = false;
    v[idx].x = v[idx].y = v[idx].z = 0.0;

    if (PROC_UNIT_GPU == comp_dev.proc_unit)
    {
        // Create aliases
        var3_t* d_v = (var3_t*)(d_y + nv);
        //nbp_t::param_t* d_p = (nbp_t::param_t*)d_par;

        copy_array_to_device((void **)&d_v[idx], (void **)&v[idx], sizeof(var3_t));
        //printf("d_v[idx]: "); redutil2::prn(&d_v[idx]);
        copy_array_to_device((void **)&d_md[idx], (void **)&h_md[idx], sizeof(nbp_t::metadata_t));
        //printf("d_md[idx]: "); redutil2::prn(&d_md[idx]);
    }
}

void nbody::rebuild_host_array()
{
    // Number of space and velocity coordinates
    const uint32_t nv = NDIM * n_obj;
    // Create aliases
    var3_t* r = (var3_t*)h_y;
    var3_t* v = (var3_t*)(h_y + nv);
    nbp_t::param_t* p = (nbp_t::param_t*)h_par;

    var_t* _y = NULL;
    nbp_t::param_t*    _p  = NULL;
    nbp_t::metadata_t* _md = NULL;

    uint32_t n_active = calc_n_active();

    ALLOCATE_HOST_ARRAY((void **)&_y,  n_active * NVPO * sizeof(var_t));
    ALLOCATE_HOST_ARRAY((void **)&_p,  n_active * sizeof(nbp_t::param_t));
    ALLOCATE_HOST_ARRAY((void **)&_md, n_active * sizeof(nbp_t::metadata_t));

    var3_t* _r = (var3_t*)_y;
    var3_t* _v = (var3_t*)(_y + NDIM * n_active);

    uint32_t k = 0;
    for (uint32_t i = 0; i < n_obj; i++)
    {
        // Active?
        if (h_md[i].active && n_active > k)
        {
            _r[k] = r[i];
            _v[k] = v[i];
            _p[k] = p[i];
            _md[k] = h_md[i];
            k++;
        }
    }

    memcpy((void *)h_y,  _y, NVPO * n_active * sizeof(var_t));
    memcpy((void *)p,    _p,  n_active * sizeof(nbp_t::param_t));
    memcpy((void *)h_md, _md, n_active * sizeof(nbp_t::metadata_t));

    set_n_obj(n_active);
    this->omd_size = n_active * sizeof(nbp_t::metadata_t);
    calc_n_types();
    if (n_obj != calc_n_active() || n_obj != (n_si + n_nsi + n_ni))
    {
        throw string("Number of active bodies does not equal to the total number of bodies.");
    }

    FREE_HOST_ARRAY((void **)&_y);
    FREE_HOST_ARRAY((void **)&_p);
    FREE_HOST_ARRAY((void **)&_md);
}

void nbody::calc_integral()
{
    static bool first_call = true;

    // Number of space and velocity coordinates
    const uint32_t nv = NDIM * n_obj;
    // Create aliases
    const var3_t* r = (var3_t*)h_y;
    const var3_t* v = (var3_t*)(h_y + nv);
    const nbp_t::param_t* p = (nbp_t::param_t*)h_par;

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

    nbp_t::param_t* p = (nbp_t::param_t*)h_par;

    uint32_t k = 0;
	for (uint32_t i = 0; i < n_obj; i++)
	{
        if ((100.0 * i) / (var_t)n_obj > k)
        {
            putchar('.'); k++;
        }
        input >> h_md[i].id >> h_md[i].body_type >> h_md[i].active >> h_md[i].mig_type >> h_md[i].mig_stop_at >> h_md[i].unused1 >> h_md[i].unused2 >> h_md[i].unused3;
        input >> p[i].density >> p[i].mass >> p[i].radius;
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
    input.read((char*)h_par, n_obj * sizeof(nbp_t::param_t));
    input.read((char*)h_y, NVPO * n_obj * sizeof(var_t));
}

void nbody::print_solution(string& path_si, string& path_sd, data_rep_t repres)
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
	file::nbp::print_solution_data(sout, n_obj, n_ppo, n_vpo, h_md, (nbp_t::param_t*)h_par, h_y, repres);
	sout.close();

    if (print_oe)
    {
        string dir = file::get_directory(path_sd);
        string fn = file::get_filename_without_ext(path_sd) + "_oe.txt";
        string path_oe = file::combine_path(dir, fn);
        sout.open(path_oe.c_str(), ios::out | ios::app);
        if (!sout)
        {
            throw string("Cannot open " + path_oe + ".");
        }

        // Number of space and velocity coordinates
        const uint32_t nv = NDIM * n_obj;
        // Create aliases
        const var3_t* r = (var3_t*)h_y;
        const var3_t* v = (var3_t*)(h_y + nv);
        const nbp_t::param_t* p = (nbp_t::param_t*)h_par;
        orbelem_t oe = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

        for (uint32_t i = 1; i < n_obj; i++)
        {
            var_t mu = K2 * (p[0].mass + p[i].mass);
            var3_t dr = { r[i].x - r[0].x, r[i].y - r[0].y, r[i].z - r[0].z };
            var3_t dv = { v[i].x - v[0].x, v[i].y - v[0].y, v[i].z - v[0].z };
            int r = tools::calc_oe(mu, &dr, &dv, &oe);
            if (r)
            {
                oe.sma = oe.ecc = oe.inc = oe.peri = oe.node = oe.mean = 0.0;
            }
            file::nbp::print_oe_record(sout, t, oe, p[i], h_md[i]);
        }
        sout.close();
    }
}

void nbody::print_dump(string& path_si, string& path_sd)
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
    file::nbp::print_solution_data(sout, n_obj, n_ppo, n_vpo, h_md, (nbp_t::param_t*)h_par, h_y, DATA_REPRESENTATION_BINARY);
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

void nbody::print_event_data(string& path, uint32_t n_event)
{
    static uint32_t int_t_w = 8;
    static uint32_t var_t_w = 25;
    string e_names[] = { "HIT_CENTRUM", "EJECTION", "COLLISION" };

    ofstream sout;
    sout.open(path.c_str(), ios::out | ios::app);
    if (sout)
    {
        sout.precision(16);
        sout.setf(ios::right);
        sout.setf(ios::scientific);

        for (uint32_t i = 0; i < n_event; i++)
        {
            sout << setw(6) << h_ed[i].id << SEP
                << setw(16) << e_names[h_ed[i].event_name] << SEP       /* type of the event                              []                  */
                << setw(var_t_w) << h_ed[i].t << SEP                    /* time of the event                              [day]               */
                << setw(var_t_w) << h_ed[i].d << SEP                    /* distance of the two bodies                     [AU]                */
                << setw(int_t_w) << h_ed[i].id0 << SEP		            /* id of the star                                 []                  */
                << setw(int_t_w) << h_ed[i].id1 << SEP		            /* id of the survivor                             []                  */
                << setw(int_t_w) << h_ed[i].id2 << SEP		            /* id of the merger                               []                  */

                                                               /* BEFORE THE EVENT */
                << setw(var_t_w) << h_ed[i].p0.mass << SEP              /* star's mass                                    [solar mass]        */
                << setw(var_t_w) << h_ed[i].p0.density << SEP           /*        density                                 [solar mass / AU^3] */
                << setw(var_t_w) << h_ed[i].p0.radius << SEP            /*        radius                                  [AU]                */
                << setw(var_t_w) << h_ed[i].r0.x << SEP                 /*        x-coordiante in barycentric system      [AU]                */
                << setw(var_t_w) << h_ed[i].r0.y << SEP                 /*        y-coordiante                            [AU]                */
                << setw(var_t_w) << h_ed[i].r0.z << SEP                 /*        z-coordiante                            [AU]                */
                << setw(var_t_w) << h_ed[i].v0.x << SEP                 /*        x-velocity                              [AU / day]          */
                << setw(var_t_w) << h_ed[i].v0.y << SEP                 /*        y-velocity                              [AU / day]          */
                << setw(var_t_w) << h_ed[i].v0.z << SEP                 /*        z-velocity                              [AU / day]          */

                << setw(var_t_w) << h_ed[i].p1.mass << SEP              /* survivor's mass                                [solar mass]        */
                << setw(var_t_w) << h_ed[i].p1.density << SEP           /*            density                             [solar mass / AU^3] */
                << setw(var_t_w) << h_ed[i].p1.radius << SEP            /*            radius                              [AU]                */
                << setw(var_t_w) << h_ed[i].r1.x << SEP                 /*            x-coordiante in barycentric system  [AU]                */
                << setw(var_t_w) << h_ed[i].r1.y << SEP                 /*            y-coordiante                        [AU]                */
                << setw(var_t_w) << h_ed[i].r1.z << SEP                 /*            z-coordiante                        [AU]                */
                << setw(var_t_w) << h_ed[i].v1.x << SEP                 /*            x-velocity                          [AU / day]          */
                << setw(var_t_w) << h_ed[i].v1.y << SEP                 /*            y-velocity                          [AU / day]          */
                << setw(var_t_w) << h_ed[i].v1.z << SEP                 /*            z-velocity                          [AU / day]          */

                << setw(var_t_w) << h_ed[i].p2.mass << SEP              /* merger's mass                                  [solar mass]        */
                << setw(var_t_w) << h_ed[i].p2.density << SEP           /*          density                               [solar mass / AU^3] */
                << setw(var_t_w) << h_ed[i].p2.radius << SEP            /*          radius                                [AU]                */
                << setw(var_t_w) << h_ed[i].r2.x << SEP                 /*          x-coordiante in barycentric system    [AU]                */
                << setw(var_t_w) << h_ed[i].r2.y << SEP                 /*          y-coordiante                          [AU]                */
                << setw(var_t_w) << h_ed[i].r2.z << SEP                 /*          z-coordiante                          [AU]                */
                << setw(var_t_w) << h_ed[i].v2.x << SEP                 /*          x-velocity                            [AU / day]          */
                << setw(var_t_w) << h_ed[i].v2.y << SEP                 /*          y-velocity                            [AU / day]          */
                << setw(var_t_w) << h_ed[i].v2.z << SEP                 /*          z-velocity                            [AU / day]          */

                                                               /* AFTER THE EVENT */
                << setw(var_t_w) << h_ed[i].ps.mass << SEP              /* survivor's mass                                [solar mass]        */
                << setw(var_t_w) << h_ed[i].ps.density << SEP           /*            density                             [solar mass / AU^3] */
                << setw(var_t_w) << h_ed[i].ps.radius << SEP            /*            radius                              [AU]                */
                << setw(var_t_w) << h_ed[i].rs.x << SEP                 /*            x-coordiante in barycentric system  [AU]                */
                << setw(var_t_w) << h_ed[i].rs.y << SEP                 /*            y-coordiante                        [AU]                */
                << setw(var_t_w) << h_ed[i].rs.z << SEP                 /*            z-coordiante                        [AU]                */
                << setw(var_t_w) << h_ed[i].vs.x << SEP                 /*            x-velocity                          [AU / day]          */
                << setw(var_t_w) << h_ed[i].vs.y << SEP                 /*            y-velocity                          [AU / day]          */
                << setw(var_t_w) << h_ed[i].vs.z << SEP                 /*            z-velocity                          [AU / day]          */
                << endl;
        }
        sout.close();
    }
    else
    {
        throw string("Cannot open " + path + ".");
    }
}


#undef NDIM
#undef NVPO
