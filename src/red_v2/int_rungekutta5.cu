#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include "ode.h"
#include "int_rungekutta5.h"

#include "macro.h"
#include "redutil2.h"

using namespace redutil2;

#if 0
/*
 * Fehlberg, E.
 * "Classical Fifth-, Sixth-, Seventh-, and Eighth-Order Runge-Kutta Formulas with Stepsize Control"
 * NASA-TR-R-287 (https://nix.nasa.gov/search.jsp?R=19680027281&qs=N%253D4294957355)
 * p. 26 Table II. RK5(6)
 */
// The Runge-Kutta matrix
var_t int_rungekutta5::a[] = 
{        /*    0            1           2           3             4          5     6  */
/* 0 */         0.0,          0.0,         0.0,         0.0,          0.0,  0.0,  0.0, 
/* 1 */     1.0/6.0,          0.0,         0.0,         0.0,          0.0,  0.0,  0.0, 
/* 2 */     4.0/75.0,   16.0/75.0,         0.0,         0.0,          0.0,  0.0,  0.0, 
/* 3 */     5.0/6.0,    -8.0/3.0,    5.0/2.0,           0.0,          0.0,  0.0,  0.0,
/* 4 */    -8.0/5.0,   144.0/25.0,      -4.0,    16.0/25.0,          0.0,   0.0,  0.0, 
/* 5 */   361.0/320.0,  -18.0/5.0, 407.0/128.0, -11.0/80.0,   55.0/128.0,   0.0,  0.0, 
/* 6 */   -11.0/640.0,        0.0,  11.0/256.0, -11.0/160.0,  11.0/256.0,   0.0,  0.0,
/* 7 */    93.0/640.0,  -18.0/5.0, 803.0/256.0, -11.0/160.0,  99.0/256.0,   0.0,  1.0
};
// weights
var_t int_rungekutta5::bh[] = { 7.0/1408, 0.0, 1125.0/2816.0, 9.0/32.0, 125.0/768.0, 0.0, 5.0/66.0, 5.0/66.0 };
// nodes
var_t int_rungekutta5::c[]  = { 0.0, 1.0/6.0, 4.0/15.0, 2.0/3.0, 4.0/5.0, 1.0, 0.0, 1.0 };
#endif

/*
 * Dormand, J. R.; Prince, P. J.
 * "New Runge-Kutta algorithms for numerical simulation in dynamical astronomy"
 * Celestial Mechanics, vol. 18, Oct. 1978, p. 223-232.
 * p. 225 Table II. Runge-Kutta 5(4)T
 */
static const var_t lambda = 1.0/60.0;

// The Runge-Kutta matrix
var_t int_rungekutta5::a[] = 
{      /*        1              2                 3                 4              5       6      */
/* 1 */         0.0,           0.0,              0.0,              0.0,           0.0,    0.0,  // -> k1
/* 2 */   1.0/8.0,             0.0,              0.0,              0.0,           0.0,    0.0, 	// -> k2
/* 3 */         0.0,     1.0/4.0,                0.0,              0.0,           0.0,    0.0, 	// -> k3
/* 4 */ 196.0/729.0,  -320.0/729.0,    448.0/729.0,                0.0,           0.0,    0.0,	// -> k4
/* 5 */ 836.0/2875.0,   64.0/575.0, -13376.0/20125.0,  21384.0/20125.0,           0.0,    0.0, 	// -> k5
/* 6 */ -73.0/48.0,            0.0,   1312.0/231.0,    -2025.0/448.0,   2875.0/2112.0,    0.0, 	// -> k6
/*-------------------------------------------------------------------------------------------------------*/
/* 7 */  17.0/192.0,           0.0,     64.0/231.0,     2187.0/8960.0,  2875.0/8448.0, 1.0/20	// -> k7
}; /* 7 x 6 matrix */
static uint16_t a_row = 7;
static uint16_t a_col = 6;
// weights
var_t int_rungekutta5::bh[] = { 17.0/192.0, 0.0, 64.0/231.0, 2187.0/8960.0, 2875.0/8448.0, 1.0/20 };
// nodes
var_t int_rungekutta5::c[]  = { 0.0, 1.0/8.0, 1.0/4.0, 4.0/9.0, 4.0/5.0, 1.0, 1.0 };

// These arrays will contain the stepsize multiplied by the constants
var_t int_rungekutta5::h_a[ sizeof(int_rungekutta5::a ) / sizeof(var_t)];
var_t int_rungekutta5::h_bh[ sizeof(int_rungekutta5::bh ) / sizeof(var_t)];

//__constant__ var_t dc_a[sizeof(int_rungekutta5::a) / sizeof(var_t)];
//__constant__ var_t dc_bh[sizeof(int_rungekutta5::bh) / sizeof(var_t)];


int_rungekutta5::int_rungekutta5(ode& f, bool adaptive, var_t tolerance, comp_dev_t comp_dev) :
	integrator(f, adaptive, tolerance, (adaptive ? 7 : 6), comp_dev)
{
	name    = "Runge-Kutta5";
	n_order = 5;

	d_a  = NULL;
	d_bh = NULL;
	check_Butcher_tableau();
	if (PROC_UNIT_GPU == comp_dev.proc_unit)
	{
		allocate_Butcher_tableau();
	}	
}

int_rungekutta5::~int_rungekutta5()
{
	if (PROC_UNIT_GPU == comp_dev.proc_unit)
	{
		deallocate_Butcher_tableau();
	}	
}

void int_rungekutta5::allocate_Butcher_tableau()
{
	ALLOCATE_DEVICE_VECTOR((void**)&d_a,  sizeof(a));
	ALLOCATE_DEVICE_VECTOR((void**)&d_bh, sizeof(bh));
}

void int_rungekutta5::deallocate_Butcher_tableau()
{
	FREE_DEVICE_VECTOR((void**)&d_a);
	FREE_DEVICE_VECTOR((void**)&d_bh);
}

void int_rungekutta5::check_Butcher_tableau()
{
	uint16_t n_c = sizeof(int_rungekutta5::c) / sizeof(var_t);
	uint16_t n_col = (sizeof(int_rungekutta5::a) / sizeof(var_t)) / n_c;

	for (uint16_t i = 0; i < n_c; i++)
	{
		var_t sum = 0.0;
		for (uint16_t j = 0; j < n_col; j++)
		{
			uint16_t k = i * n_col + j;
			sum += a[k];
		}
		if (1.0e-15 < fabs(sum - c[i]))
		{
			throw std::string("The Runge-Kutta 5 is not consistent (sum(a_ij) != c_i).");
		}
	}
}

void int_rungekutta5::calc_ytemp(uint16_t stage)
{
	if (PROC_UNIT_GPU == comp_dev.proc_unit)
	{
		var_t* coeff = d_a + stage * a_col;
		gpu_calc_lin_comb_s(ytemp, f.y, d_k, coeff, stage, f.n_var, comp_dev.id_dev, optimize);
	}
	else
	{
		var_t* coeff = h_a + stage * a_col;
		tools::calc_lin_comb_s(ytemp, f.y, h_k, coeff, stage, f.n_var);
	}
}

void int_rungekutta5::calc_y_np1()
{
	if (PROC_UNIT_GPU == comp_dev.proc_unit)
	{
		var_t* coeff = d_bh;
		gpu_calc_lin_comb_s(f.yout, f.y, d_k, coeff, 6, f.n_var, comp_dev.id_dev, optimize);
	}
	else
	{
		var_t* coeff = h_bh;
		tools::calc_lin_comb_s(f.yout, f.y, h_k, coeff, 6, f.n_var);
	}
}

void int_rungekutta5::calc_error(uint32_t n)
{
	if (PROC_UNIT_GPU == comp_dev.proc_unit)
	{
        gpu_calc_rk5_error(err, k[5], k[6], n, comp_dev.id_dev, optimize);
	}
	else
	{
		for (uint32_t i = 0; i < n; i++)
		{
			h_err[i] = fabs(k[5][i] - k[6][i]);
		}
	}
}

var_t int_rungekutta5::step()
{
	static std::string err_msg1 = "The integrator could not provide the approximation of the solution with the specified tolerance.";

	static const uint16_t n_a = sizeof(int_rungekutta5::a) / sizeof(var_t);
	static const uint16_t n_bh = sizeof(int_rungekutta5::bh) / sizeof(var_t);
	static bool first_call = true;
	static uint32_t n_var = 0;

    if (n_var != f.n_var)
	{
		optimize = true;
		n_var = f.n_var;
	}
	else
	{
		optimize = false;
	}

	uint16_t stage = 0;
	t = f.t;
	//f.calc_dy(stage, t, f.y, k[0]); // -> k1

	// The final function evaluation at the nth step is the same as the first at the (n+1)th step,
	// thus the effective number of function evaluations per step is 6.
	if (!adaptive)
	{
		// Calculate initial differentials and store them into k1
		f.calc_dy(stage, t, f.y, k[0]); // -> k1
	}
	else
	{
		if (first_call)
		{
			first_call = false;
			// Calculate initial differentials and store them into k1
			f.calc_dy(stage, t, f.y, k[0]); // -> k1
		}
		else
		{
            if (PROC_UNIT_GPU == comp_dev.proc_unit)
            {
                CUDA_SAFE_CALL(cudaMemcpy(k[0], k[6], f.n_var*sizeof(var_t), cudaMemcpyDeviceToDevice));
            }
            else
            {
    			memcpy(k[0], k[6], f.n_var*sizeof(var_t));
            }
		}
	}

	var_t max_err = 0.0;
	uint16_t iter = 0;
	do
	{
		dt_did = dt_try;
		// Compute in advance the dt_try * coefficients to save n_var multiplication per stage
		for (uint16_t i = 0; i < n_a; i++)
		{
			h_a[i] = dt_try * a[i];
		}
		for (uint16_t i = 0; i < n_bh; i++)
		{
			h_bh[i] = dt_try * bh[i];
		}
	    if (PROC_UNIT_GPU == comp_dev.proc_unit)
	    {
			copy_vector_to_device(d_a,  h_a,  sizeof(h_a) );
			copy_vector_to_device(d_bh, h_bh, sizeof(h_bh));
	    }

		for (stage = 1; stage < 6; stage++)
		{
			t = f.t + c[stage] * dt_try;
			calc_ytemp(stage);
			f.calc_dy(stage, t, ytemp, k[stage]); // -> k2, k3, k4, k5, k6
		}
		// We have stage (6) number of k vectors, approximate the solution in f.yout using the bh coeff:
		calc_y_np1();   // -> f.yout = y = ynp1 = yn + 17/192*k1 + ... + 1/20*k6

		if (adaptive)
		{
			// Here stage = 6
			t = f.t + c[stage] * dt_try;
			f.calc_dy(stage, t, f.yout, k[stage]); // -> k7
			calc_error(f.n_var);
			max_err = get_max_error(f.n_var);
			max_err *= dt_try * lambda;
			calc_dt_try(max_err);
		}
		iter++;
	} while (adaptive && max_iter > iter && dt_min < dt_try && max_err > tolerance);

	if (max_iter <= iter)
	{
		throw std::string(err_msg1 + " The number of iteration exceeded the limit.");
	}
	if (dt_min >= dt_try)
	{
		throw std::string(err_msg1 + " The stepsize is smaller than the limit.");
	}

	t = f.tout = f.t + dt_did;
	f.swap();

	update_counters(iter);

	return dt_did;
}
