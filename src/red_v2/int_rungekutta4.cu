#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include "ode.h"
#include "int_rungekutta4.h"

#include "macro.h"
#include "redutil2.h"

using namespace redutil2;

static const var_t lambda = 1.0/10.0;

/*
 * Dormand, J. R.; Prince, P. J.
 * "New Runge-Kutta algorithms for numerical simulation in dynamical astronomy"
 * Celestial Mechanics, vol. 18, Oct. 1978, p. 223-232.
 * p. 225 Table I. Runge-Kutta 4(3)T
 */
// The Runge-Kutta matrix
var_t int_rungekutta4::a[] = 
{ 
	0.0,     0.0,     0.0,     0.0,     // y = yn                  -> k1
	1.0/2.0, 0.0,     0.0,     0.0,     // y = ytmp = yn + h/2*k1  -> k2
	0.0,     1.0/2.0, 0.0,     0.0,     // y = ytmp = yn + h/2*k2  -> k3
	0.0,     0.0,     1.0,     0.0,     // y = ytmp = yn + h*k3    -> k4
/*---------------------------------------------------------------------------------------*/
	1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0  // y = yn + h/6*(k1 + 2*k1 + 2*k3 + k4)    -> k5
}; /* 5 x 4 matrix */
//static uint16_t a_row = 5;
static uint16_t a_col = 4;

// weights
var_t int_rungekutta4::bh[] = { 1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0 };
// nodes
var_t int_rungekutta4::c[]  = { 0.0, 1.0/2.0, 1.0/2.0, 1.0, 1.0 };

// These arrays will contain the stepsize multiplied by the constants
var_t int_rungekutta4::h_a[sizeof(int_rungekutta4::a) / sizeof(var_t)];
var_t int_rungekutta4::h_bh[sizeof(int_rungekutta4::bh) / sizeof(var_t)];

//__constant__ var_t dc_a[sizeof(int_rungekutta4::a) / sizeof(var_t)];
//__constant__ var_t dc_bh[sizeof(int_rungekutta4::bh) / sizeof(var_t)];


int_rungekutta4::int_rungekutta4(ode& f, bool adaptive, var_t tolerance, comp_dev_t comp_dev) :
	integrator(f, adaptive, tolerance, (adaptive ? 5 : 4), comp_dev)
{
	name    = "Runge-Kutta4";
	n_order = 4;

	d_a  = NULL;
	d_bh = NULL;
	check_Butcher_tableau();
	if (PROC_UNIT_GPU == comp_dev.proc_unit)
	{
		allocate_Butcher_tableau();
	}	
}

int_rungekutta4::~int_rungekutta4()
{
	if (PROC_UNIT_GPU == comp_dev.proc_unit)
	{
		deallocate_Butcher_tableau();
	}	
}

void int_rungekutta4::allocate_Butcher_tableau()
{
	ALLOCATE_DEVICE_ARRAY((void**)&d_a,  sizeof(a));
	ALLOCATE_DEVICE_ARRAY((void**)&d_bh, sizeof(bh));
}

void int_rungekutta4::deallocate_Butcher_tableau()
{
	FREE_DEVICE_ARRAY((void**)&d_a);
	FREE_DEVICE_ARRAY((void**)&d_bh);
}

void int_rungekutta4::check_Butcher_tableau()
{
	uint16_t n_c = sizeof(int_rungekutta4::c) / sizeof(var_t);
	uint16_t n_col = (sizeof(int_rungekutta4::a) / sizeof(var_t)) / n_c;

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
			throw std::string("The Runge-Kutta 4 is not consistent (sum(a_ij) != c_i).");
		}
	}
}

void int_rungekutta4::calc_ytemp(uint16_t stage)
{
	if (PROC_UNIT_GPU == comp_dev.proc_unit)
	{
		var_t* coeff = d_a + stage * a_col;
		gpu_calc_lin_comb_s(ytemp, f.y, d_k, coeff, stage, f.get_n_var(), comp_dev.id_dev);
	}
	else
	{
		var_t* coeff = h_a + stage * a_col;
		tools::calc_lin_comb_s(ytemp, f.y, h_k, coeff, stage, f.get_n_var());
	}
}

void int_rungekutta4::calc_y_np1()
{
	if (PROC_UNIT_GPU == comp_dev.proc_unit)
	{
		var_t* coeff = d_bh;
		gpu_calc_lin_comb_s(f.yout, f.y, d_k, coeff, 4, f.get_n_var(), comp_dev.id_dev);
	}
	else
	{
		var_t* coeff = h_bh;
		tools::calc_lin_comb_s(f.yout, f.y, h_k, coeff, 4, f.get_n_var());
	}
}

void int_rungekutta4::calc_error(uint32_t n)
{
	if (PROC_UNIT_GPU == comp_dev.proc_unit)
	{
        gpu_calc_rk4_error(err, k[3], k[4], n, comp_dev.id_dev);
    }
	else
	{
		for (uint32_t i = 0; i < n; i++)
		{
			h_err[i] = fabs(h_k[3][i] - h_k[4][i]);
		}
	}
}

var_t int_rungekutta4::step()
{
	static std::string err_msg1 = "The integrator could not provide the approximation of the solution with the specified tolerance.";

	static const uint16_t n_a = sizeof(int_rungekutta4::a) / sizeof(var_t);
	static const uint16_t n_bh = sizeof(int_rungekutta4::bh) / sizeof(var_t);
	static bool first_call = true;

	uint16_t stage = 0;
	t = f.t;
	//f.calc_dy(stage, t, f.y, k[0]); // -> k1

    // The final function evaluation at the nth step is the same as the first at the (n+1)th step,
	// thus the effective number of function evaluations per step is 4.
    if (!adaptive)
    {
        // Calculate initial differentials and store them into k
        f.calc_dy(stage, t, f.y, k[0]); // -> k1
    }
    else
    {
        if (first_call)
        {
            first_call = false;
            // Calculate initial differentials and store them into k
            f.calc_dy(stage, t, f.y, k[0]); // -> k1
        }
        else
        {
            if (PROC_UNIT_GPU == comp_dev.proc_unit)
            {
                CUDA_SAFE_CALL(cudaMemcpy(k[0], k[4], f.get_n_var()*sizeof(var_t), cudaMemcpyDeviceToDevice));
            }
            else
            {
                memcpy(k[0], k[4], f.get_n_var()*sizeof(var_t));
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
			copy_array_to_device(d_a,  h_a,  sizeof(h_a) );
			copy_array_to_device(d_bh, h_bh, sizeof(h_bh));
	    }

		for (stage = 1; stage < 4; stage++) // stage = 1, 2, 3 
		{
			t = f.t + c[stage] * dt_try;          // -> tn + h2, tn + h/2, tn + h
			calc_ytemp(stage);                    // -> ytmp = yn + h/2*k1,  ytmp = yn + h/2*k2,  ytmp = yn + h*k3
			f.calc_dy(stage, t, ytemp, k[stage]); // -> k2, k3, k4
		}
		// We have stage (4) number of k vectors, approximate the solution in f.yout using the bh coeff:
		calc_y_np1();   // -> f.yout = y = ynp1 = yn + h/6*(k1 + 2*k2 + 2*k3 + k4)

		if (adaptive)
		{
			// Here stage = 4
			t = f.t + c[stage] * dt_try;
			f.calc_dy(stage, t, f.yout, k[stage]); // -> k5			
			calc_error(f.get_n_var());
			max_err = get_max_error(f.get_n_var());
			max_err *= dt_try * lambda;
			calc_dt_try(max_err);
		}
		iter++;
	} while (adaptive && max_iter > iter && dt_min < dt_try && max_err > tolerance);

	if (max_iter <= iter)
	{
		throw std::string(err_msg1 + " The number of iteration exceeded the limit.");
	}
	if (dt_min > dt_try)
	{
		throw std::string(err_msg1 + " The stepsize is smaller than the limit.");
	}

    t = f.tout = f.t + dt_did;
	f.swap();

	update_counters(iter);

    return dt_did;
}
