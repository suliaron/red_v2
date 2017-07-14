#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include "ode.h"
#include "int_rungekutta7.h"

#include "macro.h"
#include "redutil2.h"

using namespace redutil2;

static const var_t lambda = 41.0/840.0;

/*
 * Fehlberg, E.
 * "Classical Fifth-, Sixth-, Seventh-, and Eighth-Order Runge-Kutta Formulas with Stepsize Control"
 * NASA-TR-R-287 (https://nix.nasa.gov/search.jsp?R=19680027281&qs=N%253D4294957355)
 * p. 65 Table X. RK7(8)
 */
// The Runge-Kutta matrix
var_t int_rungekutta7::a[] = 
{    /*         1             2           3            4            5               6             7           8          9         10     11    12  */
/* 1 */        0.0,          0.0,        0.0,         0.0,         0.0,            0.0,          0.0,        0.0,       0.0,       0.0,   0.0,  0.0, 
/* 2 */     2.0/27.0,        0.0,        0.0,         0.0,         0.0,            0.0,          0.0,        0.0,       0.0,       0.0,   0.0,  0.0, 
/* 3 */     1.0/36.0,   1.0/12.0,        0.0,         0.0,         0.0,            0.0,          0.0,        0.0,       0.0,       0.0,   0.0,  0.0, 
/* 4 */     1.0/24.0,        0.0,   1.0/8.0,          0.0,         0.0,            0.0,          0.0,        0.0,       0.0,       0.0,   0.0,  0.0, 
/* 5 */     5.0/12.0,        0.0, -25.0/16.0,   25.0/16.0,         0.0,            0.0,          0.0,        0.0,       0.0,       0.0,   0.0,  0.0, 
/* 6 */     1.0/20.0,        0.0,        0.0,    1.0/4.0,      1.0/5.0,            0.0,          0.0,        0.0,       0.0,       0.0,   0.0,  0.0, 
/* 7 */   -25.0/108.0,       0.0,        0.0,  125.0/108.0,  -65.0/27.0,    125.0/54.0,          0.0,        0.0,       0.0,       0.0,   0.0,  0.0, 
/* 8 */    31.0/300.0,       0.0,        0.0,          0.0,   61.0/225.0,    -2.0/9.0,    13.0/900.0,        0.0,       0.0,       0.0,   0.0,  0.0,
/* 9 */           2.0,       0.0,        0.0,  -53.0/6.0,    704.0/45.0,   -107.0/9.0,    67.0/90.0,         3.0,       0.0,       0.0,   0.0,  0.0,
/*10 */   -91.0/108.0,       0.0,        0.0,   23.0/108.0, -976.0/135.0,   311.0/54.0,  -19.0/60.0,   17.0/6.0,  -1.0/12.0,       0.0,   0.0,  0.0,
/*11 */  2383.0/4100.0,      0.0,        0.0, -341.0/164.0, 4496.0/1025.0, -301.0/82.0, 2133.0/4100.0, 45.0/82.0, 45.0/164.0, 18.0/41.0,  0.0,  0.0,
/*----------------------------------------------------------------------------------------------------------------------------------------------------*/
/*12 */     3.0/205.0,       0.0,        0.0,          0.0,           0.0,   -6.0/41.0,   -3.0/205.0,  -3.0/41.0,  3.0/41.0,   6.0/41.0,  0.0,  0.0,
/*13 */ -1777.0/4100.0,      0.0,        0.0, -341.0/164.0, 4496.0/1025.0, -289.0/82.0, 2193.0/4100.0, 51.0/82.0, 33.0/164.0, 12.0/41.0,  0.0,  1.0
};
//static uint16_t a_row = 13;
static uint16_t a_col = 12;
// weights
var_t int_rungekutta7::bh[] = { 41.0/840.0, 0.0, 0.0, 0.0, 0.0, 34.0/105.0, 9.0/35.0, 9.0/35.0, 9.0/280.0, 9.0/280.0, 41.0/840.0 };
// nodes
var_t int_rungekutta7::c[]  = { 0.0, 2.0/27.0, 1.0/9.0, 1.0/6.0, 5.0/12.0, 1.0/2.0, 5.0/6.0, 1.0/6.0, 2.0/3.0, 1.0/3.0, 1.0, 0.0, 1.0 };

// These arrays will contain the stepsize multiplied by the constants
var_t int_rungekutta7::h_a[ sizeof(int_rungekutta7::a ) / sizeof(var_t)];
var_t int_rungekutta7::h_bh[ sizeof(int_rungekutta7::bh ) / sizeof(var_t)];


int_rungekutta7::int_rungekutta7(ode& f, bool adaptive, var_t tolerance, comp_dev_t comp_dev) :
	integrator(f, adaptive, tolerance, (adaptive ? 13 : 11), comp_dev)
{
	name    = "Runge-Kutta7";
	n_order = 7;

	d_a  = NULL;
	d_bh = NULL;
	check_Butcher_tableau();
	if (PROC_UNIT_GPU == comp_dev.proc_unit)
	{
		allocate_Butcher_tableau();
	}	
}

int_rungekutta7::~int_rungekutta7()
{
	if (PROC_UNIT_GPU == comp_dev.proc_unit)
	{
		deallocate_Butcher_tableau();
	}	
}

void int_rungekutta7::allocate_Butcher_tableau()
{
	ALLOCATE_DEVICE_VECTOR((void**)&d_a,  sizeof(a));
	ALLOCATE_DEVICE_VECTOR((void**)&d_bh, sizeof(bh));
}

void int_rungekutta7::deallocate_Butcher_tableau()
{
	FREE_DEVICE_VECTOR((void**)&d_a);
	FREE_DEVICE_VECTOR((void**)&d_bh);
}

void int_rungekutta7::check_Butcher_tableau()
{
	uint16_t n_c = sizeof(int_rungekutta7::c) / sizeof(var_t);
	uint16_t n_col = (sizeof(int_rungekutta7::a) / sizeof(var_t)) / n_c;

	for (uint16_t i = 0; i < n_c; i++)
	{
		var_t sum = 0.0;
		for (uint16_t j = 0; j < n_col; j++)
		{
			uint16_t k = i * n_col + j;
			sum += a[k];
		}
		if (1.0e-14 < fabs(sum - c[i]))
		{
			throw std::string("The Runge-Kutta 7 is not consistent (sum(a_ij) != c_i).");
		}
	}
}

void int_rungekutta7::calc_ytemp(uint16_t stage)
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

void int_rungekutta7::calc_y_np1()
{
	if (PROC_UNIT_GPU == comp_dev.proc_unit)
	{
		var_t* coeff = d_bh;
		gpu_calc_lin_comb_s(f.yout, f.y, d_k, coeff, 11, f.get_n_var(), comp_dev.id_dev);
	}
	else
	{
		var_t* coeff = h_bh;
		tools::calc_lin_comb_s(f.yout, f.y, h_k, coeff, 11, f.get_n_var());
	}
}

void int_rungekutta7::calc_error(uint32_t n)
{
	if (PROC_UNIT_GPU == comp_dev.proc_unit)
	{
        gpu_calc_rk7_error(err, k[0], k[10], k[11], k[12], n, comp_dev.id_dev);
	}
	else
	{
		for (uint32_t i = 0; i < n; i++)
		{
			h_err[i] = fabs(k[0][i] + k[10][i] - k[11][i] - k[12][i]);
		}
	}
}

var_t int_rungekutta7::step()
{
	static std::string err_msg1 = "The integrator could not provide the approximation of the solution with the specified tolerance.";

	static const uint16_t n_a = sizeof(int_rungekutta7::a) / sizeof(var_t);
	static const uint16_t n_bh = sizeof(int_rungekutta7::bh) / sizeof(var_t);

	uint16_t stage = 0;
	t = f.t;
#if 0
    std::string path("C:\\Work\\red.cuda.Results\\v2.0\\Test\\nbp\\TwoBody\\20170624\\gpu_rk7.txt");
    char buffer[512];
    sprintf(buffer, "stage = %2d, f.y:", stage);
    std::string comment(buffer);
    redutil2::print_array(path, comment, f.get_n_var(), f.y, (comp_dev.proc_unit == PROC_UNIT_CPU ? MEM_LOC_HOST : MEM_LOC_DEVICE));
#endif
    // Calculate initial differentials and store them into k
	f.calc_dy(stage, t, f.y, k[0]);  // -> k1

#if 0
    sprintf(buffer, "k[%2d]:", stage);
    comment.assign(buffer);
    redutil2::print_array(path, comment, f.get_n_var(), k[stage], (comp_dev.proc_unit == PROC_UNIT_CPU ?  MEM_LOC_HOST : MEM_LOC_DEVICE));
#endif

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

		for (stage = 1; stage < 11; stage++)
		{
			t = f.t + c[stage] * dt_try;
			calc_ytemp(stage);
#if 0
            sprintf(buffer, "stage = %2d, ytemp:", stage);
            comment.assign(buffer);
            redutil2::print_array(path, comment, f.get_n_var(), ytemp, (comp_dev.proc_unit == PROC_UNIT_CPU ? MEM_LOC_HOST : MEM_LOC_DEVICE));
#endif
            f.calc_dy(stage, t, ytemp, k[stage]); // -> k2, k3, ..., k11
#if 0
            sprintf(buffer, "k[%2d]:", stage);
            comment.assign(buffer);
            redutil2::print_array(path, comment, f.get_n_var(), k[stage], (comp_dev.proc_unit == PROC_UNIT_CPU ? MEM_LOC_HOST : MEM_LOC_DEVICE));
#endif
        }
		// We have stage (11) number of k vectors, approximate the solution in f.yout using the bh coeff:
		calc_y_np1();   // -> f.yout = y = ynp1 = yn + 41/840*k1 - 34/105*k6 ... + 41/840*k11
#if 0
        sprintf(buffer, "stage = %2d, f.yout:", stage);
        comment.assign(buffer);
        redutil2::print_array(path, comment, f.get_n_var(), f.yout, (comp_dev.proc_unit == PROC_UNIT_CPU ? MEM_LOC_HOST : MEM_LOC_DEVICE));
#endif

		if (adaptive)
		{
			// Here stage = 11
			for ( ; stage < n_stage; stage++)
			{
				t = f.t + c[stage] * dt_try;
				calc_ytemp(stage);
#if 0
                sprintf(buffer, "stage = %2d, ytemp:", stage);
                comment.assign(buffer);
                redutil2::print_array(path, comment, f.get_n_var(), ytemp, (comp_dev.proc_unit == PROC_UNIT_CPU ? MEM_LOC_HOST : MEM_LOC_DEVICE));
#endif
                f.calc_dy(stage, t, ytemp, k[stage]); // -> k12, k13
#if 0
                sprintf(buffer, "k[%2d]:", stage);
                comment.assign(buffer);
                redutil2::print_array(path, comment, f.get_n_var(), k[stage], (comp_dev.proc_unit == PROC_UNIT_CPU ? MEM_LOC_HOST : MEM_LOC_DEVICE));
#endif
            }
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
	if (dt_min >= dt_try)
	{
		throw std::string(err_msg1 + " The stepsize is smaller than the limit.");
	}

	update_counters(iter);

	t = f.t + dt_did;
	f.tout = t;
	f.swap();

	return dt_did;
}

#undef LAMBDA
