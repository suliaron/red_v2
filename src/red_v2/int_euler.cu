#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include "ode.h"
#include "int_euler.h"

#include "redutil2.h"
#include "macro.h"

using namespace redutil2;

euler::euler(ode& f, comp_dev_t comp_dev) :
	integrator(f, false, 0.0, 1, comp_dev)
{
	name    = "Euler";
	n_order = 1;
}

euler::~euler()
{ }

void euler::allocate_Butcher_tableau()
{ }

void euler::deallocate_Butcher_tableau()
{ }

void euler::check_Butcher_tableau()
{ }

void euler::calc_y_np1()
{
	if (PROC_UNIT_GPU == comp_dev.proc_unit)
	{
		gpu_calc_lin_comb_s(f.yout, f.y, k[0], dt_try, f.n_var, comp_dev.id_dev, optimize);
	}
	else
	{
		tools::calc_lin_comb_s(f.yout, f.y, k[0], dt_try, f.n_var);
	}
}

var_t euler::step()
{
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
	// Calculate initial differentials and store them into k
	f.calc_dy(stage, t, f.y, k[stage]);

	calc_y_np1();

	dt_did = dt_try;
	f.tout = t = f.t + dt_did;
	f.swap();

	update_counters(1);

    return dt_did;
}
