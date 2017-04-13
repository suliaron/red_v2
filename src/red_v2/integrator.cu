#include "thrust\device_ptr.h"
#include "thrust\fill.h"
#include "thrust\extrema.h"

#include "integrator.h"
#include "ode.h"
#include "macro.h"

#include "redutil2.h"

using namespace std;
using namespace redutil2;

integrator::integrator(ode& f, bool adaptive, var_t tolerance, uint16_t n_stage, comp_dev_t comp_dev) : 
	f(f),
	adaptive(adaptive),
	tolerance(tolerance),
	n_stage(n_stage),
	comp_dev(comp_dev)
{
	initialize();

	allocate_storage(f.n_var);
	create_aliases();
}

integrator::~integrator()
{
	deallocate_storage();
}

void integrator::initialize()
{
    optimize      = true;

    t             = f.t;
	dt_try        = f.dt;
	dt_did        = 0.0;

	h_k           = NULL;
	d_k           = NULL;
	k             = NULL;
	cpy_dk        = NULL;

	h_ytemp       = NULL;
	d_ytemp       = NULL;
	ytemp         = NULL;

	h_err         = NULL;
	d_err         = NULL;
	err           = NULL;

	max_iter      = 100;
	dt_min        = 1.0e10;

	n_tried_step  = 0;
	n_passed_step = 0;
	n_failed_step = 0;
}

void integrator::allocate_storage(uint32_t n_var)
{
	allocate_host_storage(n_var);
	if (PROC_UNIT_GPU == comp_dev.proc_unit)
	{
		allocate_device_storage(n_var);
	}
}

void integrator::allocate_host_storage(uint32_t n_var)
{
	ALLOCATE_HOST_VECTOR((void**)&h_k,    n_stage*sizeof(var_t*));
	ALLOCATE_HOST_VECTOR((void**)&cpy_dk, n_stage*sizeof(var_t*));

	for (uint16_t i = 0; i < n_stage; i++)
	{
		ALLOCATE_HOST_VECTOR((void**)(h_k + i), n_var*sizeof(var_t));
	}

	ALLOCATE_HOST_VECTOR((void**)&(h_ytemp), n_var*sizeof(var_t));
	if (adaptive)
	{
		ALLOCATE_HOST_VECTOR((void**)&(h_err), n_var*sizeof(var_t));
	}
}

void integrator::allocate_device_storage(uint32_t n_var)
{
	ALLOCATE_DEVICE_VECTOR((void**)(&d_k), n_stage*sizeof(var_t*));
	for (uint16_t i = 0; i < n_stage; i++)
	{
		ALLOCATE_DEVICE_VECTOR((void**)(cpy_dk + i), n_var*sizeof(var_t));
	}
	CUDA_SAFE_CALL(cudaMemcpy(d_k, cpy_dk, n_stage * sizeof(var_t*), cudaMemcpyHostToDevice));

	ALLOCATE_DEVICE_VECTOR((void**)&(d_ytemp), n_var*sizeof(var_t));
	if (adaptive)
	{
		ALLOCATE_DEVICE_VECTOR((void**)&(d_err), n_var*sizeof(var_t));
	}
}

void integrator::deallocate_storage()
{
	//NOTE : First always release the DEVICE memory
	if (PROC_UNIT_GPU == comp_dev.proc_unit)
	{
		deallocate_device_storage();
	}
	deallocate_host_storage();
}

void integrator::deallocate_host_storage()
{
	for (uint16_t i = 0; i < n_stage; i++)
	{
		FREE_HOST_VECTOR((void **)(h_k + i));
	}
	//free(k);           k = NULL;
	free(h_k);       h_k = NULL;
	free(cpy_dk); cpy_dk = NULL;

	FREE_HOST_VECTOR((void **)&(h_ytemp));
	if (adaptive)
	{
		FREE_HOST_VECTOR((void **)&(h_err));
	}
}

void integrator::deallocate_device_storage()
{
	for (uint16_t i = 0; i < n_stage; i++)
	{
		FREE_DEVICE_VECTOR((void **)&(cpy_dk[i]));
	}
	FREE_DEVICE_VECTOR((void **)&(d_k));

	FREE_DEVICE_VECTOR((void **)&(d_ytemp));
	if (adaptive)
	{
		FREE_DEVICE_VECTOR((void **)&(d_err));
	}
}

// Date of creation: 2016.08.02.
// Last edited: 
// Status: Not tested
void integrator::create_aliases()
{
// TODO: Ezt kell tesztelni!!
	switch (comp_dev.proc_unit)
	{
	case PROC_UNIT_CPU:
		ytemp = h_ytemp;
		k = h_k;
		//memcpy(k, h_k, n_stage * sizeof(var_t*));
		if (adaptive)
		{
			err = h_err;
		}
		break;
	case PROC_UNIT_GPU:
		ytemp = d_ytemp;
		k = cpy_dk;
		//memcpy(k, cpy_dk, n_stage * sizeof(var_t*));
		if (adaptive)
		{
			err = d_err;
		}
		break;
	default:
		throw string("Parameter 'PROC_UNIT' is out of range.");
	}
}

void integrator::set_computing_device(comp_dev_t comp_dev)
{
	// If the execution is already on the requested device than nothing to do
	if (this->comp_dev.proc_unit == comp_dev.proc_unit)
	{
		return;
	}
	// TODO: implement

	//int n_body = ppd->n_bodies->get_n_total_playing();

	//switch (device)
	//{
	//case PROC_UNIT_CPU:
	//	deallocate_device_storage();
	//	break;
	//case PROC_UNIT_GPU:
	//	allocate_device_storage(n_body);
	//	break;
	//default:
	//	throw string("Parameter 'device' is out of range.");
	//}

	//this->PROC_UNIT = device;
	//create_aliases();
	//f->set_computing_device(device);
}


var_t integrator::get_max_error(uint32_t n_var)
{
	var_t max_err = 0.0;

	if (PROC_UNIT_GPU == comp_dev.proc_unit)
	{
		// Wrap raw pointer with a device_ptr
		thrust::device_ptr<var_t> d_ptr(d_err);
		// Use thrust to find the maximum element
		thrust::device_ptr<var_t> d_ptr_max = thrust::max_element(d_ptr, d_ptr + n_var);
		// Copy the max element from device memory to host memory
		cudaMemcpy((void*)&max_err, (void*)d_ptr_max.get(), sizeof(var_t), cudaMemcpyDeviceToHost);
	}
	else
	{
		for (uint32_t i = 0; i < n_var; i++)
		{
			if (max_err < fabs(h_err[i]))
			{
				max_err = fabs(h_err[i]);
			}
		}		
	}
	return (max_err);
}

void integrator::calc_dt_try(var_t max_err)
{
	if (1.0e-20 < max_err)
	{
		dt_try *= 0.9 * pow(tolerance / max_err, 1.0/(n_order));
	}
	else
	{
		dt_try *= 5.0;
	}
}

void integrator::update_counters(uint16_t iter)
{
	n_tried_step  += iter;
	n_failed_step += (iter - 1);
	n_passed_step++;
}
