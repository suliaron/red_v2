#include <string>

#include "ode.h"
#include "redutil2.h"

using namespace redutil2;

ode::ode(uint16_t n_dim, uint32_t n_obj, uint16_t n_vpo, uint16_t n_ppo, size_t omd_size, comp_dev_t comp_dev) :
	n_dim(n_dim),
	n_obj(n_obj),
	n_vpo(n_vpo),
	n_ppo(n_ppo),
    omd_size(omd_size),
    comp_dev(comp_dev)
{
	initialize();

	n_var  = n_obj * n_vpo;
	n_par  = n_obj * n_ppo;
    event_size = 0;
	allocate_storage(n_var, n_par, omd_size, event_size);
	create_aliases();
}

ode::ode(uint16_t n_dim, uint32_t n_obj, uint16_t n_vpo, uint16_t n_ppo, size_t omd_size, size_t event_size, comp_dev_t comp_dev) :
    n_dim(n_dim),
    n_obj(n_obj),
    n_vpo(n_vpo),
    n_ppo(n_ppo),
    omd_size(omd_size),
    event_size(event_size),
    comp_dev(comp_dev)
{
    initialize();

    n_var = n_obj * n_vpo;
    n_par = n_obj * n_ppo;
    allocate_storage(n_var, n_par, omd_size, event_size);
    create_aliases();
}

ode::ode(uint16_t n_dim, uint32_t n_obj, uint16_t n_vpo, uint16_t n_ppo, uint32_t n_var, uint32_t n_par, size_t omd_size, comp_dev_t comp_dev) :
	n_dim(n_dim),
	n_obj(n_obj),
	n_vpo(n_vpo),
	n_ppo(n_ppo),
	n_var(n_var),
	n_par(n_par),
    omd_size(omd_size),
	comp_dev(comp_dev)
{
	initialize();

    event_size = 0;
    allocate_storage(n_var, n_par, omd_size, event_size);
	create_aliases();
}

ode::~ode()
{
	deallocate_storage();
}

void ode::initialize()
{
    t_wc     = 0.0;

    t        = 0.0;
	tout     = 0.0;
	dt       = 0.0;

	h_y      = NULL;
	h_yout   = NULL;
	d_y	     = NULL;
	d_yout   = NULL;
	y	     = NULL;
	yout     = NULL;
	h_par    = NULL;
	d_par    = NULL;
	par      = NULL;

    h_omd    = NULL;
    d_omd    = NULL;
    omd      = NULL;

    h_event  = NULL;
    d_event  = NULL;
    event    = NULL;

	var3_t zero = {0, 0, 0};
	integral.h0 = integral.h = 0;
	integral.c0 = integral.c = zero;
	integral.R0 = integral.R = zero;
	integral.V0 = integral.V = zero;
}

void ode::set_n_obj(uint32_t n)
{
    n_obj = n;
    n_var = n_obj * n_vpo;
    n_par = n_obj * n_ppo;
}

void ode::allocate_storage(uint32_t n_var, uint32_t n_par, size_t omd_size, size_t event_size)
{
	allocate_host_storage(n_var, n_par, omd_size, event_size);
	if (PROC_UNIT_GPU == comp_dev.proc_unit)
	{
		allocate_device_storage(n_var, n_par, omd_size, event_size);
	}
}

void ode::allocate_host_storage(uint32_t n_var, uint32_t n_par, size_t omd_size, size_t event_size)
{
	ALLOCATE_HOST_ARRAY((void **)&(h_y),     n_var * sizeof(var_t));
	ALLOCATE_HOST_ARRAY((void **)&(h_yout),  n_var * sizeof(var_t));
	ALLOCATE_HOST_ARRAY((void **)&(h_par),   n_par * sizeof(var_t));
    ALLOCATE_HOST_ARRAY((void **)&(h_omd),   omd_size);
    ALLOCATE_HOST_ARRAY((void **)&(h_event), event_size);
}

void ode::allocate_device_storage(uint32_t n_var, uint32_t n_par, size_t omd_size, size_t event_size)
{
	ALLOCATE_DEVICE_ARRAY((void **)&(d_y),     n_var * sizeof(var_t));
	ALLOCATE_DEVICE_ARRAY((void **)&(d_yout),  n_var * sizeof(var_t));
	ALLOCATE_DEVICE_ARRAY((void **)&(d_par),   n_par * sizeof(var_t));
    ALLOCATE_DEVICE_ARRAY((void **)&(d_omd),   omd_size);
    ALLOCATE_DEVICE_ARRAY((void **)&(d_event), event_size);
}

void ode::deallocate_storage()
{
	//NOTE : First always release the DEVICE memory
	if (PROC_UNIT_GPU == comp_dev.proc_unit)
	{
		deallocate_device_storage();
	}
	deallocate_host_storage();
}

void ode::deallocate_host_storage()
{
	FREE_HOST_ARRAY((void **)&(h_y));
	FREE_HOST_ARRAY((void **)&(h_yout));
	FREE_HOST_ARRAY((void **)&(h_par));
    FREE_HOST_ARRAY((void **)&(h_omd));
    FREE_HOST_ARRAY((void **)&(h_event));
}

void ode::deallocate_device_storage()
{
	FREE_DEVICE_ARRAY((void **)&(d_y));
	FREE_DEVICE_ARRAY((void **)&(d_yout));
	FREE_DEVICE_ARRAY((void **)&(d_par));
    FREE_DEVICE_ARRAY((void **)&(d_omd));
    FREE_DEVICE_ARRAY((void **)&(d_event));
}

void ode::copy_vars(copy_direction_t dir)
{
	switch (dir)
	{
	case COPY_DIRECTION_TO_DEVICE:
		copy_array_to_device(d_y, h_y, n_var*sizeof(var_t));
		break;
	case COPY_DIRECTION_TO_HOST:
		copy_array_to_host(h_y, d_y, n_var*sizeof(var_t));
		break;
	default:
		throw std::string("Parameter 'dir' is out of range.");
	}
}

void ode::copy_params(copy_direction_t dir)
{
    switch (dir)
    {
    case COPY_DIRECTION_TO_DEVICE:
        copy_array_to_device(d_par, h_par, n_par * sizeof(var_t));
        break;
    case COPY_DIRECTION_TO_HOST:
        copy_array_to_host(h_par, d_par, n_par * sizeof(var_t));
        break;
    default:
        throw std::string("Parameter 'dir' is out of range.");
    }
}

void ode::copy_metadata(copy_direction_t dir)
{
    switch (dir)
    {
    case COPY_DIRECTION_TO_DEVICE:
        copy_array_to_device(d_omd, h_omd, omd_size);
        break;
    case COPY_DIRECTION_TO_HOST:
        copy_array_to_host(h_omd, d_omd, omd_size);
        break;
    default:
        throw std::string("Parameter 'dir' is out of range.");
    }
}

void ode::copy_eventdata(copy_direction_t dir)
{
    switch (dir)
    {
    case COPY_DIRECTION_TO_DEVICE:
        copy_array_to_device(d_event, h_event, event_size);
        break;
    case COPY_DIRECTION_TO_HOST:
        copy_array_to_host(h_event, d_event, event_size);
        break;
    default:
        throw std::string("Parameter 'dir' is out of range.");
    }
}

void ode::swap()
{
	std::swap(t, tout);
	std::swap(y, yout);
	std::swap(h_y, h_yout);
	std::swap(d_y, d_yout);
}

// Date of creation: 2016.08.03.
// Last edited: 
// Status: Not tested
void ode::create_aliases()
{
    switch (comp_dev.proc_unit)
    {
    case PROC_UNIT_CPU:
        y     = h_y;
        yout  = h_yout;
        par   = h_par;
        omd   = h_omd;
        event = h_event;
        break;
    case PROC_UNIT_GPU:
        y     = d_y;
        yout  = d_yout;
        par   = d_par;
        omd   = d_omd;
        event = d_event;
        break;
    default:
        throw std::string("Parameter 'proc_unit' is out of range.");
    }
}

void ode::set_comp_dev(comp_dev_t cd)
{
    // If the execution is already on the requested device than nothing to do
    if (comp_dev.proc_unit == cd.proc_unit)
    {
        return;
    }

    switch (cd.proc_unit)
    {
    case PROC_UNIT_CPU:
        copy_params(COPY_DIRECTION_TO_HOST);
        copy_vars(COPY_DIRECTION_TO_HOST);
        copy_metadata(COPY_DIRECTION_TO_HOST);
        copy_eventdata(COPY_DIRECTION_TO_HOST);
        deallocate_device_storage();
        break;
    case PROC_UNIT_GPU:
        allocate_device_storage(n_var, n_par, omd_size, event_size);
        copy_params(COPY_DIRECTION_TO_DEVICE);
        copy_vars(COPY_DIRECTION_TO_DEVICE);
        copy_metadata(COPY_DIRECTION_TO_DEVICE);
        copy_eventdata(COPY_DIRECTION_TO_DEVICE);
        break;
    default:
        throw std::string("Parameter 'cd.proc_unit' is out of range.");
    }

    comp_dev = cd;
    create_aliases();
}
