#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include "ode.h"
#include "int_hermite4.h"

#include "macro.h"
#include "redutil2.h"

using namespace redutil2;

/*
 * http://www.artcompsci.org/kali/vol/two_body_problem_2/ch11.html#rdocsect76
 */
int_hermite4::int_hermite4(ode& f, bool adaptive, var_t eta, comp_dev_t comp_dev) :
    integrator(f, adaptive, 0.0, 4, comp_dev)
{
    name    = "Hermite4";
    n_order = 4;

    this->eta = eta;
}

int_hermite4::~int_hermite4()
{ }

void int_hermite4::allocate_Butcher_tableau()
{ }

void int_hermite4::deallocate_Butcher_tableau()
{ }

void int_hermite4::check_Butcher_tableau()
{ }

void int_hermite4::predict()
{
	if (PROC_UNIT_GPU == comp_dev.proc_unit)
	{
	}
	else
	{
        // 2 FLOPS
        const var_t dt = dt_try;
        const var_t dt_2 = (1./2.) * dt_try;
        const var_t dt_3 = (1./3.) * dt_try;

        const var3_t* acc = (var3_t*)k[0];
        const var3_t* jrk = (var3_t*)k[1];
        const var3_t* r = (var3_t*)f.y;
        const var3_t* v = (var3_t*)(f.y + 3*f.n_obj);

        var3_t* pr = (var3_t*)ytemp;
        var3_t* pv = (var3_t*)(ytemp + 3*f.n_obj);

        /* Eq. (101) : I first determine trial values for the position and velocity,
         * simply by expanding both as a Taylor series, using only what we know at time,
         * which are the position, velocity, acceleration and jerk.
         */
        for (uint32_t i = 0; i < f.n_obj; i++)
        {
            pr[i].x = r[i].x + dt * (v[i].x + dt_2 * (acc[i].x + dt_3 * jrk[i].x));
            pr[i].y = r[i].y + dt * (v[i].y + dt_2 * (acc[i].y + dt_3 * jrk[i].y));
            pr[i].z = r[i].z + dt * (v[i].z + dt_2 * (acc[i].z + dt_3 * jrk[i].z));

            pv[i].x = v[i].x + dt * (acc[i].x + dt_2 * jrk[i].x);
            pv[i].y = v[i].y + dt * (acc[i].y + dt_2 * jrk[i].y);
            pv[i].z = v[i].z + dt * (acc[i].z + dt_2 * jrk[i].z);
        }
	}
}

void int_hermite4::correct()
{
    if (PROC_UNIT_GPU == comp_dev.proc_unit)
	{
	}
	else
	{
        // 2 FLOPS
        const var_t dt_2 = (1./2.) * dt_try;
        const var_t dt_6 = (1./6.) * dt_try;

        const var3_t* acc  = (var3_t*)k[0];
        const var3_t* jrk  = (var3_t*)k[1];
        const var3_t* pacc = (var3_t*)k[2];
        const var3_t* pjrk = (var3_t*)k[3];
        const var3_t* r = (var3_t*)f.y;
        const var3_t* v = (var3_t*)(f.y + 3*f.n_obj);

        var3_t* cr = (var3_t*)f.yout;
        var3_t* cv = (var3_t*)(f.yout + 3*f.n_obj);

        /* Eq. (103) : By switching the order of computation for the corrected 
         * form of the position and velocity, you are able to use the corrected
         * version of the velocity, rather than the predicted version, in determining 
         * the corrected version of the position.
         */
        for (uint32_t i = 0; i < f.n_obj; i++)
        {
            var3_t am = {acc[i].x - pacc[i].x, acc[i].y - pacc[i].y, acc[i].z - pacc[i].z};
            var3_t ap = {acc[i].x + pacc[i].x, acc[i].y + pacc[i].y, acc[i].z + pacc[i].z};
            var3_t jm = {jrk[i].x - pjrk[i].x, jrk[i].y - pjrk[i].y, jrk[i].z - pjrk[i].z};

            cv[i].x = v[i].x + dt_2 * (ap.x + dt_6 * jm.x);
            cv[i].y = v[i].y + dt_2 * (ap.y + dt_6 * jm.y);
            cv[i].z = v[i].z + dt_2 * (ap.z + dt_6 * jm.z);

            var3_t vp = {v[i].x + cv[i].x, v[i].y + cv[i].y, v[i].z + cv[i].z};

            cr[i].x = r[i].x + dt_2 * (vp.x + dt_6 * am.x);
            cr[i].y = r[i].y + dt_2 * (vp.y + dt_6 * am.y);
            cr[i].z = r[i].z + dt_2 * (vp.z + dt_6 * am.z);
        }
	}
}

var_t int_hermite4::step()
{
	t = f.t;
    dt_did = dt_try;

    f.calc_dy(0, t, f.y, k[0], k[1]);
    predict();

    f.calc_dy(0, t, ytemp, k[2], k[3]);
    correct();

    if (adaptive)
    {
        const var3_t* acc  = (var3_t*)k[2];
        const var3_t* jrk  = (var3_t*)k[3];

        var_t min_dt2 = 1.0e20;
        for (uint32_t i = 0; i < f.n_obj; i++)
        {
            var_t acc2 = SQR(acc[i].x) + SQR(acc[i].y) + SQR(acc[i].z);
            var_t jrk2 = SQR(jrk[i].x) + SQR(jrk[i].y) + SQR(jrk[i].z);
            var_t dt2 = acc2 / jrk2;
            if (min_dt2 > dt2)
            {
                min_dt2 = dt2;
            }
        }
        dt_try = eta * sqrt(min_dt2);
    }
    
    t = f.tout = f.t + dt_did;
	f.swap();

	update_counters(1);

    return dt_did;
}


int_hermite4b::int_hermite4b(ode& f, bool adaptive, var_t eta, comp_dev_t comp_dev) :
    integrator(f, adaptive, 0.0, 4, comp_dev)
{
    name    = "Hermite4b";
    n_order = 4;

    this->eta = eta;
}

int_hermite4b::~int_hermite4b()
{ }

void int_hermite4b::allocate_Butcher_tableau()
{ }

void int_hermite4b::deallocate_Butcher_tableau()
{ }

void int_hermite4b::check_Butcher_tableau()
{ }

void int_hermite4b::predict()
{
	if (PROC_UNIT_GPU == comp_dev.proc_unit)
	{
	}
	else
	{
        // 2 FLOPS
        const var_t dt = dt_try;
        const var_t dt_2 = (1./2.) * dt_try;
        const var_t dt_3 = (1./3.) * dt_try;

        const var3_t* acc = (var3_t*)k[0];
        const var3_t* jrk = (var3_t*)k[1];
        const var3_t* r = (var3_t*)f.y;
        const var3_t* v = (var3_t*)(f.y + 3*f.n_obj);

        var3_t* pr = (var3_t*)ytemp;
        var3_t* pv = (var3_t*)(ytemp + 3*f.n_obj);

        /* Eq. (101) : I first determine trial values for the position and velocity,
         * simply by expanding both as a Taylor series, using only what we know at time,
         * which are the position, velocity, acceleration and jerk.
         */
        for (uint32_t i = 0; i < f.n_obj; i++)
        {
            pr[i].x = r[i].x + dt * (v[i].x + dt_2 * (acc[i].x + dt_3 * jrk[i].x));
            pr[i].y = r[i].y + dt * (v[i].y + dt_2 * (acc[i].y + dt_3 * jrk[i].y));
            pr[i].z = r[i].z + dt * (v[i].z + dt_2 * (acc[i].z + dt_3 * jrk[i].z));

            pv[i].x = v[i].x + dt * (acc[i].x + dt_2 * jrk[i].x);
            pv[i].y = v[i].y + dt * (acc[i].y + dt_2 * jrk[i].y);
            pv[i].z = v[i].z + dt * (acc[i].z + dt_2 * jrk[i].z);
        }
	}
}

void int_hermite4b::correct()
{
    if (PROC_UNIT_GPU == comp_dev.proc_unit)
	{
	}
	else
	{
        // 2 FLOPS
        const var_t dt_2 = (1./2.) * dt_try;
        const var_t dt_6 = (1./6.) * dt_try;

        const var3_t* acc  = (var3_t*)k[0];
        const var3_t* jrk  = (var3_t*)k[1];
        const var3_t* pacc = (var3_t*)k[2];
        const var3_t* pjrk = (var3_t*)k[3];
        const var3_t* r = (var3_t*)f.y;
        const var3_t* v = (var3_t*)(f.y + 3*f.n_obj);

        var3_t* cr = (var3_t*)f.yout;
        var3_t* cv = (var3_t*)(f.yout + 3*f.n_obj);

        /* Eq. (103) : By switching the order of computation for the corrected 
         * form of the position and velocity, you are able to use the corrected
         * version of the velocity, rather than the predicted version, in determining 
         * the corrected version of the position.
         */
        for (uint32_t i = 0; i < f.n_obj; i++)
        {
            var3_t am = {acc[i].x - pacc[i].x, acc[i].y - pacc[i].y, acc[i].z - pacc[i].z};
            var3_t ap = {acc[i].x + pacc[i].x, acc[i].y + pacc[i].y, acc[i].z + pacc[i].z};
            var3_t jm = {jrk[i].x - pjrk[i].x, jrk[i].y - pjrk[i].y, jrk[i].z - pjrk[i].z};
            var3_t jp = {jrk[i].x + pjrk[i].x, jrk[i].y + pjrk[i].y, jrk[i].z + pjrk[i].z};

            cv[i].x = v[i].x + dt_2 * (ap.x - dt_6 * jm.x);
            cv[i].y = v[i].y + dt_2 * (ap.y - dt_6 * jm.y);
            cv[i].z = v[i].z + dt_2 * (ap.z - dt_6 * jm.z);

            var3_t vp = {v[i].x + cv[i].x, v[i].y + cv[i].y, v[i].z + cv[i].z};

            cr[i].x = r[i].x + dt_2 * (vp.x + dt_6 * am.x);
            cr[i].y = r[i].y + dt_2 * (vp.y + dt_6 * am.y);
            cr[i].z = r[i].z + dt_2 * (vp.z + dt_6 * am.z);
        }
	}
}

var_t int_hermite4b::step()
{
	t = f.t;
    dt_did = dt_try;

    f.calc_dy(0, t, f.y, k[0], k[1]);
    predict();

    f.calc_dy(0, t, ytemp, k[2], k[3]);
    correct();

    if (adaptive)
    {
        const var3_t* acc  = (var3_t*)k[2];
        const var3_t* jrk  = (var3_t*)k[3];

        var_t min_dt2 = 1.0e20;
        for (uint32_t i = 0; i < f.n_obj; i++)
        {
            var_t acc2 = SQR(acc[i].x) + SQR(acc[i].y) + SQR(acc[i].z);
            var_t jrk2 = SQR(jrk[i].x) + SQR(jrk[i].y) + SQR(jrk[i].z);
            var_t dt2 = acc2 / jrk2;
            if (min_dt2 > dt2)
            {
                min_dt2 = dt2;
            }
        }
        dt_try = eta * sqrt(min_dt2);
    }
    
    t = f.tout = f.t + dt_did;
	f.swap();

	update_counters(1);

    return dt_did;
}
