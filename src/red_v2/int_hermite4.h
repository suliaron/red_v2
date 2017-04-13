#pragma once

#include "integrator.h"
#include "type.h"

class ode;

class int_hermite4 : public integrator
{
public:
    int_hermite4(ode& f, bool adaptive, var_t eta, comp_dev_t comp_dev);
    ~int_hermite4();

    var_t get_eta()       { return eta; }
    void set_eta(var_t e) { eta = e;    }

	void allocate_Butcher_tableau();
	void deallocate_Butcher_tableau();
	void check_Butcher_tableau();
	var_t step();

private:
    var_t eta;

    void predict();
    void correct();
};

class int_hermite4b : public integrator
{
public:
    int_hermite4b(ode& f, bool adaptive, var_t eta, comp_dev_t comp_dev);
    ~int_hermite4b();

    var_t get_eta()       { return eta; }
    void set_eta(var_t e) { eta = e;    }

	void allocate_Butcher_tableau();
	void deallocate_Butcher_tableau();
	void check_Butcher_tableau();
	var_t step();

private:
    var_t eta;

    void predict();
    void correct();
};
