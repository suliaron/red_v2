#pragma once

#include "integrator.h"

#include "type.h"

class ode;

class euler : public integrator
{
public:
	euler(ode& f, comp_dev_t comp_dev);
	~euler();

	void allocate_Butcher_tableau();
	void deallocate_Butcher_tableau();
	void check_Butcher_tableau();
	var_t step();

private:
	void calc_y_np1();
};
