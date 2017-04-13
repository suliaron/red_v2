#pragma once

#include "integrator.h"

#include "type.h"

class ode;

class int_rungekutta2 : public integrator
{
public:
	static var_t a[], h_a[];
	static var_t b[];
	static var_t c[];

	int_rungekutta2(ode& f, comp_dev_t comp_dev);
	~int_rungekutta2();

	void allocate_Butcher_tableau();
	void deallocate_Butcher_tableau();
	void check_Butcher_tableau();
	var_t step();

private:
	var_t *d_a;

	void calc_ytemp(uint16_t stage);
	void calc_y_np1();
};
