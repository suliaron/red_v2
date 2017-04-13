#pragma once

#include "integrator.h"
#include "type.h"

class ode;

class int_rungekutta4 : public integrator
{
public:
	static var_t a[],  h_a[];
	static var_t bh[], h_bh[];
	static var_t c[];

	int_rungekutta4(ode& f, bool adaptive, var_t tolerance, comp_dev_t comp_dev);
	~int_rungekutta4();

	void allocate_Butcher_tableau();
	void deallocate_Butcher_tableau();
	void check_Butcher_tableau();
	var_t step();

private:
	var_t *d_a;
	var_t *d_bh;

	void calc_ytemp(uint16_t stage);
	void calc_y_np1();

	void calc_error(uint32_t n);
};
