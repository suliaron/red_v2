#pragma once
#include <type.h>

void cpu_calc_grav_accel(var_t t, uint32_t n_obj, const var3_t* r, const nbp_t::param_t* p, var3_t* a, bool use_sym);
void cpu_calc_grav_accel(var_t t, uint32_t n_obj, uint2_t snk, uint2_t src, const var3_t* r, const nbp_t::param_t* p, var3_t* a, bool use_sym);

