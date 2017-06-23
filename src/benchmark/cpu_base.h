#pragma once
#include <type.h>

void cpu_calc_grav_accel(uint32_t n_obj, const nbp_t::metadata_t* md, const var3_t* r, const nbp_t::param_t* p, var3_t* a, bool use_sym);
void cpu_calc_grav_accel(uint2_t snk, uint2_t src, const nbp_t::metadata_t* md, const var3_t* r, const nbp_t::param_t* p, var3_t* a, bool use_sym);

void benchmark_CPU(uint32_t n_obj, const nbp_t::metadata_t* h_md, const var_t* h_y, const var_t* h_p, var_t* h_dy, std::ofstream& o_result);
