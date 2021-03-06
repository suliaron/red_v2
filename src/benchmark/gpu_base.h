#pragma once
#include "type.h"

float gpu_calc_grav_accel_naive(uint32_t n_obj, unsigned int n_tpb, cudaEvent_t& start, cudaEvent_t& stop, const nbp_t::metadata_t* d_md, const var3_t* r, const nbp_t::param_t* p, var3_t* a);
float gpu_calc_grav_accel_naive(uint2_t snk, uint2_t src, unsigned int n_tpb, cudaEvent_t& start, cudaEvent_t& stop, const nbp_t::metadata_t* d_md, const var3_t* r, const nbp_t::param_t* p, var3_t* a);
float gpu_calc_grav_accel_tile(uint32_t n_obj, unsigned int n_tpb, cudaEvent_t& start, cudaEvent_t& stop, const nbp_t::metadata_t* d_md, const var3_t* r, const nbp_t::param_t* p, var3_t* a);
float gpu_calc_grav_accel_tile(uint2_t snk, uint2_t src, unsigned int n_tpb, cudaEvent_t& start, cudaEvent_t& stop, const nbp_t::metadata_t* d_md, const var3_t* r, const nbp_t::param_t* p, var3_t* a);

void benchmark_GPU(int id_dev, uint32_t n_obj, const nbp_t::metadata_t* d_md, const var_t* d_y, const var_t* d_p, var_t* d_dy, std::ofstream& o_result);
