#pragma once

#include <fstream>
#include <string>

#include "type.h"

typedef enum benchmark_output_name
{
    BENCHMARK_OUTPUT_NAME_RESULT,
    BENCHMARK_OUTPUT_NAME_N
} benchmark_output_name_t;

typedef enum job_name
{
    JOB_NAME_UNDEFINED,
    JOB_NAME_COMPARE,
    JOB_NAME_COMPARE_TEST_CASE,
    JOB_NAME_BENCMARK_CPU,
    JOB_NAME_BENCMARK_GPU,
    JOB_NAME_N
} job_name_t;

typedef struct option
{
    std::string     o_dir;
    std::string     base_fn;
    comp_dev_t comp_dev;
    bool       compare;
    var_t      tol;
    int        id_dev;
    uint32_t   n0;
    uint32_t   n1;
    uint32_t   dn;
    uint32_t   n_iter;
    uint32_t   n_si;
    uint32_t   n_nsi;
    uint32_t   n_ni;
    job_name_t job_name;
} option_t;

void create_filename(option_t& opt, std::string& filename);
void open_stream(std::string& o_dir, std::string& filename, std::ofstream** output, benchmark_output_name_t name);
void print(proc_unit_t proc_unit, std::string& method_name, std::string& param_name, uint2_t snk, uint2_t src, int n_body, int n_tpb, var_t Dt_CPU, var_t Dt_GPU, std::ofstream& sout, bool prn_to_scr);

void allocate_host_storage(uint32_t n_obj, var_t** h_y, var_t** h_dy, var_t** h_p, nbp_t::metadata_t** h_md);
void allocate_device_storage(uint32_t n_obj, var_t** d_y, var_t** d_a, var_t** d_p, nbp_t::metadata_t** d_md);

void deallocate_host_storage(var_t** h_y, var_t** h_dy, var_t** h_p, nbp_t::metadata_t** h_md);
void deallocate_device_storage(var_t** d_y, var_t** d_a, var_t** d_p, nbp_t::metadata_t** d_md);

void populate(uint32_t seed, uint32_t n_si, uint32_t n_nsi, uint32_t n_ni, var_t* h_y, var_t* h_p, nbp_t::metadata_t* h_md);
void populate(uint32_t seed, uint32_t n_obj, var_t* h_y, var_t* h_p, nbp_t::metadata_t* h_md);
void populate(uint32_t n_obj, var_t* h_y, var_t* h_p, nbp_t::metadata_t* h_md);

bool compare(uint32_t n, var_t tol, const var3_t* y1, const var3_t* y2);

int parse_options(int argc, const char **argv, option_t& opt, bool& verbose);
void create_default_option(option_t& opt);