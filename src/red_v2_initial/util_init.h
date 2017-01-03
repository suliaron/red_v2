#pragma once

#include <string>

#include "type.h"
#include "initial_type.h"

std::string create_name(int i, int type);

void generate_oe(oe_dist_t *oe_d, orbelem_t& oe);

void generate_pp(pp_dist_t *pp_d, nbp_t::param_t& param);
void generate_pp(pp_dist_t *pp_d, pp_disk_t::param_t& param);

//! Extract the orbital elements (a, e, i, w, O, M) from the data
/*
    \param data contains information got from the HORIZONS Web-Interface (http://ssd.jpl.nasa.gov/horizons.cgi#top)
    \param oe the orbital element structure will hold the data extracetd from data
    \return the epoch for which the orbital elements are valid
*/
var_t extract_from_horizon_output(std::string &data, orbelem_t& oe);

template <typename T>
void print_number(std::string& path, T number);
void print_oe(std::string &path, uint32_t n, var_t t, pp_disk_t::sim_data_t *sd);
void print_oe(std::string &path, uint32_t n, var_t t, nbp_t::metadata_t* md, nbp_t::param_t* p, orbelem_t* oe);
