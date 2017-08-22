#pragma once

#include <string>

#include "type.h"
#include "initial_type.h"

void generate_oe(oe_dist_t *oe_d, orbelem_t& oe);

void generate_pp(pp_dist_t *pp_d, nbp_t::param_t& param);
void generate_pp(pp_dist_t *pp_d, pp_disk_t::param_t& param);

void create_star(uint32_t n_total, uint32_t& bodyIdx, uint32_t& bodyId, var_t m, var_t R, nbp_t::metadata_t* md, nbp_t::param_t* p, var_t* y);

//! Creates the disk with the specified distribution
/*
\param n_total the total number of bodies
\param n_disk the number of bodies in the current disk
\param bt the type of bodies in the current disk
\param bodyIdx the staring index of the bodies (will be updated)
\param bodyId the staring id of the bodies (will be updated)
\param oe_d the distribution of the orbital elements of the disk
\param pp_d the distribution of the physical properties of the bodies
\param md the metadata of the bodies will be stored here
\param p the parameters of the bodies will be stored here
\param oe the actual orbital elements will be stored here
\param y the coordinates and velocities will be stored here
\return the shortest orbital period of the disk
*/
double create_disk(uint32_t n_total, uint32_t n_disk, body_type_t bt, uint32_t& bodyIdx, uint32_t& bodyId, oe_dist_t& oe_d, pp_dist_t& pp_d, nbp_t::metadata_t* md, nbp_t::param_t* p, orbelem_t* oe, var_t* y);

//! Extract the orbital elements (a, e, i, w, O, M) from the data
/*
    \param data contains information got from the HORIZONS Web-Interface (http://ssd.jpl.nasa.gov/horizons.cgi#top)
    \param oe the orbital element structure will hold the data extracetd from data
    \return the epoch for which the orbital elements are valid
*/
var_t extract_from_horizon_output(std::string &data, orbelem_t& oe);

template <typename T>
void print_number(std::string& path, T number);
void print_oe(std::string &path, uint32_t n, var_t t, nbp_t::metadata_t* md, nbp_t::param_t* p, orbelem_t* oe);
