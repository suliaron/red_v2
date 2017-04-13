#pragma once

#include "ode.h"

#include "type.h"

class rtbp3D : public ode
{
public:
	rtbp3D(uint16_t n_ppo, comp_dev_t comp_dev);
	rtbp3D(uint16_t n_ppo, var_t t, tbp_t::metadata_t *md, tbp_t::param_t* p, var_t* r, var_t* v, comp_dev_t comp_dev);
	~rtbp3D();

	void load(std::string& path);
	void load_ascii(std::ifstream& input);
	void load_ascii_record(std::ifstream& input, var_t* _t, tbp_t::metadata_t *md, tbp_t::param_t* p, var_t* r, var_t* v);
	void load_binary(std::ifstream& input);

	//! Copies N-body metadata between HOST and DEVICE memory
	/*!
		\param dir The direction of the copy
	*/
	void copy_metadata(copy_direction_t dir);

	//! Print the solution (the numerical approximation of the solution)
	/*!
		\param path_si  full file path where the info about the solution will be printed
		\param path_sd  full file path where the the solution will be printed
		\param repres indicates the data representation of the file, i.e. text or binary
	*/
	void print_solution(std::string& path_si, std::string& path_sd, data_rep_t repres);
    void print_dump(std::string& path_si, std::string& path_sd);
    //! Print the solution for each object in text format
	/*   
		\param sout print the data to this stream
	*/
	void print_solution_ascii(std::ofstream& sout);
	//! Print the solution for each object in binary format
	/*!
		\param sout print the data to this stream
	*/
	void print_solution_binary(std::ofstream& sout);

	//! Print the integral(s) of the problem
	/*!
		\param path   full file path where the integrals of the problem will be printed
	*/
	void print_integral(std::string& path);

	static void trans_to_descartes(const var4_t& u, const var4_t& u_prime, var3_t& r, var3_t& v);
	static void trans_to_parameter(const var3_t& r, const var3_t& v, var4_t& u, var4_t& u_prime);
	void trans_to_descartes_var(var_t& x, var_t& y, var_t& z, var_t& vx, var_t& vy, var_t& vz);

    void calc_dy(uint16_t stage, var_t curr_t, const var_t* y_temp, var_t* acc, var_t* jrk);
	void calc_dy(uint16_t stage, var_t curr_t, const var_t* y_temp, var_t* dy);
	void calc_integral();

//private:
	void initialize();
	void allocate_storage();
	void allocate_host_storage();
	void allocate_device_storage();
	
	void deallocate_storage();
	void deallocate_host_storage();
	void deallocate_device_storage();

	void cpu_calc_dy(uint16_t stage, var_t curr_t, const var_t* y_temp, var_t* dy);
	void gpu_calc_dy(uint16_t stage, var_t curr_t, const var_t* y_temp, var_t* dy);

	var_t h;               //! Energy of the system

	tbp_t::metadata_t *h_md, *d_md, *md;
};