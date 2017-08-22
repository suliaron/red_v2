#pragma once

#include "ode.h"
#include "type.h"

class rtbp2D : public ode
{
public:
	//! Constructs a rtbp2D object
	/*!
        \param path_si   the path of the file which conatins data about the initial conditions
        \param path_sd   the path of the file which conatins the initial conditions
		\param n_ppo     the number of parameters per object
        \param omd_size  the size of the metadata in bytes
        \param PROC_UNIT the name of the executing device
	*/
	rtbp2D(std::string& path_si, std::string& path_sd, uint16_t n_ppo, size_t omd_size, comp_dev_t comp_dev);
	//! Destructor
	~rtbp2D();

	//! Print the solution (the numerical approximation of the solution)
	/*!
		\param path_si  full file path where the info about the solution will be printed
		\param path_sd  full file path where the the solution will be printed
		\param repres indicates the data representation of the file, i.e. text or binary
	*/
	void print_solution(std::string& path_si, std::string& path_sd, data_rep_t repres);
    void print_dump(std::string& path_si, std::string& path_sd);

	//! Print the integral(s) of the problem
	/*!
		\param path   full file path where the integrals of the problem will be printed
	*/
	void print_integral(std::string& path);

    void calc_dy(uint16_t stage, var_t curr_t, const var_t* y_temp, var_t* acc, var_t* jrk);
	void calc_dy(uint16_t stage, var_t curr_t, const var_t* y_temp, var_t* dy);
	void calc_integral();

private:
	void initialize();

	//! Load data about the initial conditions
	/*!
		\param path   full file path of the file
	*/
    void load_solution_info(std::string& path);

	//! Load the initial conditions from the hard disk
	/*!
		\param path   full file path of the file
	*/
	void load_solution_data(std::string& path);

	//! Load the initial conditions from the input stream (text mode)
	/*!
		\param input   the input stream associated with the file containing the initial conditions
	*/
	void load_ascii(std::ifstream& input);

	//! Load the initial conditions from the input stream (binary mode)
	/*!
		\param input   the input stream associated with the file containing the initial conditions
	*/
	void load_binary(std::ifstream& input);

	void cpu_calc_dy(uint16_t stage, var_t curr_t, const var_t* y_temp, var_t* dy);
	void gpu_calc_dy(uint16_t stage, var_t curr_t, const var_t* y_temp, var_t* dy);

	tbp_t::metadata_t *h_md, *d_md, *md;
};
