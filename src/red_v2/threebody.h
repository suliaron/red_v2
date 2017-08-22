#pragma once

#include "ode.h"

#include "type.h"

class threebody : public ode
{
public:
    //! Constructs a threebody object
    /*!
    \param path_si   the path of the file which conatins data about the initial conditions
    \param path_sd   the path of the file which conatins the initial conditions
    \param n_ppo     the number of parameters per object
    \param omd_size  the size of the metadata in bytes
    \param PROC_UNIT the name of the executing device
    */
    threebody(uint16_t n_ppo, size_t omd_size, comp_dev_t comp_dev);
	~threebody();

	void load(std::string& path);
	void load_ascii(std::ifstream& input);
	void load_ascii_record(std::ifstream& input, var_t* t, threebody_t::metadata_t *md, threebody_t::param_t* p);
	void load_binary(std::ifstream& input);

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

    void calc_dy(uint16_t stage, var_t curr_t, const var_t* y_temp, var_t* acc, var_t* jrk);
	void calc_dy(uint16_t stage, var_t curr_t, const var_t* y_temp, var_t* dy);
	void calc_integral();

//	void trans_to_descartes(var3_t& q1, var3_t& p1, var3_t& q2, var3_t& p2, var3_t& q3, var3_t& p3, const var4_t& Q1, const var4_t& P1, const var4_t& Q2, const var4_t& P2);
//	void trans_to_threebody(const var3_t& qv1, const var3_t& pv1, const var3_t& qv2, const var3_t& pv2, const var3_t& qv3, const var3_t& pv3, var4_t& Q1, var4_t& P1, var4_t& Q2, var4_t& P2);

//private:
	void initialize();

	void cpu_calc_dy(uint16_t stage, var_t curr_t, const var_t* y_temp, var_t* dy);
	void gpu_calc_dy(uint16_t stage, var_t curr_t, const var_t* y_temp, var_t* dy);

	var_t h;               //! Energy of the system

	threebody_t::metadata_t* h_md;
	threebody_t::metadata_t* d_md;
	threebody_t::metadata_t* md;

	var_t* h_epoch;
	var_t* d_epoch;
	var_t* epoch;

	bool first_open_solution;
	bool first_open_integral;
};