#pragma once

#include <string>

#include "parameter.h"

#include "type.h"

class integrator;
class ode;

class options
{
public:
	options(int argc, const char** argv);
	~options();

	ode* create_model();
	integrator* create_integrator(ode& f);
	
	dyn_model_t dyn_model;

	bool verbose;                   //!< Print every event to the log file
	bool print_to_screen;           //!< Print every event to the standard output stream (cout) 
	bool append;                    //!< If true than only one file will hold the solution
	bool ef;                        //!< Extend the file names with command line information. Only for developer and debugger purposes.
	bool benchmark;                 //!< Run benchmark to measure the time it takes to compute the right-hand-side

	uint32_t n_change_to_cpu;       //!< The threshold value for the total number of SI bodies to change to the CPU

	comp_dev_t comp_dev;            //!< The computing device to carry out the calculations (cpu or gpu)
	gas_disk_model_t g_disk_model;

	parameter* param;

	std::string fn_prefix;
	std::string out_fn[OUTPUT_NAME_N];   //!< Array for the output filenames
	std::string in_fn[INPUT_NAME_N];     //!< Array for the input filenames
	std::string dir[DIRECTORY_NAME_N];   //!< Array for the input and output directories

private:
	void create_default();
	void parse(int argc, const char** argv);
    void get_solution_path(std::string& path_si, std::string &path_sd);

	void print_usage();
};
