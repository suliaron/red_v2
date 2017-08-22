#include <iostream>
#include <iomanip>
#include <fstream>

#include "tbp2D.h"

#include "redutil2.h"
#include "constants.h"

using namespace std;
using namespace redutil2;


tbp2D::tbp2D(string& path_si, string& path_sd, uint16_t n_ppo, size_t omd_size, comp_dev_t comp_dev) :
	ode(2, 1, 4, n_ppo, omd_size, comp_dev)
{
	name = "Singular 2D two-body problem";

	initialize();

    load_solution_info(path_si);
    load_solution_data(path_sd);

    calc_integral();
    tout = t;
}

tbp2D::~tbp2D()
{ }

void tbp2D::initialize()
{
    h_md = (tbp_t::metadata_t*)h_omd;
    d_md = (tbp_t::metadata_t*)d_omd;
      md = (tbp_t::metadata_t*)omd;
}

void tbp2D::calc_dy(uint16_t stage, var_t curr_t, const var_t* y_temp, var_t* acc, var_t* jrk)
{
	throw string("The tbp2D::calc_dy is not implemented.");
}

void tbp2D::calc_dy(uint16_t stage, var_t curr_t, const var_t* y_temp, var_t* dy)
{
	if (PROC_UNIT_CPU == comp_dev.proc_unit)
	{
		cpu_calc_dy(stage, curr_t, y_temp, dy);
	}
	else
	{
		gpu_calc_dy(stage, curr_t, y_temp, dy);
	}
}

void tbp2D::calc_integral()
{
	static bool first_call = true;
	const tbp_t::param_t* p = (tbp_t::param_t*)h_par;

	var_t r  = sqrt(SQR(h_y[0]) + SQR(h_y[1]));
	var_t v2 = SQR(h_y[2]) + SQR(h_y[3]);
	integral.h = 0.5 * v2 - p[0].mu / r;
	if (first_call)
	{
		integral.h0 = integral.h;
		first_call = false;
	}
}

void tbp2D::cpu_calc_dy(uint16_t stage, var_t curr_t, const var_t* y_temp, var_t* dy)
{
	static const var_t mu = h_par[0];

	var_t r = sqrt(SQR(y_temp[0]) + SQR(y_temp[1]));
	var_t r3 = r*r*r;

	dy[0] = y_temp[2];                 // dx1 / dt = x3
	dy[1] = y_temp[3];                 // dx2 / dt = x4

	dy[2] = -(mu / r3) * y_temp[0];    // dx3 / dt = -mu / (r^3) * x1
	dy[3] = -(mu / r3) * y_temp[1];    // dx4 / dt = -mu / (r^3) * x2
}

void tbp2D::gpu_calc_dy(uint16_t stage, var_t curr_t, const var_t* y_temp, var_t* dy)
{
	throw string("The gpu_calc_dy() is not implemented.");
}

void tbp2D::load_solution_info(string& path)
{
	ifstream input;

	cout << "Loading " << path << " ";

	data_rep_t repres = file::get_data_repres(path);
	switch (repres)
	{
	case DATA_REPRESENTATION_ASCII:
		input.open(path.c_str(), ios::in);
		if (input) 
		{
			input >> t >> dt;
		}
		else 
		{
			throw string("Cannot open " + path + ".");
		}
		break;
	case DATA_REPRESENTATION_BINARY:
		input.open(path.c_str(), ios::in | ios::binary);
		if (input) 
		{
    		input.read((char*)&t, sizeof(var_t));
	        input.read((char*)&dt, sizeof(var_t));
		}
		else 
		{
			throw string("Cannot open " + path + ".");
		}
		break;
	default:
		throw string("Parameter 'repres' is out of range.");
	}
	input.close();

	cout << " done" << endl;
}

void tbp2D::load_solution_data(string& path)
{
	ifstream input;

	cout << "Loading " << path << " ";

	data_rep_t repres = file::get_data_repres(path);
	switch (repres)
	{
	case DATA_REPRESENTATION_ASCII:
		input.open(path.c_str(), ios::in);
		if (input) 
		{
			load_ascii(input);
		}
		else 
		{
			throw string("Cannot open " + path + ".");
		}
		break;
	case DATA_REPRESENTATION_BINARY:
		input.open(path.c_str(), ios::in | ios::binary);
		if (input) 
		{
			load_binary(input);
		}
		else 
		{
			throw string("Cannot open " + path + ".");
		}
		break;
	default:
		throw string("Parameter 'repres' is out of range.");
	}
	input.close();

	cout << " done" << endl;
}

void tbp2D::load_ascii(ifstream& input)
{
	// id
	input >> h_md[0].id;
	// mu = k^2*(m1 + m2)
	input >> h_par[0];
	// x - position
	input >> h_y[0];
	// y - position
	input >> h_y[1];
	// x - velocity
	input >> h_y[2];
	// y - velocity
	input >> h_y[3];
}

void tbp2D::load_binary(ifstream& input)
{
	throw string("The load_binary() is not implemented.");
}

void tbp2D::print_solution(std::string& path_si, std::string& path_sd, data_rep_t repres)
{
	ofstream sout;

	switch (repres)
	{
	case DATA_REPRESENTATION_ASCII:
		sout.open(path_si.c_str(), ios::out | ios::app);
		break;
	case DATA_REPRESENTATION_BINARY:
		sout.open(path_si.c_str(), ios::out | ios::app | ios::binary);
		break;
	default:
		throw string("Parameter 'repres' is out of range.");
	}
	if (!sout)
	{
		throw string("Cannot open " + path_si + ".");
	}
	file::tbp::print_solution_info(sout, t, dt, repres);
	sout.close();

	switch (repres)
	{
	case DATA_REPRESENTATION_ASCII:
		sout.open(path_sd.c_str(), ios::out | ios::app);
		break;
	case DATA_REPRESENTATION_BINARY:
		sout.open(path_sd.c_str(), ios::out | ios::app | ios::binary);
		break;
	default:
		throw string("Parameter 'repres' is out of range.");
	}
	if (!sout)
	{
		throw string("Cannot open " + path_sd + ".");
	}
	file::tbp::print_solution_data(sout, n_obj, n_ppo, n_vpo, h_md, h_par, h_y, repres);
	sout.close();
}

void tbp2D::print_dump(std::string& path_si, std::string& path_sd)
{
    throw string("The tbp2D::print_dump function is not implemented.");
}

void tbp2D::print_integral(string& path)
{
	ofstream sout;

	sout.open(path.c_str(), ios::out | ios::app);
	if (sout)
	{
		sout.precision(16);
		sout.setf(ios::right);
		sout.setf(ios::scientific);

	    sout << setw(VAR_T_W) << t << SEP             /* time of the record [day] (double)           */
		     << setw(VAR_T_W) << integral.h << endl;  /* energy of the system                        */
	}
	else
	{
		throw string("Cannot open " + path + ".");
	}
	sout.close();
}
