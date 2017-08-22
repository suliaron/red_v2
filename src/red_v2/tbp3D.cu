#include <iostream>
#include <iomanip>
#include <fstream>

#include "tbp3D.h"

#include "redutil2.h"
#include "constants.h"

using namespace std;
using namespace redutil2;


tbp3D::tbp3D(uint16_t n_ppo, size_t omd_size, comp_dev_t comp_dev) :
	ode(3, 1, 6, n_ppo, omd_size, comp_dev)
{
	name = "Singular 3D two-body problem";

	initialize();
}

tbp3D::~tbp3D()
{ }

void tbp3D::initialize()
{
    h_md = (tbp_t::metadata_t*)h_omd;
    d_md = (tbp_t::metadata_t*)d_omd;
      md = (tbp_t::metadata_t*)omd;
    
    h    = 0.0;
}

void tbp3D::calc_dy(uint16_t stage, var_t curr_t, const var_t* y_temp, var_t* acc, var_t* jrk)
{
	throw string("The tbp3D::calc_dy is not implemented.");
}

void tbp3D::calc_dy(uint16_t stage, var_t curr_t, const var_t* y_temp, var_t* dy)
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

void tbp3D::calc_integral()
{
	static bool first_call = true;
	const tbp_t::param_t* p = (tbp_t::param_t*)h_par;

	var_t r  = sqrt( SQR(h_y[0]) + SQR(h_y[1]) + SQR(h_y[2]) );
	var_t v2 = SQR(h_y[3]) + SQR(h_y[4]) + SQR(h_y[5]);

	h = 0.5 * v2 - p[0].mu / r;
	if (first_call)
	{
		integral.h0 = integral.h;
		first_call = false;
	}
}

void tbp3D::cpu_calc_dy(uint16_t stage, var_t curr_t, const var_t* y_temp, var_t* dy)
{
	const tbp_t::param_t* p = (tbp_t::param_t*)h_par;

	var_t r = sqrt( SQR(y_temp[0]) + SQR(y_temp[1]) + SQR(y_temp[2]) );
	var_t r3 = r*r*r;

	dy[0] = y_temp[3];                    // dx1 / dt = x4
	dy[1] = y_temp[4];                    // dx2 / dt = x5
	dy[2] = y_temp[5];                    // dx3 / dt = x6

	dy[3] = -p[0].mu / r3 * y_temp[0];    // dx4 / dt = -(mu/r^3) * x1
	dy[4] = -p[0].mu / r3 * y_temp[1];    // dx5 / dt = -(mu/r^3) * x2
	dy[5] = -p[0].mu / r3 * y_temp[2];    // dx6 / dt = -(mu/r^3) * x3
}

void tbp3D::gpu_calc_dy(uint16_t stage, var_t curr_t, const var_t* y_temp, var_t* dy)
{
	throw string("The gpu_calc_dy() is not implemented.");
}

void tbp3D::load(string& path)
{
	ifstream input;

	cout << "Loading " << path << " ";

	data_rep_t repres = file::get_data_repres(path);
	switch (repres)
	{
	case DATA_REPRESENTATION_ASCII:
		input.open(path.c_str());
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
		input.open(path.c_str(), ios::binary);
		if (input) 
		{
			load_binary(input);
		}
		else 
		{
			throw string("Cannot open " + path + ".");
		}
		break;
	}
	input.close();

	cout << " done" << endl;
}

void tbp3D::load_ascii(ifstream& input)
{
	tbp_t::param_t* p = (tbp_t::param_t*)h_par;

	var_t _t;
	for (uint32_t i = 0; i < n_obj; i++)
	{
		load_ascii_record(input, &_t, &h_md[i], &p[i], &h_y[i], &h_y[i+3]);
	}
}

void tbp3D::load_ascii_record(ifstream& input, var_t* t, tbp_t::metadata_t *md, tbp_t::param_t* p, var_t* r, var_t* v)
{
	string name;

	// epoch
	input >> *t;
	// id
	input >> md->id;
	// mu = k^2*(m1 + m2)
	input >> p->mu;

	// position
	var3_t* _r = (var3_t*)r;
	input >> _r->x >> _r->y >> _r->z;
	// velocity
	var3_t* _v = (var3_t*)v;
	input >> _v->x >> _v->y >> _v->z;
}

void tbp3D::load_binary(ifstream& input)
{
	throw string("The load_binary() is not implemented.");
}

void tbp3D::print_solution(std::string& path_si, std::string& path_sd, data_rep_t repres)
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

	if (sout)
	{
		switch (repres)
		{
		case DATA_REPRESENTATION_ASCII:
			print_solution_ascii(sout);
			break;
		case DATA_REPRESENTATION_BINARY:
			print_solution_binary(sout);
			break;
		default:
			throw string("Parameter 'repres' is out of range.");
		}
	}
	else
	{
		throw string("Cannot open " + path_si + ".");
	}
	sout.close();
}

void tbp3D::print_dump(std::string& path_si, std::string& path_sd)
{
    throw string("The tbp3D::print_dump function is not implemented.");
}

void tbp3D::print_solution_ascii(ofstream& sout)
{
	sout.precision(16);
	sout.setf(ios::right);
	sout.setf(ios::scientific);

	for (uint32_t i = 0; i < n_obj; i++)
    {
		sout << setw(VAR_T_W) << t << SEP
		// Print the metadata for each object
        << setw(INT_T_W) << h_md[i].id << SEP;

		// Print the parameters for each object
		for (uint16_t j = 0; j < n_ppo; j++)
		{
			uint32_t param_idx = i * n_ppo + j;
			sout << setw(VAR_T_W) << h_par[param_idx] << SEP;
		}
		// Print the variables for each object
		for (uint16_t j = 0; j < n_vpo; j++)
		{
			uint32_t var_idx = i * n_vpo + j;
			sout << setw(VAR_T_W) << h_y[var_idx];
			if (j < n_vpo - 1)
			{
				sout << SEP;
			}
			else
			{
				sout << endl;
			}
		}
	}
	sout.flush();
}

void tbp3D::print_solution_binary(ofstream& sout)
{
	throw string("The print_solution_binary() is not implemented.");
}

void tbp3D::print_integral(string& path)
{
	ofstream sout;

	sout.open(path.c_str(), ios::out | ios::app);
	if (sout)
	{
		sout.precision(16);
		sout.setf(ios::right);
		sout.setf(ios::scientific);

	sout << setw(VAR_T_W) << t << SEP                       /* time of the record [day] (double)           */
		 << h << endl;
	}
	else
	{
		throw string("Cannot open " + path + ".");
	}
	sout.close();
}
