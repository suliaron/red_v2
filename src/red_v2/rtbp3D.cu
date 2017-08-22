#include <iostream>
#include <iomanip>
#include <fstream>

#include "rtbp3D.h"

#include "redutil2.h"
#include "constants.h"

using namespace std;
using namespace redutil2;


rtbp3D::rtbp3D(uint16_t n_ppo, size_t omd_size, comp_dev_t comp_dev) :
	ode(3, 1, 9, n_ppo, omd_size, comp_dev)
{
	name = "Regularized 3D two-body problem";

	initialize();
}

rtbp3D::~rtbp3D()
{ }

void rtbp3D::initialize()
{
    h_md = (tbp_t::metadata_t*)h_omd;
    d_md = (tbp_t::metadata_t*)d_omd;
      md = (tbp_t::metadata_t*)omd;

	h       = 0.0;            // energy
	h_y[8]  = 0.0;            // s_0: fictitious time (4 position, 4 velocity, 1 time)
}

void rtbp3D::trans_to_descartes_var(var_t& x, var_t& y, var_t& z, var_t& vx, var_t& vy, var_t& vz)
{
	//var_t r = SQR(h_y[0]) + SQR(h_y[1]) + SQR(h_y[2]) + SQR(h_y[3]);							// r = u1^2 + u2^2 + u3^2 + u4^2

//	x  = SQR(h_y[0]) - SQR(h_y[1]) - SQR(h_y[2]) + SQR(h_y[3]);									// x = u1^2 - u2^2 - u3^2 + u4^2
//	y  = 2.0 * (h_y[0] * h_y[1] - h_y[2] * h_y[3]);												// y = 2*(u1*u2 - u3*u4)
//	z  = 2.0 * (h_y[0] * h_y[2] - h_y[1] * h_y[3]);												// z = 2*(u1*u3 - u2*u4)
//	vx = (2.0/r) * (h_y[0] * h_y[4] - h_y[1] * h_y[5] - h_y[2] * h_y[6] + h_y[3] * h_y[7]);		// vx = 2/r * (u1*vu1 - u2*vu2 - u3*vu3 + u4*vu4)
//	vy = (2.0/r) * (h_y[1] * h_y[4] + h_y[0] * h_y[5] - h_y[3] * h_y[6] - h_y[2] * h_y[7]);		// vy = 2/r * (u2*vu1 - u1*vu2 - u4*vu3 + u3*vu4)
//	vz = (2.0/r) * (h_y[2] * h_y[4] + h_y[3] * h_y[5] + h_y[0] * h_y[6] + h_y[1] * h_y[7]);		// vz = 2/r * (u3*vu1 - u4*vu2 - u1*vu3 + u2*vu4)
}

static void trans_to_descartes(const var4_t& u, const var4_t& u_prime, var3_t& r, var3_t& v)
{
	var_t d = SQR(u.x) + SQR(u.y) + SQR(u.z) + SQR(u.w);							// r = u1^2 + u2^2 + u3^2 + u4^2

	r.x  = SQR(u.x) - SQR(u.y) - SQR(u.z) + SQR(u.w);									// x = u1^2 - u2^2 - u3^2 + u4^2
//TODO
	//y  = 2.0 * (h_y[0] * h_y[1] - h_y[2] * h_y[3]);												// y = 2*(u1*u2 - u3*u4)
	//z  = 2.0 * (h_y[0] * h_y[2] - h_y[1] * h_y[3]);												// z = 2*(u1*u3 - u2*u4)
	//vx = (2.0/r) * (h_y[0] * h_y[4] - h_y[1] * h_y[5] - h_y[2] * h_y[6] + h_y[3] * h_y[7]);		// vx = 2/r * (u1*vu1 - u2*vu2 - u3*vu3 + u4*vu4)
	//vy = (2.0/r) * (h_y[1] * h_y[4] + h_y[0] * h_y[5] - h_y[3] * h_y[6] - h_y[2] * h_y[7]);		// vy = 2/r * (u2*vu1 - u1*vu2 - u4*vu3 + u3*vu4)
	//vz = (2.0/r) * (h_y[2] * h_y[4] + h_y[3] * h_y[5] + h_y[0] * h_y[6] + h_y[1] * h_y[7]);		// vz = 2/r * (u3*vu1 - u4*vu2 - u1*vu3 + u2*vu4)
}

//void rtbp3D::trans_to_regular_var(var_t& x, var_t& y, var_t& z, var_t& vx, var_t& vy, var_t& vz)
//{
//	var_t r = SQR(h_y[0]) + SQR(h_y[1]) + SQR(h_y[2]) + SQR(h_y[3]);							// r = u1^2 + u2^2 + u3^2 + u4^2
//
//	x  = SQR(h_y[0]) - SQR(h_y[1]) - SQR(h_y[2]) + SQR(h_y[3]);									// x = u1^2 - u2^2 - u3^2 + u4^2
//	y  = 2.0 * (h_y[0] * h_y[1] - h_y[2] * h_y[3]);												// y = 2*(u1*u2 - u3*u4)
//	z  = 2.0 * (h_y[0] * h_y[2] - h_y[1] * h_y[3]);												// z = 2*(u1*u3 - u2*u4)
//	vx = (2.0/r) * (h_y[0] * h_y[4] - h_y[1] * h_y[5] - h_y[2] * h_y[6] + h_y[3] * h_y[7]);		// vx = 2/r * (u1*vu1 - u2*vu2 - u3*vu3 + u4*vu4)
//	vy = (2.0/r) * (h_y[1] * h_y[4] + h_y[0] * h_y[5] - h_y[3] * h_y[6] - h_y[2] * h_y[7]);		// vy = 2/r * (u2*vu1 - u1*vu2 - u4*vu3 + u3*vu4)
//	vz = (2.0/r) * (h_y[2] * h_y[4] + h_y[3] * h_y[5] + h_y[0] * h_y[6] + h_y[1] * h_y[7]);		// vz = 2/r * (u3*vu1 - u4*vu2 - u1*vu3 + u2*vu4)
//}

void rtbp3D::calc_dy(uint16_t stage, var_t curr_t, const var_t* y_temp, var_t* acc, var_t* jrk)
{
	throw string("The rtbp3D::calc_dy is not implemented.");
}

void rtbp3D::calc_dy(uint16_t stage, var_t curr_t, const var_t* y_temp, var_t* dy)
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

void rtbp3D::calc_integral()
{
	static bool first_call = true;

	const tbp_t::param_t* p = (tbp_t::param_t*)h_par;

	var_t r  = SQR(h_y[0]) + SQR(h_y[1]) + SQR(h_y[2]) + SQR(h_y[3]);
	var_t vx = (2.0/r) * (h_y[0] * h_y[4] - h_y[1] * h_y[5] - h_y[2] * h_y[6] + h_y[3] * h_y[7]);		// vx = 2/r * (u1*vu1 - u2*vu2 - u3*vu3 + u4*vu4)
	var_t vy = (2.0/r) * (h_y[1] * h_y[4] + h_y[0] * h_y[5] - h_y[3] * h_y[6] - h_y[2] * h_y[7]);		// vy = 2/r * (u2*vu1 - u1*vu2 - u4*vu3 + u3*vu4)
	var_t vz = (2.0/r) * (h_y[2] * h_y[4] + h_y[3] * h_y[5] + h_y[0] * h_y[6] + h_y[1] * h_y[7]);		// vz = 2/r * (u3*vu1 - u4*vu2 - u1*vu3 + u2*vu4)
	var_t v2 = SQR(vx) + SQR(vy) + SQR(vz);
	h = 0.5 * v2 - p[0].mu / r;

	if (first_call)
	{
		integral.h0 = integral.h;
		first_call = false;
	}
}

void rtbp3D::cpu_calc_dy(uint16_t stage, var_t curr_t, const var_t* y_temp, var_t* dy)
{
	dy[0] = y_temp[4];																// dy1 / ds = y5
	dy[1] = y_temp[5];																// dy2 / ds = y6
	dy[2] = y_temp[6];																// dy3 / ds = y7
	dy[3] = y_temp[7];																// dy4 / ds = y8

	dy[4] = (h / 2.0) * y_temp[0];													// dy5 / ds = h/2 * y1
	dy[5] = (h / 2.0) * y_temp[1];													// dy6 / ds = h/2 * y2
	dy[6] = (h / 2.0) * y_temp[2];													// dy7 / ds = h/2 * y3
	dy[7] = (h / 2.0) * y_temp[3];													// dy8 / ds = h/2 * y4

	dy[8] = SQR(y_temp[0]) + SQR(y_temp[1])  + SQR(y_temp[2]) + SQR(y_temp[3]);     // dy9 / ds = y1^2 + y2^2 + y3^2 + y4^2
}

void rtbp3D::gpu_calc_dy(uint16_t stage, var_t curr_t, const var_t* y_temp, var_t* dy)
{
	throw string("The gpu_calc_dy() is not implemented.");
}

void rtbp3D::load(string& path)
{
	ifstream input;

	cout << "Loading " << path << " ";

	data_rep_t repres = (file::get_extension(path) == "txt" ? DATA_REPRESENTATION_ASCII : DATA_REPRESENTATION_BINARY);
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

void rtbp3D::load_ascii(ifstream& input)
{
	tbp_t::param_t* p = (tbp_t::param_t*)h_par;

	var_t _t;
	for (uint32_t i = 0; i < n_obj; i++)
	{
		load_ascii_record(input, &_t, &h_md[i], &p[i], &h_y[i], &h_y[i+4]);
	}
}

void rtbp3D::load_ascii_record(ifstream& input, var_t* t, tbp_t::metadata_t *md, tbp_t::param_t* p, var_t* r, var_t* v)
{
	string name;

	// epoch
	input >> *t;
	// id
	input >> md->id;
	// mu = k^2*(m1 + m2)
	input >> p->mu;

	// position
	var4_t* _r = (var4_t*)r;
	input >> _r->x >> _r->y >> _r->z >> _r->w;
	// velocity
	var4_t* _v = (var4_t*)v;
	input >> _v->x >> _v->y >> _v->z >> _v->w;

}

void rtbp3D::load_binary(ifstream& input)
{
	throw string("The load_binary() is not implemented.");
}

void rtbp3D::print_solution(std::string& path_si, std::string& path_sd, data_rep_t repres)
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

void rtbp3D::print_dump(std::string& path_si, std::string& path_sd)
{
    throw string("The rtbp3D::print_dump function is not implemented.");
}

void rtbp3D::print_solution_ascii(ofstream& sout) //TODO: implement correctly
{	sout.precision(16);
	sout.setf(ios::right);
	sout.setf(ios::scientific);

	var_t x = 0.0;
	var_t y = 0.0;
	var_t z = 0.0;
	var_t vx = 0.0;
	var_t vy = 0.0;
	var_t vz = 0.0;
	trans_to_descartes_var(x, y, z, vx, vy, vz);

	for (uint32_t i = 0; i < n_obj; i++)
    {
		sout << setw(VAR_T_W) << t << SEP                       /* time of the record [day] (double)           */
		// Print the metadata for each object
        << setw(INT_T_W) << h_md[i].id << SEP;

		// Print the parameters for each object
		for (uint16_t j = 0; j < n_ppo; j++)
		{
			uint32_t param_idx = i * n_ppo + j;
			sout << setw(VAR_T_W) << h_par[param_idx] << SEP;
		}
		// Print the regularized variables for each object
		for (uint16_t j = 0; j < n_vpo; j++)
		{
			uint32_t var_idx = i * n_vpo + j;
			sout << setw(VAR_T_W) << h_y[var_idx] << SEP;
		}
		// Print the descartes non-regularized variables for each object
		sout << setw(VAR_T_W) << x << SEP << y << SEP << z << SEP
			 << setw(VAR_T_W) << vx << SEP << vy << SEP << vz << endl;
	}
	sout.flush();
}

void rtbp3D::print_solution_binary(ofstream& sout)
{
	throw string("The print_solution_binary() is not implemented.");
}

void rtbp3D::print_integral(string& path)
{
	ofstream sout;

	sout.open(path.c_str(), ios::out | ios::app);
	if (sout)
	{
		sout.precision(16);
		sout.setf(ios::right);
		sout.setf(ios::scientific);

		sout << setw(VAR_T_W) << t << SEP                       /* fictitious time of the record (double)           */
			 << setw(VAR_T_W) << h_y[8] << SEP                  /* real time of the record [day] double             */
			 << h << endl;                                      /* energy of the system                             */
	}
	else
	{
		throw string("Cannot open " + path + ".");
	}
	sout.close();
}
