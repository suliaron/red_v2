#include <iostream>
#include <iomanip>
#include <fstream>

#include "threebody.h"

#include "redutil2.h"
#include "constants.h"

#include "tools.h"

using namespace std;
using namespace redutil2;


threebody::threebody(uint16_t n_ppo, size_t omd_size, comp_dev_t comp_dev) :
	ode(3, 3, 8, n_ppo, 17, 3, omd_size, comp_dev)
{
	name = "Regularized 3D three-body problem";
	
	initialize();
}

threebody::~threebody()
{ }

void threebody::initialize()
{
    h_md = (threebody_t::metadata_t*)h_omd;
    d_md = (threebody_t::metadata_t*)d_omd;
      md = (threebody_t::metadata_t*)omd;

	h       = 0.0;

	first_open_integral = true;
	first_open_solution = true;
}

void threebody::calc_dy(uint16_t stage, var_t curr_t, const var_t* y_temp, var_t* acc, var_t* jrk)
{
	throw string("The threebody::calc_dy is not implemented.");
}

void threebody::calc_dy(uint16_t stage, var_t curr_t, const var_t* y_temp, var_t* dy)
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

void threebody::calc_integral()
{
	static bool first_call = true;
	const threebody_t::param_t* p = (threebody_t::param_t*)h_par;

	var3_t q1, q2, q3, p1, p2, p3;
	/* Eq. (4.51): */
	const var4_t Q1 = {h_y[0], h_y[1], h_y[2], h_y[3]};
	const var4_t Q2 = {h_y[4], h_y[5], h_y[6], h_y[7]};
	/* Eq. (4.52): */
	const var4_t P1 = {h_y[8], h_y[9], h_y[10], h_y[11]};
	const var4_t P2 = {h_y[12], h_y[13], h_y[14], h_y[15]};

	tools::trans_to_descartes(p[0].m, p[1].m, p[2].m, q1, p1, q2, p2, q3, p3, Q1, P1, Q2, P2);

	/* Eq. (4.9): */
	var_t R = sqrt(SQR(q1.x - q2.x) + SQR(q1.y - q2.y) + SQR(q1.z - q2.z));
	var_t R1 = sqrt(SQR(q1.x - q3.x) + SQR(q1.y - q3.y) + SQR(q1.z - q3.z));
	var_t R2 = sqrt(SQR(q2.x - q3.x) + SQR(q2.y - q3.y) + SQR(q2.z - q3.z));

	var_t k1 = (SQR(p1.x) + SQR(p1.y) + SQR(p1.z)) / (2 * p[0].m);
	var_t k2 = (SQR(p2.x) + SQR(p2.y) + SQR(p2.z)) / (2 * p[1].m);
	var_t k3 = (SQR(p3.x) + SQR(p3.y) + SQR(p3.z)) / (2 * p[2].m);

	var_t u = -((p[0].m * p[2].m) / R1 + (p[1].m * p[2].m) / R2 + (p[0].m * p[1].m) / R);

	/* Eq. (4.10): */
	h = (k1 + k2 + k3) + u;

	if (first_call)
	{
		integral.h0 = integral.h;
		first_call = false;
	}
}

void threebody::cpu_calc_dy(uint16_t stage, var_t curr_t, const var_t* y_temp, var_t* dy)
{
	const threebody_t::param_t* p = (threebody_t::param_t*)h_par;

	/*
	Q1 = y_temp[0], Q2 = y_temp[1], Q3 = y_temp[2], Q4 = y_temp[3],
	Q5 = y_temp[4], Q6 = y_temp[5], Q7 = y_temp[6], Q8 = y_temp[7],
	P1 = y_temp[8], P2 = y_temp[9], P3 = y_temp[10], P4 = y_temp[11],
	P5 = y_temp[12], P6 = y_temp[13], P7 = y_temp[14], P8 = y_temp[15]
	*/

	/* Eq. (4.46): */
	var_t R_temp1 = SQR(y_temp[0]) - SQR(y_temp[1]) - SQR(y_temp[2]) + SQR(y_temp[3]) - SQR(y_temp[4]) + SQR(y_temp[5]) + SQR(y_temp[6]) - SQR(y_temp[7]);
	var_t R_temp2 = 2 * (y_temp[0]*y_temp[1] - y_temp[2]*y_temp[3] - y_temp[4]*y_temp[5] + y_temp[6]*y_temp[7]);
	var_t R_temp3 = 2 * (y_temp[0]*y_temp[2] + y_temp[1]*y_temp[3] - y_temp[4]*y_temp[6] - y_temp[5]*y_temp[7]);
	
	/* Eq. (4.45): */
	var_t R = sqrt(SQR(R_temp1) + SQR(R_temp2) + SQR(R_temp3));
	var_t R1 = SQR(y_temp[0]) + SQR(y_temp[1]) + SQR(y_temp[2]) + SQR(y_temp[3]);
	var_t R2 = SQR(y_temp[4]) + SQR(y_temp[5]) + SQR(y_temp[6]) + SQR(y_temp[7]);

	/* Eq. (4.19): */
	var_t mu13 =  (p[0].m * p[2].m) / (p[0].m + p[2].m);
	var_t mu23 =  (p[1].m * p[2].m) / (p[1].m + p[2].m);

	threebody::calc_integral(); // check if its needed or not!!

	/* p. 25 last paragraph */
	var_t P9 = -h; // check!!!!!


	/* Eq. (4.51): */
	var4_t q1 = {y_temp[0], y_temp[1], y_temp[2], y_temp[3]};
	var4_t q5 = {y_temp[4], y_temp[5], y_temp[6], y_temp[7]};
	/* Eq. (4.52): */
	var4_t p1 = {y_temp[8], y_temp[9], y_temp[10], y_temp[11]};
	var4_t p5 = {y_temp[12], y_temp[13], y_temp[14], y_temp[15]};
	var_t sqr_p1 = SQR(y_temp[8]) + SQR(y_temp[9]) + SQR(y_temp[10]) + SQR(y_temp[11]);
	var_t sqr_p5 = SQR(y_temp[12]) + SQR(y_temp[13]) + SQR(y_temp[14]) + SQR(y_temp[15]);

	// create KS matrices
	/* Eq. (4.53): */
	matrix4_t Q1 = {{y_temp[0],  -y_temp[1],  -y_temp[2],  y_temp[3]}, {y_temp[1],  y_temp[0],  -y_temp[3],  -y_temp[2]}, {y_temp[2],  y_temp[3],  y_temp[0],  y_temp[1]}, {0.0, 0.0, 0.0, 0.0}};	
	matrix4_t Q5 = {{y_temp[4],  -y_temp[5],  -y_temp[6],  y_temp[7]}, {y_temp[5],  y_temp[4],  -y_temp[7],  -y_temp[6]}, {y_temp[6],  y_temp[7],  y_temp[4],  y_temp[5]}, {0.0, 0.0, 0.0, 0.0}};
	/* Eq. (4.54): */
	matrix4_t P1 = {{y_temp[8],  -y_temp[9],  -y_temp[10], y_temp[11]},{y_temp[9],  y_temp[8],  -y_temp[11], -y_temp[10]},{y_temp[10], y_temp[11], y_temp[8],  y_temp[9]}, {0.0, 0.0, 0.0, 0.0}};
	matrix4_t P5 = {{y_temp[12], -y_temp[13], -y_temp[14], y_temp[15]},{y_temp[13], y_temp[12], -y_temp[15], -y_temp[14]},{y_temp[14], y_temp[15], y_temp[12], y_temp[13]},{0.0, 0.0, 0.0, 0.0}};

	var_t c1 = 1/(4*p[2].m);
	var_t c2 = R2/(4*mu13);
	var_t c3 = R1/(4*mu23);

	var4_t c4 = tools::calc_matrix_vector_product(tools::calc_matrix_matrix_product(Q1,tools::calc_matrix_transpose(Q5)), p5);
	var4_t c5 = tools::calc_matrix_vector_product(tools::calc_matrix_matrix_product(Q5,tools::calc_matrix_transpose(Q1)), p1);
	var4_t c6 = tools::calc_matrix_vector_product(tools::calc_matrix_matrix_product(P1,tools::calc_matrix_transpose(Q5)), p5);
	var4_t c7 = tools::calc_matrix_vector_product(tools::calc_matrix_matrix_product(P5,tools::calc_matrix_transpose(Q1)), p1);

	var_t c8 = 2*R2*P9 - 2*p[1].m*p[2].m + sqr_p5/(4*mu23) - 2*p[0].m*p[1].m*R2/R;
	var_t c9 = 2*R1*P9 - 2*p[0].m*p[2].m + sqr_p1/(4*mu13) - 2*p[0].m*p[1].m*R1/R;
	var_t c10 = 2*p[0].m*p[1].m*R1*R2/CUBE(R);

	/* Eq. (4.50): */
	//Qdot (1-8)
	dy[0] = c1 * c4.x + c2 * p1.x;
	dy[1] = c1 * c4.y + c2 * p1.y;
	dy[2] = c1 * c4.z + c2 * p1.z;
	dy[3] = c1 * c4.w + c2 * p1.w;
	
	dy[4] = c1 * c5.x + c3 * p5.x;
	dy[5] = c1 * c5.y + c3 * p5.y;
	dy[6] = c1 * c5.z + c3 * p5.z;
	dy[7] = c1 * c5.w + c3 * p5.w;

	//Pdot (1-8)
	dy[8]  = -(c1 * c6.x + c8 * q1.x + c10 * ( R_temp1*q1.x + R_temp2*q1.y + R_temp3*q1.z ) ); 
	dy[9]  = -(c1 * c6.y + c8 * q1.y + c10 * (-R_temp1*q1.y + R_temp2*q1.x + R_temp3*q1.w ) ); 
	dy[10] = -(c1 * c6.z + c8 * q1.z + c10 * ( R_temp1*q1.z + R_temp2*q1.w - R_temp3*q1.x ) ); 
	dy[11] = -(c1 * c6.w + c8 * q1.w + c10 * ( R_temp1*q1.w - R_temp2*q1.z + R_temp3*q1.y ) );
	
	dy[12] = -(c1 * c7.x + c9 * q5.x - c10 * ( R_temp1*q5.x + R_temp2*q5.y + R_temp3*q5.z ) ); 
	dy[13] = -(c1 * c7.y + c9 * q5.y - c10 * (-R_temp1*q5.y + R_temp2*q5.x + R_temp3*q5.w ) ); 
	dy[14] = -(c1 * c7.z + c9 * q5.z - c10 * ( R_temp1*q5.z + R_temp2*q5.w - R_temp3*q5.x ) ); 
	dy[15] = -(c1 * c7.w + c9 * q5.w - c10 * ( R_temp1*q5.w - R_temp2*q5.z + R_temp3*q5.y ) );

	dy[16] = R1*R2;
}

void threebody::gpu_calc_dy(uint16_t stage, var_t curr_t, const var_t* y_temp, var_t* dy)
{
	throw string("The gpu_calc_dy() is not implemented.");
}

void threebody::load(string& path)
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

void threebody::load_ascii(ifstream& input)
{
	threebody_t::param_t* p = (threebody_t::param_t*)h_par;

	for (uint32_t i = 0; i < n_obj; i++)
	{
		load_ascii_record(input, &h_epoch[i], &h_md[i], &p[i]);
	}

	for (uint32_t i = 0; i < n_var; i++)
	{
		input >> h_y[i];
	}

	//load_ascii_record(input, &h_epoch[i], &h_md[i], &p[i], &h_y[i], &h_y[i+8]);
}

void threebody::load_ascii_record(ifstream& input, var_t* t, threebody_t::metadata_t *md, threebody_t::param_t* p)
{
	string name;

	// epoch
	input >> *t;
	// id
	input >> md->id;
	// mass
	input >> p->m;
}

void threebody::load_binary(ifstream& input)
{
	throw string("The load_binary() is not implemented.");
}

void threebody::print_solution(std::string& path_si, std::string& path_sd, data_rep_t repres)
{
	ofstream sout;

	switch (repres)
	{
	case DATA_REPRESENTATION_ASCII:
		if (first_open_solution) {
			sout.open(path_si.c_str(), ios::out);
			first_open_solution = false;
		}
		else {
			sout.open(path_si.c_str(), ios::out | ios::app);
		}
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

void threebody::print_dump(std::string& path_si, std::string& path_sd)
{
    throw string("The threebody::print_dump function is not implemented.");
}

void threebody::print_solution_ascii(ofstream& sout)
{
	sout.precision(16);
	sout.setf(ios::right);
	sout.setf(ios::scientific);

	const threebody_t::param_t* p = (threebody_t::param_t*)h_par;
	var3_t q1, q2, q3, p1, p2, p3;
	const var4_t Q1 = {h_y[0], h_y[1], h_y[2], h_y[3]};
	const var4_t Q2 = {h_y[4], h_y[5], h_y[6], h_y[7]};
	const var4_t P1 = {h_y[8], h_y[9], h_y[10], h_y[11]};
	const var4_t P2 = {h_y[12], h_y[13], h_y[14], h_y[15]};

	tools::trans_to_descartes(p[0].m, p[1].m, p[2].m, q1, p1, q2, p2, q3, p3, Q1, P1, Q2, P2);


	sout << t << SEP << setw(VAR_T_W) << h_y[n_var-1] << setw(INT_T_W) << SEP;
	for (uint32_t i = 0; i < n_var-1; i++)
    {
		sout << setw(VAR_T_W) << h_y[i];
		if (i < n_var - 2)
		{
			sout << SEP;
		}
		else
		{
			sout << setw(VAR_T_W);
		}
	}
	sout << setw(VAR_T_W) << q1.x << SEP << setw(VAR_T_W) << q1.y << SEP << setw(VAR_T_W) << q1.z <<  SEP << setw(VAR_T_W) << q2.x <<  SEP << setw(VAR_T_W) << q2.y <<  SEP << setw(VAR_T_W) << q2.z <<  SEP << setw(VAR_T_W) << q3.x <<  SEP << setw(VAR_T_W) << q3.y <<  SEP << setw(VAR_T_W) << q3.z << setw(INT_T_W) << SEP;
	sout << setw(VAR_T_W) << p1.x << SEP << setw(VAR_T_W) << p1.y << SEP << setw(VAR_T_W) << p1.z <<  SEP << setw(VAR_T_W) << p2.x <<  SEP << setw(VAR_T_W) << p2.y <<  SEP << setw(VAR_T_W) << p2.z <<  SEP << setw(VAR_T_W) << p3.x <<  SEP << setw(VAR_T_W) << p3.y <<  SEP << setw(VAR_T_W) << p3.z << endl;
	sout.flush();

	//for (uint32_t i = 0; i < n_obj; i++)
 //   {
	//	uint32_t orig_idx = h_md[i].id - 1;

	//	sout << setw(VAR_T_W) << t << SEP                       /* time of the record [day] (double)           */
	//		 << setw(     30) << obj_names[orig_idx] << SEP     /* name of the body         (string = 30 char) */ 
	//	// Print the metadata for each object
 //       << setw(INT_T_W) << h_md[i].id << SEP;

	//	// Print the parameters for each object
	//	for (uint16_t j = 0; j < n_ppo; j++)
	//	{
	//		uint32_t param_idx = i * n_ppo + j;
	//		sout << setw(VAR_T_W) << h_p[param_idx] << SEP;
	//	}
	//	// Print the variables for each object
	//	for (uint16_t j = 0; j < n_vpo; j++)
	//	{
	//		uint32_t var_idx = i * n_vpo + j;
	//		sout << setw(VAR_T_W) << h_y[var_idx];
	//		if (j < n_vpo - 1)
	//		{
	//			sout << SEP;
	//		}
	//		else
	//		{
	//			sout << endl;
	//		}
	//	}
	//}

}

void threebody::print_solution_binary(ofstream& sout)
{
	throw string("The print_solution_binary() is not implemented.");
}

void threebody::print_integral(string& path)
{
	ofstream sout;

	if (first_open_integral)
	{
		sout.open(path.c_str(), ios::out);
		first_open_integral = false;
	}
	else
	{
		sout.open(path.c_str(), ios::out | ios::app);
	}
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

