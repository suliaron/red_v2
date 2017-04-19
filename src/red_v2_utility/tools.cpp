#include <cctype>
#include <cmath>
#include <ctime>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <cstdlib>
#include <string.h>

#include "tools.h"
#include "constants.h"
#include "macro.h"
#include "type.h"

using namespace std;

namespace redutil2
{
namespace tools
{
bool is_number(const string& str)
{
   for (size_t i = 0; i < str.length(); i++)
   {
	   if (std::isdigit(str[i]) || str[i] == 'e' || str[i] == 'E' || str[i] == '.' || str[i] == '-' || str[i] == '+')
	   {
           continue;
	   }
	   else
	   {
		   return false;
	   }
   }
   return true;
}

string& ltrim(string& s)
{
	s.erase(0, s.find_first_not_of(ws));
    return s;
}

string& ltrim(string& s, const char* t)
{
	s.erase(0, s.find_first_not_of(t));
    return s;
}

string& rtrim(string& s)
{
    s.erase(s.find_last_not_of(ws) + 1);
    return s;
}

string& rtrim(string& s, const char* t)
{	
    s.erase(s.find_last_not_of(t) + 1);
    return s;
}

string& trim(string& s)
{
    return ltrim(rtrim(s));
}

string& trim(string& s, const char* t)
{
    return ltrim(rtrim(s, t));
}

string& trim_comment(string& s)
{
	size_t p = s.find_first_of(comment);
	// If not found
	if (string::npos == p)
	{
		return s;
	}
	s.erase(p);
    return s;
}

string& trim_comment(string& s, const char* c)
{
	size_t p = s.find_first_of(c);
	// If not found
	if (string::npos == p)
	{
		return s;
	}
	s.erase(p);
    return s;
}

string get_time_stamp()
{
	static char time_stamp[20];

	time_t now = time(NULL);
    strftime(time_stamp, 20, "%Y-%m-%d %H:%M:%S", localtime(&now));

	return string(time_stamp);
}

string convert_time_t(time_t t)
{
	ostringstream convert;	// stream used for the conversion
	convert << t;			// insert the textual representation of 't' in the characters in the stream
	return convert.str();
}

string convert_var_t(var_t v)
{
    ostringstream convert;	// stream used for the conversion
    convert << v;			// insert the textual representation of 't' in the characters in the stream
    return convert.str();
}


// Draw a number from a given distribution
var_t generate_random(var_t xmin, var_t xmax, var_t p(var_t))
{
	var_t x;
	var_t y;

	do
	{
		x = xmin + (var_t)rand() / RAND_MAX * (xmax - xmin);
		y = (var_t)rand() / RAND_MAX;
	}
	while (y > p(x));

	return x;
}

var_t pdf_const(var_t x)
{
	return 1;
}

void populate_data(uint32_t* n_bodies, pp_disk_t::sim_data_t *sim_data)
{
	int idx = 0;
	int id = 1;

	// Create aliases
	var3_t* r             = sim_data->h_y[0];
	var3_t* v             = sim_data->h_y[1];
	pp_disk_t::param_t* p = sim_data->h_p;
	orbelem_t* oe         = sim_data->h_oe;
	pp_disk_t::body_metadata_t* bmd = sim_data->h_body_md;
	var_t* epoch         = sim_data->h_epoch;

	int upper = n_bodies[BODY_TYPE_STAR];
	for (idx = 0; idx < upper; idx++)
	{
		bmd[idx].body_type = BODY_TYPE_STAR;
		bmd[idx].id = id++;

		p[idx].mass = 1.0;
		p[idx].radius = 1.0 * constants::SolarRadiusToAu;
		
		epoch[idx] = 0.0;

		r[idx].x = r[idx].y = r[idx].z = 0.0;
		v[idx].x = v[idx].y = v[idx].z = 0.0;
	}

	upper += n_bodies[BODY_TYPE_GIANTPLANET];
	for ( ; idx < upper; idx++)
	{
		bmd[idx].body_type = BODY_TYPE_GIANTPLANET;
		bmd[idx].id = id++;
		bmd[idx].mig_type = MIGRATION_TYPE_NO;
		bmd[idx].mig_stop_at = 0.0;

		p[idx].mass = generate_random(1.0, 10.0, pdf_const) * constants::JupiterToSolar;
		p[idx].radius = generate_random(8.0e4, 1.0e5, pdf_const) * constants::KilometerToAu;
		p[idx].density = calc_density(p[idx].mass, p[idx].radius);
		p[idx].cd = 0.0;
		
		epoch[idx] = 0.0;

		oe[idx].sma = generate_random(1.0, 100.0, pdf_const);
		oe[idx].ecc = generate_random(0.0, 0.8, pdf_const);
		oe[idx].inc = generate_random(0.0, 20.0, pdf_const) * constants::DegreeToRadian;
		oe[idx].peri = generate_random(0.0, 360.0, pdf_const) * constants::DegreeToRadian;
		oe[idx].node = generate_random(0.0, 360.0, pdf_const) * constants::DegreeToRadian;
		oe[idx].mean = generate_random(0.0, 360.0, pdf_const) * constants::DegreeToRadian;
	}

	upper += n_bodies[BODY_TYPE_ROCKYPLANET];
	for ( ; idx < upper; idx++)
	{
		bmd[idx].body_type = BODY_TYPE_ROCKYPLANET;
		bmd[idx].id = id++;
		bmd[idx].mig_type = MIGRATION_TYPE_NO;
		bmd[idx].mig_stop_at = 0.0;

		p[idx].mass = generate_random(1.0, 10.0, pdf_const) * constants::EarthToSolar;
		p[idx].radius = generate_random(5.0e3, 8.0e3, pdf_const) * constants::KilometerToAu;
		p[idx].density = calc_density(p[idx].mass, p[idx].radius);
		p[idx].cd = 0.0;
		
		epoch[idx] = 0.0;

		oe[idx].sma = generate_random(1.0, 100.0, pdf_const);
		oe[idx].ecc = generate_random(0.0, 0.8, pdf_const);
		oe[idx].inc = generate_random(0.0, 20.0, pdf_const) * constants::DegreeToRadian;
		oe[idx].peri = generate_random(0.0, 360.0, pdf_const) * constants::DegreeToRadian;
		oe[idx].node = generate_random(0.0, 360.0, pdf_const) * constants::DegreeToRadian;
		oe[idx].mean = generate_random(0.0, 360.0, pdf_const) * constants::DegreeToRadian;
	}

	upper += n_bodies[BODY_TYPE_PROTOPLANET];
	for ( ; idx < upper; idx++)
	{
		bmd[idx].body_type = BODY_TYPE_PROTOPLANET;
		bmd[idx].id = id++;
		bmd[idx].mig_type = MIGRATION_TYPE_NO;
		bmd[idx].mig_stop_at = 0.0;

		p[idx].mass = generate_random(1.0, 10.0, pdf_const) * constants::CeresToSolar;
		p[idx].radius = generate_random(1.0e3, 2.0e3, pdf_const) * constants::KilometerToAu;
		p[idx].density = calc_density(p[idx].mass, p[idx].radius);
		p[idx].cd = 0.0;
		
		epoch[idx] = 0.0;

		oe[idx].sma = generate_random(1.0, 100.0, pdf_const);
		oe[idx].ecc = generate_random(0.0, 0.8, pdf_const);
		oe[idx].inc = generate_random(0.0, 20.0, pdf_const) * constants::DegreeToRadian;
		oe[idx].peri = generate_random(0.0, 360.0, pdf_const) * constants::DegreeToRadian;
		oe[idx].node = generate_random(0.0, 360.0, pdf_const) * constants::DegreeToRadian;
		oe[idx].mean = generate_random(0.0, 360.0, pdf_const) * constants::DegreeToRadian;
	}

	upper += n_bodies[BODY_TYPE_SUPERPLANETESIMAL];
	for ( ; idx < upper; idx++)
	{
		bmd[idx].body_type = BODY_TYPE_SUPERPLANETESIMAL;
		bmd[idx].id = id++;
		bmd[idx].mig_type = MIGRATION_TYPE_NO;
		bmd[idx].mig_stop_at = 0.0;

		p[idx].mass = generate_random(1.0e-2, 1.0e-1, pdf_const) * constants::CeresToSolar;
		p[idx].radius = generate_random(1.0e1, 1.0e2, pdf_const) * constants::KilometerToAu;
		p[idx].density = generate_random(1.0, 3.0, pdf_const) * constants::GramPerCm3ToSolarPerAu3;
		p[idx].cd = 1.0;
		
		epoch[idx] = 0.0;

		oe[idx].sma = generate_random(1.0, 100.0, pdf_const);
		oe[idx].ecc = generate_random(0.0, 0.8, pdf_const);
		oe[idx].inc = generate_random(0.0, 20.0, pdf_const) * constants::DegreeToRadian;
		oe[idx].peri = generate_random(0.0, 360.0, pdf_const) * constants::DegreeToRadian;
		oe[idx].node = generate_random(0.0, 360.0, pdf_const) * constants::DegreeToRadian;
		oe[idx].mean = generate_random(0.0, 360.0, pdf_const) * constants::DegreeToRadian;
	}

	upper += n_bodies[BODY_TYPE_PLANETESIMAL];
	for ( ; idx < upper; idx++)
	{
		bmd[idx].body_type = BODY_TYPE_PLANETESIMAL;
		bmd[idx].id = id++;
		bmd[idx].mig_type = MIGRATION_TYPE_NO;
		bmd[idx].mig_stop_at = 0.0;

		p[idx].mass = generate_random(1.0e-4, 1.0e-3, pdf_const) * constants::CeresToSolar;
		p[idx].radius = generate_random(1.0e1, 1.0e2, pdf_const) * constants::KilometerToAu;
		p[idx].density = calc_density(p[idx].mass, p[idx].radius);
		p[idx].cd = generate_random(0.7, 2.0, pdf_const);
		
		epoch[idx] = 0.0;

		oe[idx].sma = generate_random(1.0, 100.0, pdf_const);
		oe[idx].ecc = generate_random(0.0, 0.8, pdf_const);
		oe[idx].inc = generate_random(0.0, 20.0, pdf_const) * constants::DegreeToRadian;
		oe[idx].peri = generate_random(0.0, 360.0, pdf_const) * constants::DegreeToRadian;
		oe[idx].node = generate_random(0.0, 360.0, pdf_const) * constants::DegreeToRadian;
		oe[idx].mean = generate_random(0.0, 360.0, pdf_const) * constants::DegreeToRadian;
	}

	upper += n_bodies[BODY_TYPE_TESTPARTICLE];
	for ( ; idx < upper; idx++)
	{
		bmd[idx].body_type = BODY_TYPE_TESTPARTICLE;
		bmd[idx].id = id++;
		bmd[idx].mig_type = MIGRATION_TYPE_NO;
		bmd[idx].mig_stop_at = 0.0;

		p[idx].mass = 0.0;
		p[idx].radius = 0.0;
		p[idx].density = 0.0;
		p[idx].cd = 0.0;
		
		epoch[idx] = 0.0;

		oe[idx].sma = generate_random(1.0, 100.0, pdf_const);
		oe[idx].ecc = generate_random(0.0, 0.8, pdf_const);
		oe[idx].inc = generate_random(0.0, 20.0, pdf_const) * constants::DegreeToRadian;
		oe[idx].peri = generate_random(0.0, 360.0, pdf_const) * constants::DegreeToRadian;
		oe[idx].node = generate_random(0.0, 360.0, pdf_const) * constants::DegreeToRadian;
		oe[idx].mean = generate_random(0.0, 360.0, pdf_const) * constants::DegreeToRadian;
	}

	// Calculate coordinates and velocities
	{
		// The mass of the central star
		var_t m0 = sim_data->h_p[0].mass;
		var3_t rVec = {0.0, 0.0, 0.0};
		var3_t vVec = {0.0, 0.0, 0.0};

		// The coordinates of the central star
		sim_data->h_y[0][0] = rVec;
		sim_data->h_y[1][0] = vVec;
		for (int i = 1; i < upper; i++)
		{
			var_t mu = K2 *(m0 + sim_data->h_p[i].mass);
			tools::calc_phase(mu, &sim_data->h_oe[i], &rVec, &vVec);
			sim_data->h_y[0][i] = rVec;
			sim_data->h_y[1][i] = vVec;
		}
	}
}

var_t get_total_mass(uint32_t n, body_type_t type, const pp_disk_t::sim_data_t *sim_data)
{
	var_t total_mass = 0.0;

	pp_disk_t::param_t* p = sim_data->h_p;
	for (int j = n - 1; j >= 0; j--)
	{
		if (type == sim_data->h_body_md[j].body_type && 0 < sim_data->h_body_md[j].id)
		{
			total_mass += p[j].mass;
		}
	}

	return total_mass;
}

var_t get_total_mass(uint32_t n, const pp_disk_t::sim_data_t *sim_data)
{
	var_t M0 = 0.0;

	pp_disk_t::param_t* p = sim_data->h_p;
	for (int j = n - 1; j >= 0; j--)
	{
		if (0 > sim_data->h_body_md[j].id)
		{
			continue;
		}
		M0 += p[j].mass;
	}

	return M0 ;
}

void calc_bc(uint32_t n, const pp_disk_t::sim_data_t *sim_data, var_t M0, var3_t* R0, var3_t* V0)
{
	const var3_t* r = sim_data->h_y[0];
	const var3_t* v = sim_data->h_y[1];
	const pp_disk_t::param_t* p = sim_data->h_p;

	R0->x = R0->y = R0->z = 0.0;
	V0->x = V0->y = V0->z = 0.0;

	for (int j = n - 1; j >= 0; j-- )
	{
		if (0 > sim_data->h_body_md[j].id)
		{
			continue;
		}
		var_t m_j = p[j].mass;
		R0->x += m_j * r[j].x;
		R0->y += m_j * r[j].y;
		R0->z += m_j * r[j].z;

		V0->x += m_j * v[j].x;
		V0->y += m_j * v[j].y;
		V0->z += m_j * v[j].z;
	}
	R0->x /= M0;	R0->y /= M0;	R0->z /= M0;
	V0->x /= M0;	V0->y /= M0;	V0->z /= M0;
}

void transform_to_bc(uint32_t n, const pp_disk_t::sim_data_t *sim_data)
{
	// Position and velocity of the system's barycenter
	var3_t R0 = {0.0, 0.0, 0.0};
	var3_t V0 = {0.0, 0.0, 0.0};

	var_t M0 = get_total_mass(n, sim_data);
	calc_bc(n, sim_data, M0, &R0, &V0);

	var3_t* r = sim_data->h_y[0];
	var3_t* v = sim_data->h_y[1];
	// Transform the bodies coordinates and velocities
	for (int j = n - 1; j >= 0; j--)
	{
		if (0 > sim_data->h_body_md[j].id)
		{
			continue;
		}
		r[j].x -= R0.x;		r[j].y -= R0.y;		r[j].z -= R0.z;
		v[j].x -= V0.x;		v[j].y -= V0.y;		v[j].z -= V0.z;
	}
}

void transform_time(uint32_t n, const pp_disk_t::sim_data_t *sim_data)
{
	for (uint32_t i = 0; i < n; i++)
	{
		sim_data->h_epoch[i] *= constants::Gauss;
	}
}

void transform_velocity(uint32_t n, const pp_disk_t::sim_data_t *sim_data)
{
	var3_t* v = sim_data->h_y[1];
	// Transform the bodies' velocities
	for (uint32_t i = 0; i < n; i++)
	{
		v[i].x /= constants::Gauss;	
		v[i].y /= constants::Gauss;	
		v[i].z /= constants::Gauss;
	}
}

var_t calc_radius(var_t m, var_t density)
{
	static var_t four_pi_over_three = 4.188790204786391;

	return pow(1.0/four_pi_over_three * m/density, 1.0/3.0);
}

var_t calc_density(var_t m, var_t R)
{
	static var_t four_pi_over_three = 4.188790204786391;

	if (R == 0.0)
	{
		return 0.0;
	}
	return m / (four_pi_over_three * CUBE(R));
}

var_t calc_mass(var_t R, var_t density)
{
	static var_t four_pi_over_three = 4.188790204786391;

	return four_pi_over_three * CUBE(R) * density;
}

void calc_position_after_collision(var_t m1, var_t m2, const var3_t* r1, const var3_t* r2, var3_t& r)
{
	const var_t M = m1 + m2;

	r.x = (m1 * r1->x + m2 * r2->x) / M;
	r.y = (m1 * r1->y + m2 * r2->y) / M;
	r.z = (m1 * r1->z + m2 * r2->z) / M;
}

void calc_velocity_after_collision(var_t m1, var_t m2, const var3_t* v1, const var3_t* v2, var3_t& v)
{
	const var_t M = m1 + m2;

	v.x = (m1 * v1->x + m2 * v2->x) / M;
	v.y = (m1 * v1->y + m2 * v2->y) / M;
	v.z = (m1 * v1->z + m2 * v2->z) / M;
}

void calc_physical_properties(const pp_disk_t::param_t &p1, const pp_disk_t::param_t &p2, pp_disk_t::param_t &p)
{
	// Calculate V = V1 + V2
	var_t volume = 4.188790204786391 * (CUBE(p1.radius) + CUBE(p2.radius));

	p.mass	  = p1.mass + p2.mass;
	p.density = p.mass / volume;
	p.radius  = calc_radius(p.mass, p.density);
	p.cd      = p1.cd;
}

var_t norm(const var4_t* r)
{
	return sqrt(SQR(r->x) + SQR(r->y) + SQR(r->z));
}

var_t norm(const var3_t* r)
{
	return sqrt(SQR(r->x) + SQR(r->y) + SQR(r->z));
}

var_t calc_dot_product(const var3_t& u, const var3_t& v)
{
    return (u.x * v.x + u.y * v.y + u.z * v.z);
}

var3_t calc_cross_product(const var3_t& u, const var3_t& v)
{
    var3_t result = {0.0, 0.0, 0.0};
    
    result.x = u.y * v.z - u.z * v.y;
    result.y = u.z * v.x - u.x * v.z;
    result.z = u.x * v.y - u.y * v.x;

    return result;
}

matrix4_t calc_matrix_matrix_product(const matrix4_t& u, const matrix4_t& v)
{
	matrix4_t result = {0.0, 0.0, 0.0, 0.0};

	result.a.x = u.a.x * v.a.x + u.b.x * v.a.y  + u.c.x * v.a.z + u.d.x * v.a.w;
	result.a.y = u.a.y * v.a.x + u.b.y * v.a.y  + u.c.y * v.a.z + u.d.y * v.a.w;
	result.a.z = u.a.z * v.a.x + u.b.z * v.a.y  + u.c.z * v.a.z + u.d.z * v.a.w;
	result.a.w = u.a.w * v.a.x + u.b.w * v.a.y  + u.c.w * v.a.z + u.d.w * v.a.w;

	result.b.x = u.a.x * v.b.x + u.b.x * v.b.y  + u.c.x * v.b.z + u.d.x * v.b.w;
	result.b.y = u.a.y * v.b.x + u.b.y * v.b.y  + u.c.y * v.b.z + u.d.y * v.b.w;
	result.b.z = u.a.z * v.b.x + u.b.z * v.b.y  + u.c.z * v.b.z + u.d.z * v.b.w;
	result.b.w = u.a.w * v.b.x + u.b.w * v.b.y  + u.c.w * v.b.z + u.d.w * v.b.w;

	result.c.x = u.a.x * v.c.x + u.b.x * v.c.y  + u.c.x * v.c.z + u.d.x * v.c.w;
	result.c.y = u.a.y * v.c.x + u.b.y * v.c.y  + u.c.y * v.c.z + u.d.y * v.c.w;
	result.c.z = u.a.z * v.c.x + u.b.z * v.c.y  + u.c.z * v.c.z + u.d.z * v.c.w;
	result.c.w = u.a.w * v.c.x + u.b.w * v.c.y  + u.c.w * v.c.z + u.d.w * v.c.w;

	result.d.x = u.a.x * v.d.x + u.b.x * v.d.y  + u.c.x * v.d.z + u.d.x * v.d.w;
	result.d.y = u.a.y * v.d.x + u.b.y * v.d.y  + u.c.y * v.d.z + u.d.y * v.d.w;
	result.d.z = u.a.z * v.d.x + u.b.z * v.d.y  + u.c.z * v.d.z + u.d.z * v.d.w;
	result.d.w = u.a.w * v.d.x + u.b.w * v.d.y  + u.c.w * v.d.z + u.d.w * v.d.w;

	return result;
}

var4_t calc_matrix_vector_product(const matrix4_t& u, const var4_t& v)
{
	var4_t result = {0.0, 0.0, 0.0, 0.0};
    
	result.x = u.a.x * v.x + u.b.x * v.y  + u.c.x * v.z + u.d.x * v.w;
	result.y = u.a.y * v.x + u.b.y * v.y  + u.c.y * v.z + u.d.y * v.w;
	result.z = u.a.z * v.x + u.b.z * v.y  + u.c.z * v.z + u.d.z * v.w;
	result.w = u.a.w * v.x + u.b.w * v.y  + u.c.w * v.z + u.d.w * v.w;

	return result;
}

matrix4_t calc_matrix_transpose(const matrix4_t& u)
{
	matrix4_t result = {0.0, 0.0, 0.0, 0.0};

	result.a.x = u.a.x;
	result.a.y = u.b.x;
	result.a.z = u.c.x;
	result.a.w = u.d.x;

	result.b.x = u.a.y;
	result.b.y = u.b.y;
	result.b.z = u.c.y;
	result.b.w = u.d.y;

	result.c.x = u.a.z;
	result.c.y = u.b.z;
	result.c.z = u.c.z;
	result.c.w = u.d.z;

	result.d.x = u.a.w;
	result.d.y = u.b.w;
	result.d.z = u.c.w;
	result.d.w = u.d.w;

	return result;
}

void trans_to_descartes(const var_t m1, const var_t m2, const var_t m3, var3_t& qv1, var3_t& pv1, var3_t& qv2, var3_t& pv2, var3_t& qv3, var3_t& pv3, const var4_t& Q1, const var4_t& P1, const var4_t& Q2, const var4_t& P2)
{
	//const threebody_t::param_t* p = (threebody_t::param_t*)h_p;

	/* Eq. (4.53): */
	matrix4_t A1 = {{2*Q1.x, -2*Q1.y, -2*Q1.z, 2*Q1.w},{2*Q1.y, 2*Q1.x, -2*Q1.w, -2*Q1.z},{2*Q1.z, 2*Q1.w, 2*Q1.x, 2*Q1.y},{0,0,0,0}};
	matrix4_t A2 = {{2*Q2.x, -2*Q2.y, -2*Q2.z, 2*Q2.w},{2*Q2.y, 2*Q2.x, -2*Q2.w, -2*Q2.z},{2*Q2.z, 2*Q2.w, 2*Q2.x, 2*Q2.y},{0,0,0,0}};

	/* Eq. (4.45): */
	var_t R1 = SQR(Q1.x) + SQR(Q1.y) + SQR(Q1.z) + SQR(Q1.w);
	var_t R2 = SQR(Q2.x) + SQR(Q2.y) + SQR(Q2.z) + SQR(Q2.w);

	/* Eq. (4.59): */
	var4_t q1 = calc_matrix_vector_product(tools::calc_matrix_transpose(A1),Q1);
	q1.x /= 2;
	q1.y /= 2;
	q1.z /= 2;
	
	var4_t q2 = calc_matrix_vector_product(tools::calc_matrix_transpose(A2),Q2);
	q2.x /= 2;
	q2.y /= 2;
	q2.z /= 2;

	/* Eq. (4.60): */
	var4_t p1 = calc_matrix_vector_product(tools::calc_matrix_transpose(A1),P1);
	p1.x /= 4*R1;
	p1.y /= 4*R1;
	p1.z /= 4*R1;

	var4_t p2 = calc_matrix_vector_product(tools::calc_matrix_transpose(A2),P2);
	p2.x /= 4*R2;
	p2.y /= 4*R2;
	p2.z /= 4*R2;

	/* Eq. (4.61): */
	// q3' , q1', q2'
	qv3.x = -(m1 * q1.x + m2 * q2.x) / (m1 + m2 + m3);
	qv3.y = -(m1 * q1.y + m2 * q2.y) / (m1 + m2 + m3);
	qv3.z = -(m1 * q1.z + m2 * q2.z) / (m1 + m2 + m3);

	qv1.x = qv3.x + q1.x;
	qv1.y = qv3.y + q1.y;
	qv1.z = qv3.z + q1.z;

	qv2.x = qv3.x + q2.x;
	qv2.y = qv3.y + q2.y;
	qv2.z = qv3.z + q2.z;

	/* Eq. (4.62): */
	// p1', p2', p3'
	pv1.x = p1.x;
	pv1.y = p1.y;
	pv1.z = p1.z;

	pv2.x = p2.x;
	pv2.y = p2.y;
	pv2.z = p2.z;

	pv3.x = -(p1.x + p2.x);
	pv3.y = -(p1.y + p2.y);
	pv3.z = -(p1.z + p2.z);
}

void trans_to_threebody(const var3_t& qv1, const var3_t& pv1, const var3_t& qv2, const var3_t& pv2, const var3_t& qv3, const var3_t& pv3, var4_t& Q1, var4_t& P1, var4_t& Q2, var4_t& P2)
{
	var3_t q1, q2;
	/* Eq. (4.56): */
	q1.x = qv1.x - qv3.x;
	q1.y = qv1.y - qv3.y;
	q1.z = qv1.z - qv3.z;
	
	q2.x = qv2.x - qv3.x;
	q2.y = qv2.y - qv3.y;
	q2.z = qv2.z - qv3.z;

	var4_t p1 = {pv1.x, pv1.y, pv1.z, 0};
	var4_t p2 = {pv2.x, pv2.y, pv2.z, 0};

	/* Eq. (4.57): */
	if (q1.x >= 0) {
		Q1.x = sqrt(0.5 * (norm(&q1)  + q1.x));
		Q1.y = 0.5 * q1.y / Q1.x;
		Q1.z = 0.5 * q1.z / Q1.x;
		Q1.w = 0;
	}
	/* Eq. (4.58): */
	else {
		Q1.y = sqrt(0.5 * (norm(&q1)  - q1.x));
		Q1.x = 0.5 * q1.y / Q1.y;
		Q1.z = 0;
		Q1.w = 0.5 * q1.z / Q1.y;	
	}

	/* Eq. (4.57): */
	if (q2.x >= 0) {
		Q2.x = sqrt(0.5 * (norm(&q2)  + q2.x));
		Q2.y = 0.5 * q2.y / Q2.x;
		Q2.z = 0.5 * q2.z / Q2.x;
		Q2.w = 0;
	}
	/* Eq. (4.58): */
	else {
		Q2.y = sqrt(0.5 * (norm(&q2)  - q2.x));
		Q2.x = 0.5 * q2.y / Q2.y;
		Q2.z = 0;
		Q2.w = 0.5 * q2.z / Q2.y;	
	}

	/* Eq. (4.53): */
	matrix4_t A1 = {{2*Q1.x, -2*Q1.y, -2*Q1.z, 2*Q1.w},{2*Q1.y, 2*Q1.x, -2*Q1.w, -2*Q1.z},{2*Q1.z, 2*Q1.w, 2*Q1.x, 2*Q1.y},{0,0,0,0}};
	matrix4_t A2 = {{2*Q2.x, -2*Q2.y, -2*Q2.z, 2*Q2.w},{2*Q2.y, 2*Q2.x, -2*Q2.w, -2*Q2.z},{2*Q2.z, 2*Q2.w, 2*Q2.x, 2*Q2.y},{0,0,0,0}};

	/* Eq. (4.30): */
	P1 = calc_matrix_vector_product(A1,p1);
	P2 = calc_matrix_vector_product(A2,p2);	
}

var_t calc_kinetic_energy(const var3_t* v)
{
	return (SQR(v->x) + SQR(v->y) + SQR(v->z)) / 2.0;
}

var_t calc_pot_energy(var_t mu, const var3_t* r)
{
    return -mu / norm(r);
}

var_t calc_energy(var_t mu, const var3_t* r, const var3_t* v)
{
	return calc_kinetic_energy(v) + calc_pot_energy(mu, r);
}

var3_t calc_angular_momentum(uint32_t n, const pp_disk_t::sim_data_t *sim_data)
{
    var3_t result = {0.0, 0.0, 0.0};
    
	var3_t* r = sim_data->h_y[0];
	var3_t* v = sim_data->h_y[1];
    pp_disk_t::param_t* p = sim_data->h_p;

    for (uint32_t i = 0; i < n; i++)
    {
		if (0 > sim_data->h_body_md[i].id)
		{
			continue;
		}
        var3_t c = calc_cross_product(r[i], v[i]);
        c.x *= p[i].mass; c.y *= p[i].mass; c.z *= p[i].mass;
		result.x += c.x; result.y += c.y; result.z += c.z;
	}

	return result;
}

var3_t calc_angular_momentum_CMU(uint32_t n, const pp_disk_t::sim_data_t *sim_data)
{
    var3_t result = {0.0, 0.0, 0.0};
    
	var3_t* r = sim_data->h_y[0];
	var3_t* v = sim_data->h_y[1];
    pp_disk_t::param_t* p = sim_data->h_p;

    for (uint32_t i = 0; i < n; i++)
    {
		if (0 > sim_data->h_body_md[i].id)
		{
			continue;
		}
		var3_t vv = {v[i].x * constants::Gauss, v[i].y * constants::Gauss, v[i].z * constants::Gauss};
        var3_t c = calc_cross_product(r[i], vv);

		c.x *= p[i].mass; c.y *= p[i].mass; c.z *= p[i].mass;
		result.x += c.x; result.y += c.y; result.z += c.z;
	}

	return result;
}

var_t calc_potential_energy(uint32_t n, const pp_disk_t::sim_data_t *sim_data)
{
	var_t result = 0.0;

	var3_t* r = sim_data->h_y[0];
	var3_t* v = sim_data->h_y[1];
    pp_disk_t::param_t* p = sim_data->h_p;

    for (uint32_t i = 0; i < n; i++)
    {
		if (0 > sim_data->h_body_md[i].id)
		{
			continue;
		}
        for (uint32_t j = 0; j < n; j++)
        {
            if (i == j)
            {
                continue;
            }
			if (0 > sim_data->h_body_md[j].id)
			{
				continue;
			}
            var_t r_ij = sqrt(SQR(r[j].x - r[i].x) + SQR(r[j].y - r[i].y) + SQR(r[j].z - r[i].z));
            result += p[i].mass * p[j].mass / r_ij;
        }
	}

    return (result / 2.0);
}

var_t calc_potential_energy_CMU(uint32_t n, const pp_disk_t::sim_data_t *sim_data)
{
	var_t result = 0.0;

	var3_t* r = sim_data->h_y[0];
	var3_t* v = sim_data->h_y[1];
    pp_disk_t::param_t* p = sim_data->h_p;

    for (uint32_t i = 0; i < n; i++)
    {
		if (0 > sim_data->h_body_md[i].id)
		{
			continue;
		}
        for (uint32_t j = 0; j < n; j++)
        {
            if (i == j)
            {
                continue;
            }
			if (0 > sim_data->h_body_md[j].id)
			{
				continue;
			}
            var_t r_ij = sqrt(SQR(r[j].x - r[i].x) + SQR(r[j].y - r[i].y) + SQR(r[j].z - r[i].z));
            result += p[i].mass * p[j].mass / r_ij;
        }
	}

	return (constants::Gauss2 * result / 2.0);
}

var_t calc_kinetic_energy(uint32_t n, const pp_disk_t::sim_data_t *sim_data)
{
	var_t result = 0.0;

	var3_t* v = sim_data->h_y[1];
    pp_disk_t::param_t* p = sim_data->h_p;

    for (uint32_t i = 0; i < n; i++)
    {
		if (0 > sim_data->h_body_md[i].id)
		{
			continue;
		}
        result += p[i].mass * (SQR(v[i].x) + SQR(v[i].y) + SQR(v[i].z));
    }

    return (result / 2.0);
}

var_t calc_kinetic_energy_CMU(uint32_t n, const pp_disk_t::sim_data_t *sim_data)
{
	var_t result = 0.0;

	var3_t* v = sim_data->h_y[1];
    pp_disk_t::param_t* p = sim_data->h_p;

    for (uint32_t i = 0; i < n; i++)
    {
		if (0 > sim_data->h_body_md[i].id)
		{
			continue;
		}
		var4_t vv = {v[i].x * constants::Gauss, v[i].y * constants::Gauss, v[i].z * constants::Gauss, 0.0};
        result += p[i].mass * (SQR(vv.x) + SQR(vv.y) + SQR(vv.z));
    }

    return (result / 2.0);
}

var_t calc_total_energy(uint32_t n, const pp_disk_t::sim_data_t *sim_data)
{
    return (calc_kinetic_energy(n, sim_data) - calc_potential_energy(n, sim_data));
}

var_t calc_total_energy_CMU(uint32_t n, const pp_disk_t::sim_data_t *sim_data)
{
	return (calc_kinetic_energy_CMU(n, sim_data) - calc_potential_energy_CMU(n, sim_data));
}

void shift_into_range(var_t lower, var_t upper, var_t &value)
{
    var_t range = upper - lower;
    while (value >= upper)
    {
        value -= range;
    }
    while (value < lower)
    {
        value += range;
    }
}

// TODO if ecc >= 0.99 the iteration does not converge!!
void kepler_equation_solver(var_t ecc, var_t mean, var_t eps, var_t* E)
{
	if (ecc == 0.0 || mean == 0.0 || mean == PI)
	{
        *E = mean;
		return;
    }
    *E = mean + ecc * (sin(mean)) / (1.0 - sin(mean + ecc) + sin(mean));
    var_t E1 = 0.0;
    var_t error;
    int step = 0;
    do
	{
        E1 = *E - (*E - ecc * sin(*E) - mean) / (1.0 - ecc * cos(*E));
        error = fabs(E1 - *E);
        *E = E1;
    } while (error > eps && step++ <= 15);
	if (15 < step)
	{
		throw string("The kepler_equation_solver() failed: solution did not converge.");
	}
}

void calc_phase(var_t mu, const orbelem_t* oe, var3_t* rVec, var3_t* vVec)
{
    var_t ecc = oe->ecc;
	var_t E = 0.0;
	kepler_equation_solver(ecc, oe->mean, 1.0e-14, &E);
    var_t v = 2.0 * atan(sqrt((1.0 + ecc) / (1.0 - ecc)) * tan(E / 2.0));

    var_t p = oe->sma * (1.0 - SQR(ecc));
    var_t r = p / (1.0 + ecc * cos(v));
    var_t kszi = r * cos(v);
    var_t eta = r * sin(v);
    var_t vKszi = -sqrt(mu / p) * sin(v);
    var_t vEta = sqrt(mu / p) * (ecc + cos(v));

    var_t cw = cos(oe->peri);
    var_t sw = sin(oe->peri);
    var_t cO = cos(oe->node);
    var_t sO = sin(oe->node);
    var_t ci = cos(oe->inc);
    var_t si = sin(oe->inc);

    var3_t P;
	P.x = cw * cO - sw * sO * ci;
	P.y = cw * sO + sw * cO * ci;
	P.z = sw * si;
    var3_t Q;
	Q.x = -sw * cO - cw * sO * ci;
	Q.y = -sw * sO + cw * cO * ci;
	Q.z = cw * si;

	rVec->x = kszi * P.x + eta * Q.x;
	rVec->y = kszi * P.y + eta * Q.y;
	rVec->z = kszi * P.z + eta * Q.z;

	vVec->x = vKszi * P.x + vEta * Q.x;
	vVec->y = vKszi * P.y + vEta * Q.y;
	vVec->z = vKszi * P.z + vEta * Q.z;
}

void calc_oe(var_t mu, const var3_t* rVec, const var3_t* vVec, orbelem_t* oe)
{
    const var_t sq2 = 1.0e-14;
    const var_t sq3 = 1.0e-14;

	var_t r_norm = norm(rVec);
	var_t v_norm = norm(vVec);

	// Calculate energy, h
    var_t h = calc_energy(mu, rVec, vVec);
    if (h >= 0.0)
    {
		throw string("The Kepler-energy is positive. calc_oe() failed.");
    }

	var4_t cVec;
    cVec.x = rVec->y * vVec->z - rVec->z * vVec->y;
    cVec.y = rVec->z * vVec->x - rVec->x * vVec->z;
    cVec.z = rVec->x * vVec->y - rVec->y * vVec->x;
	cVec.w = 0.0;
	var_t c_norm = norm(&cVec);

	var4_t lVec;
	lVec.x = -mu / r_norm * rVec->x + vVec->y * cVec.z - vVec->z * cVec.y;
	lVec.y = -mu / r_norm * rVec->y + vVec->z * cVec.x - vVec->x * cVec.z;
	lVec.z = -mu / r_norm * rVec->z + vVec->x * cVec.y - vVec->y * cVec.x;
	lVec.w = 0.0;
	var_t l_norm = norm(&lVec);

	var_t rv = rVec->x * vVec->x + rVec->y * vVec->y + rVec->z * vVec->z;

    var_t e2 = 1.0 + 2.0 * SQR(c_norm) * h / (SQR(mu));
    if (abs(e2) < sq3)
    {
        e2 = 0.0;
    }
    var_t e = sqrt(e2);
    /*
    * Calculate semi-major axis, a
    */
    var_t a = -mu / (2.0 * h);
    /*
    * Calculate inclination, incl
    */
    var_t cosi = cVec.z / c_norm;
    var_t sini = sqrt(cVec.x * cVec.x + cVec.y * cVec.y) / c_norm;
    var_t incl = acos(cosi);
    if (incl < sq2)
    {
        incl = 0.0;
    }
    /*
    * Calculate longitude of node, O
    */
    var_t node = 0.0;
    if (incl != 0.0)
    {
        var_t tmpx = -cVec.y / (c_norm * sini);
        var_t tmpy =  cVec.x / (c_norm * sini);
		node = atan2(tmpy, tmpx);
		shift_into_range(0.0, 2.0*PI, node);
    }
    /*
    * Calculate argument of pericenter, w
    */
    var_t E = 0.0;
    var_t peri = 0.0;
    if (e2 != 0.0)
    {
        var_t tmpx = (lVec.x * cos(node) + lVec.y * sin(node)) / l_norm;
        var_t tmpy = (-lVec.x * sin(node) + lVec.y * cos(node)) / (l_norm * cosi);
        peri = atan2(tmpy, tmpx);
        shift_into_range(0.0, 2.0*PI, peri);

        tmpx = 1.0 / e * (1.0 - r_norm / a);
        tmpy = rv / (sqrt(mu * a) * e);
        E = atan2(tmpy, tmpx);
        shift_into_range(0.0, 2.0*PI, E);
    }
    else
    {
        peri = 0.0;
        E = atan2(rVec->y, rVec->x);
        shift_into_range(0, 2.0*PI, E);
    }
    /*
    * Calculate mean anomaly, M
    */
    var_t M = E - e * sin(E);
    shift_into_range(0, 2.0*PI, M);

	oe->sma  = a;
	oe->ecc	 = e;
	oe->inc  = incl;
	oe->peri = peri;
	oe->node = node;
	oe->mean = M;
}

var_t calc_orbital_period(var_t mu, var_t a)
{
	return (TWOPI*sqrt(CUBE(a)/mu));
}

// Date of creation: 2016.11.02.
// Last edited: 
// Status: Tested
void calc_lin_comb(var_t* a, const var_t* const *c, const var_t* coeff, uint16_t n_vct, uint32_t n_var)
{
	for (uint32_t i = 0; i < n_var; i++)
	{
		a[i] = 0.0;
		for (uint16_t j = 0; j < n_vct; j++)
		{
			if (0.0 == coeff[j])
			{
				continue;
			}
			a[i] += coeff[j] * c[j][i];
		}
	}
}

// Date of creation: 2016.11.02.
// Last edited: 
// Status: Tested
void calc_lin_comb_s(var_t* a, const var_t* b, const var_t* c, var_t f, uint32_t n_var)
{
	for (uint32_t i = 0; i < n_var; i++)
	{
		a[i] = b[i] + f*c[i];
	}
}

// Date of creation: 2016.11.02.
// Last edited: 
// Status: Tested
void calc_lin_comb_s(var_t* a, const var_t* b, const var_t* const *c, const var_t* coeff, uint16_t n_vct, uint32_t n_var)
{
	for (uint32_t i = 0; i < n_var; i++)
	{
		var_t d = 0.0;
		for (uint16_t j = 0; j < n_vct; j++)
		{
			if (0.0 == coeff[j])
			{
				continue;
			}
			d += coeff[j] * c[j][i];
		}
		a[i] = b[i] + d;
	}
}

void print_vector(const var3_t *v)
{
	cout.precision(16);
	cout.setf(ios::right);
	cout.setf(ios::scientific);

	cout << setw(VAR_T_W) << v->x << SEP
		 << setw(VAR_T_W) << v->y << SEP
		 << setw(VAR_T_W) << v->z << endl;
}

void print_vector(const var4_t *v)
{
	cout.precision(16);
	cout.setf(ios::right);
	cout.setf(ios::scientific);

	cout << setw(VAR_T_W) << v->x 
		 << setw(VAR_T_W) << v->y
		 << setw(VAR_T_W) << v->z
		 << setw(VAR_T_W) << v->w << endl;
}

void print_parameter(const pp_disk_t::param_t *p)
{
	cout.precision(16);
	cout.setf(ios::right);
	cout.setf(ios::scientific);

	cout << setw(VAR_T_W) << p->mass 
		 << setw(VAR_T_W) << p->radius
		 << setw(VAR_T_W) << p->density
		 << setw(VAR_T_W) << p->cd << endl;
}

void print_body_metadata(const pp_disk_t::body_metadata_t *b)
{
	static const char* body_type_name[] = 
	{
		"STAR",
		"GIANTPLANET",
		"ROCKYPLANET",
		"PROTOPLANET",
		"SUPERPLANETESIMAL",
		"PLANETESIMAL",
		"TESTPARTICLE"
	};
	static const char* migration_type_name[] = 
	{
		"NO",
		"TYPE_I",
		"TYPE_II"
	};

	cout << setw(INT_T_W) << b->id << endl;
	cout << (char)(48 + b->body_type) << " (" << body_type_name[b->body_type] << ")" << endl
		 << (char)(48 + b->mig_type) << " (" << migration_type_name[b->mig_type] << ")" << endl
		 << (char)(48 + b->active) << (b->active ? " (true)" : " (false)") << endl
		 << (char)(48 + b->unused) << (b->unused ? " (true)" : " (false)") << endl;
	cout.precision(16);
	cout.setf(ios::right);
	cout.setf(ios::scientific);
	cout << setw(VAR_T_W) << b->mig_stop_at << endl;
}

namespace tbp
{
var_t calc_integral(var_t mu, var2_t r, var2_t v)
{
	var_t d  = sqrt(SQR(r.x) + SQR(r.y));
	var_t v2 = SQR(v.x) + SQR(v.y);
	var_t h = 0.5*v2 - mu / d;
	return h;
}
} /* namespace tbp */

namespace rtbp1D
{
// Calculate parametric coordinate and velocity
void transform_x2u(var_t x, var_t vx, var_t& u, var_t& v)
{
	// TODO: what if x < 0 ?!
    u = sqrt(x);         // u
    v = 0.5* u * vx;     // u'
}

// Calculate descartes coordinate and velocity
void transform_u2x(var_t u, var_t v, var_t& x, var_t& vx)
{
    x =  u*u;            // x       ÉB. p.85 Eq. (10.21)
	// TODO: what if u = 0 ?!
    vx = 2.0 * (v/u);    // vx
}
} /* namespace rtbp1D */

namespace rtbp2D
{
void transform_x2u(var2_t x, var2_t& u)
{
	var_t r = sqrt(SQR(x.x) + SQR(x.y));

	// If x_1 = 0 and x_2 = 0
	if (0.0 == x.x && 0.0 == x.y)
	{
		u.x = u.y = 0.0;
	}
	// If x_1 = 0 and x_2 != 0
	else if (0.0 == x.x && 0.0 != x.y)
	{
		if (0.0 < x.y)
		{
			u.y = sqrt(x.y/2.0);
			u.x = u.y;
		}
		else
		{
			u.y = sqrt(-x.y/2.0);
			u.x = -u.y;
		}
	}
	// If x_1 != 0 and x_2 = 0
	else if (0.0 != x.x && 0.0 == x.y)
	{
		if (0.0 < x.x)
		{
			u.x = sqrt(x.x);
			u.y = 0.0;
		}
		else
		{
			u.x = 0.0;
			u.y = -sqrt(-x.x);
		}
	}
	/*
	 * I. Quadrant: x_1 > 0 and x_2 > 0
	 * -> u_1 > u_2
	 * -> sgn(u_1) = sgn(u_2)
	 */
	else if (0.0 < x.x && 0.0 < x.y) 
	{
		u.x = sqrt((x.x + r)/2.0);
		u.y = x.y/(2.0*u.x);
	}
	// II. Quadrant
	else if (0.0 > x.x && 0.0 < x.y) 
	{
		u.y = sqrt((-x.x + r)/2.0);
		u.x = x.y/(2.0*u.y);
	}
	// III. Quadrant
	else if (0.0 > x.x && 0.0 > x.y) 
	{
		u.y = sqrt((-x.x + r)/2.0);
		u.x = x.y/(2.0*u.y);
	}
	// IV. Quadrant
	else
	{
		u.x = sqrt((x.x + r)/2.0);
		u.y = x.y/(2.0*u.x);
	}
}

void transform_xd2up(var2_t u, var2_t xd, var2_t& up)
{
	var_t r = SQR(u.x) + SQR(u.y);

	var2_t xp = {r*xd.x, r*xd.y};

	var_t f = 1.0/(2.0*r);
	up.x = f * ( u.x * xp.x + u.y * xp.y);
	up.y = f * (-u.y * xp.x + u.x * xp.y);
}

void transform_u2x(var2_t u, var2_t& x)
{
	x.x = SQR(u.x) - SQR(u.y);
	x.y = 2.0*(u.x*u.y);
}

void transform_up2xd(var2_t u, var2_t up, var2_t& xd)
{
	var2_t xp = {0, 0};

	xp.x = 2.0*(u.x*up.x - u.y*up.y);
	xp.y = 2.0*(u.y*up.x + u.x*up.y);

	var_t f = 1.0 / (SQR(u.x) + SQR(u.y));
	xd.x = f * xp.x;
	xd.y = f * xp.y;
}

var_t calc_integral(var_t mu, var2_t u, var2_t up)
{
	var_t d  = SQR( u.x) + SQR( u.y);
	var_t v2 = SQR(up.x) + SQR(up.y);
	var_t h = (2.0*v2 - mu) / d;

	return h;
}

} /* namespace rtbp2D */

namespace nbp
{
var_t get_total_mass(uint32_t n, const nbp_t::param_t* p)
{
	var_t M0 = 0.0;

	for (uint32_t i = 0; i < n; i++)
	{
		M0 += p[i].mass;
	}

	return M0 ;
}
	
var_t calc_total_energy(uint32_t n, const nbp_t::param_t* p, const var3_t* r, const var3_t* v)
{
	var_t T = calc_kinetic_energy(n, p, v);
	var_t U = calc_potential_energy(n, p, r);
	var_t h = T - U;
	return h;
}

var_t calc_kinetic_energy(uint32_t n, const nbp_t::param_t* p, const var3_t* v)
{	
	var_t result = 0.0;

    for (uint32_t i = 0; i < n; i++)
    {
        result += p[i].mass * (SQR(v[i].x) + SQR(v[i].y) + SQR(v[i].z));
    }

    return (result / 2.0);
}

var_t calc_potential_energy(uint32_t n, const nbp_t::param_t* p, const var3_t* r)
{
	var_t result = 0.0;

    for (uint32_t i = 0; i < n; i++)
    {
        for (uint32_t j = 0; j < n; j++)
        {
            if (i == j)
            {
                continue;
            }
            var_t r_ij = sqrt(SQR(r[j].x - r[i].x) + SQR(r[j].y - r[i].y) + SQR(r[j].z - r[i].z));
            result += p[i].mass * p[j].mass / r_ij;
        }
	}

	return (K2 * result / 2.0);
}

var3_t calc_angular_momentum(uint32_t n, nbp_t::param_t* p, const var3_t* r, const var3_t* v)
{
    var3_t result = {0.0, 0.0, 0.0};
    
    for (uint32_t i = 0; i < n; i++)
    {
        var3_t c = calc_cross_product(r[i], v[i]);

		c.x *= p[i].mass; c.y *= p[i].mass; c.z *= p[i].mass;
		result.x += c.x; result.y += c.y; result.z += c.z;
	}

	return result;
}

var3_t calc_position_of_bc(uint32_t n, const nbp_t::param_t* p, const var3_t* r)
{
    var3_t R0 = {0.0, 0.0, 0.0};

	for (uint32_t i = 0; i < n; i++)
	{
		var_t m_i = p[i].mass;
		R0.x += m_i * r[i].x; R0.y += m_i * r[i].y; R0.z += m_i * r[i].z;
	}
	var_t M0 = get_total_mass(n, p);
	R0.x /= M0; R0.y /= M0; R0.z /= M0;

	return R0;
}

var3_t calc_velocity_of_bc(uint32_t n, const nbp_t::param_t* p, const var3_t* v)
{
    var3_t V0 = {0.0, 0.0, 0.0};

	for (uint32_t i = 0; i < n; i++)
	{
		var_t m_i = p[i].mass;
		V0.x += m_i * v[i].x; V0.y += m_i * v[i].y; V0.z += m_i * v[i].z;
	}
	var_t M0 = get_total_mass(n, p);
	V0.x /= M0; V0.y /= M0; V0.z /= M0;

	return V0;
}

void transform_to_bc(uint32_t n, const nbp_t::param_t* p, var3_t* r, var3_t* v)
{
    var_t M0 = get_total_mass(n, p);
    var3_t R0 = calc_position_of_bc(n, p, r);
    var3_t V0 = calc_velocity_of_bc(n, p, v);

    // Transform the bodies coordinates and velocities
    for (int j = n - 1; j >= 0; j--)
    {
        r[j].x -= R0.x;		r[j].y -= R0.y;		r[j].z -= R0.z;
        v[j].x -= V0.x;		v[j].y -= V0.y;		v[j].z -= V0.z;
    }
}

} /* namespace nbp */


} /* tools */
} /* redutil2 */
