#pragma once
#include <string>

#include "type.h"

namespace redutil2
{
	namespace tools
	{
		/// Default white-space characters
		static const char* ws = " \t\n\r\f\v";
		/// Default comment character
		static const char* comment = "#";

		bool is_number(const std::string& str);

        //! Removes all leading white-space characters from the current string object. The default white spaces are: " \t\n\r\f\v"
        /*
            \param s the string to process
        */
        std::string& ltrim(std::string& s);
        //! Removes all leading white-space characters from the current string object.
        /*
            \param s the string to process
            \param t holds the white spaces
        */
        std::string& ltrim(std::string& s, const char* t);

        //! Removes all trailing white-space characters from the current string object. The default white spaces are: " \t\n\r\f\v"
        /*
            \param s the string to process
        */
        std::string& rtrim(std::string& s);
        //! Removes all trailing white-space characters from the current string object.
        /*
            \param s the string to process
            \param t holds the white spaces
        */
        std::string& rtrim(std::string& s, const char* t);

        //! Removes all leading and trailing white-space characters from the current string object. The default white spaces are: " \t\n\r\f\v"
        /*
            \param s the string to process
        */
        std::string& trim(std::string& s);
        //! Removes all leading and trailing white-space characters from the current string object.
        /*
            \param s the string to process
            \param t holds the white spaces
        */
        std::string& trim(std::string& s, const char* t);

        //! Removes all trailing characters after the first comment character from the current string object. The default comment character is "#"
        /*
            \param s the string to process
        */
        std::string& trim_comment(std::string& s);
        //! Removes all trailing characters after the first comment character from the current string object. 
        /*
            \param s the string to process
            \param c holds the comment character
        */
        std::string& trim_comment(std::string& s, const char* c);

		std::string get_time_stamp();
		std::string convert_time_t(time_t t);
        std::string convert_var_t(var_t v);

		void populate_data(uint32_t* n_bodies, pp_disk_t::sim_data_t *sim_data);

		//! Computes the total mass of the system
		var_t get_total_mass(uint32_t n, const pp_disk_t::sim_data_t *sim_data);
		//! Computes the total mass of the bodies with type in the system
		var_t get_total_mass(uint32_t n, body_type_t type, const pp_disk_t::sim_data_t *sim_data);
		void calc_bc(uint32_t n, const pp_disk_t::sim_data_t *sim_data, var_t M, var3_t* R0, var3_t* V0);
		void transform_to_bc(uint32_t n, const pp_disk_t::sim_data_t *sim_data);
		void transform_time( uint32_t n, const pp_disk_t::sim_data_t *sim_data);
		void transform_velocity(uint32_t n, const pp_disk_t::sim_data_t *sim_data);

		var_t calc_radius(var_t m, var_t density);
		var_t calc_density(var_t m, var_t R);
		var_t calc_mass(var_t R, var_t density);

		void calc_position_after_collision(var_t m1, var_t m2, const var3_t* r1, const var3_t* r2, var3_t& r);
		void calc_velocity_after_collision(var_t m1, var_t m2, const var3_t* v1, const var3_t* v2, var3_t& v);
		void calc_physical_properties(const pp_disk_t::param_t &p1, const pp_disk_t::param_t &p2, pp_disk_t::param_t &p);

		var_t norm(const var4_t* r);
		var_t norm(const var3_t* r);
		var_t calc_dot_product(const var3_t& u, const var3_t& v);
		var3_t calc_cross_product(const var3_t& u, const var3_t& v);
		var_t calc_kinetic_energy(const var3_t* v);
		var_t calc_pot_energy(var_t mu, const var3_t* r);

		matrix4_t calc_matrix_matrix_product(const matrix4_t& u, const matrix4_t& v);
		var4_t calc_matrix_vector_product(const matrix4_t& u, const var4_t& v);
		matrix4_t calc_matrix_transpose(const matrix4_t& u);

		void trans_to_descartes(const var_t m1, const var_t m2, const var_t m3, var3_t& q1, var3_t& p1, var3_t& q2, var3_t& p2, var3_t& q3, var3_t& p3, const var4_t& Q1, const var4_t& P1, const var4_t& Q2, const var4_t& P2);
		void trans_to_threebody(const var3_t& qv1, const var3_t& pv1, const var3_t& qv2, const var3_t& pv2, const var3_t& qv3, const var3_t& pv3, var4_t& Q1, var4_t& P1, var4_t& Q2, var4_t& P2);

        var_t calc_total_energy(         uint32_t n, const pp_disk_t::sim_data_t *sim_data);
        var_t calc_total_energy_CMU(     uint32_t n, const pp_disk_t::sim_data_t *sim_data);
        var3_t calc_angular_momentum(    uint32_t n, const pp_disk_t::sim_data_t *sim_data);
        var3_t calc_angular_momentum_CMU(uint32_t n, const pp_disk_t::sim_data_t *sim_data);
       
		void kepler_equation_solver(var_t ecc, var_t mean, var_t eps, var_t* E);
		void calc_phase(var_t mu, const orbelem_t* oe, var3_t* rVec, var3_t* vVec);
		void calc_oe(   var_t mu, const var3_t* rVec, const var3_t* vVec, orbelem_t* oe);
		var_t calc_orbital_period(var_t mu, var_t a);

		//! Calculate the linear combination of vectors, a[i] = sum (coeff[j] * c[j][i])
		/*
			\param a     vector which will contain the result
			\param c     vectors which will linear combined
			\param coeff vector which contains the weights (coefficients)
			\param n_vct the number of vectors to combine
			\param n_var the number of elements in the vectors
		*/
		void calc_lin_comb(var_t* a, const var_t* const *c, const var_t* coeff, uint16_t n_vct, uint32_t n_var);

		//! Calculate the special linear combination of two vectors, a[i] = b[i] + f*c[i]
		/*
			\param a     vector which will contain the result
			\param b     vector to which the linear combination will be added
			\param c     vector to add with weight f
			\param f     the weight of c
			\param n_var the number of elements in the vectors
		*/
		void calc_lin_comb_s(var_t* a, const var_t* b, const var_t* c, var_t f, uint32_t n_var);

		//! Calculate the special case of linear combination of vectors, a[i] = b[i] + sum (coeff[j] * c[j][i])
		/*
			\param a     vector which will contain the result
			\param b     vector to which the linear combination will be added
			\param c     vectors which will linear combined
			\param coeff vector which contains the weights (coefficients)
			\param n_vct the number of vectors to combine
			\param n_var the number of elements in the vectors
		*/
		void calc_lin_comb_s(var_t* a, const var_t* b, const var_t* const *c, const var_t* coeff, uint16_t n_vct, uint32_t n_var);

		void print_vector(const var3_t *v);
		void print_vector(const var4_t *v);
		void print_parameter(const pp_disk_t::param_t *p);
		void print_body_metadata(const pp_disk_t::body_metadata_t *b);

		namespace tbp
		{
		//! Calculate the energy
		/*!
			\param mu  the parameter of the problem
			\param r   the physical position vector of the object
			\param v   the physical velocity vector of the object
		*/
		var_t calc_integral(var_t mu, var2_t r, var2_t v);
		} /* tbp */

		namespace rtbp1D
		{
		//! Calculate parametric coordinate and velocity
		/*!
			\param x   the physical coordinate of the object
			\param vx  the physical velocity of the object
			\param u   the parametric position of the object
			\param v   the parametric velocity of the object
		*/
		void transform_x2u(var_t x, var_t vx, var_t& u, var_t& v);
		//! Calculate physical coordinate and velocity
		/*!
			\param u   the parametric position of the object
			\param v   the parametric velocity of the object
			\param x   the physical coordinate of the object
			\param vx  the physical velocity of the object
		*/
		void transform_u2x(var_t u, var_t v, var_t& x, var_t& vx);
		} /* namespace rtbp1D */

		namespace rtbp2D
		{
		//! Calculate parametric coordinate from physical coordinate
		/*!
			\param x   the (x1, x2) physical coordinates of the object
			\param u   the parametric (u1, u2) position of the object
		*/
		void transform_x2u(var2_t x, var2_t& u);

		//! Calculate parametric velocity from physical velocity
		/*!
			\param u   the parametric (u1, u2) position of the object
			\param xd  the (xd1, xd2) physical velocity of the object (xd = x dot)
			\param up  the parametric velocity of the object (up = u prime)
		*/
		void transform_xd2up(var2_t u, var2_t xd, var2_t& up);

		void transform_u2x(var2_t u, var2_t& x);
		void transform_up2xd(var2_t u, var2_t up, var2_t& xd);
		
		//! Calculate the energy for the regularized problem
		/*!
			\param mu  the parameter of the problem
			\param u   the parametric (u1, u2) position of the object
			\param up  the parametric velocity of the object (up = u prime)
		*/
		var_t calc_integral(var_t mu, var2_t u, var2_t up);
		} /* namespace rtbp2D */

		namespace nbp
		{
			var_t get_total_mass(uint32_t n, nbp_t::param_t* p);
			var_t calc_total_energy(uint32_t n, nbp_t::param_t* p, var3_t* r, var3_t* v);
			var_t calc_kinetic_energy(uint32_t n, nbp_t::param_t* p, var3_t* v);
			var_t calc_potential_energy(uint32_t n, nbp_t::param_t* p, var3_t* r);
			var3_t calc_angular_momentum(uint32_t n, nbp_t::param_t* p, var3_t* r, var3_t* v);
			var3_t calc_position_of_bc(uint32_t n, nbp_t::param_t* p, var3_t* r);
			var3_t calc_velocity_of_bc(uint32_t n, nbp_t::param_t* p, var3_t* v);
            void transform_to_bc(uint32_t n, nbp_t::param_t* p, var3_t* r, var3_t* v);
		} /* namespace nbp */

	} /* namespace tools */
} /* namespace redutil2 */
