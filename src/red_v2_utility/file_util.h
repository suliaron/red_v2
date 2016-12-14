#pragma once
#include <string>

#include "type.h"

namespace redutil2
{
	namespace file 
	{
		std::string combine_path(const std::string& dir, const std::string& filename);

		std::string get_filename(const std::string& path);
		std::string get_filename_without_ext(const std::string& path);
		std::string get_directory(const std::string& path);
		std::string get_extension(const std::string& path);

		data_rep_t get_data_repres(const std::string& path);

		uint32_t load_ascii_file(const std::string& path, std::string& result);
		void load_binary_file(const std::string& path, size_t n_data, var_t* data);

		void log_start(std::ostream& sout, int argc, const char** argv, const char** env, std::string params, bool print_to_screen);
		void log_message(std::ostream& sout, std::string msg, bool print_to_screen);

		void print_data_info_record_ascii_RED( std::ofstream& sout, var_t t, var_t dt, n_objects_t* n_bodies);
		void print_data_info_record_binary_RED(std::ofstream& sout, var_t t, var_t dt, n_objects_t* n_bodies);

		void print_body_record_ascii_RED( std::ofstream &sout, std::string name, pp_disk_t::param_t *p, pp_disk_t::body_metadata_t *bmd, var4_t *r, var4_t *v);
		void print_body_record_binary_RED(std::ofstream &sout, std::string name, pp_disk_t::param_t *p, pp_disk_t::body_metadata_t *bmd, var4_t *r, var4_t *v);
		void print_body_record_Emese(     std::ofstream &sout, std::string name, var_t epoch, pp_disk_t::param_t *p, pp_disk_t::body_metadata_t *bmd, var4_t *r, var4_t *v);
		void print_body_record_HIPERION(  std::ofstream &sout, std::string name, var_t epoch, pp_disk_t::param_t *p, pp_disk_t::body_metadata_t *bmd, var4_t *r, var4_t *v);

		void print_oe_record(std::ofstream &sout, orbelem_t* oe);
		void print_oe_record(std::ofstream &sout, orbelem_t* oe, pp_disk_t::param_t *p);
		void print_oe_record(std::ofstream &sout, var_t epoch, orbelem_t* oe, pp_disk_t::param_t *p, pp_disk_t::body_metadata_t *bmd);

		void load_data_info_record_ascii( std::ifstream& input, var_t& t, var_t& dt, n_objects_t** n_bodies);
		void load_data_info_record_binary(std::ifstream& input, var_t& t, var_t& dt, n_objects_t** n_bodies);
		void load_data_record_ascii( std::ifstream& input, std::string& name, pp_disk_t::param_t *p, pp_disk_t::body_metadata_t *bmd, var4_t *r, var4_t *v);
		void load_data_record_binary(std::ifstream& input, std::string& name, pp_disk_t::param_t *p, pp_disk_t::body_metadata_t *bmd, var4_t *r, var4_t *v);


		namespace tbp
		{
		void print_solution_info(std::ofstream& sout, var_t t, var_t dt, data_rep_t repres);
        //! Print the solution for each object in text format
        /*   
            \param sout print the data to this stream
        */
        void print_solution_data(std::ofstream& sout, uint32_t n_obj, uint16_t n_ppo, uint16_t n_vpo, tbp_t::metadata_t* md, var_t* p, var_t* y, data_rep_t repres);
		} /* namespace tbp */

		namespace rtbp
		{
        //! Print the solution for each object in text format
        /*   
            \param sout print the data to this stream
        */
        void print_solution_data(std::ofstream& sout, uint32_t n_obj, uint16_t n_ppo, uint16_t n_vpo, tbp_t::metadata_t* md, var_t* p, var_t* y, uint16_t n_dim, data_rep_t repres);
		} /* namespace rtbp */

        namespace nbp
        {
		void print_solution_info(std::ofstream& sout, var_t t, var_t dt, uint32_t n_obj, data_rep_t repres);
        //! Print the solution for each object in text format
        /*   
            \param sout print the data to this stream
        */
        void print_solution_data(std::ofstream& sout, uint32_t n_obj, uint16_t n_ppo, uint16_t n_vpo, nbp_t::metadata_t* md, var_t* p, var_t* y, data_rep_t repres);
        } /* namespace nbp */


	} /* file */
} /* redutil2 */
