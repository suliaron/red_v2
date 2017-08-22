#pragma once

#include "ode.h"
#include "type.h"

class nbody : public ode
{
public:
	//! Constructs a collisionless nbody object
	/*!
        \param path_si   the path of the file which conatins data about the initial conditions
        \param path_sd   the path of the file which conatins the initial conditions
        \param n_obj     the number of objects
		\param n_ppo     the number of parameters per object
        \param omd_size  the size of the metadata in bytes
        \param PROC_UNIT the name of the executing device
	*/
	nbody(std::string& path_si, std::string& path_sd, uint32_t n_obj, uint16_t n_ppo, size_t omd_size, comp_dev_t comp_dev);
    //! Constructs a collisional nbody object
    /*!
    \param path_si    the path of the file which conatins data about the initial conditions
    \param path_sd    the path of the file which conatins the initial conditions
    \param n_obj      the number of objects
    \param n_ppo      the number of parameters per object
    \param omd_size   the size of the metadata in bytes
    \param event_size the total size of the event data in bytes
    \param PROC_UNIT the name of the executing device
    */
    nbody(std::string& path_si, std::string& path_sd, uint32_t n_obj, uint16_t n_ppo, size_t omd_size, size_t event_size, comp_dev_t comp_dev);
    //! Destructor
	~nbody();

	//! Print the solution (the numerical approximation of the solution)
	/*!
		\param path_si  full file path where the info about the solution will be printed
		\param path_sd  full file path where the the solution will be printed
		\param repres   indicates the data representation of the file, i.e. text or binary
	*/
	void print_solution(std::string& path_si, std::string& path_sd, data_rep_t repres);
    void print_dump(std::string& path_si, std::string& path_sd);

	//! Print the integral(s) of the problem
	/*!
		\param path   full file path where the integrals of the problem will be printed
	*/
	void print_integral(std::string& path);
    void print_event_data(std::string& path, uint32_t n_event);

    void calc_n_types();
    uint32_t calc_n_active();

    bool get_print_oe(void)   { return print_oe; }
    void set_print_oe(bool b) { print_oe = b;    }

    uint32_t get_n_coll()     { return n_coll; }
    uint32_t get_n_hitc()     { return n_hitc; }
    uint32_t get_n_ejec()     { return n_ejec; }

    uint32_t get_n_si()       { return n_si;   }
    uint32_t get_n_nsi()      { return n_nsi;  }
    uint32_t get_n_ni()       { return n_ni;   }

    void calc_dy(uint16_t stage, var_t curr_t, const var_t* y_temp, var_t* acc, var_t* jrk);
	void calc_dy(uint16_t stage, var_t curr_t, const var_t* y_temp, var_t* dy);
	void calc_integral();
   
    uint32_t chk_event(std::string path, var_t* threshold);
    void rebuild_host_array();

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

    void cpu_calc_dy(uint16_t stage, var_t curr_t, const var_t* y_temp, var_t* acc, var_t* jrk, bool use_symm_prop);
    void cpu_calc_dy(uint16_t stage, var_t curr_t, const var_t* y_temp, var_t* dy);
    void body_body_grav_accel(const var3_t& ri, const var3_t& rj, var_t mj, var3_t& ai);
    void body_body_grav_accel(const var3_t& ri, const var3_t& rj, var_t mi, var_t mj, var3_t& ai, var3_t& aj);

    void gpu_calc_dy(uint16_t stage, var_t curr_t, const var_t* y_temp, var_t* dy);
    uint16_t calc_opt_tpb(uint16_t stage, uint2_t snk, uint2_t src, var_t curr_t, const var3_t* r, const nbp_t::param_t* p);
    //float gpu_calc_grav_accel_naive(uint16_t stage, uint2_t snk, uint2_t src, unsigned int n_tpb, var_t curr_t, const var3_t* r, const nbp_t::param_t* p, var3_t* a);
    float gpu_calc_grav_accel_tile(uint16_t stage, uint2_t snk, uint2_t src, unsigned int n_tpb, var_t curr_t, const var3_t* r, const nbp_t::param_t* p, var3_t* a);

    uint32_t cpu_chk_event(var_t* threshold);
    uint32_t gpu_chk_event(var_t* threshold);
    void handle_event(uint32_t n_event);
    void handle_collision_pair(nbp_t::event_data_t *collision);
    void handle_ejection(nbp_t::event_data_t *collision);

    bool print_oe;           //! Print the orbital elements of the bodies
    nbp_t::metadata_t *h_md, *d_md, *md;
    nbp_t::event_data_t *h_ed, *d_ed, *ed;

    uint32_t* d_n_event;     //! The number of events detected on the GPU by chk_coll() or ....

    uint32_t n_si;           //! The total number of SI type bodies
    uint32_t n_nsi;          //! The total number of NSI type bodies
    uint32_t n_ni;           //! The total number of NI type bodies

    uint16_t n_tpb_si;       //! Holds the number of threads per block for the SI type bodies
    uint16_t n_tpb_nsi;      //! Holds the number of threads per block for the NSI type bodies
    uint16_t n_tpb_ni;       //! Holds the number of threads per block for the NI type bodies

    uint32_t n_coll;         //! The total number of collisions
    uint32_t n_hitc;         //! The total number of hit centrum
    uint32_t n_ejec;         //! The total number of ejections

};
