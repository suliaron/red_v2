#pragma once
#include <string>
#include <vector>

#include "vector_types.h"

#include "type.h"

class ode
{
public:
	//! Constructs an ode object
	/*!
		\param n_dim       The space dimension of the problem 
		\param n_obj       The total number of objets in the problem
		\param n_vpo       The number of variables per object (vpo)
		\param n_ppo       The number of parameters per object (ppo)
		\param PROC_UNIT    The name of the executing device
	*/
	ode(uint16_t n_dim, uint32_t n_obj, uint16_t n_vpo, uint16_t n_ppo, comp_dev_t comp_dev);
	//! Constructs an ode object
	/*!
		\param n_dim       The space dimension of the problem 
		\param n_obj       The total number of objets in the problem
		\param n_vpo       The number of variables per object (vpo)
		\param n_ppo       The number of parameters per object (ppo)
		\param n_var       The total number of variables
		\param n_par       The total number of parameters
		\param PROC_UNIT    The name of the executing device
	*/
	ode(uint16_t n_dim, uint32_t n_obj, uint16_t n_vpo, uint16_t n_ppo, uint32_t n_var, uint32_t n_par, comp_dev_t comp_dev);
	//! Destructor
	~ode();

	void initialize();

	void allocate_storage(       uint32_t n_var, uint32_t n_par);
	void allocate_host_storage(  uint32_t n_var, uint32_t n_par);
	void allocate_device_storage(uint32_t n_var, uint32_t n_par);

	void deallocate_storage();
	void deallocate_host_storage();
	void deallocate_device_storage();

	void create_aliases();
	//! Copies ODE variables between HOST and DEVICE memory
	/*!
		\param dir The direction of the copy
	*/
	void copy_vars(copy_direction_t dir);
	//! Copies ODE parameters between HOST and DEVICE memory
	/*!
		\param dir The direction of the copy
	*/
	void copy_params(copy_direction_t dir);

	void swap();

	virtual void copy_metadata(copy_direction_t dir) = 0;

	virtual void calc_dy(uint16_t stage, var_t curr_t, const var_t* y_temp, var_t* acc, var_t* jrk) = 0;
	virtual void calc_dy(uint16_t stage, var_t curr_t, const var_t* y_temp, var_t* dy) = 0;
	virtual void calc_integral() = 0;

    virtual void print_solution(std::string& path_si, std::string& path_sd, data_rep_t repres) = 0;
    virtual void print_dump(std::string& path_si, std::string& path_sd) = 0;
    virtual void print_integral(std::string& path) = 0;

	std::string name;     //! Name of the model

    bool optimize;        //!< If true than perform optimization

    var_t t_wc;           //!< The wall clock time elapsed since the start of the simlation (measured in sec)

    var_t t;              //! Independent variable (e.g. time or fictitious time)
	var_t tout;           //! Independent variable at the end of the integration step
	var_t dt;             //! Increment/decrement of the independent variable (step-size for the integrator)

	var_t* h_y;           //! Host vector (size of n_var) of the dependent variables at t
	var_t* h_yout;        //! Host vector (size of n_var) of the dependent variables at tout
	var_t* h_p;           //! Host vector (size of n_par) of parameters

	var_t* d_y;           //! Device vector (size of n_var) of the dependent variables at t
	var_t* d_yout;        //! Device vector (size of n_var) of the dependent variables at tout
	var_t* d_p;           //! Device vector (size of n_par) of parameters

	var_t* y;             //! Alias to Host or Device vector of the dependent variables at t depeding on the execution device
	var_t* yout;          //! Alias to Host or Device vector of the dependent variables at tout depeding on the execution device
	var_t* p;             //! Alias to Host or Device vector of parameters depeding on the execution device

	uint16_t n_dim;       //! The space dimension of the problem 
	uint32_t n_obj;       //! The total number of objets in the problem

	uint16_t n_vpo;       //! The number of variables per object (vpo)
	uint16_t n_ppo;       //! The number of parameters per object (ppo)

	uint32_t n_var;       //! The total number of variables of the problem
	uint32_t n_par;       //! The total number of parameters of the problem

	comp_dev_t comp_dev;  //! The name of the executing device

	dim3 grid;            //! Defines the grid of the blocks of the current execution
	dim3 block;           //! Defines the block of the threads of the current execution
	uint16_t n_tpb;       //! Holds the number of threads per block

	integral_t integral;  //! Holds the classical integrals 
};
