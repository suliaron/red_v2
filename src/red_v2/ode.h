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
        \param omd_size    The size of the metadata in bytes
        \param PROC_UNIT   The name of the executing device
	*/
	ode(uint16_t n_dim, uint32_t n_obj, uint16_t n_vpo, uint16_t n_ppo, size_t omd_size, comp_dev_t comp_dev);
    //! Constructs an ode object
    /*!
    \param n_dim       The space dimension of the problem
    \param n_obj       The total number of objets in the problem
    \param n_vpo       The number of variables per object (vpo)
    \param n_ppo       The number of parameters per object (ppo)
    \param omd_size    The total size of the metadata array in bytes
    \param event_size  The total size of the eventdata array in bytes
    \param PROC_UNIT   The name of the executing device
    */
    ode(uint16_t n_dim, uint32_t n_obj, uint16_t n_vpo, uint16_t n_ppo, size_t omd_size, size_t event_size, comp_dev_t comp_dev);
    //! Constructs an ode object
	/*!
		\param n_dim       The space dimension of the problem 
		\param n_obj       The total number of objets in the problem
		\param n_vpo       The number of variables per object (vpo)
		\param n_ppo       The number of parameters per object (ppo)
		\param n_var       The total number of variables
        \param n_par       The total number of parameters
        \param omd_size    The size of the metadata in bytes
        \param PROC_UNIT   The name of the executing device
	*/
	ode(uint16_t n_dim, uint32_t n_obj, uint16_t n_vpo, uint16_t n_ppo, uint32_t n_var, uint32_t n_par, size_t omd_size, comp_dev_t comp_dev);
	//! Destructor
	~ode();

    //! Set the computing device to calculate the forces for the next step
    /*
    \param device specifies which device will execute the computations
    */
    void set_comp_dev(comp_dev_t cd);
    comp_dev_t get_comp_dev() { return comp_dev; }

    uint32_t get_n_obj()      { return n_obj;    }
    void set_n_obj(uint32_t n);

    uint32_t get_n_var() { return n_var; }
    uint32_t get_n_par() { return n_par; }
    
    void initialize();

	void allocate_storage(       uint32_t n_var, uint32_t n_par, size_t omd_size, size_t event_size);
	void allocate_host_storage(  uint32_t n_var, uint32_t n_par, size_t omd_size, size_t event_size);
	void allocate_device_storage(uint32_t n_var, uint32_t n_par, size_t omd_size, size_t event_size);

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
    //! Copies the metadata between HOST and DEVICE memory
    /*!
    \param dir The direction of the copy
    */
    void copy_metadata(copy_direction_t dir);
    //! Copies the event data between HOST and DEVICE memory
    /*!
    \param dir The direction of the copy
    */
    void copy_eventdata(copy_direction_t dir);

	void swap();

	virtual void calc_dy(uint16_t stage, var_t curr_t, const var_t* y_temp, var_t* acc, var_t* jrk) = 0;
	virtual void calc_dy(uint16_t stage, var_t curr_t, const var_t* y_temp, var_t* dy) = 0;
	virtual void calc_integral() = 0;

    virtual void print_solution(std::string& path_si, std::string& path_sd, data_rep_t repres) = 0;
    virtual void print_dump(std::string& path_si, std::string& path_sd) = 0;
    virtual void print_integral(std::string& path) = 0;

	std::string name;     //! Name of the model

    var_t t_wc;           //!< The wall clock time elapsed since the start of the simlation (measured in sec)

    var_t t;              //! Independent variable (e.g. time or fictitious time)
	var_t tout;           //! Independent variable at the end of the integration step
	var_t dt;             //! Increment/decrement of the independent variable (step-size for the integrator)

	var_t* h_y;           //! Host array (size of n_var) of the dependent variables at t
	var_t* h_yout;        //! Host array (size of n_var) of the dependent variables at tout
	var_t* h_par;         //! Host array (size of n_par) of parameters
    void*  h_omd;         //! Host array (size of omd_size) of metadata
    void*  h_event;       //! Host array (size of event_size) of event

	var_t* d_y;           //! Device array (size of n_var) of the dependent variables at t
	var_t* d_yout;        //! Device array (size of n_var) of the dependent variables at tout
	var_t* d_par;         //! Device array (size of n_par) of parameters
    void*  d_omd;         //! Device array (size of omd_size) of metadata
    void*  d_event;       //! Device array (size of event_size) of event

	var_t* y;             //! Alias to Host or Device array of the dependent variables at t depeding on the execution device
	var_t* yout;          //! Alias to Host or Device array of the dependent variables at tout depeding on the execution device
	var_t* par;           //! Alias to Host or Device array of parameters depeding on the execution device
    void*  omd;           //! Alias to Host or Device array of metadata depeding on the execution device
    void*  event;         //! Alias to Host or Device array of event depeding on the execution device

    integral_t integral;  //! Holds the classical integrals 

protected:
    uint32_t n_obj;       //! The total number of objects in the problem
    uint32_t n_var;       //! The total number of variables of the problem
    uint32_t n_par;       //! The total number of parameters of the problem
    size_t omd_size;      //! The total size of the metadata in bytes
    size_t event_size;    //! The total size of the event array in bytes

    uint16_t n_dim;       //! The space dimension of the problem 
    uint16_t n_vpo;       //! The number of variables per object (vpo)
	uint16_t n_ppo;       //! The number of parameters per object (ppo)

	comp_dev_t comp_dev;  //! The name of the executing device

	dim3 grid;            //! Defines the grid of the blocks of the current execution
	dim3 block;           //! Defines the block of the threads of the current execution
};
