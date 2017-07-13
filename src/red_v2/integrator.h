#pragma once

#include <string>
#include <vector>

#include "type.h"

class ode;

class integrator
{
public:
	integrator(ode& f, bool adaptive, var_t tolerance, uint16_t n_stage, comp_dev_t comp_dev);
	~integrator();

	//! Set the computing device to calculate the integration step
	/*
		\param device specifies which device will execute the computations
	*/
    void set_comp_dev(comp_dev_t cd);
    comp_dev_t get_comp_dev()      { return comp_dev;      }

	void update_counters(uint16_t iter);

	uint64_t get_n_tried_step()    { return n_tried_step;  }
	uint64_t get_n_passed_step()   { return n_passed_step; }
	uint64_t get_n_failed_step()   { return n_failed_step; }
    var_t get_dt_did()             { return dt_did;        }

	void set_max_iter(uint16_t n)  { max_iter = n;         }
	uint16_t get_max_iter()        { return max_iter;      }

	void set_dt_min(var_t dt)      { dt_min = dt;          }
	var_t get_dt_min()             { return dt_min;        }

	virtual void allocate_Butcher_tableau()   = 0;
	virtual void deallocate_Butcher_tableau() = 0;
	virtual void check_Butcher_tableau() = 0;
	virtual var_t step() = 0;

	ode& f;

    bool optimize;               //!< If true than perform optimization
	bool error_check_for_tp;	 //!< Check the error also for the test particles
	std::string name;

protected:
	var_t get_max_error(uint32_t n_var);
	void calc_dt_try(var_t max_err);

	comp_dev_t comp_dev;         //!< The computing device to carry out the calculations (cpu or gpu)

	dim3 grid;
	dim3 block;

	var_t t;					 //!< Independent variable (e.g. time or fictitious time)
	var_t dt_try;                //!< The size of the step to try (based on the previous successfull step dt_did)
	var_t dt_did;                //!< The size of the previous successfull step

	uint64_t n_tried_step;
	uint64_t n_passed_step;
	uint64_t n_failed_step;                                                                                                                                                                                         
	bool adaptive;               //!< True if the method estimates the error and accordingly adjusts the step-size	
	var_t tolerance;             //!< The maximum of the allowed local truncation error
	uint16_t n_order;            //!< The order of the embedded RK formulae
	uint16_t n_stage;            //!< The number of the method's stages

	var_t** h_k;                 //!< Differentials in the HOST memory
	var_t** d_k;                 //!< Differentials in the DEVICE memory
	var_t** k;                   //!< Alias to the differentials (either in the HOST or the DEVICE memory)
	var_t** cpy_dk;              //!< Holds a copy of the allocated DEVICE vectors address in the HOST memory

	var_t* h_ytemp;	             //!< Holds the temporary solution approximation along the step in the HOST memory
	var_t* d_ytemp;	             //!< Holds the temporary solution approximation along the step in the DEVICE memory
	var_t* ytemp;	             //!< Alias either to h_ytemp or d_ytemp (either in the HOST or the DEVICE memory)

	var_t* h_err;	             //!< Holds the leading local truncation error for each variable in HOST memory
	var_t* d_err;	             //!< Holds the leading local truncation error for each variable in DEVICE memory
	var_t* err;  	             //!< Alias either to h_err or d_err (either in the HOST or the DEVICE memory)

	uint16_t max_iter;
	var_t dt_min;

private:
	void initialize();
	void create_aliases();

	//! Allocates storage for data on the host and device memory
	void allocate_storage(       uint32_t n_var);
	void allocate_host_storage(  uint32_t n_var);
	void allocate_device_storage(uint32_t n_var);

	void deallocate_storage();
	void deallocate_host_storage();
	void deallocate_device_storage();
};
