#pragma once

#include "type.h"

namespace redutil2
{
    uint64_t GetTimeMs64();

    __host__ __device__ var4_t rotate_2D_vector(var_t theta, const var4_t& r);

	template <typename T>
	std::string number_to_string( T number, uint32_t width, bool fill);
	template <typename T>
	std::string number_to_string(T number);

	void device_query(std::ostream& sout, int id_dev, bool print_to_screen);

	int get_id_fastest_cuda_device();
	int get_n_cuda_device();
	std::string get_device_name(int id_dev);
    void set_device(int id_of_target_dev, std::ostream& sout);

	void set_kernel_launch_param(uint32_t n_data, uint16_t n_tpb, dim3& grid, dim3& block);

	void allocate_host_array(  void **ptr, size_t size,           const char *file, int line);
	void allocate_device_array(void **ptr, size_t size,           const char *file, int line);
    void free_host_array(void **ptr, const char *file, int line);
    void free_device_array(void **ptr, const char *file, int line);

	#define ALLOCATE_HOST_ARRAY(  ptr, size)    (allocate_host_array(  ptr, size, __FILE__, __LINE__))
	#define ALLOCATE_DEVICE_ARRAY(ptr, size)    (allocate_device_array(ptr, size, __FILE__, __LINE__))
	#define FREE_HOST_ARRAY(  ptr)              (free_host_array(  ptr, __FILE__, __LINE__))
	#define FREE_DEVICE_ARRAY(ptr)              (free_device_array(ptr, __FILE__, __LINE__))
	
	void copy_array_to_device(void* dst, const void *src, size_t count);
	void copy_array_to_host(  void* dst, const void *src, size_t count);
	void copy_array_d2d(      void* dst, const void *src, size_t count);

	void copy_constant_to_device(const void* dst, const void *src, size_t count);

    void prn(const var_t* a);
    void prn(const var3_t* a);
    void prn(const nbp_t::param_t* a);
    void prn(const nbp_t::metadata_t* a);

	void print_array(std::string path, std::string comment, uint32_t n, var_t *data, mem_loc_t);

	//! Calculate the special linear combination of two vectors, a[i] = b[i] + f*c[i]
	/*
		\param a     vector which will contain the result
		\param b     vector to which the linear combination will be added
		\param c     vector to add with weight f
		\param f     the weight of c
		\param n_var the number of elements in the vectors
		\param id_dev the id of the GPU to use for the computation
		\param optimize If true than perform optimization of the underlying kernel
	*/
	void gpu_calc_lin_comb_s(var_t* a, const var_t* b, const var_t* c, var_t f, uint32_t n_var, int id_dev);

	//! Calculate the special case of linear combination of vectors, a[i] = b[i] + sum (coeff[j] * c[j][i])
	/*
		\param a     vector which will contain the result
		\param b     vector to which the linear combination will be added
		\param c     vectors which will linear combined
		\param coeff vector which contains the weights (coefficients)
		\param n_vct the number of vectors to combine
		\param n_var the number of elements in the vectors
		\param id_dev the id of the GPU to use for the computation
	*/
	void gpu_calc_lin_comb_s(var_t* a, const var_t* b, const var_t* const *c, const var_t* coeff, uint16_t n_vct, uint32_t n_var, int id_dev);

    //! Calculate the error for the Runge-Kutta 4 method: a = |k4 - k5|
    /*
        \param a will hold the result
        \param k4 input vector
        \param k5 input vector
		\param n_var the number of elements in the vectors
		\param id_dev the id of the GPU to use for the computation
    */
    void gpu_calc_rk4_error(var_t* a, const var_t* k4, const var_t* k5, uint32_t n_var, int id_dev);

    //! Calculate the error for the Runge-Kutta 5 method: a = |k5 - k6|
    /*
        \param a will hold the result
        \param k5 input vector
        \param k6 input vector
		\param n_var the number of elements in the vectors
		\param id_dev the id of the GPU to use for the computation
    */
    void gpu_calc_rk5_error(var_t* a, const var_t* k5, const var_t* k6, uint32_t n_var, int id_dev);

    //! Calculate the error for the Runge-Kutta 4 method: a = |k1 + k11 - k12 - k13|
    /*
        \param a will hold the result
        \param k0 input vector
        \param k10 input vector
        \param k11 input vector
        \param k12 input vector
		\param n_var the number of elements in the vectors
		\param id_dev the id of the GPU to use for the computation
    */
	void gpu_calc_rk7_error(var_t* a, const var_t* k1, const var_t* k11, const var_t* k12, const var_t* k13, uint32_t n_var, int id_dev);
} /* redutil2 */
