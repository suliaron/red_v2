#include <algorithm>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <ostream>
#include <sstream>
#include <string>

#ifdef _WIN32
#include <chrono>
#include <Windows.h>
#else
#include <sys/time.h>
#include <ctime>
#endif

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include "util.h"
#include "file_util.h"
#include "type.h"
#include "macro.h"

using namespace std;

namespace red_kernel
{
//! Print pointer and dereferennce it
__global__
void print_ptr(const var_t* a)
{
	uint32_t tid = blockIdx.x * blockDim.x + threadIdx.x;

	if (0 == tid)
	{
		printf("%p\t%25.16le\n", a, *a);
	}
}

__global__
void print_ptr(const var_t* const *a)
{
	uint32_t tid = blockIdx.x * blockDim.x + threadIdx.x;

	if (0 == tid)
	{
		printf("%p\t%p\n", a, *a);
	}
}
	
//! Calculate the special linear combination of two vectors, a[i] = b[i] + f*c[i]
__global__
void calc_lin_comb_s(var_t* a, const var_t* b, var_t f, const var_t* c, uint32_t n)
{
	uint32_t tid = blockIdx.x * blockDim.x + threadIdx.x;
	uint32_t stride = gridDim.x * blockDim.x;

	while (n > tid)
	{
		a[tid] = b[tid] + f * c[tid];
		tid += stride;
	}
}

//! Calculate the special case of linear combination of vectors, a[i] = b[i] + sum (coeff[j] * c[j][i])
__global__
void calc_lin_comb_s(var_t* a, const var_t* b, const var_t* const *c, const var_t* coeff, uint16_t n_vct, uint32_t n_var)
{
	const int tid = blockIdx.x * blockDim.x + threadIdx.x;

	if (tid < n_var)
	{
		var_t d = 0.0;
		for (uint16_t j = 0; j < n_vct; j++)
		{
			if (0.0 == coeff[j])
			{
				continue;
			}
			d += coeff[j] * c[j][tid];
		}
		a[tid] = b[tid] + d;
	}
}

//! Calculate the error for the Runge-Kutta 4 method: error = |k4 - k5|
__global__
void calc_rk4_error(var_t* a, const var_t* k4, const var_t* k5, uint32_t n)
{
	uint32_t tid = blockIdx.x * blockDim.x + threadIdx.x;
	uint32_t stride = gridDim.x * blockDim.x;

	while (n > tid)
	{
        a[tid] = fabs(k4[tid] - k5[tid]);
		tid += stride;
	}
}

//! Calculate the error for the Runge-Kutta 5 method: error = |k5 - k6|
__global__
void calc_rk5_error(var_t* a, const var_t* k5, const var_t* k6, uint32_t n)
{
	uint32_t tid = blockIdx.x * blockDim.x + threadIdx.x;
	uint32_t stride = gridDim.x * blockDim.x;

	while (n > tid)
	{
        a[tid] = fabs(k5[tid] - k6[tid]);
		tid += stride;
	}
}

//! Calculate the error for the Runge-Kutta 7 method: error = |k1 + k11 - k12 - k13|
__global__
void calc_rk7_error(var_t* a, const var_t* k1, const var_t* k11, const var_t* k12, const var_t* k13, uint32_t n)
{
	uint32_t tid = blockIdx.x * blockDim.x + threadIdx.x;
	uint32_t stride = gridDim.x * blockDim.x;

	while (n > tid)
	{
        a[tid] = fabs(k1[tid] + k11[tid] - k12[tid] - k13[tid]);
		tid += stride;
	}
}
} /* red_kernel */

namespace redutil2
{
/*
* Returns the amount of microseconds elapsed since the UNIX epoch.
* Works on windows and linux.
*/
uint64_t GetTimeMs64()
{
#ifdef _WIN32
    FILETIME ft;
    LARGE_INTEGER li;

    /* Get the amount of 100 nano seconds intervals elapsed since January 1, 1601 (UTC)
    * and copy it to a LARGE_INTEGER structure.
    */
    GetSystemTimeAsFileTime(&ft);
    li.LowPart = ft.dwLowDateTime;
    li.HighPart = ft.dwHighDateTime;

    uint64_t ret = li.QuadPart;
    /* Convert from file time to UNIX epoch time. */
    ret -= 116444736000000000LL;
    ret /= 10;      /* From 100 nano seconds (10^-7) to 1 microsecond (10^-6) intervals */

    return ret;
#else
    // Linux
    struct timeval tv;

    gettimeofday(&tv, NULL);
    uint64_t ret = tv.tv_usec;
    /* Convert from micro seconds (10^-6) to milliseconds (10^-3) */
    //ret /= 1000;

    /* Adds the seconds (10^0) after converting them to microseconds (10^-6) */
    ret += (tv.tv_sec * 1000000);

    return ret;
#endif
}

    template <typename T>
std::string number_to_string( T number, uint32_t width, bool fill)
{
	std::ostringstream ss;

	if (fill)
	{
		if (0 < width)
		{
			ss << setw(width) << setfill('0') << number;
		}
		else
		{
			ss << setfill('0') << number;
		}
	}
	else
	{
		if (0 < width)
		{
			ss << setw(width) << number;
		}
		else
		{
			ss << number;
		}
	}

	return ss.str();
}

template std::string number_to_string<char>(                  char, uint32_t width, bool fill);
template std::string number_to_string<unsigned char>(unsigned char, uint32_t width, bool fill);
template std::string number_to_string<int>(                    int, uint32_t width, bool fill);
template std::string number_to_string<uint32_t>(  uint32_t, uint32_t width, bool fill);
template std::string number_to_string<long>(                  long, uint32_t width, bool fill);
template std::string number_to_string<unsigned long>(unsigned long, uint32_t width, bool fill);

template <typename T>
std::string number_to_string( T number )
{
	std::ostringstream ss;
	ss << number;
	return ss.str();
}

template std::string number_to_string<char>(char);
template std::string number_to_string<unsigned char>(unsigned char);
template std::string number_to_string<int>(int);
template std::string number_to_string<uint32_t>(uint32_t);
template std::string number_to_string<long>(long);
template std::string number_to_string<unsigned long>(unsigned long);
template std::string number_to_string<float>(float);
template std::string number_to_string<double>(double);

__host__ __device__
	var4_t rotate_2D_vector(var_t theta, const var4_t& r)
{
	var_t ct = cos(theta);
	var_t st = sin(theta);

	var4_t result = {ct * r.x - st * r.y, st * r.x + ct * r.y, 0.0, 0.0};
	return result;
}

// Beginning of GPU Architecture definitions
inline int _ConvertSMVer2Cores(int major, int minor)
{
    // Defines for GPU Architecture types (using the SM version to determine the # of cores per SM
    typedef struct
    {
        int SM; // 0xMm (hexidecimal notation), M = SM Major version, and m = SM minor version
        int Cores;
    } sSMtoCores;

    sSMtoCores nGpuArchCoresPerSM[] =
    {
        { 0x10,  8 }, // Tesla Generation (SM 1.0) G80 class
        { 0x11,  8 }, // Tesla Generation (SM 1.1) G8x class
        { 0x12,  8 }, // Tesla Generation (SM 1.2) G9x class
        { 0x13,  8 }, // Tesla Generation (SM 1.3) GT200 class
        { 0x20, 32 }, // Fermi Generation (SM 2.0) GF100 class
        { 0x21, 48 }, // Fermi Generation (SM 2.1) GF10x class
        { 0x30, 192}, // Kepler Generation (SM 3.0) GK10x class
        { 0x35, 192}, // Kepler Generation (SM 3.5) GK11x class
        {   -1, -1 }
    };

    int index = 0;
    while (nGpuArchCoresPerSM[index].SM != -1)
    {
        if (nGpuArchCoresPerSM[index].SM == ((major << 4) + minor))
        {
            return nGpuArchCoresPerSM[index].Cores;
        }
        index++;
    }

    // If we don't find the values, we default use the previous one to run properly
    printf("MapSMtoCores for SM %d.%d is undefined.  Default to use %d Cores/SM\n", major, minor, nGpuArchCoresPerSM[7].Cores);
    return nGpuArchCoresPerSM[7].Cores;
}
// end of GPU Architecture definitions

string get_device_name(int id_dev)
{
	cudaDeviceProp deviceProp;
    CUDA_SAFE_CALL(cudaGetDeviceProperties(&deviceProp, id_dev));

	string result(deviceProp.name);
	return result;
}

// TODO: implement
int get_id_fastest_cuda_device()
{
	return 0;
}

int get_n_cuda_device()
{
	int n_device = 0;
	CUDA_SAFE_CALL(cudaGetDeviceCount(&n_device));
	return n_device;
}

void device_query(ostream& sout, int id_dev)
{
    int deviceCount = 0;
    CUDA_SAFE_CALL(cudaGetDeviceCount(&deviceCount));

    // This function call returns 0 if there are no CUDA capable devices.
    if (deviceCount == 0)
    {
		throw string("There are no available device(s) that support CUDA.");
    }

    int dev, driverVersion = 0, runtimeVersion = 0;

    cudaDeviceProp deviceProp;
    CUDA_SAFE_CALL(cudaGetDeviceProperties(&deviceProp, id_dev));

    sout << "The code runs on a " << deviceProp.name << " device:" << endl;

    // Console log
    cudaDriverGetVersion(&driverVersion);
    cudaRuntimeGetVersion(&runtimeVersion);
	sout << "  CUDA Driver Version / Runtime Version          " << driverVersion/1000 << "." << (driverVersion%100)/10 << " / " << runtimeVersion/1000 << "." << (runtimeVersion%100)/10 << endl;
	sout << "  CUDA Capability Major/Minor version number:    " << deviceProp.major << "." << deviceProp.minor << endl;
	sout << "  Total amount of global memory:                 " << deviceProp.totalGlobalMem/1048576.0f << " MBytes" << endl;
	sout << "  " << deviceProp.multiProcessorCount <<  " Multiprocessors, " << _ConvertSMVer2Cores(deviceProp.major, deviceProp.minor) << " CUDA Cores/MP:     " << _ConvertSMVer2Cores(deviceProp.major, deviceProp.minor) * deviceProp.multiProcessorCount << " CUDA Cores" << endl;
	sout << "  GPU Clock rate:                                " << deviceProp.clockRate * 1e-3f << ".0 MHz" << endl;
	sout << "  Total amount of constant memory:               " << deviceProp.totalConstMem << " bytes" << endl;
	sout << "  Total amount of shared memory per block:       " << deviceProp.sharedMemPerBlock << " bytes" << endl;
	sout << "  Total number of registers available per block: " << deviceProp.regsPerBlock << endl;
	sout << "  Warp size:                                     " << deviceProp.warpSize << endl;
	sout << "  Maximum number of threads per multiprocessor:  " << deviceProp.maxThreadsPerMultiProcessor << endl;
	sout << "  Maximum number of threads per block:           " << deviceProp.maxThreadsPerBlock << endl;
	sout << "  Max dimension size of a thread block (x,y,z): (" << deviceProp.maxThreadsDim[0] << "," << deviceProp.maxThreadsDim[1] << "," << deviceProp.maxThreadsDim[2] << ")" << endl;
	sout << "  Max dimension size of a grid size    (x,y,z): (" <<deviceProp.maxGridSize[0] << "," << deviceProp.maxGridSize[1] << "," << deviceProp.maxGridSize[2] << ")" << endl << endl;

    std::string sProfileString = "deviceQuery, CUDA Driver = CUDART";
    char cTemp[16];

    // driver version
    sProfileString += ", CUDA Driver Version = ";
#ifdef WIN32
    sprintf_s(cTemp, 10, "%d.%d", driverVersion/1000, (driverVersion%100)/10);
#else
    sprintf(cTemp, "%d.%d", driverVersion/1000, (driverVersion%100)/10);
#endif
    sProfileString +=  cTemp;

    // Runtime version
    sProfileString += ", CUDA Runtime Version = ";
#ifdef WIN32
    sprintf_s(cTemp, 10, "%d.%d", runtimeVersion/1000, (runtimeVersion%100)/10);
#else
    sprintf(cTemp, "%d.%d", runtimeVersion/1000, (runtimeVersion%100)/10);
#endif
    sProfileString +=  cTemp;

    // Device count
    sProfileString += ", NumDevs = ";
#ifdef WIN32
    sprintf_s(cTemp, 10, "%d", deviceCount);
#else
    sprintf(cTemp, "%d", deviceCount);
#endif
    sProfileString += cTemp;

    // Print Out all device Names
    for (dev = 0; dev < deviceCount; ++dev)
    {
#ifdef _WIN32
        sprintf_s(cTemp, 13, ", Device%d = ", dev);
#else
        sprintf(cTemp, ", Device%d = ", dev);
#endif
        cudaDeviceProp deviceProp;
        cudaGetDeviceProperties(&deviceProp, dev);
        sProfileString += cTemp;
        sProfileString += deviceProp.name;
    }
    sProfileString += "\n";

	sout << sProfileString;
}

void set_kernel_launch_param(uint32_t n_data, uint16_t n_tpb, dim3& grid, dim3& block)
{
	uint32_t n_thread = min((uint32_t)n_tpb, n_data);
	uint32_t n_block = (n_data + n_thread - 1)/n_thread;

	grid.x	= n_block;
	block.x = n_thread;
}

void device_query(ostream& sout, int id_dev, bool print_to_screen)
{
	device_query(sout, id_dev);
	if (print_to_screen)
	{
		device_query(cout, id_dev);
	}
}

void allocate_host_vector(void **ptr, size_t size, const char *file, int line)
{
	*ptr = (void *)malloc(size);
	if (NULL == ptr)
	{
		throw string("Out of memory.");
	}

	// Clear memory 
	memset(*ptr, 0, size);
}

void allocate_device_vector(void **ptr, size_t size, const char *file, int line)
{
	// Allocate memory
	CUDA_SAFE_CALL(cudaMalloc(ptr, size));
	// Clear memory 
	CUDA_SAFE_CALL(cudaMemset(*ptr, 0, size));
}

void allocate_vector(void **ptr, size_t size, bool cpu, const char *file, int line)
{
	if (cpu)
	{
		allocate_host_vector(ptr, size, file, line);
	}
	else
	{
		allocate_device_vector(ptr, size, file, line);
	}
}

void free_host_vector(void **ptr, const char *file, int line)
{
	if (NULL != *ptr)
	{
		delete[] *ptr;
		*ptr = NULL;
	}
}

void free_device_vector(void **ptr, const char *file, int line)
{
	if (NULL != *ptr)
	{
		CUDA_SAFE_CALL(cudaFree(*ptr));
		*ptr = NULL;
	}
}

//void free_vector(void **ptr, bool cpu, const char *file, int line)
//{
//	if (cpu)
//	{
//		free_host_vector(ptr, file, line);
//	}
//	else
//	{
//		free_device_vector(ptr, file, line);
//	}
//}

void allocate_host_storage(pp_disk_t::sim_data_t *sd, int n)
{
	sd->h_y.resize(2);
	sd->h_yout.resize(2);

	for (int i = 0; i < 2; i++)
	{
		ALLOCATE_HOST_VECTOR((void **)&(sd->h_y[i]),    n*sizeof(var4_t));
		ALLOCATE_HOST_VECTOR((void **)&(sd->h_yout[i]), n*sizeof(var4_t));
	}
	ALLOCATE_HOST_VECTOR((void **)&(sd->h_p),           n*sizeof(pp_disk_t::param_t));
	ALLOCATE_HOST_VECTOR((void **)&(sd->h_body_md),     n*sizeof(pp_disk_t::body_metadata_t));
	ALLOCATE_HOST_VECTOR((void **)&(sd->h_epoch),       n*sizeof(var_t));

	ALLOCATE_HOST_VECTOR((void **)&(sd->h_oe),          n*sizeof(orbelem_t));
}

void allocate_device_storage(pp_disk_t::sim_data_t *sd, int n)
{
	sd->d_y.resize(2);
	sd->d_yout.resize(2);

	for (int i = 0; i < 2; i++)
	{
		ALLOCATE_DEVICE_VECTOR((void **)&(sd->d_y[i]),	  n*sizeof(var4_t));
		ALLOCATE_DEVICE_VECTOR((void **)&(sd->d_yout[i]), n*sizeof(var4_t));
	}
	ALLOCATE_DEVICE_VECTOR((void **)&(sd->d_p),			  n*sizeof(pp_disk_t::param_t));
	ALLOCATE_DEVICE_VECTOR((void **)&(sd->d_body_md),	  n*sizeof(pp_disk_t::body_metadata_t));
	ALLOCATE_DEVICE_VECTOR((void **)&(sd->d_epoch),		  n*sizeof(var_t));

    ALLOCATE_DEVICE_VECTOR((void **)&(sd->d_oe),          n*sizeof(orbelem_t));
}

void deallocate_host_storage(pp_disk_t::sim_data_t *sd)
{
	for (int i = 0; i < 2; i++)
	{
		FREE_HOST_VECTOR((void **)&(sd->h_y[i]));
		FREE_HOST_VECTOR((void **)&(sd->h_yout[i]));
	}
	FREE_HOST_VECTOR((void **)&(sd->h_p));
	FREE_HOST_VECTOR((void **)&(sd->h_body_md));
	FREE_HOST_VECTOR((void **)&(sd->h_epoch));

	FREE_HOST_VECTOR((void **)&(sd->h_oe));
}

void deallocate_device_storage(pp_disk_t::sim_data_t *sd)
{
	for (int i = 0; i < 2; i++)
	{
		FREE_DEVICE_VECTOR((void **)&(sd->d_y[i]));
		FREE_DEVICE_VECTOR((void **)&(sd->d_yout[i]));
	}
	FREE_DEVICE_VECTOR((void **)&(sd->d_p));
	FREE_DEVICE_VECTOR((void **)&(sd->d_body_md));
	FREE_DEVICE_VECTOR((void **)&(sd->d_epoch));

    FREE_DEVICE_VECTOR((void **)&(sd->d_oe));
}


void copy_vector_to_device(void* dst, const void *src, size_t count)
{
	CUDA_SAFE_CALL(cudaMemcpy(dst, src, count, cudaMemcpyHostToDevice));
}

void copy_vector_to_host(void* dst, const void *src, size_t count)
{
	CUDA_SAFE_CALL(cudaMemcpy(dst, src, count, cudaMemcpyDeviceToHost));
}

void copy_vector_d2d(void* dst, const void *src, size_t count)
{
	CUDA_SAFE_CALL(cudaMemcpy(dst, src, count, cudaMemcpyDeviceToDevice));
}

void copy_constant_to_device(const void* dst, const void *src, size_t count)
{
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(dst, src, count, 0, cudaMemcpyHostToDevice));
}


void set_device(int id_of_target_dev, ostream& sout)
{
	int n_device = get_n_cuda_device();
	if (0 == n_device)
	{
		throw string("There are no available device(s) that support CUDA.");
	}

	if (n_device > id_of_target_dev && 0 <= id_of_target_dev)
	{
		// Set the desired id of the device
		CUDA_SAFE_CALL(cudaSetDevice(id_of_target_dev));
	}
	else
	{
		throw string("The device with the requested id does not exist.");
	}
}

void print_array(string path, string comment, uint32_t n, var_t *data, mem_loc_t mem_loc)
{
	var_t* h_data = NULL;

	ostream *out = NULL;
	if (0 < path.length())
	{
		out = new ofstream(path.c_str(), ios::app);
	}
	else
	{
		out = &cout;
	}

	out->setf(ios::right);
	out->setf(ios::scientific);

	if (MEM_LOC_DEVICE == mem_loc)
	{
		h_data = new var_t[n];
		copy_vector_to_host(h_data, data, n * sizeof(var_t));
	}
	else
	{
		h_data = data;
	}
    if (!comment.empty()) *out << comment << endl;
    for (uint32_t i = 0; i < n; i++)
	{
		*out << setw(5) << i << setprecision(16) << setw(25) << h_data[i] << endl;
	}

	if (MEM_LOC_DEVICE == mem_loc)
	{
		delete[] h_data;
	}
	if (0 < path.length())
	{
		out->flush();
		delete out;
	}
}

void create_aliases(comp_dev_t comp_dev, pp_disk_t::sim_data_t *sd)
{
	switch (comp_dev.proc_unit)
	{
	case PROC_UNIT_CPU:
		for (int i = 0; i < 2; i++)
		{
			sd->y[i]    = sd->h_y[i];
			sd->yout[i] = sd->h_yout[i];
		}
		sd->p       = sd->h_p;
		sd->body_md = sd->h_body_md;
		sd->epoch   = sd->h_epoch;
        sd->oe      = sd->h_oe;
		break;
	case PROC_UNIT_GPU:
		for (int i = 0; i < 2; i++)
		{
			sd->y[i]    = sd->d_y[i];
			sd->yout[i] = sd->d_yout[i];
		}
		sd->p       = sd->d_p;
		sd->body_md = sd->d_body_md;
		sd->epoch   = sd->d_epoch;
        sd->oe      = sd->d_oe;
		break;
	default:
		throw string("Parameter 'proc_unit' is out of range.");
	}
}

// Date of creation: 2016.11.22.
// Last edited: 
// Status:
void gpu_calc_lin_comb_s(var_t* a, const var_t* b, const var_t* c, var_t f, uint32_t n_var, int id_dev, bool optimize)
{
	static uint16_t n_tpb = 256;
	static bool first_call = true;

	dim3 grid;
	dim3 block;

	if (optimize || first_call)
	{
		first_call = false;

		cudaDeviceProp prop;
		cudaEvent_t start, stop;

        CUDA_SAFE_CALL(cudaGetDeviceProperties(&prop, id_dev));
		CUDA_SAFE_CALL(cudaEventCreate(&start));
		CUDA_SAFE_CALL(cudaEventCreate(&stop));

		float min_GPU_DT = 1.0e8;
		uint16_t d_nt = prop.warpSize / 2;
		for (uint16_t nt = d_nt; nt <= prop.maxThreadsPerBlock / 2; nt += d_nt)
		{
			set_kernel_launch_param(n_var, nt, grid, block);

			CUDA_SAFE_CALL(cudaEventRecord(start));
			red_kernel::calc_lin_comb_s<<<grid, block>>>(a, b, f, c, n_var);
       		CUDA_CHECK_ERROR();
			CUDA_SAFE_CALL(cudaEventRecord(stop));
			CUDA_SAFE_CALL(cudaEventSynchronize(stop));

			float GPU_DT = 0.0f;
			CUDA_SAFE_CALL(cudaEventElapsedTime(&GPU_DT, start, stop));
			if (GPU_DT < min_GPU_DT)
			{
				min_GPU_DT = GPU_DT;
				n_tpb = nt;
			}
			//printf("%4u %10.6f [ms]\n", nt, GPU_DT);
		}
		//printf("\n%4u %10.6f [ms]\n", n_tpb, min_GPU_DT);
	}
	else
	{
		set_kernel_launch_param(n_var, n_tpb, grid, block);
		red_kernel::calc_lin_comb_s<<<grid, block>>>(a, b, f, c, n_var);
   		CUDA_CHECK_ERROR();
	}
}

void gpu_calc_lin_comb_s(var_t* a, const var_t* b, const var_t* const *c, const var_t* coeff, uint16_t n_vct, uint32_t n_var, int id_dev, bool optimize)
{
	static uint16_t n_tpb = 256;
	static bool first_call = true;

	dim3 grid;
	dim3 block;

	if (optimize || first_call)
	{
		first_call = false;

		cudaDeviceProp prop;
		cudaEvent_t start, stop;

        CUDA_SAFE_CALL(cudaGetDeviceProperties(&prop, id_dev));
		CUDA_SAFE_CALL(cudaEventCreate(&start));
		CUDA_SAFE_CALL(cudaEventCreate(&stop));

		float min_GPU_DT = 1.0e8;
		uint16_t d_nt = prop.warpSize / 2;
		for (uint16_t nt = d_nt; nt <= prop.maxThreadsPerBlock / 2; nt += d_nt)
		{
			set_kernel_launch_param(n_var, nt, grid, block);

            CUDA_SAFE_CALL(cudaEventRecord(start));
			red_kernel::calc_lin_comb_s<<<grid, block>>>(a, b, c, coeff, n_vct, n_var);
    		CUDA_CHECK_ERROR();
			CUDA_SAFE_CALL(cudaEventRecord(stop));
			CUDA_SAFE_CALL(cudaEventSynchronize(stop));

			float GPU_DT = 0.0f;
			CUDA_SAFE_CALL(cudaEventElapsedTime(&GPU_DT, start, stop));
			if (GPU_DT < min_GPU_DT)
			{
				min_GPU_DT = GPU_DT;
				n_tpb = nt;
			}
			//printf("%4u %10.6f [ms]\n", nt, GPU_DT);
		}
		//printf("\n%4u %10.6f [ms]\n", n_tpb, min_GPU_DT);
	}
	else
	{
		set_kernel_launch_param(n_var, n_tpb, grid, block);
		red_kernel::calc_lin_comb_s<<<grid, block>>>(a, b, c, coeff, n_vct, n_var);
   		CUDA_CHECK_ERROR();
	}
}

// Calculate the error for the Runge-Kutta 4 method: a = |k4 - k5|
void gpu_calc_rk4_error(var_t* a, const var_t* k4, const var_t* k5, uint32_t n_var, int id_dev, bool optimize)
{
	static uint16_t n_tpb = 256;
	static bool first_call = true;

	dim3 grid;
	dim3 block;

	if (optimize || first_call)
	{
		first_call = false;

		cudaDeviceProp prop;
		cudaEvent_t start, stop;

        CUDA_SAFE_CALL(cudaGetDeviceProperties(&prop, id_dev));
		CUDA_SAFE_CALL(cudaEventCreate(&start));
		CUDA_SAFE_CALL(cudaEventCreate(&stop));

		float min_GPU_DT = 1.0e8;
		uint16_t d_nt = prop.warpSize / 2;
		for (uint16_t nt = d_nt; nt <= prop.maxThreadsPerBlock / 2; nt += d_nt)
		{
			set_kernel_launch_param(n_var, nt, grid, block);

			CUDA_SAFE_CALL(cudaEventRecord(start));
            red_kernel::calc_rk4_error<<<grid, block>>>(a, k4, k5, n_var);
       		CUDA_CHECK_ERROR();
			CUDA_SAFE_CALL(cudaEventRecord(stop));
			CUDA_SAFE_CALL(cudaEventSynchronize(stop));

			float GPU_DT = 0.0f;
			CUDA_SAFE_CALL(cudaEventElapsedTime(&GPU_DT, start, stop));
			if (GPU_DT < min_GPU_DT)
			{
				min_GPU_DT = GPU_DT;
				n_tpb = nt;
			}
			//printf("%4u %10.6f [ms]\n", nt, GPU_DT);
		}
		//printf("\n%4u %10.6f [ms]\n", n_tpb, min_GPU_DT);
	}
	else
	{
		set_kernel_launch_param(n_var, n_tpb, grid, block);
        red_kernel::calc_rk4_error<<<grid, block>>>(a, k4, k5, n_var);
   		CUDA_CHECK_ERROR();
	}
}

// Calculate the error for the Runge-Kutta 5 method: a = |k5 - k6|
void gpu_calc_rk5_error(var_t* a, const var_t* k5, const var_t* k6, uint32_t n_var, int id_dev, bool optimize)
{
	static uint16_t n_tpb = 256;
	static bool first_call = true;

	dim3 grid;
	dim3 block;

	if (optimize || first_call)
	{
		first_call = false;

		cudaDeviceProp prop;
		cudaEvent_t start, stop;

        CUDA_SAFE_CALL(cudaGetDeviceProperties(&prop, id_dev));
		CUDA_SAFE_CALL(cudaEventCreate(&start));
		CUDA_SAFE_CALL(cudaEventCreate(&stop));

		float min_GPU_DT = 1.0e8;
		uint16_t d_nt = prop.warpSize / 2;
		for (uint16_t nt = d_nt; nt <= prop.maxThreadsPerBlock / 2; nt += d_nt)
		{
			set_kernel_launch_param(n_var, nt, grid, block);

			CUDA_SAFE_CALL(cudaEventRecord(start));
            red_kernel::calc_rk5_error<<<grid, block>>>(a, k5, k6, n_var);
       		CUDA_CHECK_ERROR();
			CUDA_SAFE_CALL(cudaEventRecord(stop));
			CUDA_SAFE_CALL(cudaEventSynchronize(stop));

			float GPU_DT = 0.0f;
			CUDA_SAFE_CALL(cudaEventElapsedTime(&GPU_DT, start, stop));
			if (GPU_DT < min_GPU_DT)
			{
				min_GPU_DT = GPU_DT;
				n_tpb = nt;
			}
			//printf("%4u %10.6f [ms]\n", nt, GPU_DT);
		}
		//printf("\n%4u %10.6f [ms]\n", n_tpb, min_GPU_DT);
	}
	else
	{
		set_kernel_launch_param(n_var, n_tpb, grid, block);
        red_kernel::calc_rk5_error<<<grid, block>>>(a, k5, k6, n_var);
   		CUDA_CHECK_ERROR();
	}
}

// Calculate the error for the Runge-Kutta 7 method: a = |k1 + k11 - k12 - k13|
void gpu_calc_rk7_error(var_t* a, const var_t* k1, const var_t* k11, const var_t* k12, const var_t* k13, uint32_t n_var, int id_dev, bool optimize)
{
	static uint16_t n_tpb = 256;
	static bool first_call = true;

	dim3 grid;
	dim3 block;

	if (optimize || first_call)
	{
		first_call = false;

		cudaDeviceProp prop;
		cudaEvent_t start, stop;

        CUDA_SAFE_CALL(cudaGetDeviceProperties(&prop, id_dev));
		CUDA_SAFE_CALL(cudaEventCreate(&start));
		CUDA_SAFE_CALL(cudaEventCreate(&stop));

		float min_GPU_DT = 1.0e8;
		uint16_t d_nt = prop.warpSize / 2;
		for (uint16_t nt = d_nt; nt <= prop.maxThreadsPerBlock / 2; nt += d_nt)
		{
			set_kernel_launch_param(n_var, nt, grid, block);

			CUDA_SAFE_CALL(cudaEventRecord(start));
            red_kernel::calc_rk7_error<<<grid, block>>>(a, k1, k11, k12, k13, n_var);
       		CUDA_CHECK_ERROR();
			CUDA_SAFE_CALL(cudaEventRecord(stop));
			CUDA_SAFE_CALL(cudaEventSynchronize(stop));

			float GPU_DT = 0.0f;
			CUDA_SAFE_CALL(cudaEventElapsedTime(&GPU_DT, start, stop));
			if (GPU_DT < min_GPU_DT)
			{
				min_GPU_DT = GPU_DT;
				n_tpb = nt;
			}
			//printf("%4u %10.6f [ms]\n", nt, GPU_DT);
		}
		//printf("\n%4u %10.6f [ms]\n", n_tpb, min_GPU_DT);
	}
	else
	{
		set_kernel_launch_param(n_var, n_tpb, grid, block);
        red_kernel::calc_rk7_error<<<grid, block>>>(a, k1, k11, k12, k13, n_var);
   		CUDA_CHECK_ERROR();
	}
}
} /* redutil2 */
