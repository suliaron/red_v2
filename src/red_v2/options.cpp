#include <iostream>
#include <fstream>

#include "options.h"
#include "integrator.h"
#include "int_euler.h"
#include "int_rungekutta2.h"
#include "int_rungekutta4.h"
#include "int_rungekutta5.h"
#include "int_rungekutta7.h"
#include "int_hermite4.h"
#include "parameter.h"
#include "tbp1D.h"
#include "tbp2D.h"
#include "tbp3D.h"
#include "rtbp1D.h"
#include "rtbp2D.h"
#include "rtbp3D.h"
#include "threebody.h"
#include "nbody.h"

#include "redutil2.h"
#include "constants.h"
#include "macro.h"

using namespace std;
using namespace redutil2;

options::options(int argc, const char** argv)
{
	static const char* integrator_type_short_name[] = {"E",	"RK2", "RK4", "RK5", "RK7", "H4", "H4b"};

	create_default();
	parse(argc, argv);
	param = new parameter(dir[DIRECTORY_NAME_IN], in_fn[INPUT_NAME_PARAMETER], print_to_screen);

	if (INTEGRATOR_EULER == param->int_type || INTEGRATOR_RUNGEKUTTA2 == param->int_type)
	{
		param->adaptive = false;
	}

	if (PROC_UNIT_GPU == comp_dev.proc_unit)
	{
		int n_device = 0;
		CUDA_SAFE_CALL(cudaGetDeviceCount(&n_device));
		if (0 == n_device)
		{
			throw string("There are no cuda capable device in the machine.");
		}
		if (0 > comp_dev.id_dev || comp_dev.id_dev > n_device)
		{
			throw string("The selected device does not exist.");
		}
		CUDA_SAFE_CALL(cudaSetDevice(comp_dev.id_dev));
	}

	if (ef)
	{
		char sep = '_';
		string config;
#ifdef _DEBUG
		config = "D";
#else
		config = "R";
#endif
		string dev = (PROC_UNIT_CPU == comp_dev.proc_unit ? "cpu" : "gpu");
		// as: adaptive step-size, fs: fix step-size
		string adapt = (param->adaptive == true ? "as" : "fs");
		
		string int_name(integrator_type_short_name[param->int_type]);
		fn_prefix += config + sep + dev + sep + adapt + sep + int_name + sep;
	}
}

options::~options() 
{ }

void options::create_default()
{
	param                              = NULL;
	dyn_model                          = DYN_MODEL_N;
	verbose                            = false;
	print_to_screen                    = false;
	append                             = false;
	ef                                 = false;
	benchmark                          = false;
	n_change_to_cpu                    = 100;
	comp_dev.proc_unit                 = PROC_UNIT_CPU;
	comp_dev.id_dev                    = -1;
	g_disk_model                       = GAS_DISK_MODEL_NONE;
	out_fn[OUTPUT_NAME_LOG]            = "log";
	out_fn[OUTPUT_NAME_INFO]           = "info";
	out_fn[OUTPUT_NAME_EVENT]          = "event";
	out_fn[OUTPUT_NAME_SOLUTION_INFO]  = "solution.info";
	out_fn[OUTPUT_NAME_SOLUTION_DATA]  = "solution.data";
	out_fn[OUTPUT_NAME_INTEGRAL]       = "integral";
	out_fn[OUTPUT_NAME_INTEGRAL_EVENT] = "integral.event";
}

void options::parse(int argc, const char** argv)
{
	int i = 1;

	if (1 >= argc)
	{
		throw string("Missing command line arguments. For help use -h.");
	}

	while (i < argc)
	{
		string p = argv[i];
		if (     p == "-m")
		{
			i++;
			string value = argv[i];
			if (     value == "tbp1D")
			{
				dyn_model = DYN_MODEL_TBP1D;
			}
			else if (value == "tbp2D")
			{
				dyn_model = DYN_MODEL_TBP2D;
			}
			else if (value == "tbp3D")
			{
				dyn_model = DYN_MODEL_TBP3D;
			}
			else if (value == "rtbp1D")
			{
				dyn_model = DYN_MODEL_RTBP1D;
			}
			else if (value == "rtbp2D")
			{
				dyn_model = DYN_MODEL_RTBP2D;
			}
			else if (value == "rtbp3D")
			{
				dyn_model = DYN_MODEL_RTBP3D;
			}
			else if (value == "threebody")
			{
				dyn_model = DYN_MODEL_THREEBODY;
			}
			else if (value == "nbody")
			{
				dyn_model = DYN_MODEL_NBODY;
			}
			else
			{
				throw string("Invalid dynamical model: " + value + ".");
			}
		}
		else if (p == "-v")
		{
			verbose = true;
		}
		else if (p == "-pts")
		{
			verbose = true;
			print_to_screen = true;
		}
		else if (p == "-a")
		{
			append = true;
		}
		else if (p == "-ef")
		{
			ef = true;
		}
		else if (p == "-b")
		{
			benchmark = true;
		}

		else if (p == "-id_dev")
		{
			i++;
			if (!tools::is_number(argv[i])) 
			{
				throw string("Invalid number at: " + p + ".");
			}
			comp_dev.id_dev = atoi(argv[i]);
		}
		else if (p == "-n_chg")
		{
			i++;
			if (!tools::is_number(argv[i])) 
			{
                throw string("Invalid number at: " + p + ".");
			}
			n_change_to_cpu = atoi(argv[i]);
		}

		else if (p == "-cpu")
		{
			comp_dev.proc_unit = PROC_UNIT_CPU;
		}
		else if (p == "-gpu")
		{
			comp_dev.proc_unit = PROC_UNIT_GPU;
		}
		else if (p == "-ga")
		{
			i++;
			in_fn[INPUT_NAME_GAS_DISK_MODEL] = argv[i];
			g_disk_model = GAS_DISK_MODEL_ANALYTIC;
		}
		else if (p == "-gf")
		{
			i++;
			in_fn[INPUT_NAME_GAS_DISK_MODEL] = argv[i];
			g_disk_model = GAS_DISK_MODEL_FARGO;
		}

		else if (p == "-i")
		{
			i++;
			in_fn[INPUT_NAME_START_FILES] = argv[i];
		}
		else if (p == "-p")
		{
			i++;
			in_fn[INPUT_NAME_PARAMETER] = argv[i];
		}

		else if (p == "-idir")
		{
			i++;
			dir[DIRECTORY_NAME_IN] = argv[i];
		}
		else if (p == "-odir")
		{
			i++;
			dir[DIRECTORY_NAME_OUT] = argv[i];
		}

		else if (p == "-h")
		{
			print_usage();
			exit(EXIT_SUCCESS);
		}
		else
		{
			throw string("Invalid switch on command line: " + p + ".");
		}
		i++;
	}
	// Check for obligatory arguments
	if (PROC_UNIT_GPU == comp_dev.proc_unit && -1 == comp_dev.id_dev)
	{
		throw string("Missing the id of the GPU device (see -id_dev).");
	}
	if (0 == in_fn[INPUT_NAME_START_FILES].length())
	{
		throw string("Missing value for -i.");
	}
	if (0 == in_fn[INPUT_NAME_PARAMETER].length())
	{
		throw string("Missing value for -p.");
	}
	if (DYN_MODEL_N == dyn_model)
	{
		throw string("Missing value for -m.");
	}
	if (0 == dir[DIRECTORY_NAME_OUT].length())
	{
		dir[DIRECTORY_NAME_OUT] = dir[DIRECTORY_NAME_IN];
	}
}

void options::get_solution_path(string& path_si, string &path_sd)
{
	string path = file::combine_path(dir[DIRECTORY_NAME_IN], in_fn[INPUT_NAME_START_FILES]);
	ifstream file(path.c_str(), ifstream::in);
	if (file)
	{
		uint32_t n = 0;
		string str;
		while (getline(file, str))
		{
			if (0 == n)
			{
                in_fn[INPUT_NAME_IC_INFO] = str;
			}
			if (1 == n)
			{
                in_fn[INPUT_NAME_IC_DATA] = str;
			}
			n++;
		} 	
		file.close();
	}
	else
	{
		throw string("Cannot open " + path + ".");
	}
    path_si = file::combine_path(dir[DIRECTORY_NAME_IN], in_fn[INPUT_NAME_IC_INFO]);
    path_sd = file::combine_path(dir[DIRECTORY_NAME_IN], in_fn[INPUT_NAME_IC_DATA]);
}

ode* options::create_model()
{
	ode *model;
    string path_si;
    string path_sd;

    get_solution_path(path_si, path_sd);

	switch (dyn_model)
	{
	case DYN_MODEL_TBP1D:
		{
		    model = new tbp1D(path_si, path_sd, 1, comp_dev);
			break;
		}
	case DYN_MODEL_TBP2D:
		{
		    model = new tbp2D(path_si, path_sd, 1, comp_dev);
			break;
		}
	case DYN_MODEL_TBP3D:
		{
		    model = new tbp3D(1, comp_dev);
			break;
		}
	case DYN_MODEL_RTBP1D:
		{
		    model = new rtbp1D(path_si, path_sd, 1, comp_dev);
			break;
		}
	case DYN_MODEL_RTBP2D:
		{
		    model = new rtbp2D(path_si, path_sd, 1, comp_dev);
			break;
		}
	case DYN_MODEL_RTBP3D:
		{
		    model = new rtbp3D(1, comp_dev);
			break;
		}
	case DYN_MODEL_THREEBODY:
		{
	    	model = new threebody(1, comp_dev);
			break;
		}
	case DYN_MODEL_NBODY:
		{
			uint32_t n_obj = 0;
            ifstream input;

            data_rep_t repres = file::get_data_repres(path_si);
            switch (repres)
            {
            case DATA_REPRESENTATION_ASCII:
                input.open(path_si.c_str(), ios::in);
                if (input)
                {
                    var_t t, dt, t_wc;

                    input >> t >> dt >> t_wc >> n_obj;
                }
                else
                {
                    throw string("Cannot open " + path_si + ".");
                }
                break;
            case DATA_REPRESENTATION_BINARY:
                input.open(path_si.c_str(), ios::in | ios::binary);
                if (input)
                {
                    var_t t, dt, t_wc;

                    input.read((char*)&t, sizeof(var_t));
                    input.read((char*)&dt, sizeof(var_t));
                    input.read((char*)&t_wc, sizeof(var_t));
                    input.read((char*)&n_obj, sizeof(uint32_t));
                }
                else
                {
                    throw string("Cannot open " + path_si + ".");
                }
                break;
            default:
                throw string("Parameter 'repres' is out of range.");
            }
            input.close();

            uint16_t n_ppo = sizeof(nbp_t::param_t) / sizeof(var_t);
			model = new nbody(path_si, path_sd, n_obj, n_ppo, comp_dev);
			break;
		}
	default:
		throw string("Invalid dynamical model.");
	}

	return model;
}

integrator* options::create_integrator(ode& f)
{
	integrator* intgr = NULL;

	switch (param->int_type)
	{
	case INTEGRATOR_EULER:
		intgr = new euler(f, comp_dev);
		break;
	case INTEGRATOR_RUNGEKUTTA2:
		intgr = new int_rungekutta2(f, comp_dev);
		break;
	case INTEGRATOR_RUNGEKUTTA4:
		intgr = new int_rungekutta4(f, param->adaptive, param->tolerance, comp_dev);
		break;
	case INTEGRATOR_RUNGEKUTTA5:
		intgr = new int_rungekutta5(f, param->adaptive, param->tolerance, comp_dev);
		break;
	case INTEGRATOR_RUNGEKUTTA7:
		intgr = new int_rungekutta7(f, param->adaptive, param->tolerance, comp_dev);
		break;
	case INTEGRATOR_HERMITE4:
        intgr = new int_hermite4(f, param->adaptive, param->eta, comp_dev);
		break;
	case INTEGRATOR_HERMITE4B:
        intgr = new int_hermite4b(f, param->adaptive, param->eta, comp_dev);
		break;
	default:
		throw string("Requested integrator is not implemented.");
	}

	// TODO: For every model it should be provieded a method to determine the minimum stepsize
	// OR use the solution provided by the Numerical Recepies
	intgr->set_dt_min(1.0e-20); // day
	intgr->set_max_iter(100);

    if (param->error_check_for_tp)
	{
		intgr->error_check_for_tp = true;
	}

	return intgr;
}

void options::print_usage()
{
	cout << "Usage: red.cuda <parameterlist>" << endl;
	cout << "Parameters:" << endl;
	cout << "     -t      | --test                         : run tests" << endl;
	cout << "     -v      | --verbose                      : verbose mode (log all event during the execution fo the code to the log file)" << endl;
	cout << "     -pts    | --print_to_screen              : verbose mode and print everything to the standard output stream too" << endl;
	cout << "     -ef     |                                : use extended file names (use only for debuging purposes)" << endl;

	cout << "     -i_dt   | --info-dt <number>             : the time interval in seconds between two subsequent information print to the screen (default value is 5 sec)" << endl;
	cout << "     -d_dt   | --dump-dt <number>             : the time interval in seconds between two subsequent data dump to the hard disk (default value is 3600 sec)" << endl;

	cout << "     -id_dev | --id_active_device <number>    : the id of the device which will execute the code (default value is 0)" << endl;
	cout << "     -n_chg  | --n_change_to_cpu <number>     : the threshold value for the total number of SI bodies to change to the CPU (default value is 100)" << endl;

	cout << "     -gpu    | --gpu                          : execute the code on the graphics processing unit (GPU) (default value is true)" << endl;
	cout << "     -cpu    | --cpu                          : execute the code on the cpu if required by the user or if no GPU is installed (default value is false)" << endl;

	cout << "     -log    | --log-filename <filename>      : the name of the file where the details of the execution of the code will be stored (default value is log.txt)" << endl;
	cout << "     -info   | --info-filename <filename>     : the name of the file where the runtime output of the code will be stored (default value is info.txt)" << endl;
	cout << "     -event  | --event-filename <filename>    : the name of the file where the details of each event will be stored (default value is event.txt)" << endl;
	cout << "     -result | --result-filename <filename>   : the name of the file where the simlation data for a time instance will be stored (default value is solution.txt)" << endl;

	cout << "     -iDir   | --inputDir <directory>         : the directory containing the input files"  << endl;
	cout << "     -oDir   | --outputDir <directory>        : the directory where the output files will be stored (if omitted the input directory will be used)" << endl;
	cout << "     -ic     | --initial_condition <filename> : the file containing the initial conditions"  << endl;
	cout << "     -p      | --parameter <filename>         : the file containing the parameters of the simulation"  << endl;

	cout << "     -ga     | --analytic_gas_disk <filename> : the file containing the parameters of an analyticaly prescribed gas disk"  << endl;
	cout << "     -gf     | --fargo_gas_disk <filename>    : the file containing the details of the gas disk resulted from FARGO simulations"  << endl;

	cout << "     -h      | --help                         : print this help" << endl;
}
