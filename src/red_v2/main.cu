#include <ctime>
#include <chrono>
#include <iostream>
#include <iomanip>      // std::setw
#include <fstream>

#include "integrator.h"
#include "ode.h"
#include "options.h"
#include "rtbp1D.h"
#include "tbp1D.h"
#include "rtbp3D.h"
#include "tbp3D.h"
#include "threebody.h"
#include "nbody.h"

#include "tools.h"
#include "type.h"
#include "constants.h"
#include "redutil2.h"

using namespace std;
using namespace redutil2;


void print_dump(options* opt, ode* f, integrator* intgr)
{
    // Download data from DEVICE ?
    if (PROC_UNIT_GPU == opt->comp_dev.proc_unit)
    {
        f->copy_vars(COPY_DIRECTION_TO_HOST);
        f->copy_params(COPY_DIRECTION_TO_HOST);
        f->copy_metadata(COPY_DIRECTION_TO_HOST);
    }

    string fn_info = opt->fn_prefix + "dump.info.dat";
    string path_si = file::combine_path(opt->dir[DIRECTORY_NAME_OUT], fn_info);

    string fn_data = opt->fn_prefix + "dump.data.dat";
    string path_sd = file::combine_path(opt->dir[DIRECTORY_NAME_OUT], fn_data);

    f->print_dump(path_si, path_sd);

    string path = file::combine_path(opt->dir[DIRECTORY_NAME_OUT], "start_files.txt");
    ofstream sout(path.c_str(), ios_base::out);
    if (sout)
    {
        sout << fn_info << endl;
        sout << fn_data << endl;
    }
    else
    {
        throw string("Cannot open " + path + ".");
    }
    sout.close();
}

void print_solution(uint32_t& n_print, options* opt, ode* f, integrator* intgr, ofstream& slog)
{
	static string ext = (DATA_REPRESENTATION_ASCII == opt->param->output_data_rep ? "txt" : "dat");
	static var_t last_download = 0.0;

	string fn_info;
	string path_si;
	string fn_data;
	string path_sd;

	// Download data from DEVICE ?
	if (PROC_UNIT_GPU == opt->comp_dev.proc_unit && 0.0 < fabs(last_download - f->t))
	{
		last_download = f->t;
		f->copy_vars(COPY_DIRECTION_TO_HOST);
		f->copy_params(COPY_DIRECTION_TO_HOST);
		f->copy_metadata(COPY_DIRECTION_TO_HOST);
	}

	if (opt->append)
	{
		fn_info = opt->fn_prefix + opt->out_fn[OUTPUT_NAME_SOLUTION_INFO] + "." + ext;
		path_si = file::combine_path(opt->dir[DIRECTORY_NAME_OUT], fn_info);

		fn_data = opt->fn_prefix + opt->out_fn[OUTPUT_NAME_SOLUTION_DATA] + "." + ext;
		path_sd = file::combine_path(opt->dir[DIRECTORY_NAME_OUT], fn_data);

		f->print_solution(path_si, path_sd, opt->param->output_data_rep);
	}
	else
	{
		string n_print_str = redutil2::number_to_string(n_print, OUTPUT_ORDINAL_NUMBER_WIDTH, true);

		fn_info = opt->fn_prefix + opt->out_fn[OUTPUT_NAME_SOLUTION_INFO] + "_" + n_print_str + "." + ext;
		path_si = file::combine_path(opt->dir[DIRECTORY_NAME_OUT], fn_info);

		fn_data = opt->fn_prefix + opt->out_fn[OUTPUT_NAME_SOLUTION_DATA] + "_" + n_print_str + "." + ext;
		path_sd = file::combine_path(opt->dir[DIRECTORY_NAME_OUT], fn_data);

		f->print_solution(path_si, path_sd, opt->param->output_data_rep);
		n_print++;

		string path = file::combine_path(opt->dir[DIRECTORY_NAME_OUT], "start_files.txt");
		ofstream sout(path.c_str(), ios_base::out);
		if (sout)
		{
			sout << fn_info << endl;
			sout << fn_data << endl;
		}
		else
		{
			throw string("Cannot open " + path + ".");
		}
		sout.close();
	}
}

void print_info(options* opt, ode* f, integrator* intgr, var_t dt, var_t total_t)
{
    static char FMT_STR[] = "%s %s %13.7e %10.4e  %12.6e %9.3e  %12llu %12llu %12llu\n";
    string fn_info;
    string path;

    fn_info = opt->fn_prefix + opt->out_fn[OUTPUT_NAME_INFO] + ".txt";
    path = file::combine_path(opt->dir[DIRECTORY_NAME_OUT], fn_info);
    FILE* fout = fopen(path.c_str(), "a");
    if (NULL == fout)
    {
        throw string("Cannot open " + path + ".");
    }

    string time_stamp = tools::get_time_stamp();
    string dev = (opt->comp_dev.proc_unit == PROC_UNIT_CPU ? "CPU" : "GPU");
    var_t dt_did = intgr->get_dt_did();
    uint64_t np = intgr->get_n_passed_step();
    uint64_t nf = intgr->get_n_failed_step();
    uint64_t nt = intgr->get_n_tried_step();

    fprintf(fout, FMT_STR, time_stamp.c_str(), dev.c_str(), f->t, dt_did, total_t, dt, np, nf, nt);
    if (opt->print_to_screen)
    {
        printf(FMT_STR, time_stamp.c_str(), dev.c_str(), f->t, dt_did, total_t, dt, np, nf, nt);
    }
    fclose(fout);
}

void run_benchmark(options* opt, ode* f, integrator* intgr, ofstream& slog)
{ }

void run_simulation(options* opt, ode* f, integrator* intgr, ofstream& slog)
{
	static string ext = (DATA_REPRESENTATION_ASCII == opt->param->output_data_rep ? "txt" : "dat");

	static string outdir = opt->dir[DIRECTORY_NAME_OUT];
	static string path_info           = file::combine_path(outdir, opt->fn_prefix + opt->out_fn[OUTPUT_NAME_INFO] + ".txt");
	static string path_event          = file::combine_path(outdir, opt->fn_prefix + opt->out_fn[OUTPUT_NAME_EVENT] + ".txt");
	static string path_integral       = file::combine_path(outdir, opt->fn_prefix + opt->out_fn[OUTPUT_NAME_INTEGRAL] + ".txt");
	static string path_integral_event = file::combine_path(outdir, opt->fn_prefix + opt->out_fn[OUTPUT_NAME_INTEGRAL_EVENT] + ".txt");

	var_t dt_ls = 0.0;  //! The elapsed time since the last data save 

    chrono::time_point<chrono::system_clock> t0, t1;

	time_t time_last_info = clock();
	time_t time_last_dump = clock();

	uint32_t n_print = 0;
	// Check if it is a continuation of a previous simulation
	size_t pos = opt->in_fn[INPUT_NAME_IC_DATA].find("solution.data_");
	if (string::npos > pos)
	{
		// 14 is the length of the 'solution.data_' string
		string str = opt->in_fn[INPUT_NAME_IC_DATA].substr(pos + 14, OUTPUT_ORDINAL_NUMBER_WIDTH);
		n_print = atoi(str.c_str()) + 1;
	}
    if (0 == n_print)
    {
        print_solution(n_print, opt, f, intgr, slog);
		f->calc_integral();
		f->print_integral(path_integral);	
    }

	/* 
	 * Main cycle
	 */
	while (f->t <= opt->param->simulation_length)
	{
        //if (DYN_MODEL_NBODY == opt->dyn_model)
        //{
        //    nbody* nb = (nbody*)f;
        //    nb->chk_coll(opt->param->threshold[THRESHOLD_RADII_ENHANCE_FACTOR]);
        //}

		// make the integration step, and measure the time it takes
        t0 = chrono::system_clock::now();
		f->dt = intgr->step();
        t1 = chrono::system_clock::now();
        chrono::duration<var_t> dt = t1 - t0;
        f->t_wc += dt.count();
        dt_ls += fabs(f->dt);

		if (opt->param->output_interval <= fabs(dt_ls))
		{
            dt_ls = 0.0;
            print_solution(n_print, opt, f, intgr, slog);
			f->calc_integral();
			f->print_integral(path_integral);	
            slog << tools::get_time_stamp() << " Data were printed" << endl;
        }

        if (opt->param->dump_dt < (clock() - time_last_dump) / (double)CLOCKS_PER_SEC)
        {
            time_last_dump = clock();
            print_dump(opt, f, intgr);
            slog << tools::get_time_stamp() << " Dump file was created" << endl;
        }

		if (opt->param->info_dt < (clock() - time_last_info) / (double)CLOCKS_PER_SEC) 
		{
			time_last_info = clock();
			print_info(opt, f, intgr, dt.count(), f->t_wc);
		}
	} /* while : main cycle */
    if (0.0 < fabs(dt_ls))
    {
        print_solution(n_print, opt, f, intgr, slog);
        f->calc_integral();
        f->print_integral(path_integral);
        slog << tools::get_time_stamp() << " Data were printed" << endl;
    }
} /* run_simulation() */

#if 0
	//matrix4_t m = {{10,8,10,7},{1,9,10,7},{8,8,4,7},{2,8,1,3}};
	//var4_t v = {4,10,1,5};
	//matrix4_t n = tools::calc_matrix_matrix_product(m,m);
	//matrix4_t l = tools::calc_matrix_transpose(n);
	//var4_t u = tools::calc_matrix_vector_product(m,v);

	//var3_t qv1 = {-8.98410183670583e-09,	-9.96565986672290e-09,	0};
	//var3_t qv2 = {0.790071598004156,   0.876376842255768,                   0};
	//var3_t qv3 = {0.790015741877597,   0.876342937103978 ,                  0};
	//var3_t pv1 = {1.35584617114928e-10,	-1.18154635090028e-10,	0};
	//var3_t pv2 = {-0.012137023259470 * 5.685826099573812e-09,   0.010261361613838 * 5.685826099573812e-09, 0};
	//var3_t pv3 = {-0.011709048488151 * 5.685826099573812e-09,   0.010519195691438 * 5.685826099573812e-09, 0};

	// DEBUG
	//pv1.x /= constants::Gauss, pv1.y /= constants::Gauss, pv1.z /= constants::Gauss;
	//pv2.x /= constants::Gauss, pv2.y /= constants::Gauss, pv2.z /= constants::Gauss;
	//pv3.x /= constants::Gauss, pv3.y /= constants::Gauss, pv3.z /= constants::Gauss;

	//var4_t Q1,Q2,P1,P2;
	//
	//tools::trans_to_threebody(qv1,pv1,qv2,pv2,qv3,pv3,Q1,P1,Q2,P2);
	//tools::trans_to_descartes(1,5.685826099573812e-09,5.685826099573812e-09,qv1,pv1,qv2,pv2,qv3,pv3,Q1,P1,Q2,P2);
#endif

int main(int argc, const char** argv, const char** env)
{
    chrono::time_point<chrono::system_clock> start = chrono::system_clock::now();
	ofstream* slog = NULL;

	try
	{
		options* opt = new options(argc, argv);
		string path_log = file::combine_path(opt->dir[DIRECTORY_NAME_OUT], opt->fn_prefix + opt->out_fn[OUTPUT_NAME_LOG]) + ".txt";
		slog = new ofstream(path_log.c_str(), ios::out | ios::app);
		if (!slog)
		{
			throw string("Cannot open " + path_log + ".");
		}
		file::log_start(*slog, argc, argv, env, opt->param->get_data(), opt->print_to_screen);

		ode *f = opt->create_model();
		integrator *intgr = opt->create_integrator(*f);

		if (opt->benchmark && DYN_MODEL_NBODY == opt->dyn_model)
		{
			run_benchmark(opt, f, intgr, *slog);
		}
		else
		{
			run_simulation(opt, f, intgr, *slog);
			if (opt->verbose && opt->print_to_screen)
			{
				printf("No. of passed steps: %20llu\n", intgr->get_n_passed_step());
				printf("No. of failed steps: %20llu\n", intgr->get_n_failed_step());
				printf("No. of  tried steps: %20llu\n", intgr->get_n_tried_step());
			}
		}
	} /* try */
	catch (const string& msg)
	{
		if (NULL != slog)
		{
			file::log_message(*slog, "ERROR: " + msg, false);
		}
		cerr << "ERROR: " << msg << endl;
	}
    chrono::time_point<chrono::system_clock> end = chrono::system_clock::now();
    chrono::duration<var_t> total_time = end - start;

	if (NULL != slog)
	{
		file::log_message(*slog, "Total time: " + tools::convert_var_t(total_time.count()) + " s.", false);
	}
	cout << "Total time: " << total_time.count() << " s." << endl;

	return (EXIT_SUCCESS);
}
