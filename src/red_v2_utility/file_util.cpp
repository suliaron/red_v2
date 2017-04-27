#include <algorithm>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <ostream>
#include <fstream>

#include "file_util.h"
#include "tools.h"
#include "util.h"

#include "constants.h"
#include "macro.h"
#include "type.h"

using namespace std;

namespace redutil2
{
    namespace file
    {
        string combine_path(const string& dir, const string& filename)
        {
	        if (0 < dir.size())
	        {
		        if (*(dir.end() - 1) != '/' && *(dir.end() - 1) != '\\')
		        {
			        return dir + '/' + filename;
		        }
		        else
		        {
			        return dir + filename;
		        }
	        }
	        else
	        {
		        return filename;
	        }
        }

        string get_filename(const string& path)
        {
	        string result;

	        if (path.size() > 0)
	        {
		        size_t pos = path.find_last_of("/\\");
		        result = path.substr(pos + 1);
	        }

	        return result;
        }

        string get_filename_without_ext(const string& path)
        {
	        string result;

	        if (path.size() > 0)
	        {
		        size_t pos = path.find_last_of("/\\");
		        result = path.substr(pos + 1);
		        pos = result.find_last_of('.');
		        result = result.substr(0, pos);
	        }

	        return result;
        }

        string get_directory(const string& path)
        {
	        string result;

	        if (0 < path.size())
	        {
		        size_t pos = path.find_last_of("/\\");
		        // If path does not contain / or \ than path does not contain any directory so return an empty string
		        if (pos == string::npos)
		        {
			        return result;
		        }
		        // Copy the directory part of path into result
		        result = path.substr(0, pos);
	        }

	        return result;
        }

        string get_extension(const string& path)
        {
	        string result;

	        if (path.size() > 0)
	        {
		        size_t pos = path.find_last_of('.');
		        result = path.substr(pos + 1);
	        }

	        return result;
        }

        data_rep_t get_data_repres(const string& path)
        {
	        data_rep_t repres;
	        string ext = file::get_extension(path);
	        if (     "txt" == ext)
	        {
		        repres = DATA_REPRESENTATION_ASCII;
	        }
	        else if ("dat" == ext)
	        {
		        repres = DATA_REPRESENTATION_BINARY;
	        }
	        else
	        {
		        throw string("The extension of the path must be either 'txt' or 'dat'.");
	        }

	        return repres;
        }

        uint32_t load_ascii_file(const string& path, string& result)
        {
	        uint32_t n_line = 0;

	        ifstream file(path.c_str(), ifstream::in);
	        if (file)
	        {
		        string str;
		        while (getline(file, str))
		        {
			        // delete everything after the comment '#' character and the '#'
			        str = tools::trim_comment(str);
			        str = tools::trim(str);
			        if (0 == str.length())
			        {
				        continue;
			        }
			        result += str;
			        result.push_back('\n');
			        n_line++;
		        } 	
		        file.close();
	        }
	        else
	        {
		        throw string("Cannot open " + path + ".");
	        }

	        return n_line;
        }

        void load_binary_file(const string& path, size_t n_data, var_t* data)
        {
	        ifstream file(path.c_str(), ios::in | ios::binary);
	        if (file)
	        {
		        file.seekg(0, file.end);     //N is the size of file in byte
		        size_t N = file.tellg();              
		        file.seekg(0, file.beg);
		        size_t size = n_data * sizeof(var_t);
		        if (size != N)
		        {
			        throw string("The file " + path + " has different number of data than expected.");
		        }
		        file.read(reinterpret_cast<char*>(data), size);
		        file.close();
	        }
	        else
	        {
		        throw string("Cannot open " + path + ".");
	        }
	        file.close();
        }

        void log_start(ostream& sout, int argc, const char** argv, const char** env, string params)
        {
	        sout << tools::get_time_stamp() << " starting " << argv[0] << endl;
	        sout << "Command line arguments: " << endl;
	        for (int i = 1; i < argc; i++)
	        {
		        sout << argv[i] << SEP;
	        }
	        sout << endl << endl;

	        while (*env)
	        {
		        string s = *env;
        #ifdef __GNUC__
		        // TODO
		        if(      s.find("HOSTNAME=") < s.length())
		        {
			        sout << s << endl;
		        }
		        else if (s.find("USER=") < s.length())
		        {
			        sout << s << endl;
		        }
		        else if (s.find("OSTYPE=") < s.length())
		        {
			        sout << s << endl;
		        }
        #else
		        if(      s.find("COMPUTERNAME=") < s.length())
		        {
			        sout << s << endl;
		        }
		        else if (s.find("USERNAME=") < s.length())
		        {
			        sout << s << endl;
		        }
		        else if (s.find("OS=") < s.length())
		        {
			        sout << s << endl;
		        }
        #endif
		        env++;
	        }
	        sout << endl;

	        sout << "Parameters:" << endl << params << endl;
        }

        void log_start(ostream& sout, int argc, const char** argv, const char** env, string params, bool print_to_screen)
        {
	        log_start(sout, argc, argv, env, params);
	        if (print_to_screen)
	        {
		        log_start(cout, argc, argv, env, params);
	        }
        }

        void log_message(ostream& sout, string msg, bool print_to_screen)
        {
	        sout << tools::get_time_stamp() << SEP << msg << endl;
	        //if (print_to_screen && sout != cout)
	        //{
	        //	cout << tools::get_time_stamp(false) << SEP << msg << endl;
	        //}
        }

        void print_data_info_record_ascii_RED(ofstream& sout, var_t t, var_t dt, n_objects_t* n_bodies)
        {
	        static uint32_t int_t_w  =  8;
	        static uint32_t var_t_w  = 25;

	        sout.precision(16);
	        sout.setf(ios::right);
	        sout.setf(ios::scientific);

	        sout << setw(var_t_w) << t << SEP
		         << setw(var_t_w) << dt << SEP;
	        for (uint32_t type = 0; type < BODY_TYPE_N; type++)
	        {
		        sout << setw(int_t_w) << n_bodies->get_n_active_by((body_type_t)type) << SEP;
	        }
        }

        void print_data_info_record_binary_RED(ofstream& sout, var_t t, var_t dt, n_objects_t* n_bodies)
        {
	        sout.write((char*)&(t), sizeof(var_t));
	        sout.write((char*)&(dt), sizeof(var_t));
	        for (uint32_t type = 0; type < BODY_TYPE_N; type++)
	        {
		        uint32_t n = n_bodies->get_n_active_by((body_type_t)type);
		        sout.write((char*)&n, sizeof(n));
	        }
        }

        void print_body_record_ascii_RED(ofstream &sout, string name, pp_disk_t::param_t *p, pp_disk_t::body_metadata_t *bmd, var4_t *r, var4_t *v)
        {
	        static int int_t_w  =  8;
	        static int var_t_w  = 25;

	        sout.precision(16);
	        sout.setf(ios::right);
	        sout.setf(ios::scientific);

	        sout << setw(     30) << name             << SEP
		         << setw(int_t_w) << bmd->id          << SEP
		         << setw(      2) << bmd->body_type   << SEP 
		         << setw(      2) << bmd->mig_type    << SEP
		         << setw(var_t_w) << bmd->mig_stop_at << SEP
		         << setw(var_t_w) << p->mass          << SEP
		         << setw(var_t_w) << p->radius        << SEP
		         << setw(var_t_w) << p->density       << SEP
		         << setw(var_t_w) << p->cd            << SEP
		         << setw(var_t_w) << r->x             << SEP
		         << setw(var_t_w) << r->y             << SEP
		         << setw(var_t_w) << r->z             << SEP
		         << setw(var_t_w) << v->x             << SEP
		         << setw(var_t_w) << v->y             << SEP
		         << setw(var_t_w) << v->z             << endl;

            sout.flush();
        }

        void print_body_record_binary_RED(ofstream &sout, string name, pp_disk_t::param_t *p, pp_disk_t::body_metadata_t *bmd, var4_t *r, var4_t *v)
        {
	        sout.write((char*)name.c_str(), 30*sizeof(char));
	        sout.write((char*)bmd, sizeof(pp_disk_t::body_metadata_t));
	        sout.write((char*)p,   sizeof(pp_disk_t::param_t));
	        sout.write((char*)r,   3*sizeof(var_t));
	        sout.write((char*)v,   3*sizeof(var_t));
        }

        void print_body_record_HIPERION(ofstream &sout, string name, var_t epoch, pp_disk_t::param_t *p, pp_disk_t::body_metadata_t *body_md, var4_t *r, var4_t *v)
        {
	        static int ids[4] = {0, 10, 20, 10000000};

	        static int int_t_w  =  8;
	        static int var_t_w  = 25;

	        sout.precision(16);
	        sout.setf(ios::right);
	        sout.setf(ios::scientific);

	        int id = 0;
	        // NOTE: this is not the beta parameter but derived from it
	        var_t _beta = 0.0;
	        switch (body_md->body_type)
	        {
	        case BODY_TYPE_STAR:
		        id = ids[0];
		        ids[0]++;
		        _beta = 0.0;
		        break;
	        case BODY_TYPE_GIANTPLANET:
		        id = ids[1];
		        ids[1]++;
		        _beta = 0.0;
		        break;
	        case BODY_TYPE_ROCKYPLANET:
		        id = ids[1];
		        ids[1]++;
		        _beta = 0.0;
		        break;
	        case BODY_TYPE_PROTOPLANET:
		        id = ids[2];
		        ids[2]++;
		        _beta = p->density;
		        break;
	        case BODY_TYPE_SUPERPLANETESIMAL:
		        break;
	        case BODY_TYPE_PLANETESIMAL:
		        id = ids[2];
		        ids[2]++;
		        _beta = p->density;
		        break;
	        case BODY_TYPE_TESTPARTICLE:
		        id = ids[3];
		        ids[3]++;
		        _beta = 1.0;
		        break;
	        default:
		        throw string("Parameter 'body_type' is out of range.");
	        }

	        var_t eps = 0.0;

	        sout << setw(int_t_w) << id      << SEP
		         << setw(var_t_w) << p->mass << SEP
		         << setw(var_t_w) << r->x    << SEP
		         << setw(var_t_w) << r->y    << SEP
		         << setw(var_t_w) << r->z    << SEP
		         << setw(var_t_w) << v->x / constants::Gauss << SEP
		         << setw(var_t_w) << v->y / constants::Gauss << SEP
		         << setw(var_t_w) << v->z / constants::Gauss << SEP
		         << setw(var_t_w) << eps     << SEP
		         << setw(var_t_w) << _beta   << endl;

            sout.flush();
        }

        void print_body_record_Emese(ofstream &sout, string name, var_t epoch, pp_disk_t::param_t *p, pp_disk_t::body_metadata_t *body_md, var4_t *r, var4_t *v)
        {
	        const char* body_type_name[] = 
	        {
		        "STAR",
		        "GIANTPLANET",
		        "ROCKYPLANET",
		        "PROTOPLANET",
		        "SUPERPLANETESIMAL",
		        "PLANETESIMAL",
		        "TESTPARTICLE",
	        };

	        static int int_t_w  = 25;
	        static int var_t_w  = 25;

	        sout.precision(16);
	        sout.setf(ios::left);
	        sout.setf(ios::scientific);

	        // NOTE: Emese start the ids from 0, red starts from 1.
	        sout << setw(int_t_w) << noshowpos << body_md->id - 1
		         << setw(     25) << name
		         << setw(     25) << noshowpos << body_type_name[body_md->body_type]
		         << setw(var_t_w) << showpos << r->x
		         << setw(var_t_w) << showpos << r->y
		         << setw(var_t_w) << showpos << r->z
		         << setw(var_t_w) << showpos << v->x
		         << setw(var_t_w) << showpos << v->y
		         << setw(var_t_w) << showpos << v->z
		         << setw(var_t_w) << showpos << p->mass
		         << setw(var_t_w) << showpos << p->radius << endl;

            sout.flush();
        }

        void print_oe_record(ofstream &sout, orbelem_t* oe)
        {
	        static int var_t_w  = 15;

	        sout.precision(6);
	        sout.setf(ios::right);
	        sout.setf(ios::scientific);

	        sout << setw(var_t_w) << oe->sma << SEP 
                 << setw(var_t_w) << oe->ecc << SEP 
                 << setw(var_t_w) << oe->inc << SEP 
                 << setw(var_t_w) << oe->peri << SEP 
                 << setw(var_t_w) << oe->node << SEP 
                 << setw(var_t_w) << oe->mean << endl;

	        sout.flush();
        }

        void print_oe_record(ofstream &sout, orbelem_t* oe, pp_disk_t::param_t *p)
        {
	        static int var_t_w  = 15;

	        sout.precision(6);
	        sout.setf(ios::right);
	        sout.setf(ios::scientific);

	        sout << setw(var_t_w) << oe->sma    << SEP 
                 << setw(var_t_w) << oe->ecc    << SEP 
                 << setw(var_t_w) << oe->inc    << SEP 
                 << setw(var_t_w) << oe->peri   << SEP 
                 << setw(var_t_w) << oe->node   << SEP 
                 << setw(var_t_w) << oe->mean   << SEP
                 << setw(var_t_w) << p->mass    << SEP
                 << setw(var_t_w) << p->radius  << SEP
                 << setw(var_t_w) << p->density << SEP
                 << setw(var_t_w) << p->cd      << endl;

	        sout.flush();
        }

        void print_oe_record(ofstream &sout, var_t epoch, orbelem_t* oe, pp_disk_t::param_t *p, pp_disk_t::body_metadata_t *bmd)
        {
	        static int var_t_w  = 15;
	        static int int_t_w  = 7;

	        sout.precision(6);
	        sout.setf(ios::right);
	        sout.setf(ios::scientific);

	        sout << setw(var_t_w) << epoch          << SEP 
		         << setw(int_t_w) << bmd->id        << SEP 
		         << setw(      2) << bmd->body_type << SEP 
                 << setw(var_t_w) << p->mass        << SEP
                 << setw(var_t_w) << p->radius      << SEP
                 << setw(var_t_w) << p->density     << SEP
                 << setw(var_t_w) << p->cd          << SEP
		         << setw(var_t_w) << oe->sma        << SEP 
                 << setw(var_t_w) << oe->ecc        << SEP 
                 << setw(var_t_w) << oe->inc        << SEP 
                 << setw(var_t_w) << oe->peri       << SEP 
                 << setw(var_t_w) << oe->node       << SEP 
                 << setw(var_t_w) << oe->mean       << endl;

	        sout.flush();
        }

        void load_data_info_record_ascii(ifstream& input, var_t& t, var_t& dt, n_objects_t** n_bodies)
        {
	        uint32_t ns, ngp, nrp, npp, nspl, npl, ntp;
	        ns = ngp = nrp = npp = nspl = npl = ntp = 0;

	        input >> t >> dt; 
	        input >> ns >> ngp >> nrp >> npp >> nspl >> npl >> ntp;

	        *n_bodies = new n_objects_t(ns, ngp, nrp, npp, nspl, npl, ntp);
        }

        void load_data_info_record_binary(ifstream& input, var_t& t, var_t& dt, n_objects_t** n_bodies)
        {
	        uint32_t ns, ngp, nrp, npp, nspl, npl, ntp;
	        ns = ngp = nrp = npp = nspl = npl = ntp = 0;

	        input.read((char*)&t, sizeof(var_t));
	        input.read((char*)&dt, sizeof(var_t));

	        input.read((char*)&ns,   sizeof(ns));
	        input.read((char*)&ngp,  sizeof(ngp));
	        input.read((char*)&nrp,  sizeof(nrp));
	        input.read((char*)&npp,  sizeof(npp));
	        input.read((char*)&nspl, sizeof(nspl));
	        input.read((char*)&npl,  sizeof(npl));
	        input.read((char*)&ntp,  sizeof(ntp));

	        *n_bodies = new n_objects_t(ns, ngp, nrp, npp, nspl, npl, ntp);
        }

        void load_data_record_ascii(ifstream& input, string& name, pp_disk_t::param_t *p, pp_disk_t::body_metadata_t *bmd, var4_t *r, var4_t *v)
        {
	        int	type = 0;
	        string	buffer;

	        // name
	        input >> buffer;
	        // The names must be less than or equal to 30 chars
	        if (buffer.length() > 30)
	        {
		        buffer = buffer.substr(0, 30);
	        }
	        name = buffer;

	        // id
	        input >> bmd->id;
	        // body type
	        input >> type;
	        bmd->body_type = static_cast<body_type_t>(type);
	        // migration type
	        input >> type;
	        bmd->mig_type = static_cast<migration_type_t>(type);
	        // migration stop at
	        input >> bmd->mig_stop_at;

	        // mass, radius density and stokes coefficient
	        input >> p->mass >> p->radius >> p->density >> p->cd;

	        // position
	        input >> r->x >> r->y >> r->z;
	        // velocity
	        input >> v->x >> v->y >> v->z;
	        r->w = v->w = 0.0;
        }

        void load_data_record_binary(ifstream& input, string& name, pp_disk_t::param_t *p, pp_disk_t::body_metadata_t *bmd, var4_t *r, var4_t *v)
        {
	        char buffer[30];
	        memset(buffer, 0, sizeof(buffer));

	        input.read(buffer,      30*sizeof(char));
	        input.read((char*)bmd,  1*sizeof(pp_disk_t::body_metadata_t));
	        input.read((char*)p,    1*sizeof(pp_disk_t::param_t));
	        input.read((char*)r,    3*sizeof(var_t));
	        input.read((char*)v,    3*sizeof(var_t));

	        name = buffer;
        }

        namespace tbp
        {
            void print_solution_info(std::ofstream& sout, var_t t, var_t dt, data_rep_t repres)
            {
	            switch (repres)
	            {
	            case DATA_REPRESENTATION_ASCII:
		            sout.precision(16);
		            sout.setf(ios::right);
		            sout.setf(ios::scientific);

		            sout << setw(VAR_T_W) << t << SEP         /* time of the record [day] (double)    */
			             << setw(VAR_T_W) << dt << endl;      /* next timestep [day]      (double)    */
		            sout.flush();
		            break;

	            case DATA_REPRESENTATION_BINARY:
		            sout.write((char*)&(t),  sizeof(var_t));  /* time of the record [day] (double)    */
		            sout.write((char*)&(dt), sizeof(var_t));  /* next timestep [day]      (double)    */
		            break;

	            default:
		            throw string("Parameter 'repres' is out of range.");
	            }
            }

            void print_solution_data(std::ofstream& sout, uint32_t n_obj, uint16_t n_ppo, uint16_t n_vpo, tbp_t::metadata_t* md, var_t* p, var_t* y, data_rep_t repres)
            {
	            switch (repres)
	            {
	            case DATA_REPRESENTATION_ASCII:
		            sout.precision(16);
		            sout.setf(ios::right);
		            sout.setf(ios::scientific);

		            for (uint32_t i = 0; i < n_obj; i++)
		            {
			            sout << setw(INT_T_W) << md[i].id << SEP;
			            // Print the parameters for each object
			            for (uint16_t j = 0; j < n_ppo; j++)
			            {
				            uint32_t param_idx = i * n_ppo + j;
				            sout << setw(VAR_T_W) << p[param_idx] << SEP;
			            }
			            // Print the variables for each object
			            for (uint16_t j = 0; j < n_vpo; j++)
			            {
				            uint32_t var_idx = i * n_vpo + j;
				            sout << setw(VAR_T_W) << y[var_idx];
				            if (j < n_vpo - 1)
				            {
					            sout << SEP;
				            }
				            else
				            {
					            sout << endl;
				            }
			            }
		            }
		            sout.flush();
		            break;

	            case DATA_REPRESENTATION_BINARY:
		            throw string("The print_result_binary() is not implemented.");
		            break;

	            default:
		            throw string("Parameter 'repres' is out of range.");
	            }
            }
        } /* namespace tbp  */

        namespace rtbp
        {
            void print_solution_data(std::ofstream& sout, uint32_t n_obj, uint16_t n_ppo, uint16_t n_vpo, tbp_t::metadata_t* md, var_t* p, var_t* y, uint16_t n_dim, data_rep_t repres)
        {
	        switch (repres)
	        {
	        case DATA_REPRESENTATION_ASCII:
		        sout.precision(16);
		        sout.setf(ios::right);
		        sout.setf(ios::scientific);

		        var_t  x1, xd1;  // 1D variables
		        var2_t x2, xd2;  // 2D variables
		        switch (n_dim)
		        {
		        case 1:
			        tools::rtbp1D::transform_u2x(y[0], y[1], x1, xd1);
			        break;

		        case 2:
			        {
			        var2_t u  = {y[0], y[1]};
			        var2_t up = {y[2], y[3]};
			        tools::rtbp2D::transform_u2x(u, x2);
			        tools::rtbp2D::transform_up2xd(u, up, xd2);
			        break;
			        }
		        case 3:
			        throw string("Parameter 'n_dim = 3' is not implemented.");
			        break;

		        default:
			        throw string("Parameter 'n_dim' is out of range.");
		        }

		        for (uint32_t i = 0; i < n_obj; i++)
		        {
			        sout << setw(INT_T_W) << md[i].id << SEP;               // 1 id of the object
			        // Print the parameters for each object
			        for (uint16_t j = 0; j < n_ppo; j++)                    // 2 parameter of the problem
			        {
				        uint32_t param_idx = i * n_ppo + j;
				        sout << setw(VAR_T_W) << p[param_idx] << SEP;
			        }
			        // Print the variables for each object
			        for (uint16_t j = 0; j < n_vpo; j++)                    // 3, 4  u1, u2 coordinate of the object
			        {                                                       // 5, 6  u1', u2' velocity of the object
				        uint32_t var_idx = i * n_vpo + j;
				        sout << setw(VAR_T_W) << y[var_idx] << SEP;
			        }
			        switch (n_dim)
			        {
			        case 1:
				        sout << setw(VAR_T_W) << y[2] << SEP                // 5 t physical time
					         << setw(VAR_T_W) << x1   << SEP                // 6 x physical coordinate of the object
					         << setw(VAR_T_W) << xd1  << endl;              // 7 x physical velocity of the object
				
				        break;

			        case 2:
				        sout << setw(VAR_T_W) << y[4] << SEP                // 7 t physical time
					         << setw(VAR_T_W) << x2.x  << SEP               // 8 x physical coordinate of the object
					         << setw(VAR_T_W) << x2.y  << SEP               // 9 y physical coordinate of the object
					         << setw(VAR_T_W) << xd2.x << SEP               //10 vx physical velocity of the object
					         << setw(VAR_T_W) << xd2.y << endl;             //11 vy physical velocity of the object
				        break;

			        case 3:
				        throw string("Parameter 'n_dim = 3' is not implemented.");
				        break;

			        default:
				        throw string("Parameter 'n_dim' is out of range.");
			        }
		        }
		        sout.flush();
		        break;

	        case DATA_REPRESENTATION_BINARY:
		        throw string("The print_result_binary() is not implemented.");
		        break;

	        default:
		        throw string("Parameter 'repres' is out of range.");
	        }
        }
        } /* namespace rtbp */

        namespace nbp
        {
            void print_solution_info(ofstream& sout, var_t t, var_t dt, var_t dt_wc, uint32_t n_obj, data_rep_t repres)
            {
	            switch (repres)
	            {
	            case DATA_REPRESENTATION_ASCII:
		            sout.precision(16);
		            sout.setf(ios::right);
		            sout.setf(ios::scientific);

		            sout << setw(VAR_T_W) << t << SEP            /* time of the record [day] (double) */
                         << setw(VAR_T_W) << dt << SEP            /* next timestep [day]      (double) */
                         << setw(VAR_T_W) << dt_wc << SEP         /* the ellapsed time since the begining of the integration [sec]      (double) */
                         << setw(INT_T_W) << n_obj << endl;       /* number of bodies         (uint32) */
		            sout.flush();
		            break;

	            case DATA_REPRESENTATION_BINARY:
		            sout.write((char*)&t,     sizeof(var_t));     /* time of the record [day] (double) */
                    sout.write((char*)&dt,    sizeof(var_t));     /* next timestep [day]      (double) */
                    sout.write((char*)&dt_wc, sizeof(var_t));     /* the ellapsed time since the begining of the integration [sec]      (double) */
                    sout.write((char*)&n_obj, sizeof(uint32_t));  /* number of bodies         (uint32) */
		            break;

	            default:
		            throw string("Parameter 'repres' is out of range.");
	            }

            }

            void print_solution_data(std::ofstream& sout, uint32_t n_obj, uint16_t n_ppo, uint16_t n_vpo, nbp_t::metadata_t* md, nbp_t::param_t* p, var_t* y, data_rep_t repres)
            {
	            switch (repres)
	            {
	            case DATA_REPRESENTATION_ASCII:
		            {
		            sout.setf(ios::right);
		            sout.setf(ios::scientific);

		            uint32_t offset = 3 * n_obj;
		            for (uint32_t i = 0; i < n_obj; i++)
		            {
                        sout << setw(INT_T_W) << md[i].id << SEP
                             << setw(2) << (int)md[i].body_type << SEP
                             << setw(2) << md[i].active << SEP
                             << setw(2) << (int)md[i].mig_type << SEP
                             << setprecision(6) << setw(FLT_T_W) << md[i].mig_stop_at << SEP
                             << setw(2) << md[i].unused1 << SEP
                             << setw(2) << md[i].unused2 << SEP
                             << setw(2) << md[i].unused3 << SEP;

                        // Print the parameters for each object
                        sout << setprecision(16) << setw(VAR_T_W) << p[i].density << SEP
                             << setw(VAR_T_W) << p[i].mass << SEP
                             << setw(VAR_T_W) << p[i].radius << SEP;

			            // Print the variables for each object
			            uint32_t j = 3 * i;
			            sout << setw(VAR_T_W) << y[j  ] << SEP            /* x-coordinate */
				             << setw(VAR_T_W) << y[j+1] << SEP            /* y-coordinate */
				             << setw(VAR_T_W) << y[j+2] << SEP            /* z-coordinate */
 				             << setw(VAR_T_W) << y[offset + j  ] << SEP   /* vx-velocity  */
 				             << setw(VAR_T_W) << y[offset + j+1] << SEP   /* vy-velocity  */
 				             << setw(VAR_T_W) << y[offset + j+2] << endl; /* vz-velocity  */
		            }
		            }
		            sout.flush();
		            break;

	            case DATA_REPRESENTATION_BINARY:
		            sout.write((char*)md,    n_obj * sizeof(nbp_t::metadata_t));
		            sout.write((char*)p,     n_obj * sizeof(nbp_t::param_t));
		            sout.write((char*)y, 6 * n_obj * sizeof(var_t));
		            break;

	            default:
		            throw string("Parameter 'repres' is out of range.");
	            }

            }
        } /* nbp */

    } /* file */
} /* redutil2 */
