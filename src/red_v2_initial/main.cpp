#include <cfloat>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>
#include <time.h>       /* time                      */

#include "distribution.h"
#include "initial_type.h"
#include "util_init.h"

#include "constants.h"
#include "type.h"
#include "redutil2.h"

using namespace std;
using namespace redutil2;

namespace ephemeris_major_planets
{
    /*
    Symbol meaning [1 au=149597870.700 km, 1 day=86400.0 s]:
    JDCT     Epoch Julian Date, Coordinate Time
    EC     Eccentricity, e
    QR     Periapsis distance, q (AU)
    IN     Inclination w.r.t xy-plane, i (degrees)
    OM     Longitude of Ascending Node, OMEGA, (degrees)
    W      Argument of Perifocus, w (degrees)
    Tp     Time of periapsis (Julian day number)
    N      Mean motion, n (degrees/day)
    MA     Mean anomaly, M (degrees)
    TA     True anomaly, nu (degrees)
    A      Semi-major axis, a (AU)
    AD     Apoapsis distance (AU)
    PR     Sidereal orbit period (day)
    */
    namespace date_20150511
    {
        string mercury_oe = "2457153.500000000, A.D. 2015-May-11 00:00:00.0000,  2.056283170656704E-01,  3.075005270000471E-01,  7.004033081870772E+00,  4.831135275416084E+01,  2.917018579009503E+01,  2.457132325703553E+06,  4.092332803985572E+00,  8.665226794919347E+01,  1.098801777958040E+02,  3.870990540148295E-01,  4.666975810296119E-01,  8.796938500441402E+01,";
        string venus_oe = "2457153.500000000, A.D. 2015-May-11 00:00:00.0000,  6.761247965610893E-03,  7.184381867467915E-01,  3.394467689974794E+00,  7.663950855912078E+01,  5.463888630600258E+01,  2.457130873225775E+06,  1.602141893419786E+00,  3.625130289912041E+01,  3.671259128576724E+01,  7.233287920706470E-01,  7.282193973945026E-01,  2.246991989152577E+02,";
        string earth_oe = "2457153.500000000, A.D. 2015-May-11 00:00:00.0000,  1.706151645376220E-02,  9.828818510541610E-01,  2.832604859699792E-03,  1.991276585462631E+02,  2.644902808333621E+02,  2.457027106783887E+06,  9.856943344121811E-01,  1.245850770307868E+02,  1.261752194211184E+02,  9.999423845001243E-01,  1.017002917946088E+00,  3.652247836188346E+02,";
        string mars_oe = "2457153.500000000, A.D. 2015-May-11 00:00:00.0000,  9.345598920376896E-02,  1.381200747391636E+00,  1.848403968432629E+00,  4.951276588949865E+01,  2.865318030175783E+02,  2.457003817487961E+06,  5.240856597101560E-01,  7.844645806929105E+01,  8.912748635050312E+01,  1.523589291796774E+00,  1.665977836201911E+00,  6.869106096112168E+02,";
        string jupiter_oe = "2457153.500000000, A.D. 2015-May-11 00:00:00.0000,  4.895926045495331E-02,  4.946812500213837E+00,  1.303927595561599E+00,  1.005220570364609E+02,  2.736235124666692E+02,  2.455634711874421E+06,  8.312310735603498E-02,  1.262463884135522E+02,  1.306085971088266E+02,  5.201472759810757E+00,  5.456133019407678E+00,  4.330925676996638E+03,";
        string saturn_oe = "2457153.500000000, A.D. 2015-May-11 00:00:00.0000,  5.370451604891303E-02,  9.040758523537329E+00,  2.486273622454862E+00,  1.136214918787423E+02,  3.402344971420399E+02,  2.452848192133344E+06,  3.338097180639689E-02,  1.437153605146956E+02,  1.471679382461214E+02,  9.553843040430950E+00,  1.006692755732457E+01,  1.078458716204937E+04,";
        string uranus_oe = "2457153.500000000, A.D. 2015-May-11 00:00:00.0000,  4.928987553618186E-02,  1.821116226013173E+01,  7.717997843738144E-01,  7.400078543435393E+01,  9.677390970582711E+01,  2.470049485736452E+06,  1.175653216407963E-02,  2.083879289018873E+02,  2.058417203584891E+02,  1.915532588905842E+01,  2.009948951798512E+01,  3.062127462211412E+04,";
        string neptune_oe = "2457153.500000000, A.D. 2015-May-11 00:00:00.0000,  8.258758489207289E-03,  2.973094026694364E+01,  1.765657591232885E+00,  1.317164628477979E+02,  2.932421741299731E+02,  2.471558726412109E+06,  6.004818034445084E-03,  2.734992366503029E+02,  2.725540643806312E+02,  2.997852567031729E+01,  3.022611107369094E+01,  5.995185831359972E+04,";
    } /* date_20150511 */

    namespace date_
    {
        string mercury_oe = "";
        string venus_oe = "";
        string earth_oe = "";
        string mars_oe = "";
        string jupiter_oe = "";
        string saturn_oe = "";
        string uranus_oe = "";
        string neptune_oe = "";
    } /* date_ */
} /* ephemeris */

static string fn_info;
static string fn_data;


static void print_data(string& path, uint32_t n_data, var_t* data)
{
    printf("Writing %s to disk.\n", path.c_str());

    ofstream sout(path.c_str(), ios::out);
    if (sout)
    {
        sout.setf(ios::right);
        sout.setf(ios::scientific);

        for (uint32_t i = 0; i < n_data; i++)
        {
            sout << setprecision(16) << setw(25) << data[i] << endl;
        }
    }
    else
    {
        throw string("Cannot open " + path + ".");
    }
    sout.close();
}

static void print_start_files(string& dir)
{
    string path = file::combine_path(dir, "start_files.txt");
    printf("Writing %s to disk.\n", path.c_str());

    ofstream sout(path.c_str(), ios::out);
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

namespace model
{
	namespace tbp1D
    {
        // parameter of the problem
        tbp_t::param_t p;
       	// Epoch for the initial condition
        var_t t0;
        // Initial conditions
        var_t* y;
        // Metadata of the object
        tbp_t::metadata_t md;
        // Initial stepsize for the integrator
        var_t dt0;
        
        void print(string& dir, string& filename)
        {
           	ofstream sout;

            string path = file::combine_path(dir, fn_info);
	        printf("Writing %s to disk.\n", path.c_str());
       		sout.open(path.c_str(), ios::out);
            if (sout)
            {
				file::tbp::print_solution_info(sout, t0, dt0, DATA_REPRESENTATION_ASCII);
            }
            else
            {
                throw string("Cannot open " + path + ".");
            }
            sout.close();

	        path = file::combine_path(dir, fn_data);
	        printf("Writing %s to disk.\n", path.c_str());
       		sout.open(path.c_str(), ios::out);
            if (sout)
            {
                file::tbp::print_solution_data(sout, 1, 1, 2, &md, (var_t*)&p, y, DATA_REPRESENTATION_ASCII);
            }
            else
            {
                throw string("Cannot open " + path + ".");
            }
            sout.close();

			print_start_files(dir);
        } /* print() */

        void create(string& dir, string& filename)
        {
	        /*
	         * The units are:
	         *     Unit name         | Unit symbol | Quantity name
	         *     -----------------------------------------------
	         *     Astronomical unit |          AU | length
	         *     Solar mass        |           S | mass
	         *     Mean solar day    |           D | time
	         */
            
            ALLOCATE_HOST_VECTOR((void**)&(y), 2 * sizeof(var_t));

			// Set the parameters of the problem
            p.mu  = constants::Gauss2 * (1.0 + 1.0);
			// Set the initial conditions at t0
            t0   = 0.0;
            y[0] = 1.0;    /* x0  */
            y[1] = 0.0;    /* vx0 */
			// Set the object metadata
            md.id = 1;
			// Set the initial stepsize for the integrator (should be modell dependent)
            dt0   = 1.0e-4;

            print(dir, filename);

            FREE_HOST_VECTOR((void **)&(y));
        }
    } /* namespace tbp1D */

	namespace tbp2D
	{
        // parameter of the problem
        tbp_t::param_t p;
       	// Epoch for the initial condition
        var_t t0;
        // Initial conditions
        var_t* y;
        // Metadata of the object
        tbp_t::metadata_t md;
        // Initial stepsize for the integrator
        var_t dt0;
        
        void print(string& dir, string& filename)
        {
           	ofstream sout;

            string path = file::combine_path(dir, fn_info);
	        printf("Writing %s to disk.\n", path.c_str());
       		sout.open(path.c_str(), ios::out);
            if (sout)
            {
                file::tbp::print_solution_info(sout, t0, dt0, DATA_REPRESENTATION_ASCII);
            }
            else
            {
                throw string("Cannot open " + path + ".");
            }
            sout.close();

	        path = file::combine_path(dir, fn_data);
	        printf("Writing %s to disk.\n", path.c_str());
       		sout.open(path.c_str(), ios::out);
            if (sout)
            {
                file::tbp::print_solution_data(sout, 1, 1, 4, &md, (var_t*)&p, y, DATA_REPRESENTATION_ASCII);
            }
            else
            {
                throw string("Cannot open " + path + ".");
            }
            sout.close();

			print_start_files(dir);
        } /* print() */

        void create(string& dir, string& filename)
        {
	        /*
	         * The units are:
	         *     Unit name         | Unit symbol | Quantity name
	         *     -----------------------------------------------
	         *     Astronomical unit |          AU | length
	         *     Solar mass        |           S | mass
	         *     Mean solar day    |           D | time
	         */
            
            ALLOCATE_HOST_VECTOR((void**)&(y), 4 * sizeof(var_t));

			// Set the parameter of the problem
            p.mu = constants::Gauss2 * (1.0 + 1.0);
			// Set the initial conditions at t0
            t0   = 0.0;
			// Set the initial orbital elements
			orbelem_t oe = {1.0, 0.1, 0.0, 0.0, 0.0, 0.0};
			var3_t r0 = {0, 0, 0};
			var3_t v0 = {0, 0, 0};
			// Calculate the initial position and velocity vectors
			tools::calc_phase(p.mu, &oe, &r0, &v0);
			y[0] = r0.x;
			y[1] = r0.y;
			y[2] = v0.x;
            y[3] = v0.y;
			// Set the object metadata
            md.id = 1;
			// Set the initial stepsize for the integrator (should be modell dependent)
            dt0   = 1.0e-4;

            print(dir, filename);

            FREE_HOST_VECTOR((void **)&(y));
        }
	} /* namespace tbp2D */

	namespace rtbp1D
    {
        // parameter of the problem
        tbp_t::param_t p;
       	// Value of the independent variable for the initial condition
        var_t s0;
        // Initial conditions
        var_t* y;
        // Metadata of the object
        tbp_t::metadata_t md;
        // Initial stepsize for the integrator
        var_t ds0;
        
        void print(string& dir, string& filename)
        {
           	ofstream sout;

            string path = file::combine_path(dir, fn_info);
	        printf("Writing %s to disk.\n", path.c_str());
       		sout.open(path.c_str(), ios::out);
            if (sout)
            {
				file::tbp::print_solution_info(sout, s0, ds0, DATA_REPRESENTATION_ASCII);
            }
            else
            {
                throw string("Cannot open " + path + ".");
            }
            sout.close();

	        path = file::combine_path(dir, fn_data);
	        printf("Writing %s to disk.\n", path.c_str());
       		sout.open(path.c_str(), ios::out);
            if (sout)
            {
                file::rtbp::print_solution_data(sout, 1, 1, 2, &md, (var_t*)&p, y, 1, DATA_REPRESENTATION_ASCII);
            }
            else
            {
                throw string("Cannot open " + path + ".");
            }
            sout.close();

			print_start_files(dir);
        } /* print() */

        void create(string& dir, string& filename)
        {
	        /*
	         * The units are:
	         *     Unit name         | Unit symbol | Quantity name
	         *     -----------------------------------------------
	         *     Astronomical unit |          AU | length
	         *     Solar mass        |           S | mass
	         *     Mean solar day    |           D | time
	         */
            
            ALLOCATE_HOST_VECTOR((void**)&(y), 3 * sizeof(var_t));

			// Set the parameters of the problem
            p.mu  = constants::Gauss2 * (1.0 + 1.0);
			// Set and compute the initial conditions at s0 (fictitious time)
            s0   = 0.0;
			var_t x0   = 1.0;
			var_t vx0  = 0.0;
			tools::rtbp1D::transform_x2u(x0, vx0, y[0], y[1]);
			var_t _x0  = 0.0;
			var_t _vx0 = 0.0;
			tools::rtbp1D::transform_u2x(y[0], y[1], _x0, _vx0);

			y[2] = 0.0;   // t0
			// Set the object metadata
            md.id = 1;
			// Set the initial stepsize for the integrator (should be model dependent)
            ds0   = 1.0e-4;

            print(dir, filename);

            FREE_HOST_VECTOR((void **)&(y));
        } /* create() */
    } /* namespace tbp1D */

	namespace rtbp2D
    {
        // parameter of the problem
        tbp_t::param_t p;
       	// Value of the independent variable for the initial condition
        var_t s0;
        // Initial conditions
        var_t* y;
        // Metadata of the object
        tbp_t::metadata_t md;
        // Initial stepsize for the integrator
        var_t ds0;
        
		void print(string& dir, string& filename)
        {
           	ofstream sout;

            string path = file::combine_path(dir, fn_info);
	        printf("Writing %s to disk.\n", path.c_str());
       		sout.open(path.c_str(), ios::out);
            if (sout)
            {
                file::tbp::print_solution_info(sout, s0, ds0, DATA_REPRESENTATION_ASCII);
            }
            else
            {
                throw string("Cannot open " + path + ".");
            }
            sout.close();

	        path = file::combine_path(dir, fn_data);
	        printf("Writing %s to disk.\n", path.c_str());
       		sout.open(path.c_str(), ios::out);
            if (sout)
            {
                file::rtbp::print_solution_data(sout, 1, 1, 4, &md, (var_t*)&p, y, 2, DATA_REPRESENTATION_ASCII);
            }
            else
            {
                throw string("Cannot open " + path + ".");
            }
            sout.close();

			print_start_files(dir);
        } /* print() */

        void create(string& dir, string& filename)
        {
	        /*
	         * The units are:
	         *     Unit name         | Unit symbol | Quantity name
	         *     -----------------------------------------------
	         *     Astronomical unit |          AU | length
	         *     Solar mass        |           S | mass
	         *     Mean solar day    |           D | time
	         */
            
            ALLOCATE_HOST_VECTOR((void**)&(y), 5 * sizeof(var_t));

			// Set the parameters of the problem
            p.mu = constants::Gauss2 * (1.0 + 1.0);
			// Set and compute the initial conditions at s0 (fictitious time)
            s0   = 0.0;
			// Set the initial orbital elements
			orbelem_t oe = {1.0, 0.8, 0.0, 0.0, 0.0, 0.0};
			var3_t r0 = {0, 0, 0};
			var3_t v0 = {0, 0, 0};
			// Calculate the initial position and velocity vectors
			tools::calc_phase(p.mu, &oe, &r0, &v0);

			var2_t x0  = {r0.x, r0.y};
			var2_t xd0 = {v0.x, v0.y};
			var2_t u0  = {0, 0};
			var2_t up0 = {0, 0};
			tools::rtbp2D::transform_x2u(x0, u0);
            tools::rtbp2D::transform_xd2up(u0, xd0, up0);

			y[0] = u0.x;
			y[1] = u0.y;
			y[2] = up0.x;
			y[3] = up0.y;
			y[4] = 0.0;   // t0
			// Set the object metadata
            md.id = 1;
			// Set the initial stepsize for the integrator (should be model dependent)
            ds0   = 1.0e-4;

            print(dir, filename);

            FREE_HOST_VECTOR((void **)&(y));
        } /* create() */
    } /* namespace tbp2D */

    /*
    * The units are:
    *     Unit name         | Unit symbol | Quantity name
    *     -----------------------------------------------
    *     Astronomical unit |          AU | length
    *     Solar mass        |           S | mass
    *     Mean solar day    |           D | time
    */
	namespace nbody
	{
        // parameters of the problem (masses)
        nbp_t::param_t* p = NULL;
       	// Epoch for the initial condition
        var_t t0 = 0.0;
        // Initial conditions
        var_t* y = NULL;
        orbelem_t* oe = NULL;
        // Metadata of the object
        nbp_t::metadata_t* md = NULL;
        // Initial stepsize for the integrator
        var_t dt0 = 0.0;

		var_t random(var_t x0, var_t x1)
		{
			return (x0 + ((var_t)rand() / RAND_MAX) * (x1 - x0));
		}

        void print(string& dir, string& filename, uint32_t n_obj)
        {
           	ofstream sout;

            string path = file::combine_path(dir, fn_info);
	        printf("Writing %s to disk.\n", path.c_str());
       		sout.open(path.c_str(), ios::out);
            if (sout)
            {
                file::nbp::print_solution_info(sout, t0, dt0, 0, n_obj, DATA_REPRESENTATION_ASCII);
            }
            else
            {
                throw string("Cannot open " + path + ".");
            }
            sout.close();

	        path = file::combine_path(dir, fn_data);
	        printf("Writing %s to disk.\n", path.c_str());
       		sout.open(path.c_str(), ios::out);
            if (sout)
            {
                file::nbp::print_solution_data(sout, n_obj, 3, 6, md, p, y, DATA_REPRESENTATION_ASCII);
            }
            else
            {
                throw string("Cannot open " + path + ".");
            }
            sout.close();

			print_start_files(dir);
        } /* print() */

		void create(string& dir, string& filename, uint32_t n_obj)
        {
			uint32_t n_var = 6 * n_obj;
			uint32_t n_par = n_obj;
            ALLOCATE_HOST_VECTOR((void**)&y,  n_var * sizeof(var_t));
			ALLOCATE_HOST_VECTOR((void**)&p,  n_par * sizeof(nbp_t::param_t));
			ALLOCATE_HOST_VECTOR((void**)&md, n_obj * sizeof(nbp_t::metadata_t));

			srand((unsigned int)time(NULL));
			// Set the parameters of the problem
			for (uint32_t i = 0; i < n_obj; i++)
			{
				p[i].mass = random(0.1, 1.0);
				// Set the object metadata
	            md[i].id = i + 1;
			}
			// Set the initial conditions at t0
            t0   = 0.0;
			uint32_t offset = 3 * n_obj;
			for (uint32_t i = 0; i < n_obj; i++)
			{
				uint32_t j = 3 * i;
				y[j  ] = random(-1.0e1, 1.0e1);
				y[j+1] = random(-1.0e1, 1.0e1);
				y[j+2] = random(-1.0e1, 1.0e1);

				y[offset + j  ] = random(-1.0e-5, 1.0e-5);
				y[offset + j+1] = random(-1.0e-5, 1.0e-5);
				y[offset + j+2] = random(-1.0e-5, 1.0e-5);
			}
			
			// Set the initial stepsize for the integrator (should be model dependent)
            dt0   = 1.0e-4;

            print(dir, filename, n_obj);

            FREE_HOST_VECTOR((void **)&y);
			FREE_HOST_VECTOR((void **)&p);
			FREE_HOST_VECTOR((void **)&md);
        }

        void ceate_disk(string& dir, string& filename, uint32_t n_obj)
        {
            uint32_t seed = (uint32_t)time(NULL);
            cout << "The seed number is " << seed << endl;
            // The pseudo-random number generator is initialized using the argument passed as seed.
            // Used by the subsequent rand() function calls
            srand(seed);

            // Epoch for the disk's state
            t0 = 0.0;
            const var_t m0 = 1.0;  //! Mass of the central star
            const var_t R0 = 1.0 * constants::SolarRadiusToAu;
            oe_dist_t oe_d;
            pp_dist_t pp_d;

            oe_d.item[ORBELEM_NAME_SMA] = new uniform_distribution(rand(), 1.0, 2.0);
            // If the distribution is NULL than the corresponding orbital element is zero
            oe_d.item[ORBELEM_NAME_ECC] = new uniform_distribution(rand(), 0.0, 0.3);
            //oe_d.item[ORBELEM_NAME_INC] = new uniform_distribution(rand(), 0.0, 0.0);
            oe_d.item[ORBELEM_NAME_PERI] = new uniform_distribution(rand(), 0.0, 2.0 *PI);
            //oe_d.item[ORBELEM_NAME_NODE] = new uniform_distribution(rand(), 0.0, 0.0);
            oe_d.item[ORBELEM_NAME_MEAN] = new uniform_distribution(rand(), 0.0, 2.0 *PI);

            pp_d.item[PP_NAME_DENSITY] = new uniform_distribution(rand(), 2.0 * constants::GramPerCm3ToSolarPerAu3, 2.0 * constants::GramPerCm3ToSolarPerAu3);
            pp_d.item[PP_NAME_MASS   ] = new uniform_distribution(rand(), 1.0e2*constants::CeresToSolar, 1.0e2*constants::CeresToSolar);

            // Increase n_obj by one to include the central body
            n_obj++;
            uint32_t n_var = 6 * n_obj;
            uint32_t n_par = (sizeof(nbp_t::param_t) / sizeof(var_t)) * n_obj;
            ALLOCATE_HOST_VECTOR((void**)&oe, n_obj * sizeof(orbelem_t));
            ALLOCATE_HOST_VECTOR((void**)&y,  n_var * sizeof(var_t));
            ALLOCATE_HOST_VECTOR((void**)&p,  n_par * sizeof(nbp_t::param_t));
            ALLOCATE_HOST_VECTOR((void**)&md, n_obj * sizeof(nbp_t::metadata_t));

            {
                nbp_t::param_t param = { 0.0, 0.0, 0.0 };
                nbp_t::metadata_t body_md = { 0, 0, 0, true, false, 0.0 };
                orbelem_t _oe = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
                var3_t rVec = { 0.0, 0.0, 0.0 };
                var3_t vVec = { 0.0, 0.0, 0.0 };

                // The id of each body must be larger than 0 in order to indicate inactive body with negative id 
                // (ie. zero is not good)
                uint32_t bodyIdx = 0;
                uint32_t bodyId = 1;
                var_t min_P = DBL_MAX;

                for (uint32_t i = 0; i < n_obj; i++, bodyIdx++, bodyId++)
                {
                    body_md.id = bodyId;
                    body_md.active = true;
                    body_md.mig_type = MIGRATION_TYPE_NO;
                    body_md.mig_stop_at = 0.0;

                    body_md.unused1 = body_md.unused2 = body_md.unused3 = false;

                    // The central star
                    if (1 == bodyId)
                    {
                        body_md.body_type = BODY_TYPE_STAR;

                        param.mass = m0;
                        param.radius = R0;
                        param.density = tools::calc_density(param.mass, param.radius);
                        uint32_t offset = 3 * i;
                        y[offset + 0] = y[offset + 1] = y[offset + 2] = 0.0;
                        offset += 3 * n_obj;
                        y[offset + 0] = y[offset + 1] = y[offset + 2] = 0.0;
                    }
                    else
                    {
                        body_md.body_type = BODY_TYPE_PROTOPLANET;

                        generate_oe(&oe_d, _oe);
                        generate_pp(&pp_d, param);
                        if (0x0 == pp_d.item[PP_NAME_MASS])
                        {
                            param.mass = tools::calc_mass(param.radius, param.density);
                        }
                        if (0x0 == pp_d.item[PP_NAME_RADIUS])
                        {
                            param.radius = tools::calc_radius(param.mass, param.density);
                        }
                        if (0x0 == pp_d.item[PP_NAME_DENSITY])
                        {
                            param.density = tools::calc_density(param.mass, param.radius);
                        }

                        var_t mu = K2 * (m0 + param.mass);
                        tools::calc_phase(mu, &_oe, &rVec, &vVec);
                        uint32_t offset = 3 * i;
                        y[offset + 0] = rVec.x; y[offset + 1] = rVec.y; y[offset + 2] = rVec.z;
                        offset += 3 * n_obj;
                        y[offset + 0] = vVec.x; y[offset + 1] = vVec.y; y[offset + 2] = vVec.z;

                        var_t P = tools::calc_orbital_period(mu, _oe.sma);
                        if (min_P > P)
                        {
                            min_P = P;
                        }
                    }

                    md[bodyIdx] = body_md;
                    p[bodyIdx]  = param;
                    oe[bodyIdx] = _oe;
                } /* for */
                // Set the initial stepsize for the integrator
                dt0 = min_P / 50.0;
            }

            // Transform the coordinates into barycentric
            var3_t* r = (var3_t*)y;
            var3_t* v = (var3_t*)(y + 3 * n_obj);
            tools::nbp::transform_to_bc(n_obj, p, r, v);

            print(dir, filename, n_obj);

            FREE_HOST_VECTOR((void **)&md);
            FREE_HOST_VECTOR((void **)&p);
            FREE_HOST_VECTOR((void **)&y);
            FREE_HOST_VECTOR((void **)&oe);
        }

		void create(string& dir, string& filename)
        {
            uint32_t n_obj = 2;
			uint32_t n_var = 6 * n_obj;
			uint32_t n_par = n_obj;
            ALLOCATE_HOST_VECTOR((void**)&y,  n_var * sizeof(var_t));
			ALLOCATE_HOST_VECTOR((void**)&p,  n_par * sizeof(nbp_t::param_t));
			ALLOCATE_HOST_VECTOR((void**)&md, n_obj * sizeof(nbp_t::metadata_t));

			srand((unsigned int)time(NULL));
			// Set the parameters of the problem
			p[0].mass = 1.0;
			p[1].mass = 1.0;
			// Set the object metadata
            md[0].id = 1;
            md[1].id = 2;

			// Set the initial conditions at t0
            t0 = 0.0;
			const var_t mu = K2*(p[0].mass + p[1].mass);
			orbelem_t oe = {1.0, 0.9, 0.0, 0.0, 0.0, 0.0};
			var3_t r = {0, 0, 0};
			var3_t v = {0, 0, 0};
			tools::calc_phase(mu, &oe, &r, &v);

			var_t s = p[1].mass/(p[0].mass + p[1].mass);
			var3_t r1 = {-s*r.x, -s*r.y, -s*r.z};
			var3_t v1 = {-s*v.x, -s*v.y, -s*v.z};

			s = p[0].mass/(p[0].mass + p[1].mass);
			var3_t r2 = {s*r.x, s*r.y, s*r.z};
			var3_t v2 = {s*v.x, s*v.y, s*v.z};

			y[0] = r1.x, y[1] = r1.y, y[2] = r1.z;
			y[3] = r2.x, y[4] = r2.y, y[5] = r2.z;

			y[6] = v1.x, y[ 7] = v1.y, y[ 8] = v1.z;
			y[9] = v2.x, y[10] = v2.y, y[11] = v2.z;

			// Set the initial stepsize for the integrator (should be model dependent)
            dt0   = 1.0e-4;

            print(dir, filename, n_obj);

            FREE_HOST_VECTOR((void **)&y);
			FREE_HOST_VECTOR((void **)&p);
			FREE_HOST_VECTOR((void **)&md);
        }

        void create_grav_focusing(string& dir, string& filename)
        {
            uint32_t n_obj = 2;
            uint32_t n_var = 6 * n_obj;
            uint32_t n_par = n_obj;
            ALLOCATE_HOST_VECTOR((void**)&y, n_var * sizeof(var_t));
            ALLOCATE_HOST_VECTOR((void**)&p, n_par * sizeof(nbp_t::param_t));
            ALLOCATE_HOST_VECTOR((void**)&md, n_obj * sizeof(nbp_t::metadata_t));

            srand((unsigned int)time(NULL));
            // Set the parameters of the problem
            p[0].mass = 1.0 * constants::EarthToSolar, p[0].density = 2.0 * constants::GramPerCm3ToSolarPerAu3;
            p[0].radius = tools::calc_radius(p[0].mass, p[0].density);
            p[1].mass = 1.0e-2 * constants::EarthToSolar, p[1].density = 2.0 * constants::GramPerCm3ToSolarPerAu3;
            p[1].radius = tools::calc_radius(p[1].mass, p[1].density);

            // Set the object metadata
            md[0].id = 1, md[0].active = true, md[0].body_type = BODY_TYPE_STAR, md[0].mig_stop_at = 0.0, md[0].mig_type = MIGRATION_TYPE_NO;
            md[1].id = 2, md[1].active = true, md[1].body_type = BODY_TYPE_STAR, md[1].mig_stop_at = 0.0, md[1].mig_type = MIGRATION_TYPE_NO;

            // Set the initial conditions at t0
            t0 = 0.0;
            var3_t r1 = { 0, 0, 0 };
            var3_t v1 = { 0, 0, 0 };

            var_t b = 1.0e-3;
            var3_t r2 = { -1.0, b, 0.0 };
            var_t v2_x = 1.0e0 * constants::KmPerSecToAuPerDay;
            var3_t v2 = { v2_x, 0.0, 0.0 };

            y[0] = r1.x, y[1] = r1.y, y[2] = r1.z;
            y[3] = r2.x, y[4] = r2.y, y[5] = r2.z;

            y[6] = v1.x, y[ 7] = v1.y, y[ 8] = v1.z;
            y[9] = v2.x, y[10] = v2.y, y[11] = v2.z;

            // Set the initial stepsize for the integrator (should be model dependent)
            dt0 = 1.0e-4;

            print(dir, filename, n_obj);

            FREE_HOST_VECTOR((void **)&y);
            FREE_HOST_VECTOR((void **)&p);
            FREE_HOST_VECTOR((void **)&md);
        }

	} /* namespace nbody */
} /* namespace model */


int parse_options(int argc, const char **argv, string &odir, string &filename, uint32_t& n_obj)
{
	int i = 1;

	while (i < argc)
	{
		string p = argv[i];

		if (     p == "-odir")
		{
			i++;
			odir = argv[i];
		}
		else if (p == "-f")
		{
			i++;
			filename = argv[i];
		}
		else if (p == "-n_obj")
		{
			i++;
			n_obj = atoi(argv[i]);
		}
		else
		{
			throw string("Invalid switch on command-line.");
		}
		i++;
	}

	return 0;
}

int main(int argc, const char **argv)
{
	string odir;
	string filename;
	uint32_t n_obj;

#if 0
    /*
     * Test the distributions
     */
    try
    {
        uint32_t seed = (uint32_t)time(NULL);
        cout << "The seed number is " << seed << endl;
        string dir = "C:\\Work\\red.cuda.Results\\v2.0\\Test\\Distribution";
        string path;

        uint32_t n_data = 10000;
        var_t* data = NULL;
        ALLOCATE_HOST_VECTOR((void**)&data, n_data * sizeof(var_t));

        //{
        //    path = file::combine_path(dir, "uniform.txt");

        //    uniform_distribution dist(seed);
        //    for (uint32_t i = 0; i < n_data; i++)
        //    {
        //        data[i] = dist.get_next();
        //    }
        //    print_data(path, n_data, data);
        //}

        //{
        //    path = file::combine_path(dir, "exponential.txt");

        //    exponential_distribution dist(seed, 1.0);
        //    for (uint32_t i = 0; i < n_data; i++)
        //    {
        //        data[i] = dist.get_next();
        //    }
        //    print_data(path, n_data, data);
        //}

        //{
        //    path = file::combine_path(dir, "rayleigh.txt");

        //    rayleigh_distribution dist(seed, 1.0);
        //    for (uint32_t i = 0; i < n_data; i++)
        //    {
        //        data[i] = dist.get_next();
        //    }
        //    print_data(path, n_data, data);
        //}

        //{
        //    path = file::combine_path(dir, "normal.txt");

        //    normal_distribution dist(seed, 0, sqrt(2.0));
        //    for (uint32_t i = 0; i < n_data; i++)
        //    {
        //        data[i] = dist.get_next();
        //    }
        //    print_data(path, n_data, data);
        //}

        {
        //    path = file::combine_path(dir, "power_law_-0.99.txt");

        //    power_law_distribution dist(seed, -0.99, 0.4, 5.0);
        //    for (uint32_t i = 0; i < n_data; i++)
        //    {
        //        data[i] = dist.get_next();
        //    }
        //    print_data(path, n_data, data);
        //}

            path = file::combine_path(dir, "log_normal.txt");

            lognormal_distribution dist(seed, 1.0e-4, 0.25, 0, 10);
            for (uint32_t i = 0; i < n_data; i++)
            {
                data[i] = dist.get_next();
            }
            print_data(path, n_data, data);
        }

        FREE_HOST_VECTOR((void **)&data);

        return (EXIT_SUCCESS);
    }
    catch (const string& msg)
    {
        cerr << "Error: " << msg << endl;
        return (EXIT_FAILURE);
    }
#endif

	n_obj = 0;
	try
	{
        if (2 > argc)
        {
            throw string("Missing command line arguments.");
        }
		parse_options(argc, argv, odir, filename, n_obj);

		fn_info = filename + ".info.txt";
		fn_data = filename + ".data.txt";
        //model::tbp1D::create(odir, filename);
        //model::tbp2D::create(odir, filename);
        //model::rtbp1D::create(odir, filename);
        //model::rtbp2D::create(odir, filename);
		//model::nbody::create(odir, filename, n_obj);
		//model::nbody::create(odir, filename);               // The two-body problem
        model::nbody::create_grav_focusing(odir, filename);
        // model::nbody::ceate_disk(odir, filename, n_obj);
	}
	catch (const string& msg)
	{
		cerr << "Error: " << msg << endl;
		return (EXIT_FAILURE);
	}

	return (EXIT_SUCCESS);
}
