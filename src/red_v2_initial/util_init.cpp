// includes, system 
#include <iomanip>      // setw()
#include <cmath>
#include <fstream>		// ofstream
#include <string>		// string
#include <sstream>		// ostringstream
#include <stdlib.h>		// atof()

// includes, project
#include "util_init.h"
#include "distribution.h"
#include "constants.h"
#include "redutil2.h"

using namespace std;
using namespace redutil2;

var_t generate_oe(orbelem_name_t name, oe_dist_t *oe_d)
{
    var_t result = 0.0;

    var_t x_min = oe_d->range[name].x;
    var_t x_max = oe_d->range[name].y;
    var_t dx = fabs(x_max - x_min);
    do
    {
        result = oe_d->item[name]->get_next();
    } while (dx > 0.0 && x_min > result && x_max < result);

    return result;
}

void generate_oe(oe_dist_t *oe_d, orbelem_t& oe)
{
    if (0x0 != oe_d->item[ORBELEM_NAME_SMA])
    {
        oe.sma = generate_oe(ORBELEM_NAME_SMA, oe_d);
    }
    else
    {
        oe.sma = 0.0;
    }
    if (0x0 != oe_d->item[ORBELEM_NAME_ECC])
    {
        oe.ecc = generate_oe(ORBELEM_NAME_ECC, oe_d);
    }
    else
    {
        oe.ecc = 0.0;
    }
    if (0x0 != oe_d->item[ORBELEM_NAME_INC])
    {
        oe.inc = generate_oe(ORBELEM_NAME_INC, oe_d);
    }
    else
    {
        oe.inc = 0.0;
    }
    if (0x0 != oe_d->item[ORBELEM_NAME_PERI])
    {
        oe.peri = generate_oe(ORBELEM_NAME_PERI, oe_d);
    }
    else
    {
        oe.peri = 0.0;
    }
    if (0x0 != oe_d->item[ORBELEM_NAME_NODE])
    {
        oe.node = generate_oe(ORBELEM_NAME_NODE, oe_d);
    }
    else
    {
        oe.node = 0.0;
    }
    if (0x0 != oe_d->item[ORBELEM_NAME_MEAN])
    {
        oe.mean = generate_oe(ORBELEM_NAME_MEAN, oe_d);
    }
    else
    {
        oe.mean = 0.0;
    }
}

var_t generate_pp(pp_name_t name, pp_dist_t* pp_d)
{
    var_t result = 0.0;

    var_t x_min = pp_d->range[name].x;
    var_t x_max = pp_d->range[name].y;
    var_t dx = fabs(x_max - x_min);
    do
    {
        result = pp_d->item[name]->get_next();
    } while (dx > 0.0 && x_min > result && x_max < result);

    return result;
}

void generate_pp(pp_dist_t *pp_d, nbp_t::param_t& param)
{
    param.mass = param.radius = param.density = 0.0;

    if (0x0 != pp_d->item[PP_NAME_MASS])
    {
        param.mass = generate_pp(PP_NAME_MASS, pp_d);
    }
    if (0x0 != pp_d->item[PP_NAME_RADIUS])
    {
        param.radius = generate_pp(PP_NAME_RADIUS, pp_d);
    }
    if (0x0 != pp_d->item[PP_NAME_DENSITY])
    {
        param.density = generate_pp(PP_NAME_DENSITY, pp_d);
    }
}

void generate_pp(pp_dist_t *pp_d, pp_disk_t::param_t& param)
{
    param.mass = param.radius = param.density = param.cd = 0.0;

    if (0x0 != pp_d->item[PP_NAME_MASS])
    {
        param.mass = generate_pp(PP_NAME_MASS, pp_d);
    }
    if (0x0 != pp_d->item[PP_NAME_RADIUS])
    {
        param.radius = generate_pp(PP_NAME_RADIUS, pp_d);
    }
    if (0x0 != pp_d->item[PP_NAME_DENSITY])
    {
        param.density = generate_pp(PP_NAME_DENSITY, pp_d);
    }
    if (0x0 != pp_d->item[PP_NAME_DRAG_COEFF])
    {
        param.cd = generate_pp(PP_NAME_DRAG_COEFF, pp_d);
    }
}

template <typename T>
void print_number(string& path, T number)
{
    printf("Writing %s to disk .", path.c_str());

    ofstream sout(path.c_str(), ios_base::out);
    if (sout)
    {
        sout << number;
        sout.close();
        printf(" done\n");
    }
    else
    {
        throw string("Cannot open " + path + ".");
    }
}
template void print_number<char>(string& path, char number);
template void print_number<unsigned char>(string& path, unsigned char number);
template void print_number<int>(string& path, int number);
template void print_number<uint32_t>(string& path, uint32_t number);
template void print_number<long>(string& path, long number);
template void print_number<unsigned long>(string& path, unsigned long number);


void print_oe(std::string &path, uint32_t n, var_t t, nbp_t::metadata_t* md, nbp_t::param_t* p, orbelem_t* oe)
{
    printf("Writing %s to disk .", path.c_str());

    ofstream sout(path.c_str(), ios_base::out);
    if (sout)
    {
        int pcd = 1;
        for (uint32_t i = 0; i < n; i++)
        {
            file::nbp::print_oe_record(sout, t, oe[i], p[i], md[i]);
            if (pcd <= (int)((((var_t)(i + 1) / (var_t)n))*100.0))
            {
                putchar('.');
                pcd++;
            }
        }
        sout.close();
        printf(" done\n");
    }
    else
    {
        throw string("Cannot open " + path + ".");
    }
}

/*
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
JDCT ,   ,                                          EC,                     QR,                     IN,                     OM,                     W,                      Tp,                     N,                      MA,                     TA,                     A,                      AD,                     PR
2457153.500000000, A.D. 2015-May-11 00:00:00.0000,  2.056283170656704E-01,  3.075005270000471E-01,  7.004033081870772E+00,  4.831135275416084E+01,  2.917018579009503E+01,  2.457132325703553E+06,  4.092332803985572E+00,  8.665226794919347E+01,  1.098801777958040E+02,  3.870990540148295E-01,  4.666975810296119E-01,  8.796938500441402E+01,
*/
var_t extract_from_horizon_output(string &data, orbelem_t& oe)
{
    int start = 0;
    int len = start + 17;
    string s = data.substr(start, len);

    var_t epoch = atof(s.c_str());

    start = 52;
    len = 21;
    s = data.substr(start, len);
    oe.ecc = atof(s.c_str());

    start = 100;
    len = 21;
    s = data.substr(start, len);
    oe.inc = atof(s.c_str()) * constants::DegreeToRadian;

    start = 124;
    len = 21;
    s = data.substr(start, len);
    oe.node = atof(s.c_str()) * constants::DegreeToRadian;

    start = 148;
    len = 21;
    s = data.substr(start, len);
    oe.peri = atof(s.c_str()) * constants::DegreeToRadian;

    start = 220;
    len = 21;
    s = data.substr(start, len);
    oe.mean = atof(s.c_str()) * constants::DegreeToRadian;

    start = 268;
    len = 21;
    s = data.substr(start, len);
    oe.sma = atof(s.c_str());

    return epoch;
}
