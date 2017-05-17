#include <stdio.h>
#include <iomanip>
#include <iostream>
#include <fstream>

#include "macro.h"
#include "type.h"
#include "redutil2.h"

using namespace std;
using namespace redutil2;

#if 0
int main()
{
    printf("sizeof(unsigned short int) = %d\n", sizeof(unsigned short int));
    printf("sizeof(usi_t) = %d\n", sizeof(usi_t));

    float f = 1.23456789e10;

    printf("f = '%*e'\n", FLT_T_W, f);
    cout << "'" << setw(FLT_T_W) << f << "'" << endl;

    printf("sizeof(nbp_t::metadata_t) = %3d\n", sizeof(nbp_t::metadata_t));
    printf("sizeof(nbp_t::param_t)    = %3d\n", sizeof(nbp_t::param_t));


    try
    {
        nbp_t::metadata_t *h_md;
        var_t* h_y = NULL;
        var_t* h_p = NULL;

        uint32_t n_obj = 11;
        uint32_t n_var = 6 * n_obj;
        uint32_t n_par = 3 * n_obj;
        ALLOCATE_HOST_VECTOR((void**)&(h_md), n_obj * sizeof(nbp_t::metadata_t));
        ALLOCATE_HOST_VECTOR((void**)&(h_y), n_var * sizeof(var_t));
        ALLOCATE_HOST_VECTOR((void**)&(h_p), n_par * sizeof(var_t));

        ifstream input;

        string path = "C:\\Work\\red.cuda.Results\\v2.0\\Test\\nbp\\02\\input.data.txt";
        input.open(path.c_str(), ios::in);
        if (input)
        {
            nbp_t::param_t* _p = (nbp_t::param_t*)h_p;
            for (uint32_t i = 0; i < n_obj; i++)
            {
                input >> h_md[i].id >> h_md[i].body_type >> h_md[i].active >> h_md[i].mig_type >> h_md[i].mig_stop_at >> h_md[i].unused1 >> h_md[i].unused2 >> h_md[i].unused3;
                input >> _p[i].density >> _p[i].mass >> _p[i].radius;
                uint32_t offset = 3 * i;
                // position
                input >> h_y[offset + 0] >> h_y[offset + 1] >> h_y[offset + 2];
                offset += 3 * n_obj;
                // velocity
                input >> h_y[offset + 0] >> h_y[offset + 1] >> h_y[offset + 2];
            }
        }
        else
        {
            printf("Cannot open %s.\n", path.c_str());
        }
        input.close();

        nbp_t::metadata_t *md = h_md;

        uint32_t offset = 3 * n_obj;
        nbp_t::param_t* p = (nbp_t::param_t*)h_p;
        var_t* y = h_y;

        cout.setf(ios::right);
        cout.setf(ios::scientific);

        for (uint32_t i = 0; i < n_obj; i++)
        {
            cout << setw(INT_T_W) << md[i].id << SEP
                << setw(2) << (int)md[i].body_type << SEP
                << setw(2) << md[i].active << SEP
                << setw(2) << (int)md[i].mig_type << SEP
                << setprecision(6) << setw(FLT_T_W) << md[i].mig_stop_at << SEP
                << setw(2) << md[i].unused1 << SEP
                << setw(2) << md[i].unused2 << SEP
                << setw(2) << md[i].unused3 << SEP;

            // Print the parameters for each object
            cout << setprecision(16) << setw(VAR_T_W) << p[i].density << SEP
                << setw(VAR_T_W) << p[i].mass << SEP
                << setw(VAR_T_W) << p[i].radius << SEP;

            // Print the variables for each object
            uint32_t j = 3 * i;
            cout << setw(VAR_T_W) << y[j] << SEP            /* x-coordinate */
                << setw(VAR_T_W) << y[j + 1] << SEP            /* y-coordinate */
                << setw(VAR_T_W) << y[j + 2] << SEP            /* z-coordinate */
                << setw(VAR_T_W) << y[offset + j] << SEP   /* vx-velocity  */
                << setw(VAR_T_W) << y[offset + j + 1] << SEP   /* vy-velocity  */
                << setw(VAR_T_W) << y[offset + j + 2] << endl; /* vz-velocity  */
        }


        FREE_HOST_VECTOR((void **)&(h_y));
        FREE_HOST_VECTOR((void **)&(h_p));
        FREE_HOST_VECTOR((void **)&(h_md));

        cout << " done" << endl;
    }
    catch (const string& msg)
    {
        cerr << "ERROR: " << msg << endl;
    }


    return 0;
}
#endif

#if 1
#include <stdio.h>
#include <time.h>
#include <sys/time.h>

int main()
{
    time_t start = time(NULL);
    int same = 0;
    int different = 0;
    int max_usec = 0;
    while (1) {
    time_t t;
    struct timeval tv;
    gettimeofday(&tv, NULL);
    t = time(NULL);
    if (t < tv.tv_sec) {
    different++;
    if (tv.tv_usec > max_usec) {
    max_usec = tv.tv_usec;
    }
    }
    else {
    same++;
    }
    if (t > start + 5) {
    break;
    }
    }
    printf("Same:      %i\n", same);
    printf("Different: %i\n", different);
    printf("Largest difference seen at %i\n", max_usec);
}

#endif

