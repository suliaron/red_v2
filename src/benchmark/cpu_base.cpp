#include <macro.h>
#include <type.h>

#define NDIM   3        // Number of space dimension
#define NVPO   6        // Number of variables per object (3 space and 3 velocity coordinates)

static inline
void body_body_grav_accel(const var3_t& ri, const var3_t& rj, var_t mj, var3_t& ai)
{
    var3_t r_ij = { 0.0, 0.0, 0.0 };

    // compute r_ij = r_j - r_i [3 FLOPS] [6 read, 3 write]
    r_ij.x = rj.x - ri.x;
    r_ij.y = rj.y - ri.y;
    r_ij.z = rj.z - ri.z;

    // compute norm square of d vector [5 FLOPS] [3 read, 1 write]
    var_t d2 = SQR(r_ij.x) + SQR(r_ij.y) + SQR(r_ij.z);
    var_t d = sqrt(d2);
    var_t s = K2 * mj / (d * d2);

    ai.x += s * r_ij.x;
    ai.y += s * r_ij.y;
    ai.z += s * r_ij.z;
}

static inline
void body_body_grav_accel_sym(const var3_t& ri, const var3_t& rj, var_t mi, var_t mj, var3_t& ai, var3_t& aj)
{
    var3_t r_ij = { 0.0, 0.0, 0.0 };

    // compute r_ij = r_j - r_i [3 FLOPS] [6 read, 3 write]
    r_ij.x = rj.x - ri.x;
    r_ij.y = rj.y - ri.y;
    r_ij.z = rj.z - ri.z;

    // compute norm square of d vector [5 FLOPS] [3 read, 1 write]
    var_t d2 = SQR(r_ij.x) + SQR(r_ij.y) + SQR(r_ij.z);
    var_t d = sqrt(d2);
    var_t d_3 = K2 / (d * d2);

    var_t s = mj * d_3;
    ai.x += s * r_ij.x;
    ai.y += s * r_ij.y;
    ai.z += s * r_ij.z;

    s = mi * d_3;
    aj.x -= s * r_ij.x;
    aj.y -= s * r_ij.y;
    aj.z -= s * r_ij.z;
}

void cpu_calc_grav_accel(var_t t, uint32_t n_obj, const var3_t* r, const nbp_t::param_t* p, var3_t* a, bool use_sym)
{
    // Number of space and velocity coordinates
    const uint32_t nv = NDIM * n_obj;

    // Create aliases
    //const var3_t* r = (var3_t*)h_y;
    //const nbp_t::param_t* p = (nbp_t::param_t*)h_p;
    //var3_t* a = (var3_t*)(h_dy + nv);

    // Clear the acceleration array: the += op can be used
    memset(a, 0, nv * sizeof(var_t));

    if (use_sym)
    {
        for (uint32_t i = 0; i < n_obj; i++)
        {
            var3_t r_ij = { 0.0, 0.0, 0.0 };
            for (uint32_t j = i + 1; j < n_obj; j++)
            {
                body_body_grav_accel_sym(r[i], r[j], p[i].mass, p[j].mass, a[i], a[j]);
                //r_ij.x = r[j].x - r[i].x;
                //r_ij.y = r[j].y - r[i].y;
                //r_ij.z = r[j].z - r[i].z;

                //var_t d2 = SQR(r_ij.x) + SQR(r_ij.y) + SQR(r_ij.z);
                //var_t d = sqrt(d2);
                //var_t d_3 = K2 / (d * d2);

                //var_t s = p[j].mass * d_3;
                //a[i].x += s * r_ij.x;
                //a[i].y += s * r_ij.y;
                //a[i].z += s * r_ij.z;

                //s = p[i].mass * d_3;
                //a[j].x -= s * r_ij.x;
                //a[j].y -= s * r_ij.y;
                //a[j].z -= s * r_ij.z;
            }
        }
    }
    else
    {
        for (uint32_t i = 0; i < n_obj; i++)
        {
            var3_t r_ij = { 0.0, 0.0, 0.0 };
            for (uint32_t j = 0; j < n_obj; j++)
            {
                if (i == j) continue;
                body_body_grav_accel(r[i], r[j], p[j].mass, a[i]);
                //r_ij.x = r[j].x - r[i].x;
                //r_ij.y = r[j].y - r[i].y;
                //r_ij.z = r[j].z - r[i].z;

                //var_t d2 = SQR(r_ij.x) + SQR(r_ij.y) + SQR(r_ij.z);
                //var_t d = sqrt(d2);
                //var_t s = K2 * p[j].mass / (d * d2);

                //a[i].x += s * r_ij.x;
                //a[i].y += s * r_ij.y;
                //a[i].z += s * r_ij.z;
            }
        }
    }
}

void cpu_calc_grav_accel(var_t t, uint32_t n_obj, uint2_t snk, uint2_t src, const var3_t* r, const nbp_t::param_t* p, var3_t* a, bool use_sym)
{
    // Number of space and velocity coordinates
    const uint32_t nv = NDIM * n_obj;

    // Create aliases
    //const var3_t* r = (var3_t*)h_y;
    //const nbp_t::param_t* p = (nbp_t::param_t*)h_p;
    //var3_t* a = (var3_t*)(h_dy + nv);

    // Clear the acceleration array: the += op can be used
    memset(a, 0, nv * sizeof(var_t));

    if (use_sym)
    {
        for (uint32_t i = snk.n1; i < snk.n2; i++)
        {
            var3_t r_ij = { 0.0, 0.0, 0.0 };
            for (uint32_t j = i + 1; j < src.n2; j++)
            {
                r_ij.x = r[j].x - r[i].x;
                r_ij.y = r[j].y - r[i].y;
                r_ij.z = r[j].z - r[i].z;

                var_t d2 = SQR(r_ij.x) + SQR(r_ij.y) + SQR(r_ij.z);
                var_t d = sqrt(d2);
                var_t d_3 = K2 / (d * d2);

                var_t s = p[j].mass * d_3;
                a[i].x += s * r_ij.x;
                a[i].y += s * r_ij.y;
                a[i].z += s * r_ij.z;

                s = p[i].mass * d_3;
                a[j].x -= s * r_ij.x;
                a[j].y -= s * r_ij.y;
                a[j].z -= s * r_ij.z;
            }
        }
    }
    else
    {
        for (uint32_t i = snk.n1; i < snk.n2; i++)
        {
            var3_t r_ij = { 0.0, 0.0, 0.0 };
            for (uint32_t j = src.n1; j < src.n2; j++)
            {
                if (i == j) continue;
                r_ij.x = r[j].x - r[i].x;
                r_ij.y = r[j].y - r[i].y;
                r_ij.z = r[j].z - r[i].z;

                var_t d2 = SQR(r_ij.x) + SQR(r_ij.y) + SQR(r_ij.z);
                var_t d = sqrt(d2);
                var_t s = K2 * p[j].mass / (d * d2);

                a[i].x += s * r_ij.x;
                a[i].y += s * r_ij.y;
                a[i].z += s * r_ij.z;
            }
        }
    }
}

#undef NDIM
#undef NVPO
