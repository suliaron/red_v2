// includes, system 
#include <algorithm>
#include <cassert>
#include <cmath>
#include <string>>

// includes, project
#include "distribution.h"

using namespace std;

mt19937::mt19937(uint32_t seed) :
    idx(0)
{
    I[0] = seed & 0xffffffff;
    for (int i = 1; i < 624; ++i)
    {
        I[i] = (1812433253 * (I[i - 1] ^ I[i - 1] >> 30) + i) & 0xffffffff;
    }
}

mt19937::~mt19937()
{ }

uint32_t mt19937::get_next()
{
    uint32_t j = (idx < 623) ? (idx + 1) : 0;
    uint32_t y = I[idx] & 0x80000000 | I[j] & 0x7fffffff;
    y = I[idx] = I[idx < 227 ? idx + 397 : idx - 227] ^ y >> 1 ^ (y & 1) * 0x9908b0df;
    idx = j;
    return y ^ (y ^= (y ^= (y ^= y >> 11) << 7 & 0x9d2c5680) << 15 & 0xefc60000) >> 18;
}


distribution_base::distribution_base(uint32_t seed, var_t x_min, var_t x_max) :
    prng(seed),
    x_min(x_min),
    x_max(x_max)
{
    if (x_min > x_max)
    {
        throw std::string("In distribution_base constructor the parameter 'x_min' is smaller than 'x_max'.");
    }
}

distribution_base::~distribution_base()
{ }


uniform_distribution::uniform_distribution(uint32_t seed) :
    distribution_base(seed, 0.0, 1.0)
{ }

uniform_distribution::uniform_distribution(uint32_t seed, var_t x_min, var_t x_max) :
    distribution_base(seed, x_min, x_max)
{ }

uniform_distribution::~uniform_distribution()
{ }

var_t uniform_distribution::get_next()
{
    return (x_min + (x_max - x_min)*(prng.get_next() / 4294967295.0));
}


exponential_distribution::exponential_distribution(uint32_t seed, var_t lambda) :
    distribution_base(seed, 0.0, 1.0),
    lambda(lambda)
{
    if (0.0 > lambda)
    {
        throw std::string("In exponential distribution the parameter 'lambda' is negative.");
    }
}

exponential_distribution::~exponential_distribution()
{ }

var_t exponential_distribution::get_next()
{
    var_t u = (prng.get_next() / 4294967295.0);
    return (-1.0 / lambda * log(u));
}


rayleigh_distribution::rayleigh_distribution(uint32_t seed, var_t sigma) :
    distribution_base(seed, 0.0, 1.0),
    sigma(sigma)
{
    if (0.0 > sigma)
    {
        throw std::string("In Rayleigh distribution the parameter 'sigma' is negative.");
    }
}

rayleigh_distribution::~rayleigh_distribution()
{ }

var_t rayleigh_distribution::get_next()
{
    var_t u = (prng.get_next() / 4294967295.0);
    return (sigma * sqrt(-2.0 * log(u)));
}


normal_distribution::normal_distribution(uint32_t seed, var_t mean, var_t variance) :
    distribution_base(seed, 0.0, 1.0),
    mean(mean),
    variance(variance)
{
    if (0.0 > variance)
    {
        throw std::string("In normal distribution the parameter 'variance' is negative.");
    }
}

normal_distribution::~normal_distribution()
{ }

var_t normal_distribution::get_next()
{
    var_t u = (prng.get_next() / 4294967295.0);
    return (mean + variance * sqrt(-2.0 * log(u)) * cos(1.4629180792671596E-9 * prng.get_next() ));
}


power_law_distribution::power_law_distribution(uint32_t seed, var_t x_min, var_t x_max, var_t power) :
    distribution_base(seed, x_min, x_max),
    power(power)
{
    if (0.0 == x_min)
    {
        throw std::string("In power-law distribution the parameter 'x_min' is zero.");
    }
    if (0.0 == x_max)
    {
        throw std::string("In power-law distribution the parameter 'x_max' is zero.");
    }
}

power_law_distribution::~power_law_distribution()
{ }

var_t power_law_distribution::get_next()
{
    var_t y_min = power == 0.0 ? 0.0 : pow(x_min, power);
    var_t y_max = power == 0.0 ? 1.0 : pow(x_max, power);
    if (y_min > y_max)
    {
        swap(y_min, y_max);
    }

    var_t d_y = y_max - y_min;
    var_t d_x = x_max - x_min;

    var_t x, y;
    var_t area_max = d_x * d_y;
    do
    {
        var_t ux = area_max * (prng.get_next() / 4294967295.0);
        var_t uy = y_min + (y_max - y_min) * (prng.get_next() / 4294967295.0);
        x = ux / d_y + x_min;
        y = uy;
    } while (y > pow(x, power));

    return x;
}


lognormal_distribution::lognormal_distribution(uint32_t seed, var_t x_min, var_t x_max, var_t mu, var_t sigma) :
    distribution_base(seed, x_min, x_max),
    mu(mu),
    sigma(sigma)
{
    if (0.0 >= x_min)
    {
        throw std::string("In log-normal distribution the parameter 'x_min' is not positive.");
    }
    if (0.0 >= sigma)
    {
        throw std::string("In log-normal distribution the parameter 'sigma' is not positive.");
    }
}

lognormal_distribution::~lognormal_distribution()
{ }

var_t lognormal_distribution::get_next()
{
    var_t y_min = 0.0;
    var_t max_loc = exp(mu - SQR(sigma));
    var_t y_max = pdf(max_loc);
    if (y_min > y_max)
    {
        swap(y_min, y_max);
    }

    var_t d_y = y_max - y_min;
    var_t d_x = x_max - x_min;

    var_t x, y;
    var_t area_max = d_x * d_y;
    do
    {
        var_t ux = area_max * (prng.get_next() / 4294967295.0);
        var_t uy = y_min + (y_max - y_min) * (prng.get_next() / 4294967295.0);
        x = ux / d_y + x_min;
        y = uy;
    } while (y > pdf(x));

    return x;
}

var_t lognormal_distribution::pdf(var_t x)
{
    return (1.0 / (sqrt(2.0*PI)*sigma * x) * exp(-SQR(log(x) - mu) / (2.0 * SQR(sigma))));
}
