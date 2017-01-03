#pragma once

#include "macro.h"
#include "type.h"

/*
* The Mersenne Twister is a pseudorandom number generator (PRNG). It is by far the most widely used
* general-purpose PRNG. Its name derives from the fact that its period length is chosen to be a Mersenne prime.
*/
class mt19937
{
public:
    mt19937(uint32_t seed);
    ~mt19937();

    uint32_t get_next();
private:
    uint32_t idx;
    uint32_t I[624];
};


class distribution_base
{
public:
    distribution_base(uint32_t seed, var_t x_min, var_t x_max);
    ~distribution_base();

    var_t get_x_min() { return x_min; }
    var_t get_x_max() { return x_max; }

    virtual var_t get_next() = 0;

protected:
    mt19937 prng;
    var_t x_min;
    var_t x_max;
};

class uniform_distribution : public distribution_base
{
public:
    uniform_distribution(uint32_t seed);
    uniform_distribution(uint32_t seed, var_t x_min, var_t x_max);
    ~uniform_distribution();

    var_t get_next();
};

class exponential_distribution : public distribution_base
{
public:
    exponential_distribution(uint32_t seed, var_t lambda);
    ~exponential_distribution();

    var_t get_next();
private:
    var_t lambda;
};

class rayleigh_distribution : public distribution_base
{
public:
    rayleigh_distribution(uint32_t seed, var_t sigma);
    ~rayleigh_distribution();

    var_t get_next();
private:
    var_t sigma;
};

class normal_distribution : public distribution_base
{
public:
    normal_distribution(uint32_t seed, var_t mean, var_t variance);
    ~normal_distribution();

    var_t get_next();
private:
    var_t mean;
    var_t variance;
};

class power_law_distribution : public distribution_base
{
public:
    power_law_distribution(uint32_t seed, var_t x_min, var_t x_max, var_t power);
    ~power_law_distribution();

    var_t get_next();
private:
    var_t power;
};

class lognormal_distribution : public distribution_base
{
public:
    lognormal_distribution(uint32_t seed, var_t x_min, var_t x_max, var_t mu, var_t sigma);
    ~lognormal_distribution();

    var_t get_next();
private:
    var_t pdf(var_t);

    var_t mu;
    var_t sigma;
};
