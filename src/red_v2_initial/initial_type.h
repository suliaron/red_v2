#pragma once
#include <string>

#include "distribution.h"
#include "type.h"

typedef enum orbelem_name
{
    ORBELEM_NAME_SMA,
    ORBELEM_NAME_ECC,
    ORBELEM_NAME_INC,
    ORBELEM_NAME_PERI,
    ORBELEM_NAME_NODE,
    ORBELEM_NAME_MEAN,
    ORBELEM_NAME_N
} orbelem_name_t;

typedef enum pp_name
{
    PP_NAME_MASS,
    PP_NAME_RADIUS,
    PP_NAME_DENSITY,
    PP_NAME_DRAG_COEFF,
    PP_NAME_N
} pp_name_t;

typedef struct oe_dist
{
    var2_t				range[ORBELEM_NAME_N];
    distribution_base*  item[ORBELEM_NAME_N];
} oe_dist_t;

typedef struct pp_dist
{
    var2_t				range[PP_NAME_N];
    distribution_base*  item[PP_NAME_N];
} pp_dist_t;
