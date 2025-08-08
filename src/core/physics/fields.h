#ifndef FIELDS_H
#define FIELDS_H

#include "core/defines.h"
#include "core/physics/inicial_conditions.h"

typedef struct EulerFields {
    f64 *density;
    f64 *momentum_x;
    f64 *momentum_y;
    f64 *momentum_z;
    f64 *energy; 
} EulerFields ;

#endif