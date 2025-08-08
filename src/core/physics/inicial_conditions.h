#ifndef INICIAL_CONDITIONS_H
#define INICIAL_CONDITIONS_H

#include "core/defines.h"

typedef f64     (*InicialScalarFunction)(f64 x, f64 y, f64 z);
typedef Vector3 (*InicialVectorFunction)(f64 x, f64 y, f64 z);

typedef struct EulerInicialConditions {
    InicialScalarFunction inicial_density;

    InicialScalarFunction inicial_momentum_x;
    InicialScalarFunction inicial_momentum_y;
    InicialScalarFunction inicial_momentum_z;
    
    InicialScalarFunction inicial_energy;

    InicialVectorFunction inicial_velocity;
} EulerInicialConditions ;

#endif