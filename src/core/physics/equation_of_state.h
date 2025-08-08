#ifndef EQUATION_OF_STATE_H
#define EQUATION_OF_STATE_H

#include "core/defines.h"

f64 ideal_gas_pressure (f64 rho, f64 e, f64 gamma) {
    return (gamma - 1.0f) * rho * e;
}

#endif