#ifndef SCALAR_FIELDS_H
#define SCALAR_FIELDS_H

#include <stdlib.h>

#include "core/math/vector.h"

typedef struct ScalarFieldsInitProfiles {
    f64 (*density_function)     (f64 x1, f64 x2, f64 x3);
    f64 (*pressure_function)    (f64 x1, f64 x2, f64 x3);
    f64 (*temperature_function) (f64 x1, f64 x2, f64 x3);
} ScalarFieldsInitProfiles;

typedef struct ScalarFields {
    f64 *density;
    f64 *pressure;
    f64 *temperature;
} ScalarFields ;

void init_scalar_fields ();
void free_scalar_fields ();

#endif