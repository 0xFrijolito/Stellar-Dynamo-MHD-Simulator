#ifndef VECTOR_FIELDS_H
#define VECTOR_FIELDS_H

#include <stdlib.h>

#include "core/math/vector.h"

typedef struct VectorFieldsInitProfiles {
    Vector3 (*velocity_function)     (f64 x1, f64 x2, f64 x3);
    Vector3 (*acceleration_function) (f64 x1, f64 x2, f64 x3);
} VectorFieldsInitProfiles ;

typedef struct VectorFields {
    Vector3 *velocity;
    Vector3 *acceleration;
} VectorFields ;

void init_vector_fields ();
void free_vector_fields ();

#endif // VECTOR_FIELDS_H