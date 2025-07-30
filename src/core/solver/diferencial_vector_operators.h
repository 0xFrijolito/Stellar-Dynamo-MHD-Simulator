#ifndef DIFERENCIAL_VECTOR_OPERATORS_H
#define DIFERENCIAL_VECTOR_OPERATORS_H

#include "core/defines.h"
#include "core/math/vector.h"

void vector_divergence (Vector3 *vector_field, f64* output,
                        u32 cells_per_dimension, 
                        f64 dx1, f64 dx2, f64 dx3);

void curl (Vector3 *vector_field, Vector3 *output, 
           f64 dx1, f64 dx2, f64 dx3);

#endif