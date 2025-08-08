#ifndef DIFERENCIAL_SCALAR_OPERATORS_H
#define DIFERENCIAL_SCALAR_OPERATORS_H

#include "core/defines.h"
#include "core/math/vector.h"

void scalar_gradient (f64 *scalar_field, Vector3 *metrics, u32 cells_per_dimension,
                      Vector3* output, f64 dx1, f64 dx2, f64 dx3);

void scalar_laplacian (f64 *scalar_field, Vector3 *metrics, u32 cells_per_dimension,
                       f64* output, f64 dx1, f64 dx2, f64 dx3);
       
#endif