#ifndef BOUNDARY_CONDITIONS_H
#define BOUNDARY_CONDITIONS_H

#include "core/defines.h"
#include "core/math/vector.h"

void scalar_dirichlet_boundary_conditions (f64 *scalar_field, f64 boundary_value, u32 cells_per_dimension);
void scalar_neumann_bondary_condition     (f64 *scalar_field, f64 boundary_value, u32 cells_per_dimension);

void vector_dirichlet_boundary_condition (Vector3 *scalar_field, Vector3 *metrics, Vector3* output);
void vector_neumann_bondary_condition    (Vector3 *scalar_field, Vector3 *metrics, Vector3* output);

#endif