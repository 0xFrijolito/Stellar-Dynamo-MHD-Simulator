#ifndef VECTOR_H
#define VECTOR_H

#include "core/defines.h"

typedef struct Vector3 {
    f64 x1;
    f64 x2;
    f64 x3;
} Vector3 ;

void vector3_rectangular_init (Vector3 *vector);
void vector3_cylindrical_init (Vector3 *vector);
void vector3_spherical_init   (Vector3 *vector);

void vector3_random (Vector3 *vector);
void vector3_zeros (Vector3 *vector);
void vector3_fill (Vector3 *vector);

void vector3_add (Vector3 v1, Vector3 v2, Vector3 *output);
void vector3_sub (Vector3 v1, Vector3 v2, Vector3 *output);
void vector3_scalar_mul (Vector3 vector, f64 value, Vector3 *output);

void vector3_dot_product   (Vector3 v1, Vector3 v2, f32 output);
void vector3_cross_product (Vector3 v1, Vector3 v2, Vector3 output);

void vector3_deep_copy (Vector3 *v1, Vector3 *v2);

#endif