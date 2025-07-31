#include "core/math/vector.h"

void vector3_zeros (Vector3 *vector) {
    vector->x1 = 0.00f;
    vector->x2 = 0.00f;
    vector->x3 = 0.00f;
}

void vector3_fill (Vector3 *vector, f64 value) {
    vector->x1 = value;
    vector->x2 = value;
    vector->x3 = value;
}

void vector3_add (Vector3 v1, Vector3 v2, Vector3 *output) {
    output->x1 = v1.x1 + v2.x1;
    output->x2 = v1.x2 + v2.x2;
    output->x3 = v1.x3 + v2.x3;
}

void vector3_sub (Vector3 v1, Vector3 v2, Vector3 *output) {
    output->x1 = v1.x1 - v2.x1;
    output->x2 = v1.x2 - v2.x2;
    output->x3 = v1.x3 - v2.x3;
}

void vector3_scalar_mul (Vector3 vector, f64 value, Vector3 *output) {
    output->x1 = vector.x1 * value;
    output->x2 = vector.x2 * value;
    output->x3 = vector.x3 * value;
}

// void vector3_dot_product   (Vector3 v1, Vector3 v2, f32 *output);
// void vector3_cross_product (Vector3 v1, Vector3 v2, Vector3 *output);

void vector3_deep_copy (Vector3 *v1, Vector3 *v2) {
    v1->x1 = v2->x1;
    v1->x2 = v2->x2;
    v1->x3 = v2->x3;
}