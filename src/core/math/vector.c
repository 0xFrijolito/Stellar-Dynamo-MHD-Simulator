#include "core/math/vector.h"

// void vector3_random (Vector3 *vector) {
//     vector->x = random_f32();
//     vector->y = random_f32();
//     vector->z = random_f32();
// }

// void vector3_zeros (Vector3 *vector) {
//     vector->x = 0.00f;
//     vector->y = 0.00f;
//     vector->z = 0.00f;
// }

// void vector3_fill (Vector3 *vector, f32 value) {
//     vector->x = value;
//     vector->y = value;
//     vector->z = value;
// }

// void vector3_add (Vector3 v1, Vector3 v2, Vector3 *output);
// void vector3_sub (Vector3 v1, Vector3 v2, Vector3 *output);
// void vector3_scalar_mul (Vector3 vector, float value, Vector3 *output);

// void vector3_dot_product   (Vector3 v1, Vector3 v2, f32 *output);
// void vector3_cross_product (Vector3 v1, Vector3 v2, Vector3 *output);

void vector3_deep_copy (Vector3 *v1, Vector3 *v2) {
    v1->x1 = v2->x1;
    v1->x1 = v2->x1;
    v1->x1 = v2->x1;
}