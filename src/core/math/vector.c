#include "core/math/vector.h"

Vector3 vector3_zeros () {
    Vector3 output_vector = {
        0.0f,
        0.0f,
        0.0f
    };

    return output_vector;
}

Vector3 vector3_fill (f64 value) {
    Vector3 output_vector = {
        value,
        value,
        value
    };

    return output_vector;
}

Vector3 vector3_add (Vector3 v1, Vector3 v2) {
    Vector3 output_vector;
    
    output_vector.x1 = v1.x1 + v2.x1;
    output_vector.x2 = v1.x2 + v2.x2;
    output_vector.x3 = v1.x3 + v2.x3;

    return output_vector;
}

Vector3 vector3_sub (Vector3 v1, Vector3 v2) {
    Vector3 output_vector;
    
    output_vector.x1 = v1.x1 - v2.x1;
    output_vector.x2 = v1.x2 - v2.x2;
    output_vector.x3 = v1.x3 - v2.x3;

    return output_vector;
}

Vector3 vector3_scalar_mul (Vector3 vector, f64 value) {
    Vector3 output_vector;

    output_vector.x1 = vector.x1 * value;
    output_vector.x2 = vector.x2 * value;
    output_vector.x3 = vector.x3 * value;

    return output_vector;
}

// void vector3_dot_product   (Vector3 v1, Vector3 v2, f32 *output);
// void vector3_cross_product (Vector3 v1, Vector3 v2, Vector3 *output);

void vector3_deep_copy (Vector3 *v1, Vector3 *v2) {
    v1->x1 = v2->x1;
    v1->x2 = v2->x2;
    v1->x3 = v2->x3;
}