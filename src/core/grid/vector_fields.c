#include "core/grid/vector_fields.h"

void free_vector_fields (VectorFields *fields) {
    free(fields->velocity);
    free(fields->acceleration);
}

Vector3 none_vectorial_function (f64 x, f64 y, f64 z) {
    Vector3 vector = {0.0f, 0.0f, 0.0f};
    return vector;
}