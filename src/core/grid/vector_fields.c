#include "core/grid/vector_fields.h"

void free_vector_fields (VectorFields *fields) {
    free(fields->velocity);
    free(fields->acceleration);
}