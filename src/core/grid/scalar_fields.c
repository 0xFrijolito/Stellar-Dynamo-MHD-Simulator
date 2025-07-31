#include "core/grid/scalar_fields.h"

void free_scalar_fields (ScalarFields *fields) {
    free(fields->density);
    free(fields->pressure);
    free(fields->temperature);
}

f64 none_scalar_function (f64 x1, f64 x2, f64 x3) { return 0.0f; }