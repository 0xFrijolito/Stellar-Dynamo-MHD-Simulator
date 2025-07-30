#include "core/grid/scalar_fields.h"

void free_scalar_fields (ScalarFields *fields) {
    free(fields->density);
    free(fields->pressure);
    free(fields->temperature);
}
