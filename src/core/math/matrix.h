#ifndef MATRIX_H
#define MATRIX_H

#include "core/defines.h"

typedef struct Matrix {
    u32 cols;
    u32 rows;

    f64 **values;
} Matrix;



#endif