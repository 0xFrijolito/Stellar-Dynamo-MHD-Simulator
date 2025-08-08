#ifndef TENSOR_H
#define TENSOR_H

#include "core/defines.h"

typedef struct {
    f64 *data;       // Puntero a los datos del tensor
    u32 *dimensions; // Dimensiones del tensor
    u32 rank;        // Rango del tensor (número de dimensiones)
} Tensor;

void init_tensor (Tensor *tensor, u32 rank, u32 *dimensions);
void free_tensor (Tensor *tensor);

Tensor tensor_zeros    (u32 rank, u32 *dimensions);
Tensor tensor_fill     (u32 rank, u32 *dimensions, f64 value);
Tensor tensor_identity (u32 size);
Tensor tensor_random   (u32 rank, u32 *dimensions, f64 min_value, f64 max_value);

// Tensor operations
Tensor tensor_add      (Tensor *t1, Tensor *t2);
Tensor tensor_subtract (Tensor *t1, Tensor *t2);

Tensor tensor_multiply (Tensor *t1, Tensor *t2);
Tensor tensor_divide   (Tensor *t1, Tensor *t2);

Tensor tensor_dot      (Tensor *t1, Tensor *t2);
Tensor tensor_cross    (Tensor *t1, Tensor *t2);

Tensor tensor_outer_product(Tensor *t1, Tensor *t2);

#endif