#include "core/solver/diferencial_vector_operators.h"

static u32 index_3D_to_1D (u32 i, u32 j, u32 k, u32 size) {
    return i * size * size + 
           j * size + k;
}

void vector_divergence (Vector3 *vector_field, f64* output,
                        u32 cells_per_dimension, f64 dx, f64 dy, f64 dz) {

    for (u32 k = 1; k < cells_per_dimension - 1; ++k) {
        for (u32 j = 1; j < cells_per_dimension - 1; ++j) {
            for (u32 i = 1; i < cells_per_dimension - 1; ++i) {

                u64 id = index_3D_to_1D(i, j, k, cells_per_dimension);
                u64 ip = index_3D_to_1D(i+1, j, k, cells_per_dimension);
                u64 im = index_3D_to_1D(i-1, j, k, cells_per_dimension);
                u64 jp = index_3D_to_1D(i, j+1, k, cells_per_dimension);
                u64 jm = index_3D_to_1D(i, j-1, k, cells_per_dimension);
                u64 kp = index_3D_to_1D(i, j, k+1, cells_per_dimension);
                u64 km = index_3D_to_1D(i, j, k-1, cells_per_dimension);

                // Aproximación en caras
                f64 fx_r = 0.5 * (vector_field[id].x1 + vector_field[ip].x1);
                f64 fx_l = 0.5 * (vector_field[im].x1 + vector_field[id].x1);
                
                f64 fy_u = 0.5 * (vector_field[id].x2 + vector_field[jp].x2);
                f64 fy_d = 0.5 * (vector_field[jm].x2 + vector_field[id].x2);
                
                f64 fz_f = 0.5 * (vector_field[id].x3 + vector_field[kp].x3);
                f64 fz_b = 0.5 * (vector_field[km].x3 + vector_field[id].x3);

                f64 div = (fx_r - fx_l) / dx
                        + (fy_u - fy_d) / dy
                        + (fz_f - fz_b) / dz;

                output[id] = div;
            }
        }
    }
}
