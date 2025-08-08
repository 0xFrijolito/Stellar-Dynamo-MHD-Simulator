#include "core/numerics/finite_elements/diferencial_vector_operators.h"

static u32 index_3D_to_1D (u32 i, u32 j, u32 k, u32 size) {
    return i * size * size + 
           j * size + k;
}

void vector_divergence (Vector3 *vector_field, f64* output,
                        u32 cells_per_dimension, f64 dx1, f64 dx2, f64 dx3) {

    for (u32 i=1 ; i<cells_per_dimension-1 ; i++) {
        for (u32 j=1 ; j<cells_per_dimension-1 ; j++) {
            for (u32 k=1 ; k<cells_per_dimension-1 ; k++) {

                u64 id = index_3D_to_1D(i  , j  , k  , cells_per_dimension);
                u64 ip = index_3D_to_1D(i+1, j  , k  , cells_per_dimension);
                u64 im = index_3D_to_1D(i-1, j  , k  , cells_per_dimension);
                u64 jp = index_3D_to_1D(i  , j+1, k  , cells_per_dimension);
                u64 jm = index_3D_to_1D(i  , j-1, k  , cells_per_dimension);
                u64 kp = index_3D_to_1D(i  , j  , k+1, cells_per_dimension);
                u64 km = index_3D_to_1D(i  , j  , k-1, cells_per_dimension);

                f64 fx_r = 0.5 * (vector_field[id].x1 + vector_field[ip].x1);
                f64 fx_l = 0.5 * (vector_field[im].x1 + vector_field[id].x1);
                
                f64 fy_u = 0.5 * (vector_field[id].x2 + vector_field[jp].x2);
                f64 fy_d = 0.5 * (vector_field[jm].x2 + vector_field[id].x2);
                
                f64 fz_f = 0.5 * (vector_field[id].x3 + vector_field[kp].x3);
                f64 fz_b = 0.5 * (vector_field[km].x3 + vector_field[id].x3);

                f64 div = (fx_r - fx_l) / (2.0f * dx1)
                        + (fy_u - fy_d) / (2.0f * dx2)
                        + (fz_f - fz_b) / (2.0f * dx3);

                output[id] = div;
            }
        }
    }
}

void vector_divergence_periodic (Vector3 *vector_field, f64* output,
                                 u32 cells_per_dimension, f64 dx1, f64 dx2, f64 dx3) {

    for (u32 i=0 ; i<cells_per_dimension-1 ; i++) {
        u32 ip = (i + 1 == cells_per_dimension) ? 0                       : i + 1;
        u32 im = (i == 0)                       ? cells_per_dimension - 1 : i - 1;
        
        for (u32 j=0 ; j<cells_per_dimension-1 ; j++) {
            u32 jp = (j + 1 == cells_per_dimension) ? 0                       : j + 1;
            u32 jm = (j == 0)                       ? cells_per_dimension - 1 : j - 1;

            for (u32 k=0 ; k<cells_per_dimension-1 ; k++) {
                u32 kp = (k + 1 == cells_per_dimension) ? 0                       : k + 1;
                u32 km = (k == 0)                       ? cells_per_dimension - 1 : k - 1;

                u64 idx  = index_3D_to_1D(i , j , k , cells_per_dimension);
                u64 ipjk = index_3D_to_1D(ip, j , k , cells_per_dimension);
                u64 imjk = index_3D_to_1D(im, j , k , cells_per_dimension);
                u64 ijpk = index_3D_to_1D(i , jp, k , cells_per_dimension);
                u64 ijmk = index_3D_to_1D(i , jm, k , cells_per_dimension);
                u64 ijkp = index_3D_to_1D(i , j , kp, cells_per_dimension);
                u64 ijkm = index_3D_to_1D(i , j , km, cells_per_dimension);

                f64 dFx = (vector_field[ipjk].x1 - vector_field[imjk].x1) / (2.0 * dx1);
                f64 dFy = (vector_field[ijpk].x2 - vector_field[ijmk].x2) / (2.0 * dx2);
                f64 dFz = (vector_field[ijkp].x3 - vector_field[ijkm].x3) / (2.0 * dx3);

                output[idx] = dFx + dFy + dFz;
            }
        }
    }
}
