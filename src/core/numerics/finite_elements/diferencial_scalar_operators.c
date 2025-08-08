#include "core/numerics/finite_elements/diferencial_scalar_operators.h"

static u32 index_3D_to_1D(u32 i, u32 j, u32 k, u32 size) {
    return i * size * size 
         + j * size + k;
}

void scalar_gradient(f64 *scalar_field, Vector3 *metrics, u32 cells_per_dimension,
                     Vector3* output, f64 dx1, f64 dx2, f64 dx3) {

    for (u32 i=1 ; i<cells_per_dimension-1 ; i++) {
        for (u32 j=1 ; j<cells_per_dimension-1 ; j++) {
            for (u32 k=1 ; k<cells_per_dimension-1 ; k++) {
                u64 id = index_3D_to_1D(i, j, k, cells_per_dimension);

                u64 ip = index_3D_to_1D(i+1, j, k, cells_per_dimension);
                u64 im = index_3D_to_1D(i-1, j, k, cells_per_dimension);
                
                u64 jp = index_3D_to_1D(i, j+1, k, cells_per_dimension);
                u64 jm = index_3D_to_1D(i, j-1, k, cells_per_dimension);
                
                u64 kp = index_3D_to_1D(i, j, k+1, cells_per_dimension);
                u64 km = index_3D_to_1D(i, j, k-1, cells_per_dimension);

                f64 dfdx = (scalar_field[ip] - scalar_field[im]) / (2.0 * dx1);
                f64 dfdy = (scalar_field[jp] - scalar_field[jm]) / (2.0 * dx2);
                f64 dfdz = (scalar_field[kp] - scalar_field[km]) / (2.0 * dx3);

                output[id].x1 = dfdx;
                output[id].x2 = dfdy;
                output[id].x3 = dfdz;
            }
        }
    }
}

void scalar_laplacian(f64 *scalar_field, Vector3 *metrics, u32 cells_per_dimension,
                      f64* output, f64 dx1, f64 dx2, f64 dx3) {

    for (u32 i=1 ; i<cells_per_dimension-1 ; i++) {
        for (u32 j=1 ; j<cells_per_dimension-1 ; j++) {
            for (u32 k=1 ; k<cells_per_dimension-1 ; k++) {
                u64 id = index_3D_to_1D(i, j, k, cells_per_dimension);

                u64 ip = index_3D_to_1D(i+1, j, k, cells_per_dimension);
                u64 im = index_3D_to_1D(i-1, j, k, cells_per_dimension);
                
                u64 jp = index_3D_to_1D(i, j+1, k, cells_per_dimension);
                u64 jm = index_3D_to_1D(i, j-1, k, cells_per_dimension);
                
                u64 kp = index_3D_to_1D(i, j, k+1, cells_per_dimension);
                u64 km = index_3D_to_1D(i, j, k-1, cells_per_dimension);

                f64 dfdx = (scalar_field[ip] - 2*scalar_field[id] + scalar_field[im]) / (dx1*dx1);
                f64 dfdy = (scalar_field[jp] - 2*scalar_field[id] + scalar_field[jm]) / (dx2*dx2);
                f64 dfdz = (scalar_field[kp] - 2*scalar_field[id] + scalar_field[km]) / (dx3*dx3);

                output[id] = dfdx + dfdy + dfdz;
            }
        }
    }
}