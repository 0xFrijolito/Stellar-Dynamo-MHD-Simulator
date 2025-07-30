#include "core/solver/diferencial_scalar_operators.h"

static u32 index_3D_to_1D (u32 i, u32 j, u32 k, u32 size) {
    return i * size * size + 
           j * size + k;
}

void scalar_gradient (f64 *scalar_field, Vector3 *metrics, Vector3* output,
                      u32 cells_per_dimension,
                      f64 dx1, f64 dx2, f64 dx3) {

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

                // Cada metrics[id].xN ya almacena el tamaño físico real de la celda en esa dirección
                f64 h1 = metrics[id].x1;
                f64 h2 = metrics[id].x2;
                f64 h3 = metrics[id].x3;

                // Calculamos cada derivada
                f64 dfdx = (scalar_field[ip] - scalar_field[im]) / (2.0 * dx1 * h1);
                f64 dfdy = (scalar_field[jp] - scalar_field[jm]) / (2.0 * dx2 * h2);
                f64 dfdz = (scalar_field[kp] - scalar_field[km]) / (2.0 * dx3 * h3);

                output[id].x1 = dfdx;
                output[id].x2 = dfdy;
                output[id].x3 = dfdz;
            }
        }
    }
}

void scalar_laplacian (f64 *scalar_field, Vector3 *metrics, f64* output,
                       u32 cells_per_dimension,
                       f64 dx1, f64 dx2, f64 dx3) {

    for (u32 k = 1; k < cells_per_dimension-1; ++k) {
        for (u32 j = 1; j < cells_per_dimension-1; ++j) {
            for (u32 i = 1; i < cells_per_dimension-1; ++i) {
                u64 id = index_3D_to_1D(i, j, k, cells_per_dimension);

                u64 ip = index_3D_to_1D(i+1, j, k, cells_per_dimension);
                u64 im = index_3D_to_1D(i-1, j, k, cells_per_dimension);
                
                u64 jp = index_3D_to_1D(i, j+1, k, cells_per_dimension);
                u64 jm = index_3D_to_1D(i, j-1, k, cells_per_dimension);
                
                u64 kp = index_3D_to_1D(i, j, k+1, cells_per_dimension);
                u64 km = index_3D_to_1D(i, j, k-1, cells_per_dimension);

                // Cada metrics[id].xN ya almacena el tamaño físico real de la celda en esa dirección
                f64 h1 = metrics[id].x1;
                f64 h2 = metrics[id].x2;
                f64 h3 = metrics[id].x3;

                // Paso físico en cada dirección
                f64 ds1 = h1 * dx1;
                f64 ds2 = h2 * dx2;
                f64 ds3 = h3 * dx3;

                f64 df1 = (scalar_field[ip] - 2*scalar_field[id] + scalar_field[im]) / (ds1 * ds1);
                f64 df2 = (scalar_field[jp] - 2*scalar_field[id] + scalar_field[jm]) / (ds2 * ds2);
                f64 df3 = (scalar_field[kp] - 2*scalar_field[id] + scalar_field[km]) / (ds3 * ds3);
            
                output[id] = df1 + df2 + df3;
            }
        }
    }
}