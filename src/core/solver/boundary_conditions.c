#include "core/solver/boundary_conditions.h"

static u32 index_3D_to_1D (u32 i, u32 j, u32 k, u32 size) {
    return i * size * size + 
           j * size + k;
}


void scalar_dirichlet_boundary_conditions (f64 *scalar_field, f64 boundary_value, u32 cells_per_dimension) {
    u32 coord_max_value = cells_per_dimension - 1;
    
    for (u32 j=0 ; j<cells_per_dimension ; j++) {
        for (u32 k=0 ; k<cells_per_dimension ; k++) {
            u32 i1 = index_3D_to_1D(0, j, k, cells_per_dimension);
            u32 i2 = index_3D_to_1D(coord_max_value, j, k, cells_per_dimension);

            scalar_field[i1] = 0.0f;
            scalar_field[i2] = 0.0f;
        }
    }

    for (u32 i=0 ; i<cells_per_dimension ; i++) {
        for (u32 k=0 ; k<cells_per_dimension ; k++) {
            u32 i1 = index_3D_to_1D(i, 0, k, cells_per_dimension);
            u32 i2 = index_3D_to_1D(i, coord_max_value, k, cells_per_dimension);

            scalar_field[i1] = 0.0f;
            scalar_field[i2] = 0.0f;
        }
    }

    for (u32 i=0 ; i<cells_per_dimension; i++) {
        for (u32 j=0 ; j<cells_per_dimension; j++) {
            u32 i1 = index_3D_to_1D(i, j, 0, cells_per_dimension);
            u32 i2 = index_3D_to_1D(i, j, coord_max_value, cells_per_dimension);

            scalar_field[i1] = 0.0f;
            scalar_field[i2] = 0.0f;
        }
    }
}

void vector_dirichlet_boundary_conditions (Vector3 *vector_field, Vector3 *boundary_value, u32 cell_per_dimension) {
    
} 
