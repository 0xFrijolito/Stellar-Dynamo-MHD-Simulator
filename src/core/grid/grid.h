#ifndef GRID_H
#define GRID_H

#include <math.h>
#include <stdio.h>

#include "core/defines.h"
#include "core/math/vector.h"

#include "core/grid/scalar_fields.h"
#include "core/grid/vector_fields.h"

typedef enum GEOMETRIC_TYPE {
    RECTANGULAR = 0,
    CYLINDRICAL = 1,
    SPHERICAL   = 2,
    GENERIC     = 3,
} GEOMETRIC_TYPE;

typedef struct Grid {
    GEOMETRIC_TYPE geometric;

    Vector3 *positions;
    Vector3 *metrics;
    Vector3 *sizes;

    ScalarFields scalar_fields;
    VectorFields vector_fields;

    u32 cells_per_dimension;

    f64 max_x1;
    f64 max_x2;
    f64 max_x3;

    f64 dx1;
    f64 dx2;
    f64 dx3;
} Grid ;

void init_grid (Grid *grid, GEOMETRIC_TYPE geometric, u64 cells_per_dimension,
                ScalarFieldsInitProfiles scalar_profiles, 
                VectorFieldsInitProfiles vector_profiles,
                f64 max_x1, f64 max_x2, f64 max_x3);
               
void init_rectangular_grid (Grid *grid, u64 cells_per_dimension, f64 max_x, f64 max_y, f64 max_z);
void init_cylindrical_grid (Grid *grid, u64 cells_per_dimension, f64 max_r, f64 max_t, f64 max_z);
void init_spherical_grid   (Grid *grid, u64 cells_per_dimension, f64 max_r, f64 max_t, f64 max_p);

void free_grid (Grid *grid);

// void set_grid_cell (Grid *grid, Vector3 position, Cell  cell);
// void get_grid_cell (Grid *grid, Vector3 position, Cell *cell);

bool export_grid (Grid *grid, u8 *filepath);
bool import_grid (Grid *grid, u8 *filepath);

#endif