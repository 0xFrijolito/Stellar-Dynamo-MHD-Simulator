#ifndef MESH_H
#define MESH_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "core/defines.h"
#include "core/math/vector.h"

typedef enum GEOMETRIC_TYPE {
    RECTANGULAR = 0,
    CYLINDRICAL = 1,
    SPHERICAL   = 2,
    GENERIC     = 3,
} GEOMETRIC_TYPE;

typedef struct Mesh {
    GEOMETRIC_TYPE geometric;

    Vector3 *positions;
    Vector3 *metrics;
    Vector3 *sizes;

    f64 max_x1;
    f64 max_x2;
    f64 max_x3;

    f64 dx1;
    f64 dx2;
    f64 dx3;

    u64 cells_per_dimension;
    u64 total_nodes;
} Mesh ;

void init_mesh (Mesh *mesh, GEOMETRIC_TYPE geometric, u64 cells_per_dimension, f64 max_x1, f64 max_x2, f64 max_x3);
void free_mesh (Mesh *mesh);
               
void init_rectangular_mesh (Mesh *mesh, u64 cells_per_dimension, f64 max_x, f64 max_y, f64 max_z);
void init_cylindrical_mesh (Mesh *mesh, u64 cells_per_dimension, f64 max_r, f64 max_t, f64 max_z);
void init_spherical_mesh   (Mesh *mesh, u64 cells_per_dimension, f64 max_r, f64 max_t, f64 max_p);


#endif