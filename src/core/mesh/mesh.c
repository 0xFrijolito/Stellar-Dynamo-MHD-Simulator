#include "core/mesh/mesh.h"

#ifndef PI
    #define PI 3.14159265358979323846
#endif

void init_mesh(Mesh *mesh, GEOMETRIC_TYPE geometric, u64 cells_per_dimension, f64 max_x1, f64 max_x2, f64 max_x3) {
    
    // Set the cells resolution
    mesh->cells_per_dimension = cells_per_dimension;

    // Allocate memory for the fields
    size_t cell_count = cells_per_dimension 
                      * cells_per_dimension 
                      * cells_per_dimension;

    mesh->positions = (Vector3*)malloc(sizeof(Vector3) * cell_count);
    mesh->metrics   = (Vector3*)malloc(sizeof(Vector3) * cell_count);
    mesh->sizes     = (Vector3*)malloc(sizeof(Vector3) * cell_count);

    switch (geometric) {
        case RECTANGULAR:
            init_rectangular_mesh(mesh, cells_per_dimension, max_x1, max_x2, max_x3);
            break ;
        case CYLINDRICAL:
            init_cylindrical_mesh(mesh, cells_per_dimension, max_x1, 2 * PI, max_x3);
            break ;
        case SPHERICAL:
            init_spherical_mesh(mesh, cells_per_dimension, max_x1, PI, PI * 2);
            break ;
        default:
            return ;
    }
}

void init_rectangular_mesh (Mesh *mesh, u64 cells_per_dimension, f64 max_x, f64 max_y, f64 max_z) {
    mesh->max_x1 = max_x;
    mesh->max_x2 = max_y;
    mesh->max_x3 = max_z;
    
    mesh->dx1 = max_x / cells_per_dimension;
    mesh->dx2 = max_y / cells_per_dimension;
    mesh->dx3 = max_z / cells_per_dimension;

    for (u32 i=0 ; i<cells_per_dimension ; i++) {
        for (u32 j=0 ; j<cells_per_dimension ; j++) {
            for (u32 k=0 ; k<cells_per_dimension ; k++) {
                f64 x = i * mesh->dx1;
                f64 y = j * mesh->dx2;
                f64 z = k * mesh->dx3;

                // Obtain the index of the cell with (i, j, k position)
                u32 idx = i * cells_per_dimension * cells_per_dimension + 
                          j * cells_per_dimension + k;

                // Set cell positions   
                mesh->positions[idx].x1 = x;
                mesh->positions[idx].x2 = y;
                mesh->positions[idx].x3 = z;

                // Set cell sizes
                mesh->sizes[idx].x1 = mesh->dx1;
                mesh->sizes[idx].x2 = mesh->dx2;
                mesh->sizes[idx].x3 = mesh->dx3;

                // Set metrics
                mesh->metrics[idx].x1 = 1.00f;
                mesh->metrics[idx].x2 = 1.00f;
                mesh->metrics[idx].x3 = 1.00f;
            } 
        }
    }
}

void init_cylindrical_mesh (Mesh *mesh, u64 cells_per_dimension, f64 max_r, f64 max_t, f64 max_z) {
    mesh->max_x1 = max_r;
    mesh->max_x2 = max_t;
    mesh->max_x3 = max_z;
    
    mesh->dx1 = max_r / cells_per_dimension;
    mesh->dx2 = max_t / cells_per_dimension;
    mesh->dx3 = max_z / cells_per_dimension;

    for (u32 i=0 ; i<cells_per_dimension ; i++) {
        for (u32 j=0 ; j<cells_per_dimension ; j++) {
            for (u32 k=0 ; k<cells_per_dimension ; k++) {
                f64 r = i * mesh->dx1;
                f64 t = j * mesh->dx2;
                f64 z = k * mesh->dx3;

                // Obtain the index of the cell with (i, j, k position)
                u32 idx = i * cells_per_dimension * cells_per_dimension + 
                          j * cells_per_dimension + k;

                // Set cell positions   
                mesh->positions[idx].x1 = r;
                mesh->positions[idx].x2 = t;
                mesh->positions[idx].x3 = z;

                // Set cell sizes
                mesh->sizes[idx].x1 = mesh->dx1;
                mesh->sizes[idx].x2 = r * mesh->dx2;
                mesh->sizes[idx].x3 = mesh->dx3;

                // Set metrics
                mesh->metrics[idx].x1 = 1.00f;
                mesh->metrics[idx].x2 = r;
                mesh->metrics[idx].x3 = 1.00f;
            }
        }
    }
}

void init_spherical_mesh   (Mesh *mesh, u64 cells_per_dimension, f64 max_r, f64 max_t, f64 max_p) {
    mesh->max_x1 = max_r;
    mesh->max_x2 = max_t;
    mesh->max_x3 = max_p;
    
    mesh->dx1 = max_r / cells_per_dimension;
    mesh->dx2 = max_t / cells_per_dimension;
    mesh->dx3 = max_p / cells_per_dimension;

    for (u32 i=0 ; i<cells_per_dimension ; i++) {
        for (u32 j=0 ; j<cells_per_dimension ; j++) {
            for (u32 k=0 ; k<cells_per_dimension ; k++) {
                f64 r = i * mesh->dx1;
                f64 t = j * mesh->dx2;
                f64 p = k * mesh->dx3;

                // Obtain the index of the cell with (i, j, k position)
                u32 idx = i * cells_per_dimension * cells_per_dimension + 
                          j * cells_per_dimension + k;
                
                // Set positions
                mesh->positions[idx].x1 = r;
                mesh->positions[idx].x2 = t;
                mesh->positions[idx].x3 = p;

                // Set cell sizes
                mesh->sizes[idx].x1 = mesh->dx1;
                mesh->sizes[idx].x2 = r * mesh->dx2;
                mesh->sizes[idx].x3 = r * sin(t) * mesh->dx3;

                // Set metrics
                mesh->metrics[idx].x1 = 1.00f;
                mesh->metrics[idx].x2 = r;
                mesh->metrics[idx].x3 = r * sin(t);
            }
        }
    }
}

void free_mesh (Mesh *mesh) {
    free (mesh->positions);
    free (mesh->metrics);
    free (mesh->sizes);
}