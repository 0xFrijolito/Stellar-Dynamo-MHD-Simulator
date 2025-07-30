#include "core/grid/grid.h"

#ifndef PI
    #define PI 3.14159265358979323846
#endif

void init_grid(Grid *grid, GEOMETRIC_TYPE geometric, u64 cells_per_dimension,
               ScalarFieldsInitProfiles scalar_profiles, 
               VectorFieldsInitProfiles vector_profiles,
               f64 max_x1, f64 max_x2, f64 max_x3) {
    
    // Set the cells resolution
    grid->cells_per_dimension = cells_per_dimension;

    // Allocate memory for the fields
    size_t cell_count = cells_per_dimension 
                      * cells_per_dimension 
                      * cells_per_dimension;

    grid->positions = (Vector3*)malloc(sizeof(Vector3) * cell_count);
    grid->metrics   = (Vector3*)malloc(sizeof(Vector3) * cell_count);
    grid->sizes     = (Vector3*)malloc(sizeof(Vector3) * cell_count);

    grid->scalar_fields.density     = (f64*)malloc(sizeof(f64) * cell_count);
    grid->scalar_fields.pressure    = (f64*)malloc(sizeof(f64) * cell_count);
    grid->scalar_fields.temperature = (f64*)malloc(sizeof(f64) * cell_count);

    grid->vector_fields.velocity     = (Vector3*)malloc(sizeof(Vector3) * cell_count);
    grid->vector_fields.acceleration = (Vector3*)malloc(sizeof(Vector3) * cell_count);

    switch (geometric) {
        case RECTANGULAR:
            init_rectangular_grid(grid, cells_per_dimension, max_x1, max_x2, max_x3);
            break ;
        case CYLINDRICAL:
            init_cylindrical_grid(grid, cells_per_dimension, max_x1, 2 * PI, max_x3);
            break ;
        case SPHERICAL:
            init_spherical_grid(grid, cells_per_dimension, max_x1, PI, PI * 2);
            break ;
        default:
            return ;
    }

    // Start the inicial values
    for (u32 i=0 ; i<cells_per_dimension ; i++) {
        for (u32 j=0 ; j<cells_per_dimension ; j++) {
            for (u32 k=0 ; k<cells_per_dimension ; k++) {
                // Obtain the index of the cell with (i, j, k position)
                u32 idx = i * cells_per_dimension * cells_per_dimension + 
                          j * cells_per_dimension + k;
                
                // Calculamos los valores de x, y, z de la celda (i, j, k)
                f64 x1 = grid->positions[idx].x1;
                f64 x2 = grid->positions[idx].x2;
                f64 x3 = grid->positions[idx].x3;

                // Set the inicial scalar parameters
                grid->scalar_fields.density[idx]     = scalar_profiles.density_function(x1, x2, x3);
                grid->scalar_fields.pressure[idx]    = scalar_profiles.pressure_function(x1, x2, x3);
                grid->scalar_fields.temperature[idx] = scalar_profiles.temperature_function(x1, x2, x3);
                
                // Set the inicial vectorial parameters
                Vector3 velocity     = vector_profiles.velocity_function(x1, x2, x3);
                Vector3 acceleration = vector_profiles.acceleration_function(x1, x2, x3);

                vector3_deep_copy(&grid->vector_fields.velocity[idx], &velocity);
                vector3_deep_copy(&grid->vector_fields.acceleration[idx], &acceleration);
            }
        }
    }
}

void init_rectangular_grid (Grid *grid, u64 cells_per_dimension, f64 max_x, f64 max_y, f64 max_z) {
    grid->max_x1 = max_x;
    grid->max_x2 = max_y;
    grid->max_x3 = max_z;
    
    grid->dx1 = max_x / cells_per_dimension;
    grid->dx2 = max_y / cells_per_dimension;
    grid->dx3 = max_z / cells_per_dimension;

    for (u32 i=0 ; i<cells_per_dimension ; i++) {
        for (u32 j=0 ; j<cells_per_dimension ; j++) {
            for (u32 k=0 ; k<cells_per_dimension ; k++) {
                f64 x = i * grid->dx1;
                f64 y = j * grid->dx2;
                f64 z = k * grid->dx3;

                // Obtain the index of the cell with (i, j, k position)
                u32 idx = i * cells_per_dimension * cells_per_dimension + 
                          j * cells_per_dimension + k;

                // Set cell positions   
                grid->positions[idx].x1 = x;
                grid->positions[idx].x2 = y;
                grid->positions[idx].x3 = z;

                // Set cell sizes
                grid->sizes[idx].x1 = grid->dx1;
                grid->sizes[idx].x2 = grid->dx2;
                grid->sizes[idx].x3 = grid->dx3;

                // Set metrics
                grid->metrics[idx].x1 = 1.00f;
                grid->metrics[idx].x2 = 1.00f;
                grid->metrics[idx].x3 = 1.00f;
            } 
        }
    }
}

void init_cylindrical_grid (Grid *grid, u64 cells_per_dimension, f64 max_r, f64 max_t, f64 max_z) {
    grid->max_x1 = max_r;
    grid->max_x2 = max_t;
    grid->max_x3 = max_z;
    
    grid->dx1 = max_r / cells_per_dimension;
    grid->dx2 = max_t / cells_per_dimension;
    grid->dx3 = max_z / cells_per_dimension;

    for (u32 i=0 ; i<cells_per_dimension ; i++) {
        for (u32 j=0 ; j<cells_per_dimension ; j++) {
            for (u32 k=0 ; k<cells_per_dimension ; k++) {
                f64 r = i * grid->dx1;
                f64 t = j * grid->dx2;
                f64 z = k * grid->dx3;

                // Obtain the index of the cell with (i, j, k position)
                u32 idx = i * cells_per_dimension * cells_per_dimension + 
                          j * cells_per_dimension + k;

                // Set cell positions   
                grid->positions[idx].x1 = r;
                grid->positions[idx].x2 = t;
                grid->positions[idx].x3 = z;

                // Set cell sizes
                grid->sizes[idx].x1 = grid->dx1;
                grid->sizes[idx].x2 = r * grid->dx2;
                grid->sizes[idx].x3 = grid->dx3;

                // Set metrics
                grid->metrics[idx].x1 = 1.00f;
                grid->metrics[idx].x2 = r;
                grid->metrics[idx].x3 = 1.00f;
            } 
        }
    }
}

void init_spherical_grid   (Grid *grid, u64 cells_per_dimension, f64 max_r, f64 max_t, f64 max_p) {
    grid->max_x1 = max_r;
    grid->max_x2 = max_t;
    grid->max_x3 = max_p;
    
    grid->dx1 = max_r / cells_per_dimension;
    grid->dx2 = max_t / cells_per_dimension;
    grid->dx3 = max_p / cells_per_dimension;

    for (u32 i=0 ; i<cells_per_dimension ; i++) {
        for (u32 j=0 ; j<cells_per_dimension ; j++) {
            for (u32 k=0 ; k<cells_per_dimension ; k++) {
                f64 r = i * grid->dx1;
                f64 t = j * grid->dx2;
                f64 p = k * grid->dx3;

                // Obtain the index of the cell with (i, j, k position)
                u32 idx = i * cells_per_dimension * cells_per_dimension + 
                          j * cells_per_dimension + k;
                
                // Set positions
                grid->positions[idx].x1 = r;
                grid->positions[idx].x2 = t;
                grid->positions[idx].x3 = p;

                // Set cell sizes
                grid->sizes[idx].x1 = grid->dx1;
                grid->sizes[idx].x2 = r * grid->dx2;
                grid->sizes[idx].x3 = r * sin(t) * grid->dx3;

                // Set metrics
                grid->metrics[idx].x1 = 1.00f;
                grid->metrics[idx].x2 = r;
                grid->metrics[idx].x3 = r * sin(t);
            } 
        }
    }
}

void free_grid (Grid *grid) {
    free_scalar_fields(&grid->scalar_fields);
    free_vector_fields(&grid->vector_fields);
    free (grid->positions);
    free (grid->metrics);
    free (grid->sizes);
}

bool export_grid (Grid *grid, u8 *filepath) {
    FILE *fp = fopen(filepath, "wb");
    if (!fp) return false;

    // Guardar metadata
    fwrite(&grid->cells_per_dimension, sizeof(u32), 1, fp);
    fwrite(&grid->max_x1, sizeof(f64), 1, fp);
    fwrite(&grid->max_x2, sizeof(f64), 1, fp);
    fwrite(&grid->max_x3, sizeof(f64), 1, fp);

    u64 total_cells = grid->cells_per_dimension
                    * grid->cells_per_dimension
                    * grid->cells_per_dimension;

    // Guardamos la metrica
    fwrite(grid->metrics, sizeof(Vector3), total_cells, fp);

    // Guardar campos escalares
    fwrite(grid->scalar_fields.temperature, sizeof(f64), total_cells, fp);
    fwrite(grid->scalar_fields.pressure,    sizeof(f64), total_cells, fp);
    fwrite(grid->scalar_fields.density,     sizeof(f64), total_cells, fp);

    // Guardar campos vectoriales
    fwrite(grid->vector_fields.velocity,     sizeof(Vector3), total_cells, fp);
    fwrite(grid->vector_fields.acceleration, sizeof(Vector3), total_cells, fp);

    fclose(fp);
    return true;
}

bool import_grid (Grid *grid, u8 *filepath) {
    FILE *fp = fopen(filepath, "rb");
    if (!fp) return false;

    // Leer metadata
    fread(&grid->cells_per_dimension, sizeof(u32), 1, fp);
    fread(&grid->max_x1, sizeof(f64), 1, fp);
    fread(&grid->max_x2, sizeof(f64), 1, fp);
    fread(&grid->max_x3, sizeof(f64), 1, fp);

    u64 total_cells = (u64)grid->cells_per_dimension
                    * grid->cells_per_dimension
                    * grid->cells_per_dimension;

    // Reservar memoria
    grid->positions = (Vector3*)malloc(sizeof(Vector3) * total_cells);
    grid->metrics   = (Vector3*)malloc(sizeof(Vector3) * total_cells);
    grid->sizes     = (Vector3*)malloc(sizeof(Vector3) * total_cells);

    grid->scalar_fields.temperature = (f64*)malloc(sizeof(f64) * total_cells);
    grid->scalar_fields.pressure    = (f64*)malloc(sizeof(f64) * total_cells);
    grid->scalar_fields.density     = (f64*)malloc(sizeof(f64) * total_cells);

    grid->vector_fields.velocity     = (Vector3*)malloc(sizeof(Vector3) * total_cells);
    grid->vector_fields.acceleration = (Vector3*)malloc(sizeof(Vector3) * total_cells);

    if (!grid->scalar_fields.temperature || !grid->scalar_fields.pressure ||
        !grid->scalar_fields.density || !grid->vector_fields.velocity ||
        !grid->vector_fields.acceleration) {
        fclose(fp);
        return false;
    }
    
    // Leer geometria
    fread(grid->positions, sizeof(f64), total_cells, fp);
    fread(grid->metrics,   sizeof(f64), total_cells, fp);
    fread(grid->sizes,     sizeof(f64), total_cells, fp);

    // Leer campos escalares
    fread(grid->scalar_fields.temperature, sizeof(f64), total_cells, fp);
    fread(grid->scalar_fields.pressure,    sizeof(f64), total_cells, fp);
    fread(grid->scalar_fields.density,     sizeof(f64), total_cells, fp);

    // Leer campos vectoriales
    fread(grid->vector_fields.velocity,     sizeof(Vector3), total_cells, fp);
    fread(grid->vector_fields.acceleration, sizeof(Vector3), total_cells, fp);

    fclose(fp);
    return true;
}