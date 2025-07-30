#include <stdio.h>
#include <stdlib.h>

#include "core/grid/grid.h"
#include "core/solver/diferencial_scalar_operators.h"
#include "core/solver/diferencial_vector_operators.h"

f64 temperature_function (f64 x, f64 y, f64 z) {
    f64 center_x = x - 10;
    f64 center_y = y - 10;
    f64 center_z = z - 10;

    return center_x*center_x + center_y*center_y + center_z*center_z;
}

f64 none_scalar_function (f64 x, f64 y, f64 z) {
    return 0.0f;
}

Vector3 none_vectorial_function (f64 x, f64 y, f64 z) {
    Vector3 vector = {0.0f, 0.0f, 0.0f};
    return vector;
}

void heat_ecuation_solver_step (Grid *grid, f64 alpha, f64 dt) {
    f64 dx = grid->max_x3 / grid->cells_per_dimension;
    f64 dy = grid->max_x2 / grid->cells_per_dimension;
    f64 dz = grid->max_x1 / grid->cells_per_dimension;

    u64 gradient_domain = (grid->cells_per_dimension)
                        * (grid->cells_per_dimension)
                        * (grid->cells_per_dimension);

    // Obtenemos la geometrica del sistema actual.
    RectangularGrid *current_geometry = grid->rectangular_grid;

    // Obtenemos los parametros a usar
    f64 *temperature = current_geometry->scalar_fields.temperature; 
    f64 *laplacian   = (f64*)malloc(sizeof(f64) * gradient_domain);

    // Calculamos el laplaciano
    scalar_laplacian(
        temperature, laplacian, 
        grid->cells_per_dimension,
        dx, dy, dz
    );

    // Actualizamos la temperura
    for (u32 k = 1; k < grid->cells_per_dimension-1; ++k) {
        for (u32 j = 1; j < grid->cells_per_dimension-1; ++j) {
            for (u32 i = 1; i < grid->cells_per_dimension-1; ++i) {

                u32 idx = i * grid->cells_per_dimension * grid->cells_per_dimension
                        + j * grid->cells_per_dimension + k;

                temperature[idx] += alpha * laplacian[idx] * dt;
            }
        }
    }

    free(laplacian);
}

int main() {    
    // Configuracion inicial de los parametros fisicos.
    ScalarFieldsInitProfiles scalar_init_profiles;
    scalar_init_profiles.temperature_function = temperature_function;
    scalar_init_profiles.density_function     = none_scalar_function;
    scalar_init_profiles.pressure_function    = none_scalar_function;

    VectorFieldsInitProfiles vector_init_profiles;
    vector_init_profiles.velocity_function     = none_vectorial_function;
    vector_init_profiles.acceleration_function = none_vectorial_function;

    // Tipo de grid a usar
    GRID_TYPE type = RECTANGULAR;

    // Tamano y resolucion de la simulacion
    f64 max_x = 20;
    f64 max_y = 20;
    f64 max_z = 20;

    u32 cells_per_dimension = 64;

    // Grid
    Grid grid;
    init_grid(
        &grid, type, 
        cells_per_dimension, 
        scalar_init_profiles,
        vector_init_profiles,
        max_x, 
        max_y, 
        max_z
    );

    FILE* fp = fopen("output.bin", "wb");
    if (!fp) return false;

    // Guardar metadata
    fwrite(&grid.rectangular_grid->cells_per_dimension, sizeof(u32), 1, fp);
    fwrite(&grid.rectangular_grid->max_x, sizeof(f64), 1, fp);
    fwrite(&grid.rectangular_grid->max_y, sizeof(f64), 1, fp);
    fwrite(&grid.rectangular_grid->max_z, sizeof(f64), 1, fp);

    u64 total_cells = grid.rectangular_grid->cells_per_dimension
                    * grid.rectangular_grid->cells_per_dimension
                    * grid.rectangular_grid->cells_per_dimension;

    f64 dt = 1e-2;
    u32 N  = grid.cells_per_dimension;
    f64 dx = grid.max_x1 / N;
    f64 dy = grid.max_x2 / N;
    f64 dz = grid.max_x3 / N;

    u32 steps = 500;        // puedes ajustar el número de pasos
    f64 alpha = 1.3;          // difusividad térmica, ajusta según tu modelo

    for (u32 step = 0; step < steps; ++step) {
        f64 time = step * dt;

        heat_ecuation_solver_step(&grid, alpha, dt);
                
        // Imprimir temperatura en cada celda
        for (u32 k = 0; k < N; ++k) {
            for (u32 j = 0; j < N; ++j) {
                for (u32 i = 0; i < N; ++i) {
                    // Calcualte the positions
                    f64 x = i * dx;
                    f64 y = j * dy;
                    f64 z = k * dz;
                    
                    // Obtain the index of the cell with (i, j, k position)
                    u32 idx = i * cells_per_dimension * cells_per_dimension 
                            + j * cells_per_dimension + k;
                    
                    f64 T = grid.rectangular_grid->scalar_fields.temperature[idx];
                }
            }
        }

        fwrite(grid.rectangular_grid->scalar_fields.temperature, sizeof(f64), total_cells, fp);
        fwrite(grid.rectangular_grid->scalar_fields.pressure,    sizeof(f64), total_cells, fp);
        fwrite(grid.rectangular_grid->scalar_fields.density,     sizeof(f64), total_cells, fp);
    }

    fclose(fp);

    return 0;
}