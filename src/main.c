#include <stdio.h>
#include <stdlib.h>

#include "core/grid/grid.h"
#include "core/solver/diferencial_scalar_operators.h"
#include "core/solver/diferencial_vector_operators.h"

#ifndef PI
    #define PI 3.14159265358979323846
#endif

f64 temperature_function (f64 x, f64 y, f64 z) {
    if (fabs(z - 10) < 1.0)
        return 100.0;
    else
        return 0.0;
}

f64 none_scalar_function (f64 x, f64 y, f64 z) {
    return 0.0f;
}

Vector3 none_vectorial_function (f64 x, f64 y, f64 z) {
    Vector3 vector = {0.0f, 0.0f, 0.0f};
    return vector;
}

void heat_ecuation_solver_step (Grid *grid, f64 alpha, f64 dt) {
    f64 dx1 = grid->dx1;
    f64 dx2 = grid->dx2;
    f64 dx3 = grid->dx3;

    u64 gradient_domain = (grid->cells_per_dimension)
                        * (grid->cells_per_dimension)
                        * (grid->cells_per_dimension);

    // Obtenemos la metrica de la geometria usada
    Vector3 *metric = grid->metrics;

    // Obtenemos los parametros a usar
    f64 *temperature = grid->scalar_fields.temperature; 
    f64 *laplacian   = (f64*)malloc(sizeof(f64) * gradient_domain);

    // Calculamos el laplaciano
    scalar_laplacian(
        temperature, metric, laplacian,
        grid->cells_per_dimension,
        dx1, dx2, dx3
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

    // Tipo de geometrica usamos
    GEOMETRIC_TYPE geometric = RECTANGULAR;

    // Tamano y resolucion de la simulacion
    f64 max_x = 20;
    f64 max_y = 20;
    f64 max_z = 20;

    u32 cells_per_dimension = 32;

    // Grid
    Grid grid;
    init_grid(
        &grid, geometric, 
        cells_per_dimension, 
        scalar_init_profiles,
        vector_init_profiles,
        max_x, 
        max_y, 
        max_z
    );

    FILE* fp = fopen("output.bin", "wb");
    if (!fp) return false;

    printf("[+] Archivo creado. \n");
    
    // Guardar metadata
    fwrite(&grid.cells_per_dimension, sizeof(u32), 1, fp);
    fwrite(&grid.max_x1, sizeof(f64), 1, fp);
    fwrite(&grid.max_x2, sizeof(f64), 1, fp);
    fwrite(&grid.max_x3, sizeof(f64), 1, fp);

    printf("[+] Metadatos guardados en el archivo.\n");

    u64 total_cells = grid.cells_per_dimension
                    * grid.cells_per_dimension
                    * grid.cells_per_dimension;

    fwrite(grid.positions, sizeof(Vector3), total_cells, fp);

    printf("[+] Posiciones guardas en el archivo. \n");

    f64 dt    = 1e-2;
    u32 steps = 500;
    f64 alpha = 2;

    for (u32 step = 0; step < steps; ++step) {
        heat_ecuation_solver_step(&grid, alpha, dt);
                
        fwrite(grid.scalar_fields.temperature, sizeof(f64), total_cells, fp);
        fwrite(grid.scalar_fields.pressure,    sizeof(f64), total_cells, fp);
        fwrite(grid.scalar_fields.density,     sizeof(f64), total_cells, fp);
    }

    fclose(fp);

    return 0;
}