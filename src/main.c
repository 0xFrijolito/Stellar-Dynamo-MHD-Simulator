#include <stdio.h>
#include <stdlib.h>

#include "core/grid/grid.h"
#include "core/solver/diferencial_scalar_operators.h"
#include "core/solver/diferencial_vector_operators.h"

#ifndef PI
    #define PI 3.14159265358979323846
#endif

f64 density_function (f64 x, f64 y, f64 z) {
    return exp(-(x-10)*(x-10)/25.0) * exp(-(y-10)*(y-10)/25.0);
}

Vector3 velocity_function(f64 x, f64 y, f64 z) {
    f64 A1 =  2.0, x1 = 5.0,  y1 = 5.0,  sigma1 = 2.0;
    f64 A2 = -2.0, x2 = 15.0, y2 = 5.0,  sigma2 = 2.0;
    f64 A3 =  2.0, x3 = 10.0, y3 = 15.0, sigma3 = 2.0;

    // Vórtice 1
    f64 dx1 = x - x1, dy1 = y - y1;
    f64 r2_1 = dx1*dx1 + dy1*dy1;
    f64 f1 = A1 * exp(-r2_1 / (sigma1 * sigma1));
    f64 vx1 = -dy1 * f1;
    f64 vy1 =  dx1 * f1;

    // Vórtice 2
    f64 dx2 = x - x2, dy2 = y - y2;
    f64 r2_2 = dx2*dx2 + dy2*dy2;
    f64 f2 = A2 * exp(-r2_2 / (sigma2 * sigma2));
    f64 vx2 = -dy2 * f2;
    f64 vy2 =  dx2 * f2;

    // Vórtice 3
    f64 dx3 = x - x3, dy3 = y - y3;
    f64 r2_3 = dx3*dx3 + dy3*dy3;
    f64 f3 = A3 * exp(-r2_3 / (sigma3 * sigma3));
    f64 vx3 = -dy3 * f3;
    f64 vy3 =  dx3 * f3;

    Vector3 v;
    v.x1 = vx1 + vx2 + vx3;
    v.x2 = vy1 + vy2 + vy3;
    v.x3 = 0.0;

    return v;
}

void continuity_equation_solver_step (Grid *grid, f64 dt) {
    f64 dx1 = grid->dx1;
    f64 dx2 = grid->dx2;
    f64 dx3 = grid->dx3;

    u64 gradient_domain = (grid->cells_per_dimension)
                        * (grid->cells_per_dimension)
                        * (grid->cells_per_dimension);

    Vector3 *field = (Vector3*)malloc(sizeof(Vector3) * gradient_domain);

    // Obtenemos los parametros a usar
    Vector3 *velocity = grid->vector_fields.velocity;
    f64 *density      = grid->scalar_fields.density;

    // Encontramos el campo vectorial a actualizar
    for (u32 k = 0; k < grid->cells_per_dimension; ++k) {
        for (u32 j = 0; j < grid->cells_per_dimension; ++j) {
            for (u32 i = 0; i < grid->cells_per_dimension; ++i) {
                // Calculamos el indice
                u32 idx = i * grid->cells_per_dimension * grid->cells_per_dimension
                        + j * grid->cells_per_dimension + k;
                
                vector3_scalar_mul(velocity[idx], density[idx], &field[idx]);
            }
        }
    }
    
    f64 *divergence = (f64*)malloc(sizeof(f64) * gradient_domain);

    vector_divergence(field, divergence, grid->cells_per_dimension, dx1, dx2, dx3);

    for (u32 k = 1; k < grid->cells_per_dimension-1; ++k) {
        for (u32 j = 1; j < grid->cells_per_dimension-1; ++j) {
            for (u32 i = 1; i < grid->cells_per_dimension-1; ++i) {
                // Calculamos el indice
                u32 idx = i * grid->cells_per_dimension * grid->cells_per_dimension
                        + j * grid->cells_per_dimension + k;
                
                density[idx] -=  divergence[idx] * dt;
            }
        }
    }

    free(divergence);
    free(field);
}

int main() {    
    // Configuracion inicial de los parametros fisicos.
    ScalarFieldsInitProfiles scalar_init_profiles;
    scalar_init_profiles.temperature_function = none_scalar_function;
    scalar_init_profiles.density_function     = density_function;
    scalar_init_profiles.pressure_function    = none_scalar_function;

    VectorFieldsInitProfiles vector_init_profiles;
    vector_init_profiles.velocity_function     = velocity_function;
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

    f64 dt    = 1e-1;
    u32 steps = 500;

    printf("[+] Iniciando la simulacion.\n");
    for (u32 step = 0; step < steps; step++) {

        continuity_equation_solver_step(&grid, dt);
                
        fwrite(grid.scalar_fields.temperature, sizeof(f64), total_cells, fp);
        fwrite(grid.scalar_fields.pressure,    sizeof(f64), total_cells, fp);
        fwrite(grid.scalar_fields.density,     sizeof(f64), total_cells, fp);
    }
    printf("[+] Fin de la simulacion.\n");

    fclose(fp);

    return 0;
}