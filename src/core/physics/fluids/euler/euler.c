#include "core/physics/fluids/euler/euler.h"

void init_euler_solver (EulerSolver *solver, f64 dt,
                        u32 cells_per_dimension,
                        f64 max_x, f64 max_y, f64 max_z, 
                        f64 gamma, 
                        EulerInicialConditions inicial_conditions) {

    // Inicializar la malla
    init_mesh (&solver->mesh, RECTANGULAR, cells_per_dimension, max_x, max_y, max_z);
    
    u64 total_nodes = cells_per_dimension
                    * cells_per_dimension
                    * cells_per_dimension;
                    
    // Asignar valores al solver
    solver->dt    = dt;
    solver->step  = 0;
    solver->gamma = gamma;
    
    // Reservar memoria para los campos de Euler
    solver->fields.density    = (f64 *) malloc(total_nodes * sizeof(f64));
    solver->fields.momentum_x = (f64 *) malloc(total_nodes * sizeof(f64));
    solver->fields.momentum_y = (f64 *) malloc(total_nodes * sizeof(f64));
    solver->fields.momentum_z = (f64 *) malloc(total_nodes * sizeof(f64));
    solver->fields.energy     = (f64 *) malloc(total_nodes * sizeof(f64));
    
    // Reservar memoria para los campos de flujo
    solver->density_flux    = (Vector3 *) malloc(total_nodes * sizeof(Vector3));
    solver->momentum_flux_u = (Vector3*) malloc(total_nodes * sizeof(Vector3));
    solver->momentum_flux_v = (Vector3*) malloc(total_nodes * sizeof(Vector3));
    solver->momentum_flux_w = (Vector3*) malloc(total_nodes * sizeof(Vector3));
    solver->energy_flux     = (Vector3*) malloc(total_nodes * sizeof(Vector3));
    
    // Reservar memoria para las divergencias.
    solver->div_density = (f64 *)malloc(total_nodes * sizeof(f64));
    solver->div_u       = (f64 *)malloc(total_nodes * sizeof(f64));
    solver->div_v       = (f64 *)malloc(total_nodes * sizeof(f64));
    solver->div_w       = (f64 *)malloc(total_nodes * sizeof(f64));
    solver->div_energy  = (f64 *)malloc(total_nodes * sizeof(f64));

    for (u64 i=0 ; i<cells_per_dimension ; i++) {
        for (u64 j=0 ; j<cells_per_dimension ; j++) {
            for (u64 k=0 ; k<cells_per_dimension ; k++) {
                // Calcular las coordenadas x, y, z de cada celda
                f64 x = i * solver->mesh.dx1;
                f64 y = j * solver->mesh.dx2;
                f64 z = k * solver->mesh.dx3;

                // Calcular el índice lineal para acceder a los arrays
                u64 idx = i * solver->mesh.cells_per_dimension * solver->mesh.cells_per_dimension +
                          j * solver->mesh.cells_per_dimension + k;

                // Inicializar los campos con las condiciones iniciales
                solver->fields.density[idx]    = inicial_conditions.inicial_density(x, y, z);
                solver->fields.momentum_x[idx] = inicial_conditions.inicial_momentum_x(x, y, z);
                solver->fields.momentum_y[idx] = inicial_conditions.inicial_momentum_y(x, y, z);
                solver->fields.momentum_z[idx] = inicial_conditions.inicial_momentum_z(x, y, z);
                solver->fields.energy[idx]     = inicial_conditions.inicial_energy(x, y, z);
            }
        }
    }
}

void free_euler_solver (EulerSolver *solver) {
    // Liberar memoria de los campos de Euler
    free(solver->fields.density);
    free(solver->fields.momentum_x);
    free(solver->fields.momentum_y);
    free(solver->fields.momentum_z);
    free(solver->fields.energy);
    free(solver->div_density);
    free(solver->div_u);
    free(solver->div_v);
    free(solver->div_w);
    free(solver->div_energy);

    // Liberar memoria de la malla
    free_mesh(&solver->mesh);
}

void euler_step_solver (EulerSolver *solver) {
    // Calcular el tamaño total de nodos
    u64 total_nodes = solver->mesh.cells_per_dimension
                    * solver->mesh.cells_per_dimension
                    * solver->mesh.cells_per_dimension;
    
    Mesh *mesh = &solver->mesh;

    f64 *density    = solver->fields.density;
    f64 *momentum_x = solver->fields.momentum_x;
    f64 *momentum_y = solver->fields.momentum_y;
    f64 *momentum_z = solver->fields.momentum_z;
    f64 *energy     = solver->fields.energy;

    Vector3 *density_flux    = solver->density_flux;
    Vector3 *momentum_flux_u = solver->momentum_flux_u;
    Vector3 *momentum_flux_v = solver->momentum_flux_v;
    Vector3 *momentum_flux_w = solver->momentum_flux_w;
    Vector3 *energy_flux     = solver->energy_flux;

    f64 *div_density = solver->div_density;
    f64 *div_u       = solver->div_u;
    f64 *div_v       = solver->div_v;
    f64 *div_w       = solver->div_w;
    f64 *div_energy  = solver->div_energy;

    for (u64 i=0 ; i<total_nodes ; i++) {
        f64 u = solver->fields.momentum_x[i] / solver->fields.density[i];
        f64 v = solver->fields.momentum_y[i] / solver->fields.density[i];
        f64 w = solver->fields.momentum_z[i] / solver->fields.density[i];

        f64 ke = 0.5 * density[i] * (u*u + v*v + w*w);
        f64 p  = (solver->gamma - 1.0) * (energy[i] - ke);

        density_flux[i].x1 = solver->fields.density[i] * u;
        density_flux[i].x2 = solver->fields.density[i] * v;
        density_flux[i].x3 = solver->fields.density[i] * w;

        momentum_flux_u[i].x1 = density[i] * u * u + p;  
        momentum_flux_u[i].x2 = density[i] * u * v;
        momentum_flux_u[i].x3 = density[i] * u * w;

        momentum_flux_v[i].x1 = density[i] * v * u;
        momentum_flux_v[i].x2 = density[i] * v * v + p;
        momentum_flux_v[i].x3 = density[i] * v * w;

        momentum_flux_w[i].x1 = density[i] * w * u;
        momentum_flux_w[i].x2 = density[i] * w * v;
        momentum_flux_w[i].x3 = density[i] * w * w + p;

        energy_flux[i].x1 = (energy[i] + p) * u;
        energy_flux[i].x2 = (energy[i] + p) * v;
        energy_flux[i].x3 = (energy[i] + p) * w;
    }

    vector_divergence_periodic(density_flux,    div_density, mesh->cells_per_dimension, mesh->dx1, mesh->dx2, mesh->dx3);
    vector_divergence_periodic(momentum_flux_u, div_u,       mesh->cells_per_dimension, mesh->dx1, mesh->dx2, mesh->dx3);
    vector_divergence_periodic(momentum_flux_v, div_v,       mesh->cells_per_dimension, mesh->dx1, mesh->dx2, mesh->dx3);
    vector_divergence_periodic(momentum_flux_w, div_w,       mesh->cells_per_dimension, mesh->dx1, mesh->dx2, mesh->dx3);
    vector_divergence_periodic(energy_flux,     div_energy,  mesh->cells_per_dimension, mesh->dx1, mesh->dx2, mesh->dx3);

    for (u64 i=0 ; i<total_nodes ; i++) {
        density[i]    -= solver->dt * div_density[i];
        momentum_x[i] -= solver->dt * div_u[i];
        momentum_y[i] -= solver->dt * div_v[i];
        momentum_z[i] -= solver->dt * div_w[i];
        energy[i]     -= solver->dt * div_energy[i];
    }

    solver->step += 1;
}

void euler_solve (EulerSolver *solver, u32 steps){
    Mesh *mesh = &solver->mesh;

    f64 *density    = solver->fields.density;
    f64 *momentum_x = solver->fields.momentum_x;
    f64 *momentum_y = solver->fields.momentum_y;
    f64 *momentum_z = solver->fields.momentum_z;
    f64 *energy     = solver->fields.energy;

    u32 cell_per_dimension = mesh->cells_per_dimension;
    u64 total_nodes = cell_per_dimension
                    * cell_per_dimension
                    * cell_per_dimension;
                    
    // always true
    FILE* fp = fopen("output.bin", "wb");
    if (!fp) return ;

    if (1) {
        fwrite(&mesh->cells_per_dimension, sizeof(u32), 1, fp);
        fwrite(&mesh->max_x1, sizeof(f64), 1, fp);
        fwrite(&mesh->max_x2, sizeof(f64), 1, fp);
        fwrite(&mesh->max_x3, sizeof(f64), 1, fp);

    }
    
    for (u32 i=0 ; i<steps ; i++) {
        if (1) {
            fwrite(density   , sizeof(f64), total_nodes, fp);
            fwrite(momentum_x, sizeof(f64), total_nodes, fp);
            fwrite(momentum_y, sizeof(f64), total_nodes, fp);
            fwrite(momentum_z, sizeof(f64), total_nodes, fp);
            fwrite(energy    , sizeof(f64), total_nodes, fp);

        }
        euler_step_solver(solver);
    }

    if (1) {
        fwrite(density   , sizeof(f64), total_nodes, fp);
        fwrite(momentum_x, sizeof(f64), total_nodes, fp);
        fwrite(momentum_y, sizeof(f64), total_nodes, fp);
        fwrite(momentum_z, sizeof(f64), total_nodes, fp);
        fwrite(energy    , sizeof(f64), total_nodes, fp);
    }

    fclose(fp);
}