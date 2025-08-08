#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "core/physics/fluids/euler/euler.h"
#include "core/physics/inicial_conditions.h"

#ifndef M_PI
    #define M_PI 3.1415
#endif

// ---------- Parámetros del caso (ajusta a gusto) ----------
static const f64 Lx = 10.0;
static const f64 Ly = 10.0;
static const f64 Lz = 10.0;

static const f64 rho1 = 1.0;   // banda inferior
static const f64 rho2 = 2.0;   // banda superior
static const f64 u1   =  0.5;  // velocidad en y < Ly/2
static const f64 u2   = -0.5;  // velocidad en y > Ly/2
static const f64 p0   = 2.5;   // presión uniforme inicial
static const f64 eps  = 1e-2;  // amplitud de perturbación en v
static const f64 delta= 0.1;   // grosor relativo de la interfase (tanh smoothing) en unidades de Ly (0.1 => suave)

// ---------- IC helpers ----------
static inline f64 smooth_step(f64 y){
    // Transición suave cerca de Ly/2 con espesor ~ delta*Ly
    const f64 s = (y - 0.5*Ly) / (delta*Ly);
    // map: tanh(s) in [-1,1] -> [0,1]
    return 0.5 * (1.0 + tanh(s));
}

static inline f64 rho_y(f64 y){
    // mezcla entre rho1 (abajo) y rho2 (arriba)
    f64 s = smooth_step(y);
    return (1.0 - s)*rho1 + s*rho2;
}

static inline f64 u_y(f64 y){
    // mezcla entre u1 (abajo) y u2 (arriba)
    f64 s = smooth_step(y);
    return (1.0 - s)*u1 + s*u2;
}

static inline f64 v_perturb(f64 x, f64 z){
    // perturbación pequeña para disparar K-H (periodic)
    // puedes hacerla 2D (solo x) o 3D (x y z). Aquí 2D en x.
    return eps * sin(2.0*M_PI * x / Lx);
}

// ---------- IC requeridas por tu API ----------
f64 inicial_density (f64 x, f64 y, f64 z) {
    (void)z;
    return rho_y(y);
}

f64 inicial_energy  (f64 x, f64 y, f64 z) {
    // E = p/(γ-1) + 1/2 ρ |u|^2
    // Usa el MISMO gamma que pases al solver (abajo uso 1.4)
    const f64 rho = inicial_density(x,y,z);
    const f64 u = u_y(y);
    const f64 v = v_perturb(x,z);
    const f64 w = 0.0;
    const f64 gamma = 1.4;
    const f64 ke = 0.5 * rho * (u*u + v*v + w*w);
    return p0/(gamma - 1.0) + ke;
}

f64 inicial_momentum_x (f64 x, f64 y, f64 z) {
    const f64 rho = inicial_density(x,y,z);
    return rho * u_y(y);
}

f64 inicial_momentum_y (f64 x, f64 y, f64 z) {
    const f64 rho = inicial_density(x,y,z);
    return rho * v_perturb(x,z);
}

f64 inicial_momentum_z (f64 x, f64 y, f64 z) {
    (void)x; (void)y; (void)z;
    return 0.0;
}

// ---------- Writer VTK legacy 3D (STRUCTURED_POINTS) ----------
// Escribe un campo escalar 3D (aquí densidad)
static int write_vtk_density(const char* filename,
                             const EulerSolver* solver)
{
    const u32 N = solver->mesh.cells_per_dimension;
    const f64 dx = solver->mesh.dx1;
    const f64 dy = solver->mesh.dx2;
    const f64 dz = solver->mesh.dx3;
    const f64* rho = solver->fields.density;

    FILE* f = fopen(filename, "w");
    if (!f) return -1;

    fprintf(f, "# vtk DataFile Version 3.0\n");
    fprintf(f, "Euler density step %u\n", (unsigned)solver->step);
    fprintf(f, "ASCII\n");
    fprintf(f, "DATASET STRUCTURED_POINTS\n");
    fprintf(f, "DIMENSIONS %u %u %u\n", (unsigned)N, (unsigned)N, (unsigned)N);
    fprintf(f, "ORIGIN 0 0 0\n");
    fprintf(f, "SPACING %.17g %.17g %.17g\n", dx, dy, dz);
    fprintf(f, "POINT_DATA %u\n", (unsigned)(N*N*N));
    fprintf(f, "SCALARS density double 1\n");
    fprintf(f, "LOOKUP_TABLE default\n");

    // VTK STRUCTURED_POINTS espera orden k variando más lento o más rápido?
    // Con tu index: idx = i*N*N + j*N + k  (k es el más contiguo)
    // Para ser consistente, emitimos en el MISMO orden (i,j,k).
    for (u32 i=0;i<N;++i){
        for (u32 j=0;j<N;++j){
            for (u32 k=0;k<N;++k){
                const u64 idx = (u64)i*N*N + (u64)j*N + (u64)k;
                fprintf(f, "%.17g\n", rho[idx]);
            }
        }
    }

    fclose(f);
    return 0;
}

int main () {
    // Resolución y dominio
    const u32 N = 32;       // súbelo a 128 si tu PC aguanta
    const f64 max_x = Lx;
    const f64 max_y = Ly;
    const f64 max_z = Lz;

    // IMPORTANTE: γ != 1.0
    const f64 gamma = 1.4;

    // Paso de tiempo (simple). Luego podemos hacer CFL adaptativo.
    const f64 dt = 5e-3;

    EulerInicialConditions inicial_conditions;
    inicial_conditions.inicial_density    = inicial_density;
    inicial_conditions.inicial_momentum_x = inicial_momentum_x;
    inicial_conditions.inicial_momentum_y = inicial_momentum_y;
    inicial_conditions.inicial_momentum_z = inicial_momentum_z;
    inicial_conditions.inicial_energy     = inicial_energy;

    EulerSolver solver;
    init_euler_solver(&solver, dt, N, max_x, max_y, max_z, gamma, inicial_conditions);

    // Simulación: escribe cada 'every' pasos
    const u32 steps = 600;     // duración de la simulación
    const u32 every = 2;       // frecuencia de escritura

    // Frame inicial
    {
        char path[256];
        snprintf(path, sizeof(path), "outputs/density_%05u.vtk", (unsigned)solver.step);
        if (write_vtk_density(path, &solver) != 0) {
            fprintf(stderr, "Error escribiendo %s\n", path);
        } else {
            printf("Escrito %s\n", path);
        }
    }

    for (u32 s=0; s<steps; ++s) {
        euler_step_solver(&solver);

        if (solver.step % every == 0) {
            char path[256];
            snprintf(path, sizeof(path), "outputs/density_%05u.vtk", (unsigned)solver.step);
            if (write_vtk_density(path, &solver) != 0) {
                fprintf(stderr, "Error escribiendo %s\n", path);
            } else {
                printf("Escrito %s\n", path);
            }
        }
    }

    free_euler_solver(&solver);
    return 0;
}
