#ifndef EULER_H
#define EULER_H

#include <stdio.h>
#include <stdlib.h>

#include "core/defines.h"
#include "core/mesh/mesh.h"

#include "core/math/vector.h"
#include "core/math/matrix.h"

#include "core/physics/fields.h"
#include "core/physics/inicial_conditions.h"

#include "core/numerics/finite_elements/diferencial_vector_operators.h"
#include "core/numerics/integrator.h"

typedef struct EulerSolver {
    Mesh        mesh;
    EulerFields fields;
    EulerFields fields_flux;

    Vector3 *density_flux;
    Vector3 *momentum_flux_u;
    Vector3 *momentum_flux_v;
    Vector3 *momentum_flux_w;
    Vector3 *energy_flux;

    f64 *div_density;
    f64 *div_u;
    f64 *div_v;
    f64 *div_w;
    f64 *div_energy;

    f64 gamma;

    f64 dt;
    u32 step;

    u8 is_recording;
} EulerSolver ;

void init_euler_solver (EulerSolver *solver, f64 dt,
                        u32 cells_per_dimension,
                        f64 max_x, f64 max_y, f64 max_z, 
                        f64 gamma, EulerInicialConditions inicial_conditions);

void free_euler_solver (EulerSolver *solver);

void euler_step_solver (EulerSolver *solver);
void euler_solve       (EulerSolver *solver, u32 steps);

#endif