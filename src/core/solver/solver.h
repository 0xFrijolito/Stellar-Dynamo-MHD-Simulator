#ifndef SOLVER_H
#define SOLVER_H

#include <stdio.h>

#include "core/defines.h"
#include "core/math/vector.h"
#include "core/grid/grid.h"
#include "core/grid/cell.h"

typedef struct Solver {
    void *grid;
    
    f64 dt;
    f64 time;
    u32 step;

    u8 recording;
    FILE *fd;
} Solver ;

void init_solver(Solver *solver, Grid *grid, f64 dt, u8 recording, const u8 *record_path);
void free_solver(Solver *solver);

void solver_step (Solver *solver);
void solver_run  (Solver *solver, u32 total_steps);

#endif