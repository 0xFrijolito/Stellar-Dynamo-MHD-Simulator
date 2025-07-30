#ifndef CELL_H
#define CELL_H

#include "core/defines.h"
#include "core/math/vector.h"

typedef struct Cell {
    Vector3 size;

    Vector3 position;
    Vector3 velocity;
    Vector3 acceleration;

    f64 density;
    f64 pressure;
    f64 temperature;
    f64 internal_energy;

    f64 gamma;

    Vector3 electric_field;
    Vector3 magnetic_field;
    Vector3 current_density;

    f64 electric_energy;
    f64 magnetic_energy;
    f64 total_energy;

    f64 viscosity;
    f64 magnetic_diffusivity;

    bool is_boundary;
} Cell ;

#endif