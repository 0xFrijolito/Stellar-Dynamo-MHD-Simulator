#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include "core/defines.h"

void euler_integrator ();

// Runge-Kutta integrator for solving ODEs
void runge_kutta_integrator (f64 (*func)(f64 x, f64 y), f64 x, f64 y, f64 dy);

void leapfrog_integrator ();

#endif