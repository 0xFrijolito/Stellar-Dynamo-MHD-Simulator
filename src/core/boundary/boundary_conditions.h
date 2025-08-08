#ifndef BOUNDARY_CONDITIONS_H
#define BOUNDARY_CONDITIONS_H

typedef enum BoundaryConditionType {
    BC_PERIODIC  = 0,
    BC_NEUMANN   = 1,
    BC_DIRICHLET = 2,
    BC_WALL      = 3
} BoundaryConditionType ;



#endif