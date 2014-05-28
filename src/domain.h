/******************************************************************************
 * FILE        : domain.h
 * AUTHOR      : Joshua Miller 
 * DESCRIPTION : Structures and functions for the domain
 ******************************************************************************/


#ifndef DOMAIN_H
#define DOMAIN_H

#include <stdlib.h>
#include <stdio.h>

#include "particles.h"

#define EPS .05 /* kcal/mol */

typedef enum{
    REFLECTING,
    OUTFLOW,
    PERIODIC,
} boundary_t;

typedef struct{
    double dim[3];
    double boundary[3];
    double v0[3];
    int npart;
    particle *p;
} domain;

extern double t;
extern double dt;

domain *domain_new(double x, double y, double z);
int domain_populate(domain *d, int n);
int domain_set_v0(domain *dm, double x, double y, double z);
int domain_set_boundary(domain *dm, int id, boundary_t b);
int update_positions(domain *dm);
int update_forces_velocities(domain *dm);
int print_checkpoint(char* basepath, domain *dm);
double dist(domain *dm, particle *p1, particle *p2, vec r);
double random(double min, double max);

#endif 
