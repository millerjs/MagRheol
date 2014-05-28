/******************************************************************************
 * FILE        : particles.h
 * AUTHOR      : Joshua Miller 
 * DESCRIPTION : STructures and functions for particles
 ******************************************************************************/


#ifndef PARTICLES_H
#define PARTICLES_H

#include <stdlib.h>
#include <stdio.h>

/* #include "domain.h" */

typedef double vec[3];

typedef struct{
    vec r;
    vec v;
    vec F;
    double m;
    double sigma;
} particle;

    
int particle_new(particle *p);
double norm(vec r);
void add(vec v1, vec v2, vec res);
void scale(vec v, double a);
void print_vec(vec v);


#endif
