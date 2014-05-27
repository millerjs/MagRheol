/******************************************************************************
 * FILE        : particles.h
 * AUTHOR      : Joshua Miller 
 * DESCRIPTION : STructures and functions for particles
 ******************************************************************************/


#ifndef PARTICLES_H
#define PARTICLES_H

#include <stdlib.h>
#include <stdio.h>

typedef struct{
    double r[3];
    double v[3];
    double F[3];
    double m;
} particle;

particle *particle_new(double x, double y, double z);

#endif
