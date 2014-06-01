/******************************************************************************
 * FILE        : particles.c
 * AUTHOR      : Joshua Miller 
 * DESCRIPTION : Structures and functions for particles
 ******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>

#include "particles.h"

const int VECSIZE = 3*sizeof(double);

int particle_new(particle *p)
{
    for (int i = 0; i < 3; i ++){
        p->r[i] = 0;
        p->v[i] = 0;
        p->F[i] = 0;
        p->sigma = 3.7;
    }
    return 0;
}

void add(vec v1, vec v2, vec res)
{
    for (int i = 0; i < 3; i++){
        res[i] = v1[i] + v2[i];
    }
}

void scale(vec v, double a)
{
    for (int i = 0; i < 3; i++){
        v[i] *= a;
    }
}

void print_vec(vec v)
{
    fprintf(stderr, "[%.3f, %.3f, %.3f]\n", v[0], v[1], v[2]);
}
