/******************************************************************************
 * FILE        : particles.c
 * AUTHOR      : Joshua Miller 
 * DESCRIPTION : Structures and functions for particles
 ******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

#include "particles.h"

particle *particle_new(double x, double y, double z)
{
    particle *new = malloc(sizeof(particle));
    new->r[0] = x;
    new->r[1] = y;
    new->r[2] = z;
    return new;
}



