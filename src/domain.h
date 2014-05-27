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

typedef struct{
    double dim[3];
    double v0[3];
    int npart;
    particle *p;
} domain;

extern double t;
extern double dt;

domain *domain_new(double x, double y, double z);
int domain_populate(domain *d, int n);

#endif 
