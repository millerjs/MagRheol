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
#define SIGMA 3.4 /* kcal/mol */

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
    double *F;
    double *r;
    double *oldr;
    double *temp;
    double *v;
} domain;

domain *domain_new(double x, double y, double z);
int domain_populate(domain *d, int n);
int domain_set_v0(domain *dm, double x, double y, double z);
int domain_set_boundary(domain *dm, int id, boundary_t b);
int update_positions(domain *dm, int a, int b);
int print_checkpoint(char* basepath, domain *dm);
double dist(domain *dm, int m, int i, vec r);
double randomd(double min, double max);
void check_boundary(domain *dm, int m);

#endif 
