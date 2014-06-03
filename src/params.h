/******************************************************************************
 * FILE        : magrehol.c
 * AUTHOR      : Joshua Miller 
 * DESCRIPTION : Main function for magnetorehological simulations
 ******************************************************************************/

#ifndef PARAMS_H
#define PARAMS_H

#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>

#include "libjosh/libjosh.h"
#include "libjosh/threadpool.h"
#include "particles.h"
#include "domain.h"

#define MU_0 1.25663706e-6
#define PI 3.14159265359
#define viscosity .25
#define EPS .05 /* kcal/mol */
#define SIGMA 1.5
#define R 1.5
#define MU_S .223

extern int step;

/* Parameters */
extern double maxt;
extern int checkpoint_interval;
extern double X;
extern double Y;
extern double Z;
extern int npart;

extern double t;
extern double dt;
extern double ratio;
extern double MU;
extern double H;
extern double Q;

int parse_config(char *path);

#endif
