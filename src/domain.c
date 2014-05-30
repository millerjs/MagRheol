/******************************************************************************
 * FILE        : domain.c
 * AUTHOR      : Joshua Miller 
 * DESCRIPTION : Structures and functions for the domain
 ******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <math.h>

#include "domain.h"
#include "libjosh/libjosh.h"

double t = 0;
double dt = .001;

domain *domain_new(double x, double y, double z)
{
    srand(time(NULL));
    domain *new = (domain*)malloc(sizeof(domain));
    new->dim[0] = x;
    new->dim[1] = y;
    new->dim[2] = z;
    for (int i = 0; i < 3; i ++){
        new->v0[i] = 0;
        new->boundary[i] = REFLECTING;
    }
    LOG("Created domain %.2f x %.2f x %.2f", x, y, z);
    new->npart = 0;
    return new;
}

int domain_set_boundary(domain *dm, int id, boundary_t b){
    return dm->boundary[id] = b;
}

double randomd(double min, double max)
{
    return (rand()/(double)RAND_MAX)*max + min;
}

int domain_populate(domain *dm, int n)
{
    dm->npart = n;
    LOG("Populating domain [%p] with [%d] particles", dm, n);
    ERROR_IF(!(dm->r = malloc(n*VECSIZE)), "Unable to allocate positions");
    ERROR_IF(!(dm->oldr = malloc(n*VECSIZE)), "Unable to allocate old positions");
    ERROR_IF(!(dm->v = malloc(n*VECSIZE)), "Unable to allocate old positions");
    ERROR_IF(!(dm->F = malloc(n*VECSIZE)), "Unable to allocate old positions");
    ERROR_IF(!(dm->temp = malloc(n*VECSIZE)), "Unable to allocate buffer space");

    double r = pow(dm->npart, 1/3.);
    vec cell = {dm->dim[0]/r, dm->dim[1]/r, dm->dim[2]/r};
    LOG("Cube discretization count: %f", r);
    for (int j = 0; j < 3; j ++)
        LOG("Particle discretization distance: %f", cell[j]);

    int m = 0;
    for (int i = 0; i < ceil(r); i++){
        for (int j = 0; j < ceil(r); j++){
            for (int k = 0; k < ceil(r); k++){

                dm->oldr[3*m+0] = (i+.5)*cell[0];
                dm->oldr[3*m+1] = (j+.5)*cell[1];
                dm->oldr[3*m+2] = (k+.5)*cell[2];

                dt = .001;
                for (int d = 0; d < 3; d ++){
                    dm->v[3*m+d] = dm->v0[d] + randomd(-100, 100);
                    dm->r[3*m+d] = dm->oldr[3*m+d] + dm->v[3*m+d]*dt;
                }

                m++;
                if (m >= dm->npart) break;
            }
            if (m >= dm->npart) break;
        }
        if (m >= dm->npart) break;
    }
    LOG("Populated domain with %d particles", m);
    return m;
}

int domain_set_v0(domain *dm, double x, double y, double z)
{
    LOG("Setting initial velocity to [%.3f, %.3f, %.3f]", x, y, z);
    dm->v0[0] = x;
    dm->v0[1] = y;
    dm->v0[2] = z;
    return 0;
}

double dist(domain *dm, int m, int n, vec r)
{
    double temp = 0;
    for (int d = 0; d < 3; d ++){
        r[d] = dm->r[3*m+d] - dm->r[3*n+d];
        r[d] = r[d] - (dm->boundary[d] == PERIODIC)*rint(r[d]/dm->dim[d]);
        temp += r[d]*r[d];    
    }
    return sqrt(temp);
}

void force_Drag(domain *dm, int m)
{
    for (int j = 1; j < 3; j++){
        if (dm->r[3*m+j] < SIGMA || dm->r[3*+j] > dm->dim[j] - SIGMA){
            dm->v[3*m] = 0;
            dm->F[3*m] = 0;
        }
    }
}

void force_Wall(domain *dm, int m)
{
    for (int d = 0; d < 3; d++){
        double r = MIN(2*dm->r[3*m+d], dm->dim[d] - 2*dm->r[3*m+d]);
        r = 1.0;
        double t6 = pow(SIGMA/r, 6);
        double t12 = t6*t6;
        double f = 4*EPS*(12/r*t12 - 6/r*t6);
        for (int d = 0; d < 3; d++)
            dm->F[3*m+d] += MAX(MIN(f, 1e3), -1e3);
    }
}

void force_LJ(domain *dm, int m)
{
    for (int i = 0; i < dm->npart; i++){
        if (i!=m){
            vec res;
            double d = MAX(dist(dm, m, i, res), 1e-4);
            double t6 = pow(SIGMA/d, 6);
            double t12 = t6*t6;
            double f = 4*EPS*(12/d*t12 - 6/d*t6);
            for (int d = 0; d < 3; d++)
                dm->F[3*m+d] += MAX(MIN(res[d]*f, 1e3), -1e3);
        }
    }
}

int calculate_force(domain *dm, int m)
{
    for (int d = 0; d < 3; d++)
        dm->F[3*m+d] = 0;
    force_LJ(dm, m);
    /* force_Wall(dm, m); */
    return 0;
}

void check_boundary(domain *dm, int m)
{
    for (int d = 0; d < 3; d ++){
        if (dm->boundary[d] == PERIODIC){

            if (dm->r[3*m+d] > dm->dim[d]){
                dm->r[3*m+d] -= dm->dim[d];
                dm->oldr[3*m+d] -= dm->dim[d];
            } else if (dm->r[3*m+d] < 0) {
                dm->r[3*m+d] += dm->dim[d];
                dm->oldr[3*m+d] += dm->dim[d];        
            }
        } else if (dm->boundary[d] == REFLECTING) {
            if (dm->r[3*m+d] > dm->dim[d]){
                dm->r[3*m+d] = 2*dm->dim[d] - dm->r[3*m+d];
                dm->oldr[3*m+d] = 2*dm->dim[d] - dm->oldr[3*m+d];
            } else if (dm->r[3*m+d] < 0) {
                dm->r[3*m+d] = dm->dim[d];
                dm->oldr[3*m+d] = dm->dim[d];        
            }
        }
    }
}

int update_positions(domain *dm, int a, int b)
{
    for (int m = a; m < b; m++){
        calculate_force(dm, m);
        for (int d = 0; d < 3; d++){
            dm->temp[3*m+d] = 2*dm->r[3*m+d] - dm->oldr[3*m+d] + dm->F[3*m+d]*10*dt*dt;
            dm->v[3*m+d] = dm->temp[3*m+d] - dm->oldr[3*m+d]/(2*dt);
        }
    }
    double *temp = dm->oldr;
    dm->oldr = dm->r;
    dm->r = dm->temp;
    dm->temp = temp;
    for (int m = a; m < b; m++)
        check_boundary(dm, m);
    return 0;
}



int checkpoint_count = 0;
int print_checkpoint(char *basepath, domain *dm){
    char path[1028];
    sprintf(path, "%s/chk_%05d.dat", basepath, checkpoint_count++);
    FILE *chkpnt = fopen(path, "w");
    WARN_IF(!chkpnt, "Unable to open checkpoint file [%s]", path);
    for (int m = 0; m < dm->npart; m++){
        for (int d = 0; d < 3; d++)
            fprintf(chkpnt, "%f\t", dm->r[3*m+d]);
        for (int d = 0; d < 3; d++)
            fprintf(chkpnt, "%f\t", dm->F[3*m+d]);
        fprintf(chkpnt, "\n");
    }
    LOG("Wrote to checkpoint file [%s]\t%f", path, t);
    fclose(chkpnt);
    return 0;
}
