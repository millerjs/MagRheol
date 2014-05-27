/******************************************************************************
 * FILE        : domain.c
 * AUTHOR      : Joshua Miller 
 * DESCRIPTION : Structures and functions for the domain
 ******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>

#include "domain.h"
#include "libjosh/libjosh.h"

double t = 0;
double dt = 0;

domain *domain_new(double x, double y, double z)
{
    srand(time(NULL));
    domain *new = (domain*)malloc(sizeof(domain));
    new->dim[0] = x;
    new->dim[1] = y;
    new->dim[2] = z;
    memset(new->v0, 0, 3*sizeof(double));
    LOG("Created domain %.2f x %.2f x %.2f", x, y, z);
    new->npart = 0;
    new->p = NULL;
    return new;
}

int domain_populate(domain *dm, int n)
{
    dm->npart = n;
    LOG("Populating domain [%p] with [%d] particles", dm, n);
    ERROR_IF(dm->p, "This domain has already been populated");
    ERROR_IF(!(dm->p = malloc(n*sizeof(particle))), "Unable to allocate particles");
    for (int i = 0; i < n; i++){
        memset(dm->p[i].F, 0, 3*sizeof(double));
        for (int j = 0; j < 3; j ++){
            dm->p[i].r[j] = rand()*dm->dim[j]/RAND_MAX;
            dm->p[i].v[j] = dm->v0[j] + rand()/(double)RAND_MAX;
        }
        dm->p[i].m = 1;
    }
    return n;
}

int domain_set_v0(domain *dm, double x, double y, double z)
{
    LOG("Setting initial velocity to [%.3f, %.3f, %.3f]", x, y, z);
    dm->v0[0] = x;
    dm->v0[1] = y;
    dm->v0[2] = z;
    return 0;
}

int calculate_force(domain *dm, int i)
{
    return 0;
}

int update_positions(domain *dm)
{
    for (int i = 0; i < dm->npart; i++){
        for (int j = 0; j < 3; j++)
            dm->p[i].r[j] += dm->p[i].v[j]*dt + dm->p[i].F[j]*dt*dt;
    }
    return 0;
}

int update_forces_velocities(domain *dm)
{
    for (int i = 0; i < dm->npart; i++){
        for (int j = 0; j < 3; j++){
            double F = dm->p[i].F[j];
            calculate_force(dm, i);
            dm->p[i].v[j] += 1/2.*(dm->p[i].m)*(F + dm->p[i].F[j])*dt;
        }
    }
    return 0;
}


int checkpoint_count = 0;
int print_checkpoint(char *basepath, domain *dm){
    char path[1028];
    sprintf(path, "%s/chk_%05d.dat", basepath, checkpoint_count++);
    FILE *chkpnt = fopen(path, "w");
    LOG("Writing checkpoint to checkpoint file [%s]", path);
    WARN_IF(!chkpnt, "Unable to open checkpoint file [%s]", path);
    for (int i = 0; i < dm->npart; i++){
        fprintf(chkpnt, "%f\t%f\t%f\n", dm->p[i].r[0], 
                dm->p[i].r[1], dm->p[i].r[2]);
    }
    LOG("Wrote to checkpoint file [%s]", path);
    fclose(chkpnt);
    return 0;
}
