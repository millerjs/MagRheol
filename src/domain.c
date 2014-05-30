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
double dt = 0;

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
    new->p = NULL;
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
    ERROR_IF(dm->p, "This domain has already been populated");
    ERROR_IF(!(dm->p = malloc(n*sizeof(particle))), "Unable to allocate particles");

    double r = pow(dm->npart, 1/3.);
    vec cell = {dm->dim[0]/r, dm->dim[1]/r, dm->dim[2]/r};
    LOG("Cube discretization count: %f", r);
    for (int j = 0; j < 3; j ++)
        LOG("Particle discretization distance: %f", cell[j]);

    int m = 0;
    for (int i = 0; i < ceil(r); i++){
        for (int j = 0; j < ceil(r); j++){
            for (int k = 0; k < ceil(r); k++){

                particle_new(dm->p+m);
                dm->p[m].r[0] = (i+.5)*cell[0] + randomd(-cell[0]/2, cell[0]/2);
                dm->p[m].r[1] = (j+.5)*cell[1] + randomd(-cell[1]/2, cell[1]/2);
                dm->p[m].r[2] = (k+.5)*cell[2] + randomd(-cell[2]/2, cell[2]/2);
                
                for (int d = 0; d < 3; d++)
                    dm->p[m].v[d] = dm->v0[d] + randomd(-1., 1.);
                dm->p[m].sigma = 3.4; // Angstrom
                dm->p[m].m = 1;
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

double dist(domain *dm, particle *p1, particle *p2, vec r)
{
    double temp = 0;
    for (int i = 0; i < 3; i ++){
        r[i] = p1->r[i] - p2->r[i];
        r[i] = r[i] - (dm->boundary[i] == PERIODIC)*rint(r[i]/dm->dim[i]);
        temp += r[i]*r[i];    
    }
    return sqrt(temp);
}

void force_Drag(domain *dm, particle *p, vec res)
{
    for (int j = 1; j < 3; j++){
        if (p->r[j] < p->sigma || p->r[j] > dm->dim[j] - p->sigma){
            /* res[j] = -p->m*p->v[j]; */
            p->v[0] = 0;
            p->F[0] = 0;
        }
    }
}

void force_LJ(domain *dm, particle *p1, particle *p2, vec res)
{
    if (p1==p2){
        scale(res, 0);
        return;
    }
    double d = dist(dm, p1, p2, res);
    norm(res);
    double t6 = pow(p1->sigma/d, 6);
    double t12 = t6*t6;
    double f = 4*EPS*(12/d*t12 - 6/d*t6);
    double maxf = 1e6;
    WARN_IF(f > maxf, "TIME-STEP TOO LARGE.");
    scale(res, MIN(f, maxf));
}

int calculate_force(domain *dm, int i)
{
    vec F = {0,0,0};
    vec temp;
    for (int j = 0; j < dm->npart; j++){
        force_LJ(dm, dm->p+j, dm->p+i, temp);
        add(temp, F, F);
    }
    force_Drag(dm, dm->p+i, temp);
    memcpy(dm->p[i].F, F, 3*sizeof(double));
    return 0;
}

int update_positions(domain *dm, int a, int b)
{
    for (int i = a; i < b; i++){
        for (int j = 0; j < 3; j++){
            double dx =  dm->p[i].v[j]*dt + dm->p[i].F[j]*dt*dt/(2*dm->p[i].m);
            if (dm->p[i].r[j]+dx < 0 || dm->p[i].r[j]+dx > dm->dim[j]){
                if (REFLECTING == dm->boundary[j]){
                    dm->p[i].v[j] *= -1;
                } else if (PERIODIC == dm->boundary[j]) {
                    dm->p[i].r[j] += dx;
                    double diff = ((long)(dm->p[i].r[j]/dm->dim[j])+1)*dm->dim[j];
                    if (dm->p[i].r[j] < 0)
                        dm->p[i].r[j] += diff;
                    else 
                        dm->p[i].r[j] -= diff;
                }
            } else {
                dm->p[i].r[j] += dx;
            }
        }
    }
    return 0;
}

int update_forces_velocities(domain *dm, int a, int b)
{
    for (int i = a; i < b; i++){
        for (int j = 0; j < 3; j++){
            double F = dm->p[i].F[j];
            calculate_force(dm, i);
            dm->p[i].v[j] += 1/2.*(1/dm->p[i].m)*(F + dm->p[i].F[j])*dt;
        }
    }
    return 0;
}


int checkpoint_count = 0;
int print_checkpoint(char *basepath, domain *dm){
    char path[1028];
    sprintf(path, "%s/chk_%05d.dat", basepath, checkpoint_count++);
    FILE *chkpnt = fopen(path, "w");
    WARN_IF(!chkpnt, "Unable to open checkpoint file [%s]", path);
    for (int i = 0; i < dm->npart; i++){
        for (int j = 0; j < 3; j++)
            fprintf(chkpnt, "%f\t", dm->p[i].r[j]);
        for (int j = 0; j < 3; j++)
            fprintf(chkpnt, "%f\t", dm->p[i].v[j]);
        fprintf(chkpnt, "\n");
    }
    LOG("Wrote to checkpoint file [%s]\t%f", path, t);
    fclose(chkpnt);
    return 0;
}
