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
#include "params.h"

#define FCONV 1.0e1

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

double norm(double *a)
{
    double norm = 0;
    for (int d = 0; d < 3; d++){
        norm += a[d]*a[d];
    }
    return sqrt(norm);
}

void normalize(double *a, double m)
{
    double n = norm(a);
    for (int d = 0; d < 3; d++)
        a[d] *= m/n;
}

void allocate_dm(domain *dm, int n)
{
    dm->npart = n;
    LOG("Populating domain [%p] with [%d] particles", dm, n);
    ERROR_IF(!(dm->r = malloc(n*VECSIZE)), 
             "Unable to allocate positions");
    ERROR_IF(!(dm->oldr = malloc(n*VECSIZE)), 
             "Unable to allocate old positions");
    ERROR_IF(!(dm->oldmu = malloc(n*VECSIZE)), 
             "Unable to allocate old positions");
    ERROR_IF(!(dm->v = malloc(n*VECSIZE)), 
             "Unable to allocate old positions");
    ERROR_IF(!(dm->F = malloc(n*VECSIZE)), 
             "Unable to allocate old positions");
    ERROR_IF(!(dm->T = malloc(n*VECSIZE)), 
             "Unable to allocate old positions");
    ERROR_IF(!(dm->temp = malloc(n*VECSIZE)), 
             "Unable to allocate buffer space");
    ERROR_IF(!(dm->tempmu = malloc(n*VECSIZE)), 
             "Unable to allocate buffer space");
    ERROR_IF(!(dm->mu = malloc(n*VECSIZE)), 
             "Unable to allocate buffer space");
    ERROR_IF(!(dm->E = malloc(n*sizeof(double))), 
             "Unable to allocate buffer space: E");
    ERROR_IF(!(dm->magnetic = malloc(n*sizeof(unsigned char))), 
             "Unable to allocate buffer space: mag");
    ERROR_IF(!(dm->pr = malloc(3*sizeof(double))), 
             "Unable to allocate buffer space: mag");
    ERROR_IF(!(dm->oldpr = malloc(3*sizeof(double))), 
             "Unable to allocate buffer space: mag");
    ERROR_IF(!(dm->pF = malloc(3*sizeof(double))), 
             "Unable to allocate buffer space: mag");
}

int domain_populate(domain *dm, int n)
{
    allocate_dm(dm, n);

    double r = pow(dm->npart, 1/3.);
    vec cell = {dm->dim[0]/r, dm->dim[1]/r, dm->dim[2]/r};
    LOG("Cube discretization count: %f", r);
    for (int j = 0; j < 3; j ++)
        LOG("Particle discretization distance: %f", cell[j]);

    dm->pr[0] = dm->dim[0] + 5;
    dm->pr[1] = dm->dim[1] / 2.;
    dm->pr[2] = dm->dim[2] / 2;
    dm->pR = 10;
    for (int d = 0; d < 3; d ++){
        dm->oldpr[d] = dm->pr[d];
        dm->F[d] = 0;
    }

    int m = 0;
    for (double i = 0; i < ceil(r); i++){
        for (int j = 0; j < ceil(r); j++){
            for (int k = 0; k < ceil(r); k++){

                if (randomd(0., 1.) < ratio){

                    dm->magnetic[m] = 1;
                    for (int d = 0; d < 3; d++){
                        dm->mu[3*m+d] = randomd(-1., 2.);
                    }
                    
                    normalize(dm->mu+3*m, MU);
                    for (int d = 0; d < 3; d++){
                        dm->oldmu[3*m+d] = dm->mu[3*m+d];
                    }

                } else {
                    for (int d = 0; d < 3; d++)
                        dm->mu[3*m+d] = 0.;
                }

                dm->oldr[3*m+0] = (i+.5)*cell[0];
                dm->oldr[3*m+1] = (j+.5)*cell[1];
                dm->oldr[3*m+2] = (k+.5)*cell[2];

                for (int d = 0; d < 3; d ++){
                    dm->v[3*m+d] = randomd(-5, 6.) + dm->v0[d];
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
        r[d] = r[d] - (dm->boundary[d] == PERIODIC)*
            dm->dim[d]*rint(r[d]/dm->dim[d]);
        temp += r[d]*r[d];    
    }
    return sqrt(temp);
}


double dot(double *a, double *b)
{
    double res = 0;
    for (int d = 0; d < 3; d++)
        res += a[d]*b[d];
    return res;
}

void cross(double *a, double *b, vec c)
{
    c[0] = a[2]*b[3] - a[3]*b[2];
    c[1] = a[3]*b[1] - a[1]*b[3];
    c[2] = a[1]*b[2] - a[2]*b[1];
}

void force_DipoleDipole(domain *dm, int m)
{
    if (!dm->magnetic[m])
        return;
    vec r = {0,0,0};
    double f;

    for (int j = 0; j < dm->npart; j++){
        double rmj = dist(dm, m, j, r);
        for (int d = 0; d < 3; d ++){
            if (j!=m){
                double mj = dot(dm->mu+3*m, dm->mu+3*j);
                double mr = dot(dm->mu+3*m, r);
                double jr = dot(dm->mu+3*j, r);
                f = mj*r[d]/rmj;
                f -= 5*mr*jr*r[d]/pow(rmj,3);
                f += (mr*dm->mu[3*j+d] + jr*dm->mu[3*m+d])/rmj;
                f *= 3/(4*PI*MU_0*pow(rmj,4))*1e2;
                dm->F[3*m+d] += f;
            }
        }
    }
}

void force_DLVO(domain *dm, int m)
{
    for (int i = 0; i < dm->npart; i++){
        if (i!=m){
            vec res;
            double r = dist(dm, m, i, res);
            double C = 3*MU_S*MU_S/(4*PI*MU_0*pow(2*R, 4));
            double f = C*exp(-40*(r - 2*R));
            dm->E[m] += C*exp(-40*r)*exp(40*2*R)*(40*r+1)/(r*r);
            for (int d = 0; d < 3; d++)
                dm->F[3*m+d] += res[d]/r*f;
        }
    }
}

void force_DLVO_Projectile(domain *dm)
{
    for (int m = 0; m < dm->npart; m++){
        vec res;
        double r = 0;
        for (int d = 0; d < 3; d ++){
            res[d] = dm->r[3*m+d] - dm->pr[d];
            r += res[d]*res[d];    
        }
        double C = 3*MU_S*MU_S/(4*PI*MU_0*pow(2*R, 4));
        double f = C*exp(-40*(r - R - dm->pR));
        for (int d = 0; d < 3; d++)
            dm->pF[d] += res[d]/r*f;
    }
}

void force_DLVO_Proj_Part(domain *dm, int m)
{
    vec res;
    double r = 0;
    for (int d = 0; d < 3; d ++){
        res[d] = dm->r[3*m+d] - dm->pr[d];
        r += res[d]*res[d];    
    }
    double C = 3*MU_S*MU_S/(4*PI*MU_0*pow(2*R, 4));
    double f = C*exp(-40*(r - R - dm->pR));
    for (int d = 0; d < 3; d++){
        dm->F[3*m+d] -= res[d]/r*f*1e3;
    }
}


void force_drag(domain *dm, int m)
{
    for (int d = 0; d < 3; d++)
        dm->F[3*m+d] -= dm->v[3*m+d]*1.5*R;
}

void check_boundary(domain *dm, int m)
{
    for (int d = 0; d < 3; d ++){
        if (dm->boundary[d] == PERIODIC){
            if (dm->r[3*m+d] > dm->dim[d]){
                dm->r[3*m+d] -= dm->dim[d];
                dm->oldr[3*m+d] = dm->oldr[3*m+d] - dm->dim[d];
            } else if (dm->r[3*m+d] < 0) {
                dm->r[3*m+d] += dm->dim[d];
                dm->oldr[3*m+d] = dm->oldr[3*m+d] + dm->dim[d];
            }
        } else if (dm->boundary[d] == REFLECTING) {
            if (dm->r[3*m+d] > dm->dim[d]){
                dm->r[3*m+d] = 2*dm->dim[d] - dm->r[3*m+d];
                dm->oldr[3*m+d] = 2*dm->dim[d] - dm->oldr[3*m+d];
            }
            if (dm->r[3*m+d] < 0) {
                dm->r[3*m+d] *= -1;
                dm->oldr[3*m+d] *= -1;
            }
        }
    }
}

void torque_DipoleDipole(domain *dm, int m)
{
    if (!(dm->mu[3*m] || dm->mu[3*m+1] || dm->mu[3*m+2]))
        return;

    vec r, mxj, mxr;
    double t;

    for (int j = 0; j < dm->npart; j++){
        if (j!=m){
            double rmj = dist(dm, m, j, r);
            cross(dm->mu+3*m, dm->mu+3*j, mxj);
            cross(dm->mu+3*m, r, mxr);
            double jr = dot(dm->mu+3*j, r);
            for (int d = 0; d < 3; d ++){
                t = mxj[d];
                t -= 3/(rmj*rmj)*jr*mxr[d];
                dm->T[3*m+d] += t;
            }
        }
    }
}

void torque_H(domain *dm, int m)
{
    vec M = {0,-H, 0};
    vec t = {0, 0, 0};
    cross(dm->mu+3*m, M, t);
    for (int d = 0; d < 3; d++)
        dm->T[3*m+d] += t[d]*1e3;
}

void calculate_torque(domain *dm, int m)
{
    for (int d = 0; d < 3; d++)
        dm->T[3*m+d] = 0;
    torque_DipoleDipole(dm, m);
    torque_H(dm, m);
    return;
}

int calculate_force(domain *dm, int m)
{
    dm->E[m] = 0;
    for (int d = 0; d < 3; d++)
        dm->F[3*m+d] = 0;
    force_DLVO(dm, m);
    force_DipoleDipole(dm, m);
    force_DLVO_Proj_Part(dm, m);
    force_drag(dm, m);
    return 0;
}

int update_angles(domain *dm, int a, int b)
{
    for (int m = a; m < b; m++){
        if (dm->magnetic[m]){
            calculate_torque(dm, m);
            for (int d = 0; d < 3; d++){
                dm->tempmu[3*m+d] = 2*dm->mu[3*m+d] - dm->oldmu[3*m+d]
                    + dm->T[3*m+d]*10*dt*dt;
            }
            normalize(dm->tempmu+3*m, MU);
        }
    }
    return 0;
}

int update_positions(domain *dm, int a, int b)
{
    /* update the particles */
    for (int m = a; m < b; m++){
        calculate_force(dm, m);
        for (int d = 0; d < 3; d++){
            dm->temp[3*m+d] = 2*dm->r[3*m+d] - dm->oldr[3*m+d]
                + dm->F[3*m+d]*FCONV*dt*dt;
            dm->v[3*m+d] = (dm->temp[3*m+d] - dm->oldr[3*m+d])/(2*dt);
            dm->E[m] += dot(dm->v+3*m, dm->v+3*m)*.5;
        }
    }    

    return 0;
}

int checkpoint_count = 0;
int print_checkpoint(char *basepath, domain *dm){
    char path[1028];
    checkpoint_count++;
    sprintf(path, "%s/chk_%05d.dat", basepath, checkpoint_count);
    FILE *chkpnt = fopen(path, "w");
    sprintf(path, "%s/mag_%05d.dat", basepath, checkpoint_count);
    FILE *magchkpnt = fopen(path, "w");
    WARN_IF(!chkpnt, "Unable to open checkpoint file [%s]", path);
    sprintf(path, "%s/proj_%05d.dat", basepath, checkpoint_count);

    for (int m = 0; m < dm->npart; m++){
        for (int d = 0; d < 3; d++)
            fprintf(chkpnt, "%f\t", dm->r[3*m+d]);
        for (int d = 0; d < 3; d++)
            fprintf(chkpnt, "%f\t", dm->mu[3*m+d]);
        fprintf(chkpnt, "%d\t", dm->magnetic[m]);
        fprintf(chkpnt, "\n");
        if (dm->magnetic[m]){
            for (int d = 0; d < 3; d++)
                fprintf(magchkpnt, "%f\t", dm->r[3*m+d]);
            for (int d = 0; d < 3; d++)
                fprintf(magchkpnt, "%f\t", fabs(dm->mu[3*m+d]));
            fprintf(magchkpnt, "\n");
        }
    }

    FILE *proj = fopen(path, "w");
    WARN_IF(!chkpnt, "Unable to open checkpoint file [%s]", path);
    /* for (int m = 0; m < 5; m ++){ */
        for (int d = 0; d < 3; d++)
            fprintf(proj, "%f\t", dm->pr[d]);
        fprintf(proj, "%f\n", norm(dm->pF));
    /* } */

    fclose(chkpnt);
    fclose(magchkpnt);
    fclose(proj);
    return 0;
}
