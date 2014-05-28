/******************************************************************************
 * FILE        : magrehol.c
 * AUTHOR      : Joshua Miller 
 * DESCRIPTION : Main function for magnetorehological simulations
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>

#include "libjosh/libjosh.h"
#include "particles.h"
#include "domain.h"

int step = 0;


/* Parameters */
int     npart                =  100;
double  maxt                 =  1000.;
int     checkpoint_interval  =  100;


void equillibrate(domain *dm)
{
    for (int i = 0; i < 10; i++){
        update_positions(dm);
        update_forces_velocities(dm);
        for (int k = 0; k < dm->npart; k++){
            for (int j = 0; j < 3; j++){
                dm->p[k].v[j] = 0;
            }
        }
    }
}

void evolveSystem(domain *dm)
{
    /* Run the simulation */
    update_positions(dm);
    update_forces_velocities(dm);
    if (!(step % checkpoint_interval))
        print_checkpoint("checkpoints", dm);
    t += dt;
    step += 1;
}

void setup(domain *dm)
{
    /* Establish the domain */
    domain_populate(dm, npart);
    domain_set_v0(dm, -.5, 0, 0);
    domain_set_boundary(dm, 0, PERIODIC);
    dt = .05;
}

int main(int argc, char *argv[])
{

    /* Open a log file */
    open_log_file("magrheol.log");
    domain *dm = domain_new(100, 50, 50);
    setup(dm);

    equillibrate(dm);

    /* Evolve the system */
    while (t < maxt)
        evolveSystem(dm);


    /* Close the log file */
    close_log_file();
    
    return 0;
}
