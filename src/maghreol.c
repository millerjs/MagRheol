/******************************************************************************
 * FILE        : magrehol.c
 * AUTHOR      : Joshua Miller 
 * DESCRIPTION : Main function for magnetorehological simulations
 ******************************************************************************/

#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>

#include "libjosh/libjosh.h"
#include "libjosh/threadpool.h"
#include "particles.h"
#include "domain.h"

int step = 0;


/* Parameters */
int     npart                =  512;
double  maxt                 =  1000.;
int     checkpoint_interval  =  100;
double  L                    =  50;

void equillibrate(domain *dm)
{
    LOG("Starting equilibration.");
    for (int i = 0; i < 1; i++){
        update_positions(dm, 0, dm->npart);
        update_forces_velocities(dm, 0, dm->npart);
    }
    LOG("Resampling velocities.");
    for (int k = 0; k < dm->npart; k++){
        for (int j = 0; j < 3; j++){
            dm->p[k].v[j] = dm->v0[j] + randomd(-.1, .1);
        }
    }
    LOG("Equilibration complete.");
}

void evolveSystem(domain *dm)
{
    /* Run the simulation */
    update_positions(dm, 0, dm->npart);
    update_forces_velocities(dm, 0, dm->npart);
    if (!(step % checkpoint_interval))
        print_checkpoint("checkpoints", dm);
    t += dt;
    step += 1;
}

void *evolveThreaded(void *args)
{
    thread_t *thread = (thread_t*) args;

    while (!thread->pool->start)
        pthread_yield();

    domain *dm = thread->pool->dm;
    int n = dm->npart;


    while (t < maxt){

        int m = n/thread->pool->size;
        int a = thread->id*m;

        update_positions(dm, a, a+m);
        pthread_barrier_wait(&thread->pool->barrier1);
        update_forces_velocities(dm, a, a+m);
        pthread_barrier_wait(&thread->pool->barrier2);

        if (thread->id == 0){
            LOG("Control thread at === %f === ", t);
            if (!(step % checkpoint_interval))
                print_checkpoint("checkpoints", dm);
            t += dt;
            step += 1;
        }

        pthread_barrier_wait(&thread->pool->barrier3);

    }
    
    return NULL;
}

void evolveController(domain *dm, threadpool_t *pool)
{
    threadpool_start(pool);
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
    domain *dm = domain_new(2*L, L, L);
    setup(dm);

    threadpool_t *pool = threadpool_create(&evolveThreaded, 16);
    pool->dm = dm;

    equillibrate(dm);
    evolveController(dm, pool);


    /* Clean up */
    threadpool_join(pool);
    close_log_file();
    
    return 0;
}
