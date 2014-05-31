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
#include "params.h"

void equillibrate(domain *dm)
{

    LOG("Starting equilibration.");
    int NEQUIL = 0;
    for (int i = 0; i < NEQUIL; i++){
        update_positions(dm, 0, dm->npart);
    }
    LOG("Equilibration steps %d complete", NEQUIL);
}

void *evolveThreaded(void *args)
{
    thread_t *thread = (thread_t*) args;

    while (!thread->pool->start)
        pthread_yield();

    domain *dm = thread->pool->dm;
    int n = dm->npart;

    int m = n/thread->pool->size;
    int a = thread->id*m;
    
    while (t < maxt){

        /* Update positions and regroup */
        update_positions(dm, a, a+m);
        pthread_barrier_wait(&thread->pool->barrier2);

        /* Let thread0 handle IO and timestep */
        if (thread->id == 0){
            double *temp = dm->oldr;
            dm->oldr = dm->r;
            dm->r = dm->temp;
            dm->temp = temp;

            for (int m = 0; m < dm->npart; m++)
                check_boundary(dm, m);

            if (!(step % checkpoint_interval))
                print_checkpoint("checkpoints", dm);
            t += dt;
            step += 1;
        }

        pthread_barrier_wait(&thread->pool->barrier1);

    }
    
    return NULL;
}

void setup(domain *dm)
{
    /* Establish the domain */
    domain_populate(dm, npart);
    domain_set_v0(dm, -10, 0, 0);
    domain_set_boundary(dm, 0, PERIODIC);
    domain_set_boundary(dm, 1, REFLECTING);
    domain_set_boundary(dm, 2, REFLECTING);
}

int main(int argc, char *argv[])
{
    /* Open a log file */
    open_log_file("magrheol.log");

    if(argc < 2)
        LOG("No config specified. Running with defaults.");
    else
        parse_config(argv[1]);

    domain *dm = domain_new(X,Y,Z);
    setup(dm);

    threadpool_t *pool = threadpool_create(&evolveThreaded, 16);
    pool->dm = dm;

    equillibrate(dm);
    threadpool_start(pool);

    /* Clean up */
    threadpool_join(pool);

    print_checkpoint("checkpoints", dm);

    close_log_file();
    
    return 0;
}
