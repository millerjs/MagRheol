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
#include "math.h"

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
    int b = a + m;
    if (thread->id == thread->pool->size -1)
        b = n;

    while (t < maxt){

        /* Update positions and regroup */
        update_positions(dm, a, b);
        update_angles(dm, a, b);
        pthread_barrier_wait(&thread->pool->barrier2);

        /* Let thread0 handle IO and timestep */
        if (thread->id == 0){

            double *temp = dm->oldr;
            dm->oldr = dm->r;
            dm->r = dm->temp;
            dm->temp = temp;

            temp = dm->oldmu;
            dm->oldmu = dm->mu;
            dm->mu = dm->temp;
            dm->temp = temp;
            
            vec mu = {0,0,0};
            vec v = {0,0,0};

            int nmagnetic = 0;
            for (int m = 0; m < dm->npart; m ++){
                check_boundary(dm, m);
                normalize(dm->mu+3*m, MU);
                if (dm->magnetic[m]){
                    nmagnetic ++;
                    for (int d = 0; d < 3; d++){
                        v[d] += dm->v[3*m+d];
                        mu[d] += dm->mu[3*m+d];
                    }                    
                }
            }

            if (!(step % checkpoint_interval)){
                print_checkpoint("checkpoints", dm);
                fprintf(stdout, "%04d  %5.3f\t%5.3f  %5.3f  %5.3f \t %5.3f  %5.3f  %5.3f\n",
                        step, t, 
                        mu[0]/nmagnetic, mu[1]/nmagnetic, mu[2]/nmagnetic,
                        v[0]/dm->npart, v[1]/dm->npart, v[2]/dm->npart);
            }
            
            t += dt;
            step += 1;
        }

        if (thread->id == thread->pool->size-1){
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
