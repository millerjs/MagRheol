/******************************************************************************/
/*                                                                            */
/* FILE    : threads.c                                                        */
/* DESCRIP : header file for threads                                          */
/* AUTHOR  : Joshua Miller                                                    */
/* PROJECT : Project 4                                                        */
/* CLASS   : Parallel Computing - Winter 2014                                 */
/*                                                                            */
/******************************************************************************/


#include <pthread.h>
#include <stdlib.h>

#include "cll.h"
#include "threadpool.h"
#include "libjosh.h"

void *threadpool_join(threadpool_t *pool)
{
    void *ret = Malloc(void*, pool->size);
    for (int i = 0; i < pool->size; i++){
        pthread_join(pool->threads[i].self, ret+i);
    }
    return ret;
}

void threadpool_free(threadpool_t *pool)
{
    free(pool->threads);
    free(pool);
    return;
}

threadpool_t *threadpool_create(loop_t loop, int n)
{
    threadpool_t *pool = Malloc(threadpool_t, n);
    thread_t *threads = Malloc(thread_t, n);

    pool->threads = threads;
    pool->size = n;
    pthread_barrier_init(&pool->barrier1, NULL, n);
    pthread_barrier_init(&pool->barrier2, NULL, n);
    pthread_barrier_init(&pool->barrier3, NULL, n);
    
    threadpool_stop(pool);
    
    for (int i = 0; i < n; i++){
        pool->threads[i].id = i;
        pool->threads[i].pool = pool;

        int rs = pthread_create(&pool->threads[i].self, NULL, loop, &pool->threads[i]);
        ERROR_IF(rs, "unable to create thread\n");
    }

    return pool;

}

void block_on_start(threadpool_t *pool)
{
    while (!(pool->start))
        usleep(1);
}

void threadpool_start(threadpool_t *pool)
{
    pool->stop = 0;
    pool->start = 1;
}

void threadpool_stop(threadpool_t *pool)
{
    pool->stop = 1;
    pool->start = 0;
}
