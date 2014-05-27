/******************************************************************************/
/*                                                                            */
/* FILE    : threads.h                                                        */
/* DESCRIP : header file for threads                                          */
/* AUTHOR  : Joshua Miller                                                    */
/* PROJECT : Project 3                                                        */
/* CLASS   : Parallel Computing - Winter 2014                                 */
/*                                                                            */
/******************************************************************************/


#ifndef THREADS_H_
#define THREADS_H_

#include <pthread.h>

typedef void *(*loop_t)(void*);

typedef struct threadpool_t threadpool_t;

typedef struct thread_t {
    pthread_t self;
    int id;
    threadpool_t *pool;
} thread_t;

typedef struct threadpool_t{
    thread_t *threads;
    int size;
    int start;
    int stop;
} threadpool_t;

/* Creates and spawns a pool of n threads */
threadpool_t *threadpool_create(loop_t loop, int n);

/* Joins all threads in array threads */
void * threadpool_join(threadpool_t *threads);

/* Frees the thread pool and it's threads */
void threadpool_free(threadpool_t *pool);

/* Function to allow each thread to call it in order to wait for the
   thread's start flag */
void block_on_start();

/* Signal to all threads waiting on block_on_start() to continue */
void threadpool_start(threadpool_t *pool);

/* Signal to all threads looping on (!pool->stop) to stop */
void threadpool_stop(threadpool_t *pool);


#endif

