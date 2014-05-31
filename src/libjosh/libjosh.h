/******************************************************************************/
/*                                                                            */
/* FILE        : util.h                                                       */
/* AUTHOR      : Joshua Miller                                                */
/* DESCRIPTION : General utility functions and macros                         */
/*                                                                            */
/******************************************************************************/

#ifndef _UTIL_H
#define _UTIL_H

#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <unistd.h>
#include <stdarg.h>
#include <string.h>
#include <pthread.h>
#include <setjmp.h>
#include <execinfo.h>

#define MAX_DESCRIPTOR 1028
#define RET_FAILURE -1
#define RET_SUCCESS 0

#define __red__ "\033[1;31m"
#define __lgr__ "\033[1;32m"
#define __gry__ "\033[1;30m"
#define __nrm__ "\033[0m"

extern FILE* logfile;

#define STR(_x)   #_x
#define XSTR(_x)  STR(_x)

#ifdef DEBUG
#undef DEBUG
#endif

#define ERR_NOMEM "memory does not exist"

jmp_buf __ex_loc__;
jmp_buf __seg_loc__;
jmp_buf __to_loc__;

int __timeout__;
int __catching_exit__;
int __catching_segfault__;

char __fname__[MAX_DESCRIPTOR];
char __tname__[MAX_DESCRIPTOR];

extern char *colors[];
int uerr;

#define MAX(a,b) \
    ({ __typeof__ (a) _a = (a); \
    __typeof__ (b) _b = (b); \
    _a > _b ? _a : _b; })

#define MIN(a,b) \
    ({ __typeof__ (a) _a = (a); \
    __typeof__ (b) _b = (b); \
    _a > _b ? _b : _a; })

#define INFO_DUMP(msg)                                  \
    do {                                                \
        fprintf(logfile, "%s:%d: %s: %s: ",              \
                __FILE__,                               \
                __LINE__,                               \
                __func__,                               \
                msg);                                   \
    } while(0)                                          


#define LOG(fmt, ...)                                   \
    do {                                                \
        if (!logfile) break;                            \
        char buf[51];                                   \
        snprintf(buf, 30, "[%d][%d][%s] ",              \
                 ((unsigned int)pthread_self())%63,     \
                 __LINE__,                              \
                 __func__);                             \
        fprintf(logfile, "%s%30s | %s",                 \
                __gry__, buf, __nrm__);                 \
        fprintf(logfile, fmt, ##__VA_ARGS__);           \
        fprintf(logfile, "%s\n", __nrm__);              \
        fflush(logfile);                                \
    } while (0)

int open_log_file(char* path);
int close_log_file();

#ifndef _DEBUG
#define DEBUG(fmt, ...)
#else
#define DEBUG(fmt, ...)                                 \
    do {                                                \
        char buf[51];                                   \
        int idx = ((unsigned int)pthread_self())%14;    \
        snprintf(buf, 30, "[%d][%d][%s] ",              \
                 ((unsigned int)pthread_self())%63,     \
                 __LINE__,                              \
                 __func__);                             \
        fprintf(stderr, "%s%30s | %s",                  \
                __gry__, buf, colors[idx]);             \
        fprintf(stderr, fmt, ##__VA_ARGS__);            \
        fprintf(stderr, "%s\n", __nrm__);               \
    } while (0)
#endif

#define ERROR(fmt, ...)                         \
    do {                                        \
        INFO_DUMP(__red__ "error" __nrm__);     \
        fprintf(logfile, fmt, ##__VA_ARGS__);    \
        fprintf(logfile, "\n");                  \
        fprintf(stderr, fmt, ##__VA_ARGS__);    \
        fprintf(stderr, "\n");                  \
        exit(EXIT_FAILURE);                     \
    } while(0)


#define WARN(fmt, ...)                          \
    do {                                        \
        INFO_DUMP("warning");                   \
        fprintf(logfile, fmt, ##__VA_ARGS__);    \
        fprintf(logfile, "\n");                  \
    } while(0)


#define WARN_IF(f, fmt, ...)                    \
    do {                                        \
        if ((f) != 0) {                         \
            WARN(fmt, ##__VA_ARGS__);           \
        }                                       \
    } while(0)


#define ERROR_IF(f, fmt, ...)                   \
    do {                                        \
        if ((f) != 0) {                         \
            ERROR(fmt, ##__VA_ARGS__);          \
        }                                       \
    } while(0)

#define SUCCESS(fmt, ...)                       \
    do {                                        \
        INFO_DUMP(__lgr__ "success" __nrm__);   \
        fprintf(logfile, fmt, ##__VA_ARGS__);    \
        fprintf(logfile, "\n");                  \
    } while(0)

#define SUCCESS_IF(condition, fmt, ...)         \
    do {                                        \
        if ((condition) != 0){                  \
            SUCCESS(fmt, ##__VA_ARGS__);        \
        }                                       \
    } while(0)


#define TRY do{ jmp_buf ex_buf__; if( !setjmp(ex_buf__) ){
#define CATCH } else {
#define ETRY } }while(0)
#define THROW longjmp(ex_buf__, 1)

#undef exit

#define exit(i)                                                 \
    do                                                          \
        {                                                       \
            if (__catching_exit__){                             \
                INFO_DUMP(__red__ "caught exit" __nrm__);       \
                fprintf(logfile, "with status %d\n", i);         \
                uerr = i;                                       \
                longjmp(__ex_loc__, 1);                         \
            } else {                                            \
                (exit)(i);                                      \
            }                                                   \
        }                                                       \
    while (0)                                                   \


#define CATCH_EXIT(f)                           \
    do {                                        \
        __catching_exit__ = 1;                  \
        if (!setjmp(__ex_loc__))                \
            (f);                                \
        __catching_exit__ = 0;                  \
    } while (0)                                 \

#define CATCH_ALL(f)                                    \
    do {                                                \
        signal(SIGSEGV, sig_handler);                   \
        snprintf(__fname__, MAX_DESCRIPTOR, STR(f));    \
        __catching_exit__ = 1;                          \
        __catching_segfault__ = 1;                      \
        uerr = 0;                                       \
        if (!setjmp(__seg_loc__) &&                     \
            !setjmp(__ex_loc__))                        \
            (f);                                        \
        if (uerr){                                      \
            INFO_DUMP("returned from catch_all");       \
            fprintf(logfile, " %s\n",                    \
                    STR(f));                            \
        }                                               \
        __catching_exit__ = 0;                          \
    } while (0)                                         \

#define TIMEOUT(d, f)                                   \
    do {                                                \
        signal(SIGALRM, sig_handler);                   \
        snprintf(__tname__, MAX_DESCRIPTOR, STR(f));    \
        __timeout__ = 1;                                \
        alarm(d);                                       \
        uerr = 0;                                       \
        if (!setjmp(__to_loc__))                        \
            (f);                                        \
        if (uerr){                                      \
            INFO_DUMP("returned from timeout");         \
            fprintf(logfile, " %s\n",                    \
                    STR(f));                            \
        }                                               \
        __timeout__ = 0;                                \
    } while (0)                                         \


#define RANDOM(min, max)                                       \
    (min + rand() / (RAND_MAX / (max - min)))                  \


#define DINT(v, ...)                            \
    do {                                        \
        INFO_DUMP("int");                       \
        fprintf(logfile, STR(v)": %d", v);       \
        fprintf(logfile, "\n");                  \
    } while(0)

void *_malloc_(size_t size);
void garbageCollect();

#define Malloc(typet, count)                                    \
    __malloc__(count*sizeof(typet), __FILE__, __LINE__);


#define TEST(f, sub)                                    \
    do {                                                \
        uerr = 0;                                       \
        print_test_result((f || uerr), sub, STR(f));    \
    } while (0)                                         \


#define RETURN_IF(f, err)                       \
    do {                                        \
        if ((f) != 0) {                         \
            return err;                         \
        }                                       \
    } while(0)



/* Functions */

void sig_handler(int signo);
void print_test_result(int failed, const char* test, const char* subtest);
void *__malloc__(size_t size, char *file, const int line);

    
#endif



