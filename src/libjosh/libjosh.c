/******************************************************************************/
/*                                                                            */
/* FILE    : util.c                                                           */
/* AUTHOR  : Joshua Miller                                                    */
/* PROJECT : Project 2                                                        */
/* CLASS   : Parallel Computing - Winter 2014                                 */
/*                                                                            */
/******************************************************************************/

#include <signal.h>
#include <pthread.h>

#include "libjosh.h"

int __timeout__           = 0;
int __catching_exit__     = 0;
int __catching_segfault__ = 0;

char *colors[] = {
    "\033[0;34m",
    "\033[0;32m",
    "\033[0;36m",
    "\033[0;31m",
    "\033[0;35m",
    "\033[0;33m",
    "\033[0;37m",
    "\033[1;30m",
    "\033[1;34m",
    "\033[1;32m",
    "\033[1;36m",
    "\033[1;31m",
    "\033[1;35m",
    "\033[1;33m",
};

void **trash = NULL;
FILE *logfile = NULL;

int uerr = 0;

int open_log_file(char *path)
{
    ERROR_IF(!(logfile = fopen(path, "w")), "Unable to open log file: %s", path);
    LOG("Log file successfully opened");
    LOG("########################################");
    LOG("# Magrheol: created by Josh Miller for #");
    LOG("# computational chemistry Spring 2014. #");
    LOG("########################################");
    return 0;
}

int close_log_file()
{
    LOG("Closing log file...");
    if (logfile) fclose(logfile);
    return 0;
}


/* Obtain a backtrace and print it to stdout. */
void print_trace (void)
{
    void *array[10];
    size_t size;
    char **strings;
    size_t i;
     
    size = backtrace (array, 10);
    strings = backtrace_symbols (array, size);
     
    printf ("Obtained %zd stack frames.\n", size);
     
    for (i = 0; i < size; i++)
        printf ("%s\n", strings[i]);
     
    free (strings);
}

void sig_handler(int signo)
{
    if (signo == SIGSEGV){
        fprintf(stderr, "%sSEGFAULT%s: %s: printing stack frame ...\n",
                __red__, __nrm__, __fname__);
        print_trace();
        if (__catching_segfault__){
            fprintf(stderr, "\treturning to calling function prior to SEGFAULT ... \n");
            uerr = SIGSEGV;
            longjmp(__seg_loc__, 1); 
        } else {
            exit(SIGSEGV);
        }
    }
    if (signo == SIGALRM){
        if (__timeout__){
            fprintf(stderr, "%sTIMEOUT%s: %s\n", __red__, __nrm__, __tname__);
            uerr = SIGALRM;
            longjmp(__to_loc__, 1); 
        }
    }
}

void print_test_result(int failed, const char* test, const char* subtest)
{
    if (failed){
        printf("TEST: [%30s] [%35s] [%s FAILED %s]\n", 
               test, subtest, __red__, __nrm__); 
        return;
    } 
    printf("TEST: [%30s] [%35s] [%s PASSED%s ]\n", 
           test, subtest, __lgr__, __nrm__); 
}


void *__malloc__(size_t size, char *file, const int line)
{
    void *ret = malloc(size);
    if (!ret){
        fprintf(stderr, "%s: %d: error: unable to allocate memory\n", file, line);
        exit(1);
    }
    return ret;
}
