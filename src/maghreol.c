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


int main(int argc, char *argv[])
{
    /* Open a log file */
    open_log_file("magrheol.log");
    
    /* Establish the domain */
    domain *dm = domain_new(10, 10, 10);
    domain_populate(dm, 10);
    domain_set_v0(dm, -1, 0, 0);
    dt = 1;

    /* Run the simulation */
    for (t = 0; t < 10; t += dt){
        update_positions(dm);
        update_forces_velocities(dm);
        print_checkpoint("checkpoints", dm);
    }

    /* Close the log file */
    close_log_file();
    
    return 0;
}
