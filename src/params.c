/******************************************************************************
 * FILE        : magrehol.c
 * AUTHOR      : Joshua Miller 
 * DESCRIPTION : Main function for magnetorehological simulations
 ******************************************************************************/

#define _GNU_SOURCE

#include "params.h"
#include "string.h"

int step = 0;

/* Parameters */
double maxt             =  .2;
int checkpoint_interval =  1;
double X                =  10;
double Y                =  10;
double Z                =  10;
int npart               =  46;
double part_dens        = -1;
double t                =  0;
double dt               =  .01;
double ratio            = .03;
double H                = 1.0;
double MU               = 10.;
double Q = 1;


int parse_config(char *path)
{
    FILE *config = fopen(path, "r");
    ERROR_IF(!config, "Unable to open config file [%s]\n", path);
    LOG("Reading from parameter file [%s]\n", path);
    char *line = NULL;
    size_t n = 0;
    while (getline(&line, &n, config) > 0){
        if (line[0] == '#'){
        } else if (strlen(line) > 1) {
            int chn = strchr(line, '=') - line;
            if (chn > 0){
                char * name = strndup(line, chn);
                char * val = line+chn+1;
                if (!strcmp("maxt", name))
                    LOG("Setting maxt to %f", maxt = atof(val));
                else if (!strcmp("X", name))
                    LOG("Setting X to %f", X = atof(val));
                else if (!strcmp("Y", name))
                    LOG("Setting Y to %f", Y = atof(val));
                else if (!strcmp("Z", name))
                    LOG("Setting Z to %f", Z = atof(val));
                else if (!strcmp("dt", name))
                    LOG("Setting dt to %f", dt = atof(val));
                else if (!strcmp("H", name))
                    LOG("Setting H to %f", H = atof(val));
                else if (!strcmp("MU", name))
                    LOG("Setting MU to %f", MU = atof(val));
                else if (!strcmp("ratio", name))
                    LOG("Setting ratio to %f", ratio = atof(val));
                else if (!strcmp("part_dens", name))
                    LOG("Setting part_dens to %f", part_dens = atof(val));
                else if (!strcmp("npart", name))
                    LOG("Setting naprt to %d", npart = atoi(val));
                else if (!strcmp("checkpoint_interval", name))
                    LOG("Setting checkpoint_interval to %d", 
                        checkpoint_interval = atof(val));

                else
                    fprintf(stderr, "Unknown parameter: %s\n", name);
                
            }
        }
    }
    free(line);

    /* 
     * if (part_dens > 0)
     *     npart = (int) (part_dens*X*Y*Z);
     */

    return 0;
}
