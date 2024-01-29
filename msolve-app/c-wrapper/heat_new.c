// Copyright (c) 2017, Institut National de Recherche en Informatique et en Automatique (Inria)
//               2017, EDF
//               2018-2021, Institut National de Recherche en Informatique et en Automatique (Inria)
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// * Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the following disclaimer.
//
// * Redistributions in binary form must reproduce the above copyright
//   notice, this list of conditions and the following disclaimer in the
//   documentation and/or other materials provided with the distribution.
//
// * Neither the name of the copyright holder nor the names of its
//   contributors may be used to endorse or promote products derived from
//   this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
// IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
// TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
// PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
// TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include <melissa_api.h>

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/time.h>
#include <time.h>

#include <math.h>
#include <mpi.h>

#define MAX_LINE_LENGTH 100

// Time indication
void get_current_time(int time_step, int max_time_step) {
    char buffer[26];
    int millisec;
    struct tm* tm_info;
    struct timeval tv;

    gettimeofday(&tv, NULL);

    millisec = lrint(tv.tv_usec/1000.0); // Round to nearest millisec
    if (millisec>=1000) { // Allow for rounding up to nearest second
        millisec -=1000;
        tv.tv_sec++;
    }

    tm_info = localtime(&tv.tv_sec);

    strftime(buffer, 26, "%Y:%m:%d %H:%M:%S", tm_info);

    char output[37];
    sprintf(output, "%s.%03d\n", buffer, millisec);
    printf("Advancement: %d/%d, %s.%03d\n", time_step, max_time_step, buffer, millisec);
    fflush(stdout);
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    
    printf("Simulation input number: %d\n",argc-1);
    for (int i=0; i+1<argc; ++i){
    	printf("Input parameters values: %s \n",argv[i + 1]);
    }
    printf("\n"); 


    // The initial temperature is stored in params[0]
    // The four next optional parameters are the boundary temperatures
    double params[5] = {0};

    for(int i = 0; i + 1 < argc - 3; ++i) {
        params[i] = strtod(argv[i + 4], NULL);
        if(!isfinite(params[i])) {
            fprintf(
                stderr,
                "argument %d '%s' is not a finite floating-point number\n",
                i + 1, argv[i + 1]);
            return EXIT_FAILURE;
        }
    }

    // create dedicated MPI communicators for each simulation instance
    int me = -1;
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    int world_size = -1;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    MPI_Comm comm_app = MPI_COMM_NULL;
    if(world_size == 1)
    {
        MPI_Comm_dup(MPI_COMM_WORLD, &comm_app);
    }
    else
    {
        int* p_appnum = NULL;
        int statinfo = -1;
        MPI_Comm_get_attr(MPI_COMM_WORLD, MPI_APPNUM, &p_appnum, &statinfo);
        assert(statinfo);
        assert(p_appnum);
        MPI_Comm_split(MPI_COMM_WORLD, *p_appnum, me, &comm_app);
    }

    MPI_Comm_rank(comm_app, &me);
    int num_procs = -1;
    MPI_Comm_size(comm_app, &num_procs);
    int fcomm = MPI_Comm_c2f(comm_app);


    struct timeval begin, end; // timer
    gettimeofday(&begin, 0);
    int nx = 3543;
    int ny = 3;
    int num_timesteps = 10;

    // partition work over MPI processes of this simulation
    // i1: first global cell indices assigned to this process
    int i1 = -1;
    // in: last global cell index assigned to the current process
    int in = -1;
    int num_cells_global = nx * ny;

    // initialization
    int num_cells = in - i1 + 1;
    double* u = malloc(num_cells * sizeof(double));

    // melissa_init is the first Melissa function to call, and it is called only
    // once by each process in comm_app. It mainly contacts the server.
    const char field_name[] = "temperature";
    melissa_init(field_name, num_cells, comm_app);
    
    
    const char* command = "/home/catzarakis/collision-melissa/msolve-app/app/BumperCollisionSimulation";

    char app[MAX_LINE_LENGTH];
    int chars_written = snprintf(app, sizeof(app), "%s %s", command, argv[1]);
    if (chars_written < 0 || chars_written >= sizeof(app)) {
        fprintf(stderr, "Error constructing concatenated string.\n");
        return EXIT_FAILURE;
    printf("Calling .Net App:\n%s\n", app);
    FILE* file = popen(app, "r");
    if (file == NULL) {
        printf("Failed to open pipe.\n");
        return -1;
    }
     printf(">.Net App returned!\n");

    // main loop
    for(int timestep = 0; timestep < num_timesteps; ++timestep) {

        //Read app output from file
	    int timeStep = 0;
        int dof = 0;
        char line[MAX_LINE_LENGTH];
        // send u to the Melissa server
        melissa_send(field_name, u);

        if (me == 0 && (timestep+1)) {
            get_current_time(timestep + 1, num_timesteps);
        }
    }

    // melissa_finalize closes the connection with the server.
    // No Melissa function should be called after melissa_finalize.
    melissa_finalize();
    
    free(u);
    
    gettimeofday(&end, 0);
    long seconds = end.tv_sec - begin.tv_sec;
    long microseconds = end.tv_usec - begin.tv_usec;
    double elapsed = seconds + microseconds*1e-6;
    if (me == 0){ printf("Total elapsed time: %f sec \n",elapsed); }

    MPI_Finalize();
}
