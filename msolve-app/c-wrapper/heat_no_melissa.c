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

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/time.h>

#include <math.h>
#include <mpi.h>

// Fortran interfaces
void read_file(int*, int*, double*, double*, double*);

void load(int*, int*, int*, int*, int*);

void init(
    double*, int*, int*, double*, double*, int*, double*, double*, double*);

void filling_A(double*, double*, double*, double*, int*, int*, double*);

void filling_F(
    int*, int*, double*, double*, double*, double*, double*, double*, double*,
    int*, int*, double*, double*, double*);

void conjgrad(
    double*, double*, double*, int*, int*, double*, int*, int*, int*, int*,
    int*, int*, int*);

void finalize(
    double*, double*, int*, int*, int*, int*, double*, double*, int*, int*);


int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    
    printf("Simulation input number: %d\n",argc-1);
    for (int i=0; i+1<argc; ++i){
    	printf("Input parameters values: %s \n",argv[i + 1]);
    }
    printf("\n"); 

    if(argc < 5 || argc > 9) {
        fprintf(
            stderr, "usage: %s <initial temperature> [boundary temperatures]\n",
            argv[0]);
        return EXIT_FAILURE;
    }

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

    // Neighbour ranks
    int next = me + 1;
    int previous = me - 1;

    if(me == num_procs - 1) {
        next = MPI_PROC_NULL;
    }
    if(me == 0) {
        previous = MPI_PROC_NULL;
    }

    struct timeval begin, end; // timer
    gettimeofday(&begin, 0);
    int nx = strtol(argv[1], (char **)NULL, 10); // x axis grid subdivisions
    int ny = strtol(argv[2], (char **)NULL, 10); // y axis grid subdivisions
    double lx = 10.0; // domain size along x axis
    double ly = 10.0; // domain size along y axis
    double d = 1.0; // diffusion coefficient
    int num_time_steps = strtol(argv[3], (char **)NULL, 10);
    double simulation_time = 1;
    double dx = lx / (nx + 1); // x axis step
    double dy = ly / (ny + 1); // y axis step
    double epsilon = 0.0001; // conjugated gradient precision

    // partition work over MPI processes of this simulation
    // i1: first global cell indices assigned to this process
    int i1 = -1;
    // in: last global cell index assigned to the current process
    int in = -1;
    int num_cells_global = nx * ny;
    load(&me, &num_cells_global, &num_procs, &i1, &in);

    // initialization
    int num_cells = in - i1 + 1;
    double* u = malloc(num_cells * sizeof(double));
    double* f = malloc(num_cells * sizeof(double));
    init(u, &i1, &in, &dx, &dy, &nx, &lx, &ly, params);
    // variables for data gathering on rank 0
    double* u_global = malloc(num_cells_global * sizeof(double)); // global solution
    int* displacements = malloc(num_procs * sizeof(int)); // data gathering offset
    int* rcount = malloc(num_procs * sizeof(int)); // data gathering offset
    int g1 = -1, gn = -1;
    for (int i=0; i<num_procs; i++) {
        load(&i, &num_cells_global, &num_procs, &g1, &gn);
        rcount[i] = (gn - g1) + 1;
        displacements[i] = g1 - 1;
    }   
    // initialize the tridiagonal matrix A:
    double a[3] = {0};
    double dt = simulation_time / num_time_steps;
    filling_A(&d, &dx, &dy, &dt, &nx, &ny, a);

    const char field_name[] = "temperature";

    // main loop
    for(int n = 0; n < num_time_steps; ++n) {
        double t = 1.0 * (n + 1) / num_time_steps * simulation_time;
        // fill right-hand side vector F before each iteration
        filling_F(
            &nx, &ny, u, &d, &dx, &dy, &dt, &t, f, &i1, &in, &lx, &ly, params);
        // solve Au = F
        conjgrad(
            a, f, u, &nx, &ny, &epsilon, &i1, &in, &num_procs, &me, &next,
            &previous, &fcomm);

	    // writing result to file
        char file_name[35];
        char *envvar = "MELISSA_SIMU_ID";
        if (!getenv(envvar)){
            // standard offline data generation with one solution per rank and per timestep
            snprintf (file_name, sizeof(file_name), "solution_%d_%d.dat", n, me);
            FILE *outfile = fopen (file_name, "w");
            if (outfile == NULL) {
                printf("Error opening file!\n");
                exit(1);
            }
            for (int i=0; i<num_cells; i++) {
                fprintf (outfile, "%d %f\n", i1+i, u[i]);
            }
            fclose (outfile);
        }
        else {
            // melissa combined data generation with one solution per timestep
            //gather solution on rank 0
            MPI_Gatherv(u, num_cells, MPI_DOUBLE, u_global, rcount, displacements, MPI_DOUBLE, 0, comm_app);
            //write solution on rank 0
            if (me == 0) {
                snprintf (file_name, sizeof(file_name), "solution_%d.dat", n);
                FILE *outfile = fopen (file_name, "w");
                if (outfile == NULL) {
                    printf("Error opening file!\n");
                    exit(1);
                }
                for (int i=0; i<num_cells_global; i++) {
                    fprintf (outfile, "%d %f\n", i + 1, u_global[i]);
                }
                fclose (outfile);
            }
        }

	if (me == 0 && (n+1) % (num_time_steps / 10) == 0){ printf("Advancement: %d /%d \n", n+1, num_time_steps); }
    }
    
    free(u);
    free(f);
    free (u_global);
    free(displacements);
    free(rcount);
    
    gettimeofday(&end, 0);
    long seconds = end.tv_sec - begin.tv_sec;
    long microseconds = end.tv_usec - begin.tv_usec;
    double elapsed = seconds + microseconds*1e-6;
    if (me == 0){ printf("Total elapsed time: %f sec \n",elapsed); }

    MPI_Finalize();
}
