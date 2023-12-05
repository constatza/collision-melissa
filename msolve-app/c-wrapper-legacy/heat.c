#include <melissa/api.h>

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <math.h>
#include <mpi.h>

#define MAX_LINE_LENGTH 1000
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

    if(argc < 2 || argc > 6) {
        fprintf(
            stderr, "usage: %s <initial temperature> [boundary temperatures]\n",
            argv[0]);
        return EXIT_FAILURE;
    }

    // The initial temperature is stored in params[0]
    // The four next optional parameters are the boundary temperatures
    double params[5] = {0};

    for(int i = 0; i + 1 < argc; ++i) {
        params[i] = strtod(argv[i + 1], NULL);

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

    int nx = 3; // x axis grid subdivisions
    int ny = 5134; // y axis grid subdivisions
    double lx = 10.0; // domain size along x axis
    double ly = 10.0; // domain size along y axis
    double d = 1.0; // diffusion coefficient
    int num_time_steps = 10;
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
	// int num_cells = 3;
    double* u = malloc(num_cells * sizeof(double));
    double* f = malloc(num_cells * sizeof(double));
    init(u, &i1, &in, &dx, &dy, &nx, &lx, &ly, params);
    // initialize the tridiagonal matrix A:
    double a[3] = {0};
    double dt = simulation_time / num_time_steps;
    filling_A(&d, &dx, &dy, &dt, &nx, &ny, a);

    // melissa_init is the first Melissa function to call, and it is called only
    // once by each process in comm_app. It mainly contacts the server.
    const char field_name[] = "temperature";
    melissa_init(field_name, num_cells, comm_app);

    // main loop
     //for(int n = 0; n < num_time_steps; ++n) {
        //printf("\n time step: %d\n", n);
        //Call .net application
        char app[] = "../../BumperCollitionSimulation/BumperCollitionSimulation/bin/Debug/net6.0/BumperCollitionSimulation ";
		//char app[] = "../../carbumperExample/bin/Debug/net6.0/carbumperExample ";
        strcat(app, argv[1]);
        strcat(app," 2>&1");
        
         printf("Calling .Net App:\n%s\n", app);
         FILE* file = popen(app, "r");
         printf(">.Net App returned!\n");

        //Read app output from file
		int timeStep = 0;
        int dof = 0;
        char line[MAX_LINE_LENGTH];
         while(fgets(line, MAX_LINE_LENGTH, file)){
             u[dof] = strtod(line, NULL);
			 printf("MSolve solution: %d\t%s\n", timeStep, line);
             dof++;
             if (dof >= num_cells_global)
             {
                printf("Sending %d dofs to melissa...\n", dof);
                melissa_send(field_name, u);
                dof = 0;
				timeStep++ ;
				if (timeStep>num_time_steps)
					break ;
            }
         }
		 //melissa_send(field_name, u);
		 printf("closing file\n");
         //pclose(file); uncomment this later
         //Get your exit code...
		int status=pclose(file);
		if(WIFEXITED(status)) {
		//If you need to do something when the pipe exited, this is the time.
			status=WEXITSTATUS(status);
			printf("process exited with status %d", status);
		}
		else if(WIFSIGNALED(status)) {
			//If you need to add something if the pipe process was terminated, do it here.
			status=WTERMSIG(status);
			printf("process TERMINATED with status %d", status);
		}
		else if(WIFSTOPPED(status)) {
			//If you need to act upon the process stopping, do it here.
			status=WSTOPSIG(status);
			printf("process stopped with status %d", status);
		}
		else {
			printf("something else happened...spookie action at a distance???");
		}
     //}

    // melissa_finalize closes the connection with the server.
    // No Melissa function should be called after melissa_finalize.
    printf("finished the analyses\n");
	melissa_finalize();
	printf("Melissa finalized\n");
    free(u);
    free(f);

    MPI_Finalize();
	printf("MPI finalized\n");
}
