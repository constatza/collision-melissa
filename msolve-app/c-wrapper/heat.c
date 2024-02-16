#include <melissa_api.h>

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/wait.h>
#include <math.h>
#include <mpi.h>

#define MAX_LINE_LENGTH 100

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    // create dedicated MPI communicators for each simulation instance
    int rank = -1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int world_size = -1;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    printf("world size: %d\n", world_size);

    MPI_Comm comm_app = MPI_COMM_NULL;
    MPI_Comm_dup(MPI_COMM_WORLD, &comm_app);


    int num_timesteps = 20;
    int num_cells = 15402;

    // initialization
    double* u = (double*) calloc(num_cells, sizeof(double));

    // melissa_init is the first Melissa function to call, and it is called only
    //// once by each process in comm_app. It mainly contacts the server.
    const char* field_name = "displacements";
    melissa_init(field_name, num_cells, comm_app);

    // Call .net application
    const char* command = argv[1];

    char app[MAX_LINE_LENGTH];
    int chars_written = snprintf(app, sizeof(app), "%s %s", command, argv[2]);
    if (chars_written < 0 || chars_written >= sizeof(app)) {
        fprintf(stderr, "Error constructing concatenated string.\n");
        return EXIT_FAILURE;
    }

    printf("Calling .Net App: %s\nrank: %d/%d\n", app, rank+1, world_size);

    FILE* file = popen(app, "r");

    if (file == NULL) {
        printf("Failed to open pipe.\n");
		perror("Executable not found.");
        return 1;
    } else {
        printf(".Net App ran.\n");
    }


    //Read app output from file
    int dof = 0;
    char line[MAX_LINE_LENGTH];
    // main loop
    for (int time_step = 1; time_step < num_timesteps + 1; ++time_step) {
        while(fgets(line, sizeof(line), file) != NULL){
            u[dof] = strtod(line, NULL);
            dof++;
            if (dof == num_cells){
                melissa_send(field_name, u);
                printf("Sent %d dofs to melissa, timestep:%d/%d\n", dof, time_step, num_timesteps);
                dof = 0;
                break;
            }
         }
     }
    //Get your exit code...
    int status = pclose(file);
    if(WIFEXITED(status)) {
    //If you need to do something when the pipe exited, this is the time.
        status=WEXITSTATUS(status);
        printf("process exited with status %d\n", status);
    } else if(WIFSIGNALED(status)) {
			//If you need to add something if the pipe process was terminated, do it here.
			status=WTERMSIG(status);
			printf("process TERMINATED with status %d\n", status);
    } else if(WIFSTOPPED(status)) {
			//If you need to act upon the process stopping, do it here.
			status=WSTOPSIG(status);
			printf("process stopped with status %d\n", status);
    } else {
			printf("process uknown error... spookie action at a distance???\n");
    }

    printf("Rank: %d, finished.\n", rank);
    // melissa_finalize closes the connection with the server.
    // No Melissa function should be called after melissa_finalize.
  	melissa_finalize();
  	printf("Melissa finalized.\n");

    free(u);

    MPI_Finalize();
  	printf("MPI finalized.\n");
}
