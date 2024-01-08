
#include <stdio.h>
#include <unistd.h>

#define MAXLENGTH 20


int main(int argc, char** argv)
{
    if (chdir("../msolve-app/app/")==0){
        printf("Directory changed successfully.\n");
    } else {
        perror("chdir");
    }

    const char app[] = "./BumperCollisionSimulation 2000000";
    FILE* process = popen(app, "r");


}