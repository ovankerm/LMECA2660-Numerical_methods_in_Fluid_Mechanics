#include <stdlib.h>
#include <stdio.h>
#include "thomas.h"

int main(int argc, char **argv){
    printf("Hello World\n");
    double x[2] = {1, 2};
    double q[2] = {1, 2};
    solve_Ac_thomas(1, 1, 1, 1, x, q);
}