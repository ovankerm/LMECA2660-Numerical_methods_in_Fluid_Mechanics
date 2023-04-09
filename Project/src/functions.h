#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "cfd.h"

problem_struct *create_problem(int Nx);
data_sim *init_data(int nx);

void free_problem(problem_struct *problem);
void free_data_sim(data_sim *sim_data);

// void print_mesh(problem_struct *problem);
// void problem_to_file(problem_struct *problem);
// void T_to_file(problem_struct *problem);
// void U_to_file(problem_struct *problem);
// void V_to_file(problem_struct *problem, const char* name);
// void U_star_to_file(problem_struct *problem);
// void V_star_to_file(problem_struct *problem);
// void P_to_file(problem_struct *problem);

// void compute_BC(problem_struct *problem);
// void compute_v_star(problem_struct *problem);
// void compute_v(problem_struct *problem);
// void compute_T(problem_struct *problem);
// void compute_H_n(problem_struct *problem);
// void add_phi_P(problem_struct *problem);
// void copy_H_n(problem_struct *problem);

// void iterate(problem_struct *problem);


// void test_poisson();

#endif
