#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "poisson.h"
#include "struct.h"

problem_struct *create_problem(int Nx);
data_sim *init_data(int nx);

void free_problem(problem_struct *problem);
void free_data_sim(data_sim *sim_data);

// void print_mesh(problem_struct *problem);
void problem_to_file(problem_struct *problem);
void T_to_file(problem_struct *problem);
void U_to_file(problem_struct *problem);
void V_to_file(problem_struct *problem);
// void U_star_to_file(problem_struct *problem);
// void V_star_to_file(problem_struct *problem);
// void P_to_file(problem_struct *problem);

void iterate(problem_struct *problem);

void compute_BC(data_sim *data);
void add_phi_P(data_sim *data);
void compute_T(data_sim *data, double dt);
void compute_v(data_sim *data, double dt);
void compute_v_star(data_sim *data, double dt);
void compute_H_n(data_sim *data);
void copy_H_n(data_sim *data);

// void iterate(problem_struct *problem);


// void test_poisson();

#endif
