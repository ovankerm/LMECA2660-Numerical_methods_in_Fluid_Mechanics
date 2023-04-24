#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "poisson.h"
#include "struct.h"

problem_struct *create_problem(int nx, int mixer);
data_sim *init_data(int nx, int mixer);

void get_indices(data_sim *data);

void free_problem(problem_struct *problem);
void free_data_sim(data_sim *sim_data);

void problem_to_file(problem_struct *problem);
void T_to_file(problem_struct *problem);
void U_to_file(problem_struct *problem);
void V_to_file(problem_struct *problem);
void P_to_file(problem_struct *problem);
void V_star_to_file(problem_struct *problem);
void phi_to_file(problem_struct *problem);
void H_to_file(problem_struct *problem);
void H_T_to_file(problem_struct *problem);
void omega_to_file(problem_struct *problem);

void iterate(problem_struct *problem);

void compute_BC(data_sim *data);
void add_phi_P(data_sim *data);
void compute_T(data_sim *data, double dt, double t);
void compute_v(data_sim *data, double dt);
void compute_v_star(data_sim *data, double dt, double t);
void compute_H_n(data_sim *data);
void copy_H_n(data_sim *data);

int in_mixer(double x, double y, double t, double *us, double *vs);
double compute_Ts(data_sim *data);
double compute_Tavg(data_sim *data);
double compute_sigT(data_sim *data, double Tavg);

#endif
