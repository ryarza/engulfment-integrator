// Copyright 2021 Ricardo Yarza.
//
// This file is part of orbit-integrator
//
// orbit-integrator is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// orbit-integrator is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with orbit-integrator.  If not, see <https://www.gnu.org/licenses/>.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "io.h"
#include "hdf5.h"
#include "utils.h"

/// Removes whitespace from a string. Based on Cholla.
char *trim (char * s)
{
  /* Initialize start, end pointers */
  char *s1 = s, *s2 = &s[strlen (s) - 1];

  /* Trim and delimit right side */
  while ( (isspace (*s2)) && (s2 >= s1) )
    s2--;
  *(s2+1) = '\0';

  /* Trim left side */
  while ( (isspace (*s1)) && (s1 < s2) )
    s1++;

  /* Copy finished string */
  strcpy (s, s1);
  return s;
}

void parse_params (char *param_file, struct parameters * parms)
{
  char *s, buff[MAXLEN-1];
  FILE *fp = fopen (param_file, "r");
  if (fp == NULL){
    perror("Couldn't load parameter file");
    exit(-1);
  }

  /* Read next line */
  while ((s = fgets (buff, sizeof buff, fp)) != NULL)
  {
    /* Skip blank lines and comments */
    if (buff[0] == '\n' || buff[0] == '#') continue;

    /* Parse name/value pair from line */
    char name[MAXLEN], value[MAXLEN];
    s = strtok (buff, "=");
    if (s==NULL)
      continue;
    else
      strncpy (name, s, MAXLEN-1);
    s = strtok (NULL, "=");
    if (s==NULL)
      continue;
    else
      strncpy (value, s, MAXLEN-1);
    trim (value);
    trim (name);
//    printf("%s %s\n", name, value);
    /* Copy into correct entry in parameters struct */
    if (strcmp(name, "Msb")==0)
      parms->Msb = atof(value);
    else if (strcmp(name, "Rsb")==0)
      parms->Rsb = atof(value);
    else if (strcmp(name, "f_hydro")==0)
      parms->f_hydro = atoi(value);
    else if (strcmp(name, "f_buoyancy")==0)
      parms->f_buoyancy = atoi(value);
    else if (strcmp(name, "f_grav")==0)
      parms->f_grav = atoi(value);
    else if (strcmp(name, "r0")==0)
      parms->r0 = atof(value);
    else if (strcmp(name, "mesa_profile_path")==0)
      strncpy (parms->mesa_profile_path, value, MAXLEN);
    else if (strcmp(name, "drag_data_path")==0)
      strncpy (parms->drag_data_path, value, MAXLEN);
    else if (strcmp(name, "use_drag_coefficients")==0)
      parms->use_drag_coefficients = atoi(value);
    else
      printf ("WARNING: %s/%s: Unknown parameter/value pair!\n", name, value);
  }

  /* Close file */
  fclose (fp);

}

int read_stellar_profile(struct parameters P, struct star *S){

  hid_t file_id, dataset_id;
  herr_t status;
  hid_t attr;

  // Open file
  file_id = H5Fopen(P.mesa_profile_path, H5F_ACC_RDONLY, H5P_DEFAULT);

  // Get number of cells
  attr = H5Aopen(file_id, "ncells", H5P_DEFAULT);
  status  = H5Aread(attr, H5T_NATIVE_INT, &S->ncells);
  status = H5Aclose(attr);

  // Allocate variables along profile
  S->r = (double *) malloc(S->ncells * sizeof(double));
  S->rho = (double *) malloc(S->ncells * sizeof(double));
  S->drhodr = (double *) malloc(S->ncells * sizeof(double));
  S->menc = (double *) malloc(S->ncells * sizeof(double));
  S->cs = (double *) malloc(S->ncells * sizeof(double));
  S->hrho = (double *) malloc(S->ncells * sizeof(double));

  // Load each variable
  dataset_id = H5Dopen2(file_id, "/r", H5P_DEFAULT);
  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, S->r);
  if ( status != 0) printf("Couldn't read r from profile\n");

  dataset_id = H5Dopen2(file_id, "/rho", H5P_DEFAULT);
  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, S->rho);

  dataset_id = H5Dopen2(file_id, "/drhodr", H5P_DEFAULT);
  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, S->drhodr);

  dataset_id = H5Dopen2(file_id, "/menc", H5P_DEFAULT);
  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, S->menc);

  dataset_id = H5Dopen2(file_id, "/cs", H5P_DEFAULT);
  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, S->cs);

  dataset_id = H5Dopen2(file_id, "/hrho", H5P_DEFAULT);
  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, S->hrho);

  // Close the dataset
  status = H5Dclose(dataset_id);

  // Close the file
  status = H5Fclose(file_id);

  return 0;

}

int write_scalars(struct units U, struct parameters P){

  FILE *fp;
  fp = fopen("scalars.txt", "w+");
  if ( fp == NULL ){
    perror("Couldn't open file to write scalars");
    exit(-1);
    }

  fprintf(fp, "%22s %22s %22s %22s %22s %22s %22s %22s %22s %100s\n", "l_unit", "t_unit", "m_unit", "Msb", "Rsb", "f_hydro", "f_grav", "f_buoyancy", "use_drag_coefficients", "profile_path");
  fprintf(fp, "%22.15e %22.15e %22.15e %22.15e %22.15e %22i %22i %22i %22i %100s\n", U.l_unit, U.t_unit, U.m_unit, P.Msb, P.Rsb, P.f_hydro, P.f_grav, P.f_buoyancy, P.use_drag_coefficients, P.mesa_profile_path);

  fclose(fp);

  return 0;
}

int write_output(struct units U, struct rhs_terms *rhs, double t, double *y){

  FILE *fp;
  double menc = gsl_spline_eval(spline_menc, y[0], acc_menc);
  double rho = gsl_spline_eval(spline_rho, y[0], acc_rho);
  double hrho = gsl_spline_eval(spline_hrho, y[0], acc_hrho);
  double cs = gsl_spline_eval(spline_cs, y[0], acc_cs);
  double v_theta = y[0] * y[3];
  double cg = 1.;
  double cp = 0.25;
  double mach = fabs(v_theta / cs);

  double q, eps_rho, rp_over_ra;
  compute_flow_parameters(rhs, y, menc, hrho, &q, &eps_rho, &rp_over_ra);

  if ( rhs->use_drag_coefficients == 1 )
  {
    cg = trilinear_interp(rhs->cg_data, rhs->log_q_vals, rhs->log_eps_rho_vals, rhs->log_rp_over_ra_vals, rhs->drag_data_dims, log10(q), log10(eps_rho), log10(rp_over_ra));
    cp = trilinear_interp(rhs->cp_data, rhs->log_q_vals, rhs->log_eps_rho_vals, rhs->log_rp_over_ra_vals, rhs->drag_data_dims, log10(q), log10(eps_rho), log10(rp_over_ra));
  }

  double forceg = 0.;
  double forcep = 0.;
  if ( rhs->f_grav && mach > 0.25 ) forceg = f_grav(rhs, rho, v_theta, q, eps_rho, rp_over_ra);
  if ( rhs->f_hydro) forcep = f_hydro(rhs, rho, v_theta, q, eps_rho, rp_over_ra);
  fp = fopen("output.txt", "a+");
  fprintf(fp, "%22.15e %22.15e %22.15e %22.15e %22.15e %22.15e %22.15e %22.15e %22.15e %22.15e %22.15e %22.15e %22.15e %22.15e\n", t, y[0], y[1], y[2], y[3], menc, rho, q, eps_rho, rp_over_ra, cp, cg, forceg, forcep);

  fclose(fp);

  return 0;
}

void h5safecall(herr_t status){
  if ( status != 0 ) exit(-1);
}
