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
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>
#include "io.h"
#include "star.h"
#include "model.h"
#include "units.h"
#include "utils.h"
#include "global.h"

#ifndef M_PI
// Make sure M_PI is defined
#include <gsl/gsl_math.h>
#endif

int main(int argc, char *argv[]){

  // Read arguments
  char *param_file;
  if (argc != 2){
    printf("Usage: %s <parameter file>\n", argv[0]);
    exit(-1);
  }
  else param_file = argv[1];

  struct parameters P;
  struct star S;
  struct units U;
  struct rhs_terms rhs;

  printf("Reading input file %s...\n", param_file);
  parse_params (param_file, &P);
  printf("Parameters:\n");
  printf("  Msb                   = %.10e g\n", P.Msb);
  printf("  Rsb                   = %.10e cm\n", P.Rsb);
  printf("  r0                    = %.10e L_unit\n", P.r0);
  printf("  f_hydro               = %d\n", P.f_hydro);
  printf("  f_grav                = %d\n", P.f_grav);
  printf("  f_buoyancy            = %d\n", P.f_buoyancy);
  printf("  use_drag_coefficients = %d\n", P.f_buoyancy);
  read_stellar_profile(P, &S);
  set_units(&P, &S, &U);
  write_scalars(U, P);
  init_stellar_profile_splines(S);
  printf("Loaded stellar profile:\n");
  printf("  Number of cells: %d\n", S.ncells);
  printf("  Star mass: %.2e MSun\n", S.menc[S.ncells-1] * U.m_unit / U.MSUN_CGS);
  printf("  Star radius: %.2e RSun\n", S.r[S.ncells-1] * U.l_unit / U.RSUN_CGS );

  // Initialize model (setup RHS, GSL system and driver)
  init_model_params(P, S, U, &rhs);
  gsl_odeiv2_system sys = {func, NULL, 4, &rhs};
  gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd, 1.e-10, 1.e-13, 0.);
//  gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd, 1.e-10, 0., 1.e-8);

//  gsl_odeiv2_system sys = {func, jac, 4, &rhs};
//  gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_bsimp, 1.e-3, 0., 1.e-13);
//  gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_bsimp, 1.e-10, 0., 1.e-15);
//  gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_bsimp, 1.e-10, 0., 5.e-7);
//  gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_bsimp, 1.e-10, 0., 1.e-4);

  if ( P.use_drag_coefficients == 1 ){
    printf("Using drag coefficients:\n");
    printf("  %i log(         q) values from %.2e to %.2e\n", rhs.drag_data_dims[0], rhs.log_q_vals[0], rhs.log_q_vals[-1]);
    printf("  %i log(   eps_rho) values from %.2e to %.2e\n", rhs.drag_data_dims[1], rhs.log_eps_rho_vals[0], rhs.log_eps_rho_vals[-1]);
    printf("  %i log(rp_over_ra) values from %.2e to %.2e\n", rhs.drag_data_dims[2], rhs.log_rp_over_ra_vals[0], rhs.log_rp_over_ra_vals[-1]);
  }

  // Initial conditions: Keplerian rotation at r0
  double menc = gsl_spline_eval(spline_menc, P.r0, acc_menc);
  double thetadot = sqrt(U.G * menc / P.r0) / P.r0;
  double y[4] = {P.r0, 0.0, 0., thetadot};

  double t = 0.0;
  double ti = 0.;
  double period, vr, vtheta, v, a;

  FILE *fp;
  fp = fopen("output.txt", "w+");
  // Write headers
  fprintf(fp, "%22s %22s %22s %22s %22s %22s %22s %22s %22s %22s %22s %22s %22s %22s\n", "time", "r", "theta", "rdot", "thetadot", "menc", "rho", "q", "eps_rho", "rp_over_ra", "cp", "cg", "f_grav", "f_ram");
  fclose(fp);

  int status = GSL_SUCCESS;
  printf("Starting integration...\n");
  while ( status == GSL_SUCCESS && y[0] > P.Rsb ){

    menc = gsl_spline_eval(spline_menc, y[0], acc_menc);
    
    vr = y[2];
    vtheta = y[0] * y[3];
    v = sqrt(vr * vr + vtheta * vtheta);
    a = U.G * menc * y[0] / ( 2 * U.G * menc - y[0] * v * v );
    
    period = sqrt(pow(a, 3.) / U.G / menc);
    write_output(U, &rhs, t, y);

    ti += 1.e-2 * fmin( fabs( y[0] / y[2] ), period);

    status = gsl_odeiv2_driver_apply (d, &t, ti, y);

  }

  printf ("Integration stopped at r = %.2e L_unit, with return value %i\n", y[0], status);

  gsl_odeiv2_driver_free (d);

  return 0;
}
