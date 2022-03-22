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

  // Initial conditions: Keplerian rotation at r0
  double ecc = 0.8;
  double a = P.r0 / ( ecc + 1. );
  double menc = gsl_spline_eval(spline_menc, P.r0, acc_menc);
  double rho = gsl_spline_eval(spline_rho, P.r0, acc_rho);
  double thetadot = sqrt(U.G * menc * (1 - ecc) / P.r0) / P.r0;
  double y[4] = {P.r0, 0.0, 0., thetadot};
  double v_tot_sq = y[2] * y[2] + pow(y[0] * y[3], 2.);
  double energy_0 = - U.G * menc / y[0] + 0.5 * v_tot_sq;

  double t = 0.0;
  double period = sqrt( 4. * pow(M_PI, 2.) * pow(a, 3.) / U.G / menc);
  double nperiods = 1.e4;
  double ti = nperiods * period;//1.e6 * 2. * M_PI * sqrt( U.G * pow(y[0], 3.) / menc );

  int status = GSL_SUCCESS;
  printf("Starting integration...\n");

  menc = gsl_spline_eval(spline_menc, y[0], acc_menc);
  rho = gsl_spline_eval(spline_rho, y[0], acc_rho);

  status = gsl_odeiv2_driver_apply (d, &t, ti, y);

  printf ("Integration stopped at r = %.2e L_unit, with return value %i\n", y[0], status);

  double menc_f = gsl_spline_eval(spline_menc, y[0], acc_menc);
  v_tot_sq = y[2] * y[2] + pow(y[0] * y[3], 2.);
  double energy_f = - U.G * menc_f / y[0] + 0.5 * v_tot_sq;

  double eps_r = fabs(y[0] / P.r0 - 1.);
  double eps_th = fabs(y[1] / nperiods / 2. / M_PI - 1);
  double eps_energy = fabs(energy_f / energy_0 - 1.);

  // Test

  printf("Relative error in radius: %.10e\n", eps_r);
  printf("Relative error in theta : %.10e\n", eps_th);
  printf("Relative error in energy: %.10e\n", eps_energy);

  gsl_odeiv2_driver_free (d);

  if ( status != 0 ) exit(-1);
  if ( eps_r > 1.e-10 ) exit(-1);
  if ( eps_th > 1.e-10 ) exit(-1);
  if ( eps_energy > 1.e-10 ) exit(-1);

  return 0;
}
