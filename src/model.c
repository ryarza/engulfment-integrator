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

#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_errno.h>
#include "hdf5.h"
#include "model.h"
#include "utils.h"

int jac (double t, const double y[], double *dfdy, double dfdt[], void *params)
{
  (void)(t); /* avoid unused parameter warning */
  struct rhs_terms * rhs = (struct rhs_terms *)params;

  double menc, rho, drhodr;
  int status1 = gsl_spline_eval_e(spline_menc, y[0], acc_menc, &menc);
  int status2 = gsl_spline_eval_e(spline_rho, y[0], acc_rho, &rho);
  int status3 = gsl_spline_eval_e(spline_drhodr, y[0], acc_drhodr, &drhodr);
  if ( ( status1 != GSL_SUCCESS ) || ( status2 != GSL_SUCCESS ) || ( status3 != GSL_SUCCESS ) ){
    printf("Interpolation failed at radius %.10e. Exit codes for menc, rho, drhodr: %i %i %i\n", y[0], status1, status2, status3);
    exit(-1);
  }

  gsl_matrix_view dfdy_mat = gsl_matrix_view_array (dfdy, 4, 4);
  gsl_matrix * m = &dfdy_mat.matrix;

  gsl_matrix_set (m, 0, 0, 0.);
  gsl_matrix_set (m, 0, 1, 0.);
  gsl_matrix_set (m, 0, 2, 1.);
  gsl_matrix_set (m, 0, 3, 0.);

  gsl_matrix_set (m, 1, 0, 0.);
  gsl_matrix_set (m, 1, 1, 0.);
  gsl_matrix_set (m, 1, 2, 0.);
  gsl_matrix_set (m, 1, 3, 1.);

  gsl_matrix_set (m, 2, 0, 2 * rhs->G * menc / pow(y[0], 3.) - 4. * rhs->G * M_PI * rho + pow(y[3], 2.));
//  gsl_matrix_set (m, 2, 0, 2 * rhs->G * menc / pow(y[0], 3.) + pow(y[3], 2.));
  gsl_matrix_set (m, 2, 1, 0.);
  gsl_matrix_set (m, 2, 2, 0.);
  gsl_matrix_set (m, 2, 3, 2 * y[0] * y[3]);

/*
  Old terms
  double term1 = 32. * pow(rhs->G, 2.) * pow(M_PI, 2.) * menc * pow(y[0], 3.) * pow(rho, 2.);
  double term2 = 4. * pow(rhs->G, 2.) * M_PI * pow(menc, 2.) * ( - 3. * rho + y[0] * drhodr );
  double term3 = pow(y[0], 2.) * pow(y[3], 3.) * ( 2 * rhs->Msb * y[2] + M_PI * pow(rhs->Rsb, 2.) * pow(y[0], 2.) * y[3] * ( rho + y[0] * drhodr ) );
  double num = term1 + term2 + term3;
  double den = rhs->Msb * pow(y[0], 4.) * pow(y[3], 2.);
  double j30 = num/den;
*/

  // New (hopefully still correct) terms that are more numerically stable
  double term1 = - 12. * pow(rhs->G, 2.) * M_PI * pow(menc, 2) * rho / ( rhs->Msb * pow(y[0], 4.) * pow(y[3], 2.) );
  double term2 = 32. * pow(rhs->G * M_PI * rho, 2.) * menc / ( rhs->Msb * y[0] * pow(y[3], 2.) );
  double term3 = 2. * y[2] * y[3] / pow(y[0], 2.);
  double term4 = M_PI * pow(rhs->Rsb * y[3], 2.) * rho / rhs->Msb;
  double term5 = 4. * M_PI * pow(rhs->G * menc, 2.) * drhodr / ( rhs->Msb * pow(y[0], 3.) * pow(y[3], 2.) );
  double term6 = M_PI * pow(rhs->Rsb * y[3], 2.) * y[0] * drhodr / rhs->Msb;
  double j30 = term1 + term2 + term3 + term4 + term5 + term6;

  gsl_matrix_set (m, 3, 0, j30);

  // Just Keplerian term
//  gsl_matrix_set (m, 3, 0, 2 * y[2] * y[3] / pow(y[0], 2));
  gsl_matrix_set (m, 3, 1, 0.);
  gsl_matrix_set (m, 3, 2, - 2. * y[3] / y[0]);

  term1 = -2. * y[2] / y[0];
  term2 = - 8. * pow(rhs->G, 2.) * M_PI * pow(menc, 2.) * rho / ( rhs->Msb * pow(y[0], 3.) * pow(y[3], 3.) );
  term3 = 2. * M_PI * pow(rhs->Rsb, 2.) * y[0] * rho * y[3] / rhs->Msb;
  double j33 = term1 + term2 + term3;
  gsl_matrix_set (m, 3, 3, j33);

  // Just Keplerian term
//  gsl_matrix_set (m, 3, 3, -2. * y[2] / y[0]);

  // No explicit time dependence
  dfdt[0] = 0.0;
  dfdt[1] = 0.0;
  dfdt[2] = 0.0;
  dfdt[3] = 0.0;

  return GSL_SUCCESS;

}

double f_grav(struct rhs_terms *rhs, double rho, double v_theta, double q, double eps_rho, double rp_over_ra){

  double cg = 1.;
  if ( rhs->use_drag_coefficients == 1 ){
    cg = trilinear_interp(rhs->cg_data, rhs->log_q_vals, rhs->log_eps_rho_vals, rhs->log_rp_over_ra_vals, rhs->drag_data_dims, log10(q), log10(eps_rho), log10(rp_over_ra));
  }

  return - 4. * M_PI * pow(rhs->G, 2.) * sign(v_theta) * cg * rho * pow(rhs->Msb, 2.) / fabs(v_theta * v_theta);

}

double f_hydro(struct rhs_terms *rhs, double rho, double v_theta, double q, double eps_rho, double rp_over_ra){

  double cp = 0.25;
  if ( rhs->use_drag_coefficients == 1){
    cp = trilinear_interp(rhs->cp_data, rhs->log_q_vals, rhs->log_eps_rho_vals, rhs->log_rp_over_ra_vals, rhs->drag_data_dims, log10(q), log10(eps_rho), log10(rp_over_ra));
  }

  return - sign(v_theta) * cp * rho * v_theta * v_theta * rhs->sigma_hydro;

}

int compute_flow_parameters(struct rhs_terms *rhs, const double *y, double menc, double hrho, double *q, double *eps_rho, double *rp_over_ra)
{
  double cs;
  gsl_spline_eval_e(spline_cs, y[0], acc_cs, &cs);
  double v_theta = y[0] * y[3];
//  double ra = 2. * rhs->G * rhs->Msb / ( pow(v_theta, 2.) + pow(cs, 2.) );
  double ra = 2. * rhs->G * rhs->Msb / ( pow(v_theta, 2.) );
  *q = rhs->Msb / menc;
  *eps_rho = ra / hrho;
  *rp_over_ra = rhs->Rsb / ra;
  return 0;
  
}

int func (double t, const double y[], double f[], void *params)
{

  (void)(t); /* avoid unused parameter warning */

  struct rhs_terms * rhs = (struct rhs_terms *)params;

  double menc, rho, cs, hrho, q, eps_rho, rp_over_ra;
  int status1 = gsl_spline_eval_e(spline_menc, y[0], acc_menc, &menc);
  int status2 = gsl_spline_eval_e(spline_rho, y[0], acc_rho, &rho);
  int status3 = gsl_spline_eval_e(spline_cs, y[0], acc_cs, &cs);
  int status4 = gsl_spline_eval_e(spline_hrho, y[0], acc_hrho, &hrho);
//  printf("%.2e %.2e %.2e %.2e %.2e\n", y[0], menc, rho, cs, hrho);
  if ( ( status1 != GSL_SUCCESS ) || ( status2 != GSL_SUCCESS ) || ( status3 != GSL_SUCCESS ) || ( status4 != GSL_SUCCESS ) ){
    printf("Interpolation failed at radius %.10e. Exit codes for menc, rho, cs, hrho splines: %i %i %i %i\n", y[0], status1, status2, status3, status4);
    exit(-1);
  }

  compute_flow_parameters(rhs, y, menc, hrho, &q, &eps_rho, &rp_over_ra);

  double f_r = 0.;
  double f_theta = 0.;
  double v_theta = y[0] * y[3];
  double mach = fabs(v_theta / cs);

  // Hydrodynamical drag
  if ( rhs->f_hydro == 1 ) f_theta += f_hydro(rhs, rho, v_theta, q, eps_rho, rp_over_ra);

  // Gravitational drag
  if ( rhs->f_grav == 1 && mach > 0.25 ) f_theta += f_grav(rhs, rho, v_theta, q, eps_rho, rp_over_ra);

  // Buoyancy
  if ( rhs->f_buoyancy == 1 ) f_r += rho * 4. * M_PI * pow(rhs->Rsb, 3.) / 3.;

//  printf("q = %.2e, eps_rho = %.2e, Rp/Ra = %.2e, Cp = %.2e, Cg = %.2e\n", rhs->Msb / menc, eps_rho, rhs->Rsb / ra, cp, cg);

  //  y[0] = r, y[1] = theta, y[2] = vr, y[3] = d(theta)/dt
  // dr/dt = vr
  f[0] = y[2];
//  printf("f[0] = %.2e\n", f[0]);
  // d(theta)/dt = d(theta)/dt duh
  f[1] = y[3];
//  printf("f[1] = %.2e\n", f[1]);
  // dvr/dt = - G Menc / r^2 + rotational support
  f[2] = - rhs->G * menc / y[0] / y[0] + v_theta * v_theta / y[0] + f_r / rhs->Msb;
//  printf("f[2] = %.2e\n", f[2]);
  // d2(theta)/dt2 = coordinate terms + acceleration from forces
  f[3] = -2. * y[2] * y[3] / y[0] + f_theta / rhs->Msb / y[0];
//  printf("f[3] = %.2e\n", f[3]);

//  printf(" %.10e %.10e\n", y[0], - 4. * M_PI * pow(rhs->G, 2.) * sign(v_theta) * cg * rho * pow(rhs->Msb, 2.) / fabs(v_theta * v_theta));
//  double edot = rhs->G * menc * y[2] / pow(y[0], 2.) + y[2] * ( y[0] * ( - 4. * rhs->G * M_PI * rho + pow(y[3], 2.) ) + f[2] ) + pow(y[0], 2.) * y[3] * f[3];
//  edot *= rhs->Msb;
//  printf(" %.20e %.20e\n", y[0], edot);

  return GSL_SUCCESS;

}

int init_model_params(struct parameters P, struct star S, struct units U, struct rhs_terms *rhs){

  rhs->sigma_hydro = M_PI * P.Rsb * P.Rsb;
  rhs->Msb = P.Msb;
  rhs->Rsb = P.Rsb;
  rhs->G = U.G;
  rhs->f_hydro = P.f_hydro;
  rhs->f_grav = P.f_grav;
  rhs->f_buoyancy = P.f_buoyancy;

// Numerical drag coefficients stuff
  rhs->use_drag_coefficients = P.use_drag_coefficients;

  if ( P.use_drag_coefficients == 1 ){

  hid_t file_id, dataset_id;

  // Open file
  file_id = H5Fopen(P.drag_data_path, H5F_ACC_RDONLY, H5P_DEFAULT);

  // Get number of cells
  dataset_id = H5Dopen2(file_id, "/npoints", H5P_DEFAULT);
  h5safecall(H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &rhs->drag_data_dims));
  h5safecall(H5Dclose(dataset_id));

  rhs->log_q_vals = (double *) malloc(rhs->drag_data_dims[0] * sizeof(double));
  rhs->log_eps_rho_vals = (double *) malloc(rhs->drag_data_dims[1] * sizeof(double));
  rhs->log_rp_over_ra_vals = (double *) malloc(rhs->drag_data_dims[2] * sizeof(double));

  // log10(q) values
  dataset_id = H5Dopen2(file_id, "/log10_qs", H5P_DEFAULT);
  h5safecall(H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, rhs->log_q_vals));
  h5safecall(H5Dclose(dataset_id));

/*
  printf("q values\n");
  for ( int i = 0; i < rhs->drag_data_dims[0]; i++ ) printf("%.2e\n", pow(10, rhs->log_q_vals[i]));
*/

  // log10(eps_rho) values
  dataset_id = H5Dopen2(file_id, "/log10_eps_rhos", H5P_DEFAULT);
  h5safecall(H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, rhs->log_eps_rho_vals));
  h5safecall(H5Dclose(dataset_id));

/*
  printf("eps_rho values\n");
  for ( int i = 0; i < rhs->drag_data_dims[1]; i++ ) printf("%.2e\n", pow(10, rhs->log_eps_rho_vals[i]));
*/

  // log10(rp_over_ra) values
  dataset_id = H5Dopen2(file_id, "/log10_rp_over_ras", H5P_DEFAULT);
  h5safecall(H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, rhs->log_rp_over_ra_vals));
  h5safecall(H5Dclose(dataset_id));

/*
  printf("rp_over_ra values\n");
  for ( int i = 0; i < rhs->drag_data_dims[2]; i++ ) printf("%.2e\n", pow(10, rhs->log_rp_over_ra_vals[i]));
*/

  // Load gravitational drag
//  printf("Allocating gravitational drag\n");
  rhs->cg_data = (double***)malloc(rhs->drag_data_dims[0] * sizeof(double**));
  for ( int i = 0; i < rhs->drag_data_dims[0]; i++ ){
    rhs->cg_data[i] = (double**)malloc(rhs->drag_data_dims[1] * sizeof(double*));
    for ( int j = 0; j < rhs->drag_data_dims[1]; j++ ){
      rhs->cg_data[i][j] = (double*)malloc(rhs->drag_data_dims[2] * sizeof(double));
    }
  }

//  for ( int i = 0; i < 3; i++ ) printf("%i\n", rhs->drag_data_dims[i]);

  double drag_data_temp[rhs->drag_data_dims[0]][rhs->drag_data_dims[1]][rhs->drag_data_dims[2]];

  dataset_id = H5Dopen2(file_id, "/cg", H5P_DEFAULT);
  h5safecall(H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &drag_data_temp));
  h5safecall(H5Dclose(dataset_id));

  for ( int i = 0; i < rhs->drag_data_dims[0]; i++ ){
    for ( int j = 0; j < rhs->drag_data_dims[1]; j++ ){
      for ( int k = 0; k < rhs->drag_data_dims[2]; k++ ){
        rhs->cg_data[i][j][k] = drag_data_temp[i][j][k];
//        printf("%.2e %.2e %.2e %.2e\n", pow(10,rhs->log_q_vals[i]), pow(10,rhs->log_eps_rho_vals[j]), pow(10,rhs->log_rp_over_ra_vals[k]), rhs->cg_data[i][j][k]);
      }
    }
  }

//  for ( int i = 0; i < 3; i++ ) printf("%i\n", rhs->drag_data_dims[i]);


//  printf("Allocating ram pressure drag\n");
  // Load ram pressure drag
  rhs->cp_data = (double***)malloc(rhs->drag_data_dims[0] * sizeof(double**));
  for ( int i = 0; i < rhs->drag_data_dims[0]; i++ ){
    rhs->cp_data[i] = (double**)malloc(rhs->drag_data_dims[1] * sizeof(double*));
    for ( int j = 0; j < rhs->drag_data_dims[1]; j++ ){
      rhs->cp_data[i][j] = (double*)malloc(rhs->drag_data_dims[2] * sizeof(double));
    }
  }

  dataset_id = H5Dopen2(file_id, "/cp", H5P_DEFAULT);
  h5safecall(H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, drag_data_temp));
  h5safecall(H5Dclose(dataset_id));

  for ( int i = 0; i < rhs->drag_data_dims[0]; i++ ){
    for ( int j = 0; j < rhs->drag_data_dims[1]; j++ ){
      for ( int k = 0; k < rhs->drag_data_dims[2]; k++ ){
        rhs->cp_data[i][j][k] = drag_data_temp[i][j][k];
      }
    }
  }

  // Close the file
  h5safecall(H5Fclose(file_id));

  }

  return 0;
  
}
