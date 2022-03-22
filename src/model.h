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

#ifndef MODEL_HEADER
#define MODEL_HEADER

#include "global.h"

/// Parameters needed when evaluating the right-hand-side of the GSL ODE system.
struct rhs_terms{

/// Whether to include hydrodynamical drag.
  int f_hydro;
/// Whether to include gravitational drag.
  int f_grav;
/// Whether to include gravitational drag.
  int f_buoyancy;
/// Geometrical cross section of the embedded object.
  double sigma_hydro;
/// Radius of the embedded object.
  double Rsb;
/// Mass of the embedded object.
  double Msb;
/// Gravitational constant in simulation units.
  double G;
/// Whether to use numerical drag coefficients.
  int use_drag_coefficients;
/// Gravitational drag coefficient data
  double ***cg_data;
/// Ram pressure drag coefficient data
  double ***cp_data;
/// Dimensions of the grid in which the drag coefficient data is sampled
  int drag_data_dims[3];
/// Mass ratio values at which the drag coefficient data is sampled
  double *log_q_vals;
/// Density gradient values at which the drag coefficient data is sampled
  double *log_eps_rho_vals;
/// Rp / Ra values at which the drag coefficient data is sampled
  double *log_rp_over_ra_vals;


};

/// Jacobian of the system of equations. Useful for ODE integrators that require it. Untested!
int jac (double t, const double y[], double *dfdy, double dfdt[], void *params);

/// Evaluates the vector right-hand-side of the GSL ODE system.
int func (double t, const double y[], double f[], void *params);

/// Initializes the right-hand-side parameter struct and copies the required parameters from other structs. Members of this struct are never modified and hold the value they had in the struct they were read from.
int init_model_params(struct parameters P, struct star S, struct units U, struct rhs_terms *rhs);

/// Ram pressure drag force
double f_hydro(struct rhs_terms *rhs, double rho, double v_theta, double menc, double eps_rho, double ra);

/// Gravitational drag force
double f_grav(struct rhs_terms *rhs, double rho, double v_theta, double menc, double eps_rho, double ra);

/// Relevant dimensionless flow parameters (mass ratio, density gradient, and ratio between geometrical and gravitational radii)
int compute_flow_parameters(struct rhs_terms *rhs, const double *y, double menc, double hrho, double *q, double *eps_rho, double *rp_over_ra);

#endif
