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

#include "star.h"

/// Globally accessible GSL interpolation accelerator for the density spline.
gsl_interp_accel *acc_rho;
/// Globally accessible GSL interpolation accelerator for the density gradient spline.
gsl_interp_accel *acc_drhodr;
/// Globally accessible GSL interpolation accelerator for the mass spline.
gsl_interp_accel *acc_menc;
/// Globally accessible GSL interpolation accelerator for the sound speed spline.
gsl_interp_accel *acc_cs;
/// Globally accessible GSL interpolation accelerator for the density scale height spline.
gsl_interp_accel *acc_hrho;
/// Globally accessible GSL density spline.
gsl_spline *spline_rho;
/// Globally accessible GSL density gradient spline.
gsl_spline *spline_drhodr;
/// Globally accessible GSL enclosed mass spline.
gsl_spline *spline_menc;
/// Globally accessible GSL sound speed spline.
gsl_spline *spline_cs;
/// Globally accessible GSL density scale height spline.
gsl_spline *spline_hrho;

int init_stellar_profile_splines(struct star S){

  const gsl_interp_type *t = gsl_interp_steffen;

  spline_rho = gsl_spline_alloc(t, S.ncells);
  gsl_spline_init(spline_rho, S.r, S.rho, S.ncells);

  spline_drhodr = gsl_spline_alloc(t, S.ncells);
  gsl_spline_init(spline_drhodr, S.r, S.drhodr, S.ncells);

  spline_menc = gsl_spline_alloc(t, S.ncells);
  gsl_spline_init(spline_menc, S.r, S.menc, S.ncells);

  spline_cs = gsl_spline_alloc(t, S.ncells);
  gsl_spline_init(spline_cs, S.r, S.cs, S.ncells);

  spline_hrho = gsl_spline_alloc(t, S.ncells);
  gsl_spline_init(spline_hrho, S.r, S.hrho, S.ncells);

  return 0;
}
