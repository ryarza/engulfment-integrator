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

#ifndef UTILS_HEADER
#define UTILS_HEADER

/// Trilinear interpolation.
double trilinear_interp(double ***f_vals, double *x_vals, double *y_vals, double *z_vals, int *dims, double x, double y, double z);

/// Sign function.
double sign(double x);

#endif
