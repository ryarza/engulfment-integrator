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

#ifndef IO_HEADER
#define IO_HEADER

#include "global.h"
#include "model.h"

/// Reads parameters from a text file. Based on Cholla.
void parse_params (char *param_file, struct parameters * parms);

/// Reads stellar profile from an HDF5 file. The path is contained in the parameter struct, and the properties of the profile are stored in the star struct.
int read_stellar_profile(struct parameters P, struct star *S);

/// One-time output operation that saves the simulation units and parameters
int write_scalars(struct units U, struct parameters P);

/// Writes values of important variables for one step into the output text file
int write_output(struct units U, struct rhs_terms *rhs, double t, double *y);

#endif
