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

#include "utils.h"
#include <stdio.h>
#include <stdlib.h>

/// Indices in array that enclose a given real number.
int bounding_indices(double value, double *array, int length, int *indices){

  double L = array[length - 1] - array[0];
  int test_idx = (int) ( ( value - array[0] ) * ( length - 1 ) / L );

  if ( value > array[length - 1] ){
    indices[0] = length - 1;
    indices[1] = length - 1;
  }
  else if ( value < array[0] ){
    indices[0] = 0;
    indices[1] = 0;
  }
  else{
    if (array[test_idx] > value){
      indices[0] = test_idx - 1;
      indices[1] = test_idx;
    }
    else{
      indices[0] = test_idx;
      indices[1] = test_idx + 1;
    }
  }

  return 0;

}

double trilinear_interp(double ***f_vals, double *x_vals, double *y_vals, double *z_vals, int *dims, double x, double y, double z){

  double xhat, yhat, zhat;
  int x_bounds[2], y_bounds[2], z_bounds[2];
  double c3[2][2][2];
  double c00, c01, c10, c11;
  double c0, c1;

  bounding_indices(x, x_vals, dims[0], x_bounds);
  bounding_indices(y, y_vals, dims[1], y_bounds);
  bounding_indices(z, z_vals, dims[2], z_bounds);

  for ( int i = 0; i < 2; i++ ){
    for ( int j = 0; j < 2; j++ ){
      for ( int k = 0; k < 2; k++ ){
        c3[i][j][k] = f_vals[ x_bounds[i] ][ y_bounds[j] ][ z_bounds[k] ];
      }
    }
  }

  xhat = 0.;
  yhat = 0.;
  zhat = 0.;

  if ( x_bounds[0] != x_bounds[1] ) xhat = ( x - x_vals[x_bounds[0]] ) / ( x_vals[x_bounds[1]] - x_vals[x_bounds[0]] );
  if ( y_bounds[0] != y_bounds[1] ) yhat = ( y - y_vals[y_bounds[0]] ) / ( y_vals[y_bounds[1]] - y_vals[y_bounds[0]] );
  if ( z_bounds[0] != z_bounds[1] ) zhat = ( z - z_vals[z_bounds[0]] ) / ( z_vals[z_bounds[1]] - z_vals[z_bounds[0]] );

  c00 = c3[0][0][0] * ( 1 - xhat ) + c3[1][0][0] * xhat;
  c01 = c3[0][0][1] * ( 1 - xhat ) + c3[1][0][1] * xhat;
  c10 = c3[0][1][0] * ( 1 - xhat ) + c3[1][1][0] * xhat;
  c11 = c3[0][1][1] * ( 1 - xhat ) + c3[1][1][1] * xhat;

  c0 = c00 * ( 1 - yhat ) + c10 * yhat;
  c1 = c01 * ( 1 - yhat ) + c11 * yhat;

  return c0 * ( 1 - zhat ) + c1 * zhat;

}

double sign(double x){
  if ( x >= 0.){
    return 1.;
  }
  else{
    return -1.;
  }
}
