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
#include <math.h>
#include <hdf5.h>

#include "../src/io.h"
#include "../src/utils.h"

double get_rand(double min, double max){
    return min + (max - min) * rand() / RAND_MAX;
}

// Test interpolation by loading grids and comparing against Python
double f(double x, double y, double z){
  return x * x - y + 2. * z;
}


int main(){

  double *x_vals, *y_vals, *z_vals;
  int dims[3];
  double ***data;

  hid_t file_id, dataset_id;
  herr_t status;
  hid_t attr;

  // Open file
  file_id = H5Fopen("interp_data.h5", H5F_ACC_RDONLY, H5P_DEFAULT);

  printf("Opened file!\n");
  // Get number of cells
  dataset_id = H5Dopen2(file_id, "/npoints", H5P_DEFAULT);
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &dims);
  status = H5Dclose(dataset_id);

  printf("Read number of cells\n");

  x_vals = (double *) malloc(dims[0] * sizeof(double));
  y_vals = (double *) malloc(dims[1] * sizeof(double));
  z_vals = (double *) malloc(dims[2] * sizeof(double));
  double data_temp[dims[0]][dims[1]][dims[2]];
  // Load ram pressure drag
  data = (double***)malloc(dims[0] * sizeof(double**));
  for ( int i = 0; i < dims[0]; i++ ){
    data[i] = (double**)malloc(dims[1] * sizeof(double*));
    for ( int j = 0; j < dims[1]; j++ ){
      data[i][j] = (double*)malloc(dims[2] * sizeof(double));
    }
  }

  dataset_id = H5Dopen2(file_id, "/x_vals", H5P_DEFAULT);
  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, x_vals);
  status = H5Dclose(dataset_id);
  printf("Read x vals\n");

  dataset_id = H5Dopen2(file_id, "/y_vals", H5P_DEFAULT);
  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, y_vals);
  status = H5Dclose(dataset_id);
  printf("Read y vals\n");

  dataset_id = H5Dopen2(file_id, "/z_vals", H5P_DEFAULT);
  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, z_vals);
  status = H5Dclose(dataset_id);
  printf("Read z vals\n");

  dataset_id = H5Dopen2(file_id, "/data", H5P_DEFAULT);
  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data_temp);
  status = H5Dclose(dataset_id);
  printf("Read data\n");

  for ( int i = 0; i < dims[0]; i++ ){
    for ( int j = 0; j < dims[1]; j++ ){
      for ( int k = 0; k < dims[2]; k++ ){
        data[i][j][k] = data_temp[i][j][k];
      }
    }
  }

// Just choose a bunch of random points and make sure they agree with the analytical function to some percentage
  double x, y, z, value;
  for ( int i =0; i < 1000; i++ ){
    x = get_rand(0., 1.);
    y = get_rand(0., 1.);
    z = get_rand(0., 1.);
    value = trilinear_interp(data, x_vals, y_vals, z_vals, dims, x, y, z);
    if ( fabs( value / f(x, y, z) - 1.) > 1.e-2 ) return -1;
  }
  return 0;
}
