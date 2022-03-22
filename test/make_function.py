import numpy as np
import h5py

def f(x, y, z):
  return x * x - y + 2 * z

n_x = 100
n_y = 100
n_z = 100

x = np.linspace(0, 1, n_x)
y = np.linspace(0, 1, n_y)
z = np.linspace(0, 1, n_z)

f_vals = np.empty( ( len(x), len(y), len(z) ) )

for x_idx, x_val in enumerate(x):
  for y_idx, y_val in enumerate(y):
    for z_idx, z_val in enumerate(z):
      f_vals[x_idx, y_idx, z_idx] = f(x_val, y_val, z_val)

f = h5py.File('interp_data.h5', 'w')

f.create_dataset('npoints', data = np.array([n_x, n_y, n_z]))
f.create_dataset('x_vals', data = x)
f.create_dataset('y_vals', data = y)
f.create_dataset('z_vals', data = z)
f.create_dataset('data', data = f_vals)
f.close()
