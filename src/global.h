#ifndef GLOBAL_HEADER
#define GLOBAL_HEADER

#ifndef M_PI
// Make sure M_PI is defined
#include <gsl/gsl_math.h>
#endif
#include <gsl/gsl_spline.h>
#include "math.h"
#include "hdf5.h"

/// Maximum string length.
#define MAXLEN 500

/// Exits if HDF5 operation was not successful.
void h5safecall(herr_t status);

/// Properties of the stellar profile.
struct star{

/// Number of cells in the profile.
  int ncells;
/// Array with radius values.
  double *r;
/// Array with density values as a function of radius.
  double *rho;
/// Array with density gradient values as a function of radius.
  double *drhodr;
/// Array with enclosed mass as a function of radius.
  double *menc;
/// Array with sound speed as a function of radius.
  double *cs;
/// Array with density scale height as a function of radius.
  double *hrho;

};

/// Input parameters read from text file.
struct parameters{

/// Whether to include hydrodynamical drag.
  int f_hydro;
/// Whether to include gravitational drag.
  int f_grav;
/// Whether to include buoyancy.
  int f_buoyancy;
/// Mass of the embedded object
  double Msb;
/// Radius of the embedded object
  double Rsb;
/// Derived quantity: average density of the embedded object
  double rho_av_sb;
/// Initial orbital separation
  double r0;
/// Path to MESA profile.
  char mesa_profile_path[MAXLEN];

// Drag coefficient stuff

/// Whether to use numerical drag coefficients
  int use_drag_coefficients;
/// Path to drag coefficient data
  char drag_data_path[MAXLEN];

};

/// Values in CGS for units in each dimensions and physical constants.
struct units{

/// Mass unit.
  double m_unit;
/// Length unit.
  double l_unit;
/// Time unit.
  double t_unit;
/// Speed unit.
  double v_unit;
/// Density unit.
  double rho_unit;

/// Gravitational constant in CGS.
  double G_CGS;
/// Solar mass in CGS.
  double MSUN_CGS;
/// Solar radius in CGS.
  double RSUN_CGS;

/// Gravitational constant (in simulation units).
  double G;

};

/// Globally accessible GSL interpolation accelerator for the density spline.
extern gsl_interp_accel *acc_rho;
/// Globally accessible GSL interpolation accelerator for the density gradient spline.
extern gsl_interp_accel *acc_drhodr;
/// Globally accessible GSL interpolation accelerator for the mass spline.
extern gsl_interp_accel *acc_menc;
/// Globally accessible GSL interpolation accelerator for the sound speed spline.
extern gsl_interp_accel *acc_cs;
/// Globally accessible GSL interpolation accelerator for the density scale height spline.
extern gsl_interp_accel *acc_hrho;
/// Globally accessible GSL density spline.
extern gsl_spline *spline_rho;
/// Globally accessible GSL density gradient spline.
extern gsl_spline *spline_drhodr;
/// Globally accessible GSL enclosed mass spline.
extern gsl_spline *spline_menc;
/// Globally accessible GSL sound speed spline.
extern gsl_spline *spline_cs;
/// Globally accessible GSL density scale height spline.
extern gsl_spline *spline_hrho;

#endif
