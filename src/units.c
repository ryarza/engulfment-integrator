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

#include "math.h"
#include <stdio.h>
#include "units.h"

int set_units(struct parameters * P, struct star * S, struct units * U){

  U->G_CGS = 6.67430e-8;
  U->RSUN_CGS = 6.957e10;
  U->MSUN_CGS = 1.3271244e26 / U->G_CGS;

  // Mass unit = mass of the star
  U->m_unit = S->menc[S->ncells-1];
  // Length unit = Radius of the star
  U->l_unit = S->r[S->ncells-1];
  // Time unit: orbital dynamical time at the radius of the star. Comes with perk that G = 1 in these units.
  U->t_unit = sqrt(pow(U->l_unit, 3.) / U->G_CGS / U->m_unit);

  U->v_unit = U->l_unit / U->t_unit;
  U->rho_unit = U->m_unit / pow(U->l_unit, 3.);
  U->G = U->G_CGS * U->m_unit * pow(U->t_unit, 2.) / pow(U->l_unit, 3.);

  printf("Units:\n");
  printf("  Length: %.10e cm\n", U->l_unit);
  printf("  Time  : %.10e s\n", U->t_unit);
  printf("  Mass  : %.10e g\n", U->m_unit);
  printf("Value of G in simulation units: %.10e L_unit^3 / M_unit / T_unit^2\n", U->G);

  for (int i = 0; i < S->ncells; i++){
    S->r[i] /= U->l_unit;
    S->rho[i] /= U->rho_unit;
    S->drhodr[i] /= ( U->rho_unit / U->l_unit );
    S->menc[i] /= U->m_unit;
    S->cs[i] /= ( U->l_unit / U->t_unit );
    S->hrho[i] /= U->l_unit;
  }

  P->Msb /= U->m_unit;
  P->Rsb /= U->l_unit;

// Compute derived quantities
  P->rho_av_sb = P->Msb / ( 4. * M_PI * pow(P->Rsb, 3.) / 3. );

  return 0;

}
