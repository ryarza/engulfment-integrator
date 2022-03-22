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

#ifndef UNITS_HEADER
#define UNITS_HEADER

#include "global.h"
#include "io.h"
#include "model.h"
#include "star.h"

/// Defines the units and stores them in the units struct. Also converts parameters and the stellar profile into these units.
int set_units(struct parameters * P, struct star * S, struct units * U);

#endif
