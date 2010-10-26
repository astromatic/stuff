/*
*				cosmo.h
*
* Include file for cosmo.c.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	Stuff
*
*	Copyright:		(C) 1999-2010 Emmanuel Bertin -- IAP/CNRS/UPMC
*
*	License:		GNU General Public License
*
*	Stuff is free software: you can redistribute it and/or modify
*	it under the terms of the GNU General Public License as published by
*	the Free Software Foundation, either version 3 of the License, or
*	(at your option) any later version.
*	Stuff is distributed in the hope that it will be useful,
*	but WITHOUT ANY WARRANTY; without even the implied warranty of
*	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*	GNU General Public License for more details.
*	You should have received a copy of the GNU General Public License
*	along with Stuff. If not, see <http://www.gnu.org/licenses/>.
*
*	Last modified:		26/10/2010
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef	_LF_H_
#include "lf.h"
#endif

/*--------------------------------- constants -------------------------------*/

#define		Z_MIN	0.001	/* Minimum redshift spanned by cosmology */
#define		Z_MAX	20.0	/* Maximum redshift spanned by cosmology */
#define		Z_NSTEP	10000	/* Number of integration steps */

/*------------------------------ global variables ---------------------------*/

double	H, H0, OmegaM, OmegaL, deltaMH;

/*-------------------------------- protos -----------------------------------*/

double	cosmo_dconf(double z),
	cosmo_dlc(double z),
	cosmo_dlum(double z),
	cosmo_dmodul(double z),
	cosmo_dmoduldif(double z, lfstruct *lf, double mappmax),
	cosmo_dvol(double z),
	cosmo_zmax(lfstruct *lf, double mappmax);

