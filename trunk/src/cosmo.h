/*
 				cosmo.h

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	Stuff
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	Include file for cosmo.c
*
*	Last modify:	10/06/2005
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

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

