/*
*				galaxies.h
*
* Include file for galaxies.c.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	Stuff
*
*	Copyright:		(C) 1999-2016 IAP/CNRS/UPMC
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
*	Last modified:		28/11/2016
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef _SED_H_
#include "sed.h"
#endif

#ifndef _LF_H_
#include "lf.h"
#endif

#ifndef _GALAXIES_H_
#define _GALAXIES_H_

/*------------------------------- constants ---------------------------------*/

#define	GAL_MAXNTYPE	8	/* Maximum number of distinct galaxy types */
#define	GAL_MAXNSED	8	/* Maximum number of galaxy SED components */

#define	GAL_BULGEQMEAN	0.65	/* Mean of intrinsic bulge flattening distr. */
#define	GAL_BULGEQSIG	0.18	/* Sigma of bulge flattening distribution */

#define	GAL_MAXINCLIN	85.0	/* Max. allowed galaxy inclination (deg) */
#define	GAL_INCLIN_STEP	0.01	/* Inclin. step angle for integrations (deg) */

/*--------------------------------- flags -----------------------------------*/

/*-------------------------- structure definitions --------------------------*/

typedef struct
  {
  double	hubtype;/* Hubble type -6.0 <= T <= 10.0 */
  lfstruct	*lf;	/* Type-dependent luminosity function */
  double	bt;	/* Bulge-to-Total (B/T) ratio in the ref. band */
  double	bdm;	/* -2.5*log10(bt) */
  double	ddm;	/* -2.5*log10(1-bt) */
  double	beta;	/* Tully-Fisher exponent in the ref. band */
  double	disk_rstar;	/* Disk scale-length */
  double	disk_sigmalambda;
  double	disk_extinct;	/* Extinction */
  double	disk_taufact;	/* Extinction factor in the ref. band */
  double	disk_tau_a;
  double	disk_tau_b;
  double	bulge_rstar;	/* Bulge effective radius */
  sedstruct	*bsed;	/* Bulge SED */
  sedstruct	*dsed;	/* Disk SED */
  sedstruct	*tau_i;	/* Disk extinction curve */
  } galtypestruct;

typedef struct
  {
  double	bposang;	/* Bulge position angle (in degrees) */
  double	bsize;		/* Bulge apparent effective radius (in arcsec)*/
  double	bflat;		/* Bulge apparent flattening (A/R) */
  double	dposang;	/* Disk apparent position angle (in degrees) */
  double	dsize;		/* Disk apparent scale-length (in arcsec) */
  double	dflat;		/* Disk apparent flattening (A/R) */
  double	*bt;		/* Apparent bulge/total in each passband */
  double	mabs;		/* Absolute magnitude in ref. passband */
  double	refmag;		/* Apparent magnitude in ref. passband */
  double	*mag;		/* Apparent magnitude in each passband */
  } galstruct;

/*-------------------------------- protos -----------------------------------*/

galtypestruct	*galtype_init(double hubtype, double bt, double extinct,
			lfstruct *lf,
			sedstruct *bsed, sedstruct *dsed,
			sedstruct *tau_i, sedstruct *refpb);
galstruct	*gal_init(galtypestruct *galtype, double z,
			double mabsmax, sedstruct *refpb,
			sedstruct **pb, int npb);

double		gal_bulgeflat(double cosi),
		gal_bulgesize(galtypestruct *galtype, double mabs),
		gal_cosi(void),
		gal_disksize(galtypestruct *galtype, double mabs),
		gal_mag(galtypestruct *galtype, double mabs, double z,
			sedstruct *pb, sedstruct *dsed, double *btout);

void		galtype_end(galtypestruct *galtype),
		gal_end(galstruct *gal),
		gal_shear(galstruct *gal, double kappa, double *gamma, int npb);

#endif

