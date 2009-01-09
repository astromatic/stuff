/*
 				galaxies.h

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	Stuff
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	Include file for galaxies.c
*
*	Last modify:	07/06/2005
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

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

#define	GAL_BULGEMKNEE	(-20.0)	/* Abs. mag at knee in Binggeli et al. 84 */
#define	GAL_BULGERKNEE	1581.0	/* r_e (in h-1.pc) at knee in Binggeli et al */
#define	GAL_BULGEQMEAN	0.65	/* Mean of intrinsic bulge flattening distr. */
#define	GAL_BULGEQSIG	0.18	/* Sigma of bulge flattening distribution */

#define	GAL_MAXINCLIN	85.0	/* Max. allowed galaxy inclination (deg) */

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
  double	disk_taufact;	/* Extinction factor in the ref. band */
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
  double	*mag;		/* Apparent magnitude in each passband */
  } galstruct;

/*-------------------------------- protos -----------------------------------*/

galtypestruct	*galtype_init(double hubtype, double bt, double extinct,
			lfstruct *lf,
			sedstruct *bsed, sedstruct *dsed,
			sedstruct *tau_i, sedstruct *refpb);
galstruct	*gal_init(galtypestruct *galtype, double z, double mabsmax,
			sedstruct **pb, int npb);

double		gal_bulgeflat(double cosi),
		gal_bulgesize(double mabs),
		gal_disksize(galtypestruct *galtype, double mabs),
		gal_cosi(void);

void		galtype_end(galtypestruct *galtype),
		gal_end(galstruct *gal);

#endif
