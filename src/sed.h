/*
*				sed.h
*
* Include file for sed.c.
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

#ifndef _SED_H_
#define _SED_H_

/*------------------------------- constants ---------------------------------*/

#define	SED_MAXNPB	128	/* Maximum number of bandpasses */
#define	SED_NDATA	1000	/* Default size of SED array */
#define	SED_MAXNCOMP	16	/* Maximum number of passband components */
#define	SED_DLAMBDA    	1e-6	/* Relative precision in lambda */
#define	REF_WAVELENGTH	(5556*ANGSTROEM)
				/* Reference wavelength for mag 0-point */
#define	REF_PHOTENERGY	(PLANCK*C/REF_WAVELENGTH)
				/* Photon energy at REF_WAVELENGTH */
#define	REF_ENERGY	(3640*JANSKY*C/(REF_WAVELENGTH*REF_WAVELENGTH))
#define	REF_PHOTRATE	(REF_ENERGY/REF_PHOTENERGY)
				/* phot/m2/s/m from an A0 0 mag star
				at REF_WAVELENGTH (Bessel 1979) */

/*--------------------------------- flags -----------------------------------*/

#define	SED_FLAMBDADATA	0	/* Data are originally f_lambda or passband */
#define	SED_FNUDATA    	1	/* Data are originally in f_nu units */

/*-------------------------- structure definitions --------------------------*/

typedef struct
  {
  char		name[MAXCHAR];	/* Bandpass name(s) */
  double	*wave;		/* Wavelengths */
  double	*data;		/* Response curve */
  double	wavemin,wavemax;/* Boundaries of the SED */
  int		ndata;		/* Total number of steps */
    int		fnu_flag;	/* Data in f_nu units? */
  } sedstruct;


/*-------------------------------- protos -----------------------------------*/

sedstruct	*sed_dup(sedstruct *sed),
		*sed_load(char *datadir_name, char *sed_name);

double		sed_calib(sedstruct *sed, sedstruct *pb),
		sed_kcor(sedstruct *sed, sedstruct *pb, double z),
		sed_mul(sedstruct *sed1, double wavefact1,
			sedstruct *sed2, double wavefact2,
			sedstruct **sedo);

void		pb_calib(sedstruct **pb, sedstruct **pbcalibsed, int npb,
			sedstruct *refpb, sedstruct *refcalibsed,
			sedstruct *backsed),
		sed_end(sedstruct *sed),
		sed_extinc(sedstruct *sed, sedstruct *tau, double taufact,
			sedstruct **sedo);

#endif
