/*
 				stars.h

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	Stuff
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	Include file for stars.c
*
*	Last modify:	05/06/2005
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifndef _SED_H_
#include "sed.h"
#endif

#ifndef _LF_H_
#include "lf.h"
#endif

#ifndef _STARS_H_
#define _STARS_H_

/*------------------------------- constants ---------------------------------*/

#define	STAR_MAXNSED	8	/* Maximum number of galaxy SED components */

/*--------------------------------- flags -----------------------------------*/

/*-------------------------- structure definitions --------------------------*/

typedef struct
  {
  lfstruct	*lf;	/* Type-dependent luminosity function */
  sedstruct	*sed;	/* Spectral Energy Distribution */
  } starstruct;

/*-------------------------------- protos -----------------------------------*/

starstruct	*star_init(sedstruct *sed);

void		star_end(starstruct *star);

#endif
