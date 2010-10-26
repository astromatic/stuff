/*
*				stars.h
*
* Include filoe for stars.c.
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
