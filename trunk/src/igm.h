/*
*				igm.h
*
* Include file for igm.c.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	Stuff
*
*	Copyright:		(C) 2013 Emmanuel Bertin -- IAP/CNRS/UPMC
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
*	Last modified:		27/05/2013
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef _SED_H_
#include "sed.h"
#endif

#ifndef _IGM_H_
#define _IGM_H_

/*--------------------------------- constants -------------------------------*/

#define IGM_NLYMAN		17	/* Number of Lyman lines (see igm.c) */

#define	IGM_LYMAN_POWER		3.46	/* exponent from Madau et al. 1995 */

#define	IGM_METAL_POWER		1.68	/* exponent from Madau et al. 1995 */
#define	IGM_METAL_COEFF		0.0017	/* coefficient from Madau et al. 1995 */

#define IGM_WAVE_START		100.0	/* Starting lambda for IGM absorption */
#define IGM_WAVE_STEP		20.0	/* Wavelength step for IGM absorption */
#define IGM_WAVE_END		100000.0/* End wavelength for IGM absorption */

/*---------------------------------- enum -----------------------------------*/

typedef enum {IGM_NONE, IGM_MADAU_AVERAGE}	igmenum;

/*-------------------------------- protos -----------------------------------*/

sedstruct	*sed_igmmadauextinct(double z);

#endif
