/**
* @file         stars.c
* @brief        Manage star field generation
* @date         24/09/2016
* @copyright
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*       This file part of:      Stuff
*
*       Copyright:              (C) 1999-2016 IAP/CNRS/UPMC
*
*       License:                GNU General Public License
*
*       Stuff is free software: you can redistribute it and/or modify
*       it under the terms of the GNU General Public License as published by
*       the Free Software Foundation, either version 3 of the License, or
*       (at your option) any later version.
*       Stuff is distributed in the hope that it will be useful,
*       but WITHOUT ANY WARRANTY; without even the implied warranty of
*       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*       GNU General Public License for more details.
*       You should have received a copy of the GNU General Public License
*       along with Stuff. If not, see <http://www.gnu.org/licenses/>.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifdef HAVE_CONFIG_H
#include	"config.h"
#endif

#ifdef HAVE_MATHIMF_H
#include <mathimf.h>
#else
#include <math.h>
#endif

#include <stdio.h>
#include <stdlib.h>

#include "define.h"
#include "globals.h"
#include "prefs.h"
#include "random.h"
#include "sed.h"
#include "stars.h"


/******* stars_add *******************************************************
PROTO	void stars_add(FILE **catfiles, int ncatfiles,
	               sedstruct **starseds, int nstarseds)
PURPOSE	Generate a simple star field.
INPUT	-.
OUTPUT	-.
 ***/
void	stars_add(FILE **catfiles, int ncatfiles,
		sedstruct **starseds, int nstarseds) {

   char			str[MAXCHAR];
   double		xrange,yrange, dnstars;
   int			nstars, step, gridflag;
#ifndef USE_THREADS
   static objstruct	obj={0};
   int			i, gridindex, ngrid;
#endif

/* Initialize the random generator */
  NFPRINTF(OUTPUT, "Initializing star random generator...")
  init_random(prefs.starseed);

  return;
}



