/**
* @file         celsys.c
* @brief        Convert local angular to equatorial coordinates and vice versa
* @date         23/09/2016
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

#include <stdlib.h>

#include "define.h"
#include "celsys.h"

/******* celsys_init *******************************************************
PROTO	celsysstruct *celsys_init(double *pos)
PURPOSE	Initialize Equatorial <=> Celestial coordinate system transforms.
INPUT	Equatorial coordinates of the celestial center coordinates.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	23/09/2016
 ***/
celsysstruct	*celsys_init(double *pos) {
   celsysstruct	*celsys;
   double	*mat,
		a0,d0,ap,dp,ap2,y;
   int		s,lng,lat;


  if (!(celsys = calloc(1, sizeof(celsysstruct))))
    return NULL;

// Pole
  ap = pos[0] * DEG;
  dp = pos[1] * DEG;
// Origin
  a0 = fmod(ap + 90.0, 360.0) * DEG;
  d0 = 0.0;

// First compute in the output referential the longitude of the south pole
  y = sin(ap - a0);
  ap2 = asin(cos(d0)*y) ;

// Equatorial <=> Celestial System transformation parameters
  mat = celsys->mat;
  mat[0] = ap;
  mat[1] = ap2;
  mat[2] = cos(dp);
  mat[3] = sin(dp);

  return celsys;
}


/******* celsys_end *******************************************************
PROTO	void celsys_end(celsysstruct *celsys)
PURPOSE	Initialize Equatorial <=> Celestial coordinate system transforms.
INPUT	Equatorial coordinates of the celestial center coordinates.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	23/09/2016
 ***/
void	celsys_end(celsysstruct *celsys) {

  free(celsys);
  return;
}


/******* celsys_to_eq *********************************************************
PROTO	int celsys_to_eq(celsysstruct *celsys, double *pos)
PURPOSE	Convert arbitrary celestial coordinates to equatorial.
INPUT	WCS structure,
	Coordinate vector.
OUTPUT	RETURN_OK if mapping successful, RETURN_ERROR otherwise.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	23/09/2016
 ***/
int	celsys_to_eq(celsysstruct *celsys, double *pos) {

   double	*mat,
		a2,d2,sd2,cd2cp,sd,x,y;
   int		lng, lat;

  mat = celsys->mat;
  a2 = pos[0]*DEG - mat[1];
  d2 = pos[1]*DEG;
/* A bit of spherical trigonometry... */
/* Compute the latitude... */
  sd2 = sin(d2);
  cd2cp = cos(d2)*mat[2];
  sd = sd2*mat[3]-cd2cp*cos(a2);
/* ...and the longitude */
  y = cd2cp*sin(a2);
  x = sd2 - sd*mat[3];
  pos[0] = fmod((atan2(y,x) + mat[0])/DEG+360.0, 360.0);
  pos[1] = asin(sd)/DEG;

  return RETURN_OK;
}


/******* eq_to_celsys *********************************************************
PROTO	int eq_to_celsys(celsysstruct *celsys, double *pos)
PURPOSE	Convert equatorial to arbitrary celestial coordinates.
INPUT	WCS structure,
	Coordinate vector.
OUTPUT	RETURN_OK if mapping successful, RETURN_ERROR otherwise.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	23/09/2016
 ***/
int	eq_to_celsys(celsysstruct *celsys, double *pos) {

   double	*mat,
		a,d,sd2,cdcp,sd,x,y;
   int		lng, lat;

  mat = celsys->mat;
  a = pos[0]*DEG - mat[0];
  d = pos[1]*DEG;
/* A bit of spherical trigonometry... */
/* Compute the latitude... */
  sd = sin(d);
  cdcp = cos(d)*mat[2];
  sd2 = sd*mat[3]+cdcp*cos(a);
/* ...and the longitude */
  y = cdcp*sin(a);
  x = sd2*mat[3]-sd;
  pos[0] = fmod((atan2(y,x) + mat[1])/DEG+360.0, 360.0);
  pos[1] = asin(sd2)/DEG;

  return RETURN_OK;
}


