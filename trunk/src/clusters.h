/*
*				clusters.h
*
* Include file for clusters.c.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	Stuff
*
*	Copyright:		(C) 2005-2010 Emmanuel Bertin -- IAP/CNRS/UPMC
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

#ifndef _GALAXIES_H_
#include "galaxies.h"
#endif

#ifndef	_CLUSTERS_H_
#define	_CLUSTERS_H_

/*----------------------------- Internal constants --------------------------*/

#define		CLUSTER_LISTNINC	64      /* default list increment */

/*-------------------------- structure definitions --------------------------*/

typedef struct
  {
  double	x,y;		/* Pixel coordinates of cluster's center */
  double	z;		/* Cluster redshift */
  double	mabs;		/* Absolute cluster mag. in the ref. passband */
  double	beta;		/* Exponent of the generalized King profile */
  double	rc;		/* Cluster core radius (m) */
  double	rmax;		/* Cluster max (cut-off) radius (m) */
  double	sigma_v;	/* Radial velocity dispersion */
  double	lf_eqvol;	/* Equivalent field LF volume */
  int		real_ngal;	/* Number of galaxies in realization */
  double	real_mabs;	/* Absolute cluster mag. in realization */
  } clusterstruct;


/*-------------------------------- protos -----------------------------------*/

double		cluster_normdens(clusterstruct *cluster,
			galtypestruct **galtype, int ngaltype);
int		cluster_read(char *listname, clusterstruct ***clusters,
			int *ncluster);

void	cluster_rndgalpos(clusterstruct *cluster,
			double *x, double *y, double *z),
	cluster_sortz(clusterstruct **clusters, int ncluster),
	cluster_write(char *listname, clusterstruct **clusters, int ncluster);

#endif
