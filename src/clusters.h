/*
 				cluster.h

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	Stuff
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	Include file for cluster.c
*
*	Last modify:	07/06/2005
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

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
