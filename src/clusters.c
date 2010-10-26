/*
*				clusters.c
*
* Manage galaxy clusters.
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

#ifdef HAVE_CONFIG_H
#include	"config.h"
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "define.h"
#include "globals.h"
#include "clusters.h"
#include "cosmo.h"
#include "galaxies.h"
#include "prefs.h"
#include "random.h"

int     compclusterz(const void *sample1, const void *sample2);

/**************************** cluster_normdens ******************************/
/*
Do the cluster normalization in luminosity density by computing an
equivalent LF volume.
*/
double	cluster_normdens(clusterstruct *cluster, galtypestruct **galtype,
		int ngaltype)
  {
   double 	lum;
   int		g;

/* Compute the luminosity density from the composite luminosity functions */
  lum = 0.0;
  for (g=0; g<ngaltype; g++)
    lum += lf_intlumschechter(galtype[g]->lf, galtype[g]->lf->mabsmin,
		galtype[g]->lf->mabsmax);

/* Compute the density factor of the cluster */
  cluster->lf_eqvol = (lum > 0.0? exp(-0.921034*cluster->mabs) / lum : 0.0);

  return cluster->lf_eqvol;
  }


/**************************** cluster_rndgalpos *****************************/
/*
Return a random galaxy position in a cluster.
*/
void	cluster_rndgalpos(clusterstruct *cluster,
			double *x, double *y, double *z)
  {
   double	rc2,rmax2,beta, r,theta, da;

/* x and y */
  rc2 = cluster->rc*cluster->rc;
  rmax2 = cluster->rmax*cluster->rmax;
  beta = cluster->beta;
  r = (cluster->beta == 1.0)?
		  sqrt(rc2*(pow(rmax2/rc2+1.0, random_double()) - 1.0))
		: sqrt(rc2*(pow(1.0+(pow(rmax2/rc2+1.0, 1.0-beta) - 1.0)
			*random_double(),  1.0/(1.0 - beta))-1.0));
  theta = 2*PI*random_double();
/* Redshift */
  *z = cluster->z + random_gauss(cluster->sigma_v/C);
  da = cosmo_dlum(*z)/((1.0+*z)*(1.0+*z));
  *x = r*cos(theta) / da / ARCSEC;
  *y = r*sin(theta) / da / ARCSEC;
  return;
  }


/****** cluster_read *********************************************************
PROTO   int cluster_read(char *listname, clusterstruct ***clusters,
			int *ncluster)
PURPOSE Read a list of clusters
INPUT   List file name,
	Pointer (output) to an array of pointers to clusters,
	Pointer (output) to the number of clusters listed.
OUTPUT  -.
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 07/06/2005
*/
int	cluster_read(char *listname, clusterstruct ***clusters, int *ncluster)
  {
   clusterstruct	**pcluster,
			*cluster;
   FILE			*file;
   char			str[MAXCHAR],msg[80];
   int			c,i, cmax;

  if ((file = fopen(listname,"r")) == NULL)
    error(EXIT_FAILURE,"*ERROR*: can't read ", listname);

  c = cmax = 0;
  pcluster = NULL;		/* Avoid gcc -Wall warnings */
  for (i=0; fgets(str, MAXCHAR, file);)
    {
/*-- Examine current input line (discard empty and comment lines) */
    if (!*str || strchr("#\t\n",*str))
      continue;
    if (!(++i%1))
      {
      sprintf(msg, "Reading cluster list... (%d objects)", i);
      NFPRINTF(OUTPUT, msg);
      }

    if (!c)
      {
      cmax = CLUSTER_LISTNINC;
      QMALLOC(pcluster, clusterstruct *, cmax);
      }
    else if (c >= cmax)
      {
      cmax += CLUSTER_LISTNINC;
      QREALLOC(pcluster, clusterstruct *, cmax);
      }
    QCALLOC(cluster, clusterstruct, 1);
    pcluster[c++] = cluster;
    sscanf(str, "%lf %lf %lf %lf %lf %lf %lf %lf",
	&cluster->x, &cluster->y, &cluster->z, &cluster->mabs, &cluster->beta,
	&cluster->rc, &cluster->rmax, &cluster->sigma_v);
    cluster->rc *= MPC/H;
    cluster->rmax *= MPC/H;
    cluster->sigma_v *= KM/H;
    cluster->real_ngal = 0;
    cluster->real_mabs = 99.0;
    }

  fclose(file);

/* Save memory */
  QREALLOC(pcluster, clusterstruct *, c);
  *clusters = pcluster;
  *ncluster = c;

  return c;
  }


/****** cluster_write *********************************************************
PROTO   void cluster_write(char *listname, clusterstruct **clusters,
		int ncluster)
PURPOSE Write a list of clusters
INPUT   List file name,
	Pointer to an array of pointers to clusters,
	Number of clusters listed.
OUTPUT  -.
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 07/06/2005
*/
void	cluster_write(char *listname, clusterstruct **clusters, int ncluster)
  {
   clusterstruct	*cluster;
   FILE			*file;
   double		da;
   int			c;

  if ((file = fopen(listname,"w")) == NULL)
    error(EXIT_FAILURE,"*ERROR*: can't write to ", listname);

  for (c=0; c<ncluster; c++)
    {
    cluster = clusters[c];
    da = cosmo_dlum(cluster->z)/((1.0+cluster->z)*(1.0+cluster->z));
    fprintf(file,
	"%10.3f %10.3f %8.5f %+9.4f %4.2f %9g %9g %9g %+9.4f %9d"
	" %9g %9g\n",
	cluster->x, cluster->y, cluster->z,
	cluster->mabs, cluster->beta, cluster->rc*H/MPC, cluster->rmax*H/MPC,
	cluster->sigma_v*H/KM,
	cluster->real_mabs, cluster->real_ngal,
	cluster->rc/da/ARCSEC, cluster->rmax/da/ARCSEC);
    }

  fclose(file);

  return;
  }


/****** cluster_sortz *********************************************************
PROTO   void cluster_sortz(clusterstruct **clusters, int ncluster)
PURPOSE Sort a list of cluster by decreasing z.
INPUT   Array of cluster pointers
	Number of clusters.
OUTPUT  -.
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 07/06/2005
*/
void cluster_sortz(clusterstruct **clusters, int ncluster)
  {
  qsort(clusters, ncluster, sizeof(clusterstruct *), compclusterz);

  return;
  }


/****** compclustz ************************************************************
PROTO   int compclusterz(const void *cluster1, const void *cluster2)
PURPOSE Provide a cluster z comparison function for qsort().
INPUT   pointer to 1st cluster,
        pointer to 2nd cluster.
OUTPUT  <0 if z1>z2, >0 if z1<z2, 0 otherwise .
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 07/06/2005
*/
int     compclusterz(const void *sample1, const void *sample2)
  {
   double	dz;

  dz = ((clusterstruct *)sample2)->z - ((clusterstruct *)sample1)->z;

  return dz>0.0 ? 1: (dz<0.0? -1 : 0);
  }


