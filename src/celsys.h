/**
* @file         celsys.h
* @brief        Include file for celsys.c.
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

#ifndef _CELSYS_H_
#define _CELSYS_H_

//----------------------------- Internal constants ----------------------------
//--------------------------- structure definitions ---------------------------
typedef struct {
  double	mat[4];
} celsysstruct;


//------------------------------- functions -----------------------------------
celsysstruct	*celsys_init(double *pos);

int		celsys_to_eq(celsysstruct *celsys, double *pos),
		eq_to_celsys(celsysstruct *celsys, double *pos);

void		celsys_end(celsysstruct *celsys);

#endif

