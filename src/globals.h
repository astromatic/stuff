/*
*				globals.h
*
* Global declarations
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	Stuff
*
*	Copyright:		(C) 1999-2021 IAP/CNRS/SorbonneU
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
*	Last modified:		08/01/2021
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#include	"types.h"

/*----------------------- miscellaneous variables ---------------------------*/
extern char		gstr[MAXCHAR];

/*------------------------------- functions ---------------------------------*/
extern	void	error(int, char *, char *),
		makeit(void),
		swapbytes(void *ptr, int nb, int n),
		warning(char *msg1, char *msg2);

extern double	counter_seconds(void);

