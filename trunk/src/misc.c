/*
*				misc.c
*
* Miscellaneous functions.
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


#ifdef HAVE_CONFIG_H
#include	"config.h"
#endif

#include	<ctype.h>
#include	<stdio.h>
#include	<stdlib.h>

#include	"define.h"
#include	"globals.h"


/********************************* error ************************************/
/*
I hope it will never be used!
*/
void	error(int num, char *msg1, char *msg2)
  {
  fprintf(stderr, "\n> %s%s\n\n",msg1,msg2);
  exit(num);
  }


/********************************* warning **********************************/
/*
Print a warning message on screen.
*/
void	warning(char *msg1, char *msg2)
  {
  fprintf(OUTPUT, "\n> WARNING: %s%s\n\n",msg1,msg2);
  return;
  }


/******************************* swapbytes **********************************/
/*
Swap bytes for doubles, longs and shorts (for DEC machines or PC for inst.).
*/
void    swapbytes(void *ptr, int nb, int n)
  {
   char	*cp,
	c;
   int	j;

  cp = (char *)ptr;

  if (nb&4)
    {
    for (j=n; j--; cp+=4)
      {
      c = cp[3];
      cp[3] = cp[0];
      cp[0] = c;
      c = cp[2];
      cp[2] = cp[1];
      cp[1] = c;
      }
    return;
    }

  if (nb&2)
    {
    for (j=n; j--; cp+=2)
      {
      c = cp[1];
      cp[1] = cp[0];
      cp[0] = c;
      }
    return;
    }

  if (nb&1)
    return;

  if (nb&8)
    {
    for (j=n; j--; cp+=8)
      {
      c = cp[7];
      cp[7] = cp[0];
      cp[0] = c;
      c = cp[6];
      cp[6] = cp[1];
      cp[1] = c;
      c = cp[5];
      cp[5] = cp[2];
      cp[2] = c;
      c = cp[4];
      cp[4] = cp[3];
      cp[3] = c;
      }
    return;
    }

  error(EXIT_FAILURE, "*Internal Error*: Unknown size in ", "swapbytes()");

  return;
  }


