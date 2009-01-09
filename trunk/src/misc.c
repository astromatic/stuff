 /*
 				misc.c

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	Stuff
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	Miscellaneous functions.
*
*	Last modify:	26/01/2005
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

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
void	swapbytes(void *ptr, int nb, int n)
  {
   char	*cp;
   int	j;

  cp = (char *)ptr;

  if (nb&4)
    {
    for (j=n; j--; cp+=4)
      {
      cp[0] ^= (cp[3]^=(cp[0]^=cp[3]));
      cp[1] ^= (cp[2]^=(cp[1]^=cp[2]));
      }
    return;
    }

  if (nb&2)
    {
    for (j=n; j--; cp+=2)
      cp[0] ^= (cp[1]^=(cp[0]^=cp[1]));
    return;
    }

  if (nb&1)
    return;

  if (nb&8)
    {
    for (j=n; j--; cp+=8)
      {
      cp[0] ^= (cp[7]^=(cp[0]^=cp[7]));
      cp[1] ^= (cp[6]^=(cp[1]^=cp[6]));
      cp[2] ^= (cp[5]^=(cp[2]^=cp[5]));
      cp[3] ^= (cp[4]^=(cp[3]^=cp[4]));
      }
    return;
    }

  error(EXIT_FAILURE, "*Internal Error*: Unknown size in ", "swapbytes()");

  return;
  }

