 /*
 				globals.h

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	Stuff
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	Global declarations.
*
*	Last modify:	23/05/2010
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#include	"types.h"

/*----------------------- miscellaneous variables ---------------------------*/
char		gstr[MAXCHAR];
int		bswapflag;

/*------------------------------- functions ---------------------------------*/
extern	void	error(int, char *, char *),
		makeit(void),
		swapbytes(void *ptr, int nb, int n),
		warning(char *msg1, char *msg2);

extern double	counter_seconds(void);

