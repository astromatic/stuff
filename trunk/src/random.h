/*
 				random.h

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	SkyStuff
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	Definitions related to the generation of random numbers
*
*	Last modify:	29/10/97
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

/*--------------------------------- constants -------------------------------*/

#ifndef	RAND_MAX
#define	RAND_MAX	0x7fffffffL	/* default dynamic range of rand() */
#endif

/*-------------------------------- protos -----------------------------------*/

double	random_double(void),
	random_gauss(double sigma),
	random_poisson(double xm);

int	random_int(void);

void	init_random(int seed);


