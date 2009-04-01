 /*
 				main.c

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	Stuff
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	Parsing of the command line.
*
*	Last modify:	01/04/2009
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifdef HAVE_CONFIG_H
#include	"config.h"
#endif

#include	<ctype.h>
#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>

#ifdef HAVE_PLPLOT
#include <plplot.h>
#endif

#include	"define.h"
#include	"globals.h"
#include	"prefs.h"

#define		SYNTAX \
BANNER " [-c <config_file>][-<keyword> <value>]\n" \
"> to dump a default configuration file: " BANNER " -d \n" \
"> to dump a default extended configuration file: " BANNER " -dd \n"

extern const char	notokstr[];

/********************************** main ************************************/
int	main(int argc, char *argv[])

  {
   static char	prefsname[MAXCHAR];
   char		**argkey, **argval;
   int		a, narg, opt,opt2;

  if (argc<1)
    {
    fprintf(OUTPUT, "\n         %s  version %s (%s)\n", BANNER,MYVERSION,DATE);
    fprintf(OUTPUT, "\nby %s\n", COPYRIGHT);
    fprintf(OUTPUT, "visit %s\n", WEBSITE);
    error(EXIT_SUCCESS, "SYNTAX: ", SYNTAX);
    }

#ifdef HAVE_PLPLOT
  if (argc>2)
    plParseOpts(&argc, argv, PL_PARSE_SKIP);
#endif

  QMALLOC(argkey, char *, argc);
  QMALLOC(argval, char *, argc);

/*default parameters */
  strcpy(prefsname, "stuff.conf");
  narg = 0;

  for (a=1; a<argc; a++)
    {
    if (*(argv[a]) == '-')
      {
      opt = (int)argv[a][1];
      if (strlen(argv[a])<4 || opt == '-')
        {
        opt2 = (int)tolower((int)argv[a][2]);
        if (opt == '-')
          {
          opt = opt2;
          opt2 = (int)tolower((int)argv[a][3]);
          }
        switch(opt)
          {
          case 'c':
            if (a<(argc-1))
              strcpy(prefsname, argv[++a]);
            break;
          case 'd':
            dumpprefs(opt2=='d' ? 1 : 0);
            exit(EXIT_SUCCESS);
            break;
          case 'v':
            printf("%s version %s (%s)\n", BANNER,MYVERSION,DATE);
            exit(EXIT_SUCCESS);
            break;
          case 'h':
            fprintf(OUTPUT, "\nSYNTAX: %s", SYNTAX);
#ifdef HAVE_PLPLOT
            fprintf(OUTPUT, "\nPLPLOT-specific options:\n");
            plParseOpts(&argc, argv, PL_PARSE_SKIP);
#endif
            break;
          default:
            error(EXIT_SUCCESS, "SYNTAX: ", SYNTAX);
          }
        }
      else
        {
        argkey[narg] = &argv[a][1];
        argval[narg++] = argv[++a];
        }       
      }
    else
      error(EXIT_SUCCESS,"SYNTAX: ", SYNTAX);
    }
  readprefs(prefsname, argkey, argval, narg);
  useprefs();

  free(argkey);
  free(argval);

  makeit();

  endprefs();

  NFPRINTF(OUTPUT, "");
  NPRINTF(OUTPUT, "> All done (in %d s)\n", prefs.time_diff);

  exit(EXIT_SUCCESS);
  }

