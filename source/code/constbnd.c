#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include <assert.h>

/* Code which reads in 'constbnd.dat',  which gives full-fat data on
constellation boundaries,  and emits a compacted subset dedicated to
memory and CPU-efficient computation of the constellation
corresponding to a given RA/dec.  See 'constbnd.txt' for a
description of these boundaries and a few alternative methods for
solving the point-in-constellation problem.

   This code needs to work only once,  to create the compacted boundary
list found in 'conbound.c' (which means this file is now of
intellectual/documentation use only.)

   This program _only_ generates the compacted subset.  That list
appears in,  and is used in,  'conbound.c' for the actual function
to determine which constellation a given RA/dec is in.

   The input constellation boundary data are also available at

https://cdsarc.cds.unistra.fr/ftp/VI/49/

   and similar output data are available at

https://cdsarc.cds.unistra.fr/ftp/VI/42/

   An article describing essentially the same algorithm is at

https://iopscience.iop.org/article/10.1086/132034/pdf

   For determining the constellation in which a given RA/dec falls,  we
only need a set of "horizontal" (east/west) boundary lines. Each such
line is stored twice (once for the constellation to the north, once for
the constellation to the south).  We only keep the former.  Determining
which constellation we're in requires us to look north of our desired
point until we hit a boundary.  Put more algorithmically,  we :

(1) Binary-search the list to find the border segment to the north of
the given SPD.
(2) Check to see if that segment's RA range covers the given RA.  If it
does,  we've found the constellation.
(3) If it doesn't,  look at the preceding (more northerly) segment and
go back to (2).

For some areas in Ursa Minor,  there is no segment to the north,  and
we'll just loop our way to the top of the list and say "we must be in
UMi".  (In hindsight,  I could have added a dummy segment at dec +91
running from RA 0 to 24.  But I didn't,  and it's a tiny difference.)
Some small bits of added logic are also required to handle the "prime
meridian" at RA=0/24h.

   The resulting list is sorted in order of decreasing south polar
distance (we use SPD instead of dec so that all quantities will be
positive), and duplicates removed.

   There are a few nuances in how the boundaries are extracted.
Consider the following (imaginary) constellation Joe :

                   |
  ---A--------B    |
     |        |    |
     H        C----D--E----
     |                |
     |   Joe          |
     |                |
 ----F-------+--------G----+
             |             |
   Alice     |    Bob

   For this constellation,  we would extract segments AB, CD,  and DE.
(The horizontal segments F+ and +G would be applied to constellations
Alice and Bob.)  CD and DE can be combined to make a single segment.

   Less obviously,  CD can be extended over to H;  i.e.,  if you're to
the south of HE,  you are within Joe.  This extension isn't necessary,
but (slightly) reduces the time required to get an identification.

   Next,  there is the problem of memory-efficient storage.  The borders
all have SPDs that fall on integer arcminutes and RAs that fall on
integer RA seconds.  Thus,  the SPD of a segment can fit in 14 bits and
the western (minimum) RA in 17 bits;  we can put both into a 32-bit
integer.  The difference between the east RA and western RA for a segment
can fit in a 16-bit integer.  As shown in 'conbound.c',  this lets us fit
everything into seven bytes per segment.  That will almost uniformly be
rounded up to eight bytes for alignment,  so the table requires 324*8 =
2592 bytes.  We could probably do even better,  but the code would get
much uglier.    */

typedef struct
   {
   int16_t spd;      /* in arcminutes */
   int32_t min_ra, max_ra;   /* in seconds */
   char constell_idx;      /* from 0 to 87 */
   } constbnd_t;

typedef struct
   {
   int32_t x, y;
   } point_t;

constbnd_t *bounds;
int n_bounds = 0;

const char *constell_names =
        "AndAntApsAqlAqrAraAriAurBooCaeCamCapCarCasCenCepCetChaCirCMaCMiCnc"
        "ColComCrACrBCrtCruCrvCVnCygDelDorDraEquEriForGemGruHerHorHyaHyiInd"
        "LacLeoLepLibLMiLupLynLyrMenMicMonMusNorOctOphOriPavPegPerPhePicPsA"
        "PscPupPyxRetSclScoSctSerSexSgeSgrTauTelTrATriTucUMaUMiVelVirVolVul";

static void dump_lines( const int n_pts, const point_t *p, const int constell_idx)
{
   int i, j;

   for( i = 0; i < n_pts - 1; i++)
      if( p[i].y == p[i + 1].y && p[i].x < p[i + 1].x)
         {
         int32_t max_ra = 49 * 3600, min_ra = -max_ra;
         int32_t spd0 = p[i].y, ra0 = p[i].x;

         for( j = 0; j < n_pts - 1; j++)
            if( p[j].x == p[j + 1].x)   /* vertical line */
               {
               int32_t ra1 = p[j].x;

               while( ra1 + 43200 < ra0)
                  ra1 += 86400;
               while( ra1 - 43200 > ra0)
                  ra1 -= 86400;
               if( p[j].y < spd0 && p[j + 1].y >= spd0)     /* north heading */
                  if( ra1 <= ra0 && ra1 > min_ra)
                     min_ra = ra1;
               if( p[j].y >= spd0 && p[j + 1].y < spd0)     /* south heading */
                  if( ra1 > ra0 && ra1 < max_ra)
                     max_ra = ra1;
               }
         if( min_ra >= 86400)
            {
            min_ra -= 86400;
            max_ra -= 86400;
            }
         if( min_ra < 0)
            {
            min_ra += 86400;
            max_ra += 86400;
            }
         if( spd0 == 300)      /* dec = -85; Octans; full circle */
            {                  /* A few special cases here... */
            min_ra = 0;
            max_ra = 24 * 60 * 60;
            }
         if( spd0 == 450)      /* dec = -82.5; Octans; 7h40m to 3h30m (=27h30m) */
            {
            min_ra = (7 * 60 + 40) * 60;
            max_ra = (27 * 60 + 30) * 60;
            }
         bounds[n_bounds].spd = spd0;
         bounds[n_bounds].min_ra = min_ra;
         bounds[n_bounds].max_ra = max_ra;
         bounds[n_bounds].constell_idx = (char)constell_idx;
         if( max_ra - min_ra > 65000)     /* long spans require splitting */
            {
            bounds[n_bounds + 1] = bounds[n_bounds];
            bounds[n_bounds + 1].min_ra = bounds[n_bounds].max_ra
                        = (min_ra + max_ra) / 2;
            n_bounds++;
            }
         n_bounds++;
         }
}

int main( const int argc, const char **argv)
{
   FILE *ifile = fopen( "constbnd.dat", "rb");
   char buff[100];
   int i, j;
   point_t *p = (point_t *)calloc( 1000, sizeof( point_t));

   bounds = (constbnd_t *)calloc( 2000, sizeof( constbnd_t));
   assert( ifile);
   while( fgets( buff, sizeof( buff), ifile))
      {
      char constell[4];
      int n_pts = 1, shift_24 = 0, constell_idx = -1;

      memcpy( constell, buff + 24, 3);
      constell[3] = '\0';
      for( i = 0; i < 88; i++)
         if( !memcmp( constell_names + i * 3, constell, 3))
            constell_idx = i;
      assert( constell_idx >= 0);
      p[0].x = (int32_t)( atof( buff) * 3600. + 0.5);
      p[0].y = (int32_t)( atof( buff + 12) * 60. + 5400.5);
      while( fgets( buff, sizeof( buff), ifile) && !memcmp( buff + 24, constell, 3))
         {
         int remove_pt = 0;
         point_t *tptr = p + n_pts;

         tptr->x = (int32_t)( atof( buff) * 3600. + 0.5);
         tptr->y = (int32_t)( atof( buff + 12) * 60. + 5400.5);
         if( tptr->x - tptr[-1].x > 43200)
            tptr->x -= 86400;
         if( tptr->x - tptr[-1].x < -43200)
            tptr->x += 86400;
         if( tptr->x < 0)
            shift_24 = 1;
         if( n_pts >= 2)
            {
            if( tptr[-2].x == tptr[-1].x && tptr[-1].x == tptr->x)
               remove_pt = 1;       /* all collinear horizontally */
            if( tptr[-2].y == tptr[-1].y && tptr[-1].y == tptr->y)
               remove_pt = 1;       /* all collinear vertically */
            }
         if( remove_pt)
            tptr[-1] = tptr[0];
         else
            n_pts++;
         }
      if( memcmp( constell, "Vul", 3))
         fseek( ifile, -strlen( buff), SEEK_CUR);
      printf( "Constell %s\n", constell);
      if( shift_24)
         for( i = 0; i < n_pts; i++)
            p[i].x += 86400;
      for( i = 0; i < n_pts; i++)
         printf( "%9.5f %+09.5f\n",
                  (double)p[i].x / 3600.,
                  (double)p[i].y / 60. - 90.);
      dump_lines( n_pts, p, constell_idx);
      }
   fclose( ifile);
   free( p);
   printf( "%d bounds found\n", n_bounds);
            /* Sort by decreasing dec.  No need to be efficient here. */
   for( i = 0; i < n_bounds - 1; i++)
      if( i >= 0)
         if( bounds[i + 1].spd > bounds[i].spd
                  || (bounds[i + 1].spd == bounds[i].spd && bounds[i + 1].min_ra > bounds[i].min_ra))
            {
            const constbnd_t temp = bounds[i + 1];

            bounds[i + 1] = bounds[i];
            bounds[i] = temp;
            i -= 2;
            }
   for( i = j = 1; i < n_bounds; i++)     /* remove duplicates */
      if( memcmp( bounds + i, bounds + i - 1, sizeof( bounds[0])))
         bounds[j++] = bounds[i];
   n_bounds = j;
   printf( "%d bounds found after removing duplicates\n", n_bounds);
           /* Output boundary list in human-readable form */
   for( i = 0; i < n_bounds; i++)
      printf( "Line at dec %+09.5f, RA %9.5f to %9.5f : %.3s\n",
                  (double)bounds[i].spd / 60. - 90.,
                  (double)bounds[i].min_ra / 3600.,
                  (double)bounds[i].max_ra / 3600.,
                  constell_names + bounds[i].constell_idx * 3);
           /* Output boundary list as a C array;  see conbound.c */
   for( i = 0; i < n_bounds; i++)
         printf( "    { 0x%08lx, 0x%04lx, %2d },   /* %.3s */\n",
                  (unsigned long)bounds[i].min_ra | ((unsigned long)bounds[i].spd << 17),
                  (unsigned long)(bounds[i].max_ra - bounds[i].min_ra),
                  (int)bounds[i].constell_idx,
                  constell_names + 3 * bounds[i].constell_idx);
   free( bounds);
   return( 0);
}
