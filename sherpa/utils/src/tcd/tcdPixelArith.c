/*_C_INSERT_SAO_COPYRIGHT_HERE_(1998-2007)_*/
/*_C_INSERT_GPL_LICENSE_HERE_*/

/*H*****************************************************************
 * FILE NAME:  misc/tcdPixelArith.c
 *
 * DEVELOPMENT: tools
 *
 * DESCRIPTION:
 *
 * This file contains routines needed to convert from pixel positions
 * (ie array of pixel locations) to an offset in the data array.
 * 
 *
 * REVISION HISTORY:
 *
 * Ref. No.        Date
 ----------       -----
 preClearCase     30March1998
                  22Jul1998 -- reworked  tcdPixelToOffset 

 *
 H***************************************************************** */


#include "tcd.h"
#include "tcd_private.h"

/* 
   +------------------------------------------------
   +
   + Convert pixel location to array offset
   +
   +------------------------------------------------
   */
int tcdPixelToOffset(
		     long  nAxes,  /* i: number of data axes */
		     long *lAxes,  /* i: length of data axes */
		     long *origin, /* i: iorigin of data axes */
		     long *pixel,  /* i: pixel to convert */
		     long *offset  /* o: offset into array */
		     )
{

  long ii;


  if ( origin )
    {
      *offset = pixel[nAxes-1] + origin[nAxes-1];

      /* work backwards */
      for (ii=nAxes-2; ii >= 0; ii-- )
	{
	  *offset = *offset * lAxes[ii] + pixel[ii] + origin[ii];
	}
    }
  else
    {
      *offset = pixel[nAxes-1];
      
      for (ii=nAxes-2; ii>=0; ii-- )
	{
	  *offset = *offset * lAxes[ii] + pixel[ii];
	}
    }

  return( tcdSUCCESS );
}





/*
  +------------------------------------------------
  +
  + Routine to calculate pixel location given the 
  + array index(offset).
  +
  +------------------------------------------------
  */
int tcdOffsetToPixel(
		     long nAxes,   /* i: number of data axes */
		     long *lAxes,  /* i: length of data axes */
		     long *origin, /* i: origin of data */
		     long offset,  /* i: offset into data array */
		     long *pixel   /* o: returned pixel location */
		     )
{

  long rn;     /* residual for n'th iteration */
  long sf = 1; /* product of lenghts          */
  long ii;     /* loop varialbe               */


  /* determine size of data in n-1 dimensions */
  for ( ii=0;ii<nAxes-1; ii++) sf *= lAxes[ii];

  /* initialize iterations */
  rn = offset;

  /* determine pixel location in reverse order */
  for ( ii=nAxes-1; ii>0; ii-- )
    {
      pixel[ii] = rn / sf ;
      if (origin) pixel[ii] -= origin[ii];
      rn = rn % sf;
      sf /= lAxes[ii-1];
    }

  /* pixel in 0th dimension is just the residual */
  pixel[0]=rn;
  if (origin) pixel[0] -= origin[0];
 

  return(tcdSUCCESS);

}


