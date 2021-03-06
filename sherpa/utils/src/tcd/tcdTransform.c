/*_C_INSERT_SAO_COPYRIGHT_HERE_(1998-2007)_*/
/*_C_INSERT_GPL_LICENSE_HERE_*/

/*H*****************************************************************
 * FILE NAME:  transform/tcdTransform.c
 *
 * DEVELOPMENT: tools
 *
 * DESCRIPTION:
 *
 * Compute an N-dimensional transform of input data array
 * 
 *
 * REVISION HISTORY:
 *
 * Ref. No.        Date
 ----------       -----
 preClearCase     30March1998

 *
 H***************************************************************** */

#include <string.h>
#include "tcd.h"
#include "tcd_private.h"
 

#include "fftw3.h"

/*
  +----------------------------------------------------
  +
  + Compute various N-D transforms
  +
  +----------------------------------------------------
  */
int tcdTransform(
		 tcdTRANSFORM  tType,  /* i: which transform to compute     */
		 float        *params, /* i: transform parameters(direction)*/
		 tcdComplex   *data,   /* i/o: data to transform- in place  */
		 long          nAxes,  /* i: number of axes                 */
		 long         *lAxes,  /* i: length of axes                 */
		 long         *dOrigin /* i: origin of axes                 */
		 )
{
  long nTotal, ii;
  int status;
  int *axe_len;   /* FFTW wants ints for this parameter */
  fftwf_plan plan;


  /* error checking */
  status = tcdCheckData( data, nAxes, lAxes );
  if ( status != tcdSUCCESS ) return( status );

  if ( params == NULL ) return( tcdERROR_NULLPTR);

  /* do transforms */
  switch ( tType )
    {
    case tcdFFT:
       axe_len = (int*)calloc(nAxes,sizeof(int));
      for (ii=0;ii<nAxes;ii++) axe_len[ii] = lAxes[nAxes-ii-1];

      if (params[0] == tcdFORWARD)
	plan=fftwf_plan_dft(nAxes, axe_len, (void *)data, (void *)data,       
			       FFTW_FORWARD, FFTW_ESTIMATE);
      else
	plan=fftwf_plan_dft(nAxes, axe_len, (void *)data, (void *)data,       
				FFTW_BACKWARD, FFTW_ESTIMATE);

      free(axe_len);

      if (plan == NULL)
	{
	  return(tcdERROR);
	}
      fftwf_execute( plan); 

            
      /* Normalize */
      if ( params[0] == (float )tcdFORWARD )
	{
	  nTotal=1;
	  for (ii=0;ii<nAxes  ; ii++) nTotal *= lAxes[ii];
	  for (ii=0; ii<nTotal; ii++) 
	    {
	      data[ii].r /= nTotal ; 
	      data[ii].i /= nTotal; 
	    }
	}

      /* destroy plan */
      fftwf_destroy_plan(plan);
      break;

    default:
      return( tcdERROR_NOTIMPLEMENTED );
    }

  return( tcdSUCCESS );

}
  /* double precision */
int tcdTransformD(
		 tcdTRANSFORM  tType,  /* i: which transform to compute     */
		 double        *params, /* i: transform parameters(direction)*/
		 tcdDComplex   *data,   /* i/o: data to transform- in place  */
		 long          nAxes,  /* i: number of axes                 */
		 long         *lAxes,  /* i: length of axes                 */
		 long         *dOrigin /* i: origin of axes                 */
		 )
{
  long nTotal, ii;
  int status;
  int *axe_len;
  fftw_plan plan;


  /* error checking */
  status = tcdCheckData( data, nAxes, lAxes );
  if ( status != tcdSUCCESS ) return( status );

  if ( params == NULL ) return( tcdERROR_NULLPTR);

  /* do transforms */
  switch ( tType )
    {
    case tcdFFT:
      axe_len = (int*)calloc(nAxes,sizeof(int));
      for (ii=0;ii<nAxes;ii++) axe_len[ii] = lAxes[nAxes-ii-1];

      if (params[0] == tcdFORWARD)
	plan=fftw_plan_dft(nAxes, axe_len, (void *)data, (void *)data,       
			       FFTW_FORWARD, FFTW_ESTIMATE);
      else
	plan=fftw_plan_dft(nAxes, axe_len,(void *)data, (void *)data,       
				FFTW_BACKWARD, FFTW_ESTIMATE);

      free(axe_len);

      if (plan == NULL)
	{
	  return(tcdERROR);
	}
      fftw_execute( plan); 

            
      /* Normalize */
      if ( params[0] == (float )tcdFORWARD )
	{
	  nTotal=1;
	  for (ii=0;ii<nAxes  ; ii++) nTotal *= lAxes[ii];
	  for (ii=0; ii<nTotal; ii++) 
	    {
	      data[ii].r /= nTotal ; 
	      data[ii].i /= nTotal; 
	    }
	}

      /* destroy plan */
      fftw_destroy_plan(plan);
      break;

    default:
      return( tcdERROR_NOTIMPLEMENTED );
    }

  return( tcdSUCCESS );

}

