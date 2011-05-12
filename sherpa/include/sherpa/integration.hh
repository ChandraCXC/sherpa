//_C++_INSERT_SAO_COPYRIGHT_HERE_(2007)_
//_C++_INSERT_GPL_LICENSE_HERE_
#ifndef __sherpa_integration_hh__
#define __sherpa_integration_hh__

#include <sstream>

extern "C" {
  typedef double (*integrand_1d)( double x, void* params );
  typedef double (*integrand_Nd)( unsigned int ndim, const double* x,
				  void* params );
  typedef int (*integrand_1d_vec)( double* x, int len, void* params );
}


#if defined(_INTEGRATIONMODULE) || !defined(Py_PYTHON_H)

namespace sherpa { namespace integration {

int integrate_1d( integrand_1d fct, void* params,
		  double xlo, double xhi,
		  unsigned int maxeval, double epsabs, double epsrel,
		  double& result, double& abserr );


int integrate_Nd( integrand_Nd fct, void* params,
		  unsigned int ndim, const double* xlo, const double* xhi,
		  unsigned int maxeval, double epsabs, double epsrel,
		  double& result, double& abserr );


int py_integrate_1d( integrand_1d_vec fct, void* params,
		     const double xlo, const double xhi,
		     unsigned int maxeval, double epsabs, double epsrel,
		     double &result, double& abserr, int errflag,
		     std::ostringstream& err);


}  }  /* namespace integration, namespace sherpa */

#endif  /* defined(_INTEGRATIONMODULE) || !defined(Py_PYTHON_H) */


#if !defined(_INTEGRATIONMODULE) && defined(Py_PYTHON_H)


typedef int (*_integrate_1d)( integrand_1d fct,
			      void* params, double xlo, double xhi,
			      unsigned int maxeval,
			      double epsabs, double epsrel,
			      double& result, double& abserr );

typedef int (*_integrate_Nd)( integrand_Nd fct,
			      void* params, unsigned int ndim,
			      const double* xlo, const double* xhi,
			      unsigned int maxeval,
			      double epsabs, double epsrel,
			      double& result, double& abserr );

typedef int (*_py_integrate_1d)( integrand_1d_vec fct, void* params,
				 const double xlo, const double xhi,
				 unsigned int maxeval,
				 double epsabs, double epsrel,
				 double &result, double& abserr, int errflag,
				 std::ostringstream& err );

static void **Integration_API;

#define integrate_1d ((_integrate_1d)Integration_API[0])
#define integrate_Nd ((_integrate_Nd)Integration_API[1])
#define py_integrate_1d ((_py_integrate_1d)Integration_API[2])

static int
import_integration(void)
{

  PyObject *m = NULL;
  PyObject *api_cobject = NULL;
  int rv = -1;

  if ( NULL == 
       ( m = PyImport_ImportModule( (char*)"sherpa.utils.integration" ) ) )
    goto error;

  if ( NULL == ( api_cobject = PyObject_GetAttrString( m, (char*)"_C_API" ) ) )
    goto error;

  if ( NULL ==
       ( Integration_API = (void**)PyCObject_AsVoidPtr( api_cobject ) ) )
    goto error;

  rv = 0;

 error:
  Py_XDECREF(m);
  Py_XDECREF(api_cobject);

  return rv;

}


#endif  /* !defined(_INTEGRATIONMODULE) && defined(Py_PYTHON_H) */


#endif  /* __sherpa_integration_hh__ */
