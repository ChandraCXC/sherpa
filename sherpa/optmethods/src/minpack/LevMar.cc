#ifdef testLevMar

//_C++_INSERT_SAO_COPYRIGHT_HERE_(2007)_
//_C++_INSERT_GPL_LICENSE_HERE_

#include "LevMar.hh"

using namespace minpack;


#include "sherpa/fcmp.hh"
#include "sherpa/functor.hh"
#include "tests/tstoptfct.hh"

template< typename Real >
void print_pars( const char* name, int nfev, Real fval, Real answer,
                int n, const std::vector< Real >& x,
		 const std::vector< Real >& err,
                Real tol=1.0e4*std::sqrt( std::numeric_limits< Real >::epsilon()
 ),
                const char* prefix="lmdif_" ) {


  std::cout << prefix << name << '\t';
  if ( 0 == sao_fcmp( fval, answer, std::sqrt(tol) ) )
    std::cout << nfev << '\t';
  else
    std::cout << -nfev << '\t';
  std::cout << answer << '\t';
  std::cout << fval << '\t';
  std::cout << x[0];
  for ( int ii = 1; ii < n; ++ii )
    std::cout << ',' << x[ii];
  std::cout << '\t';
  std::cout << err[0];
  for ( int ii = 1; ii < n; ++ii )
    std::cout << ',' << err[ii];
  std::cout << '\n';

}

template< typename Init, typename Fct >
void justdoit( Init init, Fct fct, int npars, 
	       std::vector<double>& pars, std::vector<double>& lo,
	       std::vector<double>& hi, std::vector<double>& covarerr,
	       double tol, const char* header ) {

  try {

      int mfcts;
      double answer;
      init( npars, mfcts, answer, &pars[0], &lo[0], &hi[0] );

      minpack::LevMar< Fct, void* > lm( fct, NULL, mfcts );

      int maxnfev=128*npars, nfev;
      double fmin=0.0;

      int nprint = 0;
      double epsfcn = 1.0e-8, factor=100.0;

      lm( npars, tol, tol, tol, maxnfev, epsfcn, factor, nprint, lo, hi, pars,
	  nfev, fmin, covarerr );

      // lm.minimize( &pars[0], tol, maxnfev, nfev, fmin );
      
      print_pars( header, nfev, fmin, answer, npars, pars, covarerr );
      
  } catch( const sherpa::OptErr& oe ) {
    
    std::cerr << oe << '\n';
    
  }

  return;

}

int main( int argc, char* argv[] ) {

  int npars=16;
  if ( argc == 2 )
    npars = atoi( argv[1] );

  if ( npars % 2 || npars < 2 ) {
    printf( "The minimum value for the free parameter must be an even "
            "and it is greater then 2\n" );
    return EXIT_FAILURE;
  }

  std::cout << "#\n# sizeof(double) = " << sizeof( double ) << "\n";
  std::cout << "# A negative value for the nfev signifies that the "
    "optimization method did not converge\n#\n";
  std::cout << "name\tnfev\tanswer\tfval\tpars\terr\nS\tN\tN\tN\tN\tN\n";

  double tol = std::sqrt( std::numeric_limits< double >::epsilon() );

  // The factor of 8 for the parameter is because the Chebyquad has 9
  // free parameters and npars has to be at least 2.
  // The factor of 32 for fvec is because the number of functions may
  // be greater then the number of free parameters.
  std::vector< double > pars( 8 * npars ), lo( 8 * npars ), hi( 8 * npars ),
    covarerr( 8 * npars );

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::Rosenbrock<double,void*> );
    
    justdoit( sherpa::fct_ptr( tstoptfct::RosenbrockInit<double> ),
	      fct, npars, pars, lo, hi, covarerr, tol, "Rosenbrock" );
  }

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::FreudensteinRoth<double,void*> );

    justdoit( sherpa::fct_ptr( tstoptfct::FreudensteinRothInit<double> ),
	      fct, npars, pars, lo, hi, covarerr, tol, "FreudensteinRoth" );
  }

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::PowellBadlyScaled<double,void*> );

    justdoit( sherpa::fct_ptr( tstoptfct::PowellBadlyScaledInit<double> ),
	      fct, npars, pars, lo, hi, covarerr, tol, "PowellBadlyScaled" );
  }

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::BrownBadlyScaled<double,void*> );

    justdoit( sherpa::fct_ptr( tstoptfct::BrownBadlyScaledInit<double> ),
	      fct, 2, pars, lo, hi, covarerr, tol, "BrownBadlyScaled" );
  }

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::Beale<double,void*> );

    justdoit( sherpa::fct_ptr( tstoptfct::BealeInit<double> ),
	      fct, npars, pars, lo, hi, covarerr, tol, "Beale" );
  }

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::JennrichSampson<double,void*> );

    justdoit( sherpa::fct_ptr( tstoptfct::JennrichSampsonInit<double> ),
	      fct, npars, pars, lo, hi, covarerr, tol, "JennrichSampson" );
  }

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::HelicalValley<double,void*> );

    justdoit( sherpa::fct_ptr( tstoptfct::HelicalValleyInit<double> ),
	      fct, 3*npars, pars, lo, hi, covarerr, tol, "HelicalValley" );
  }

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::Bard<double,void*> );

    justdoit( sherpa::fct_ptr( tstoptfct::BardInit<double> ),
	      fct, 3*npars, pars, lo, hi, covarerr, tol, "Bard" );
  }

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::Gaussian<double,void*> );

    justdoit( sherpa::fct_ptr( tstoptfct::GaussianInit<double> ),
	      fct, 3, pars, lo, hi, covarerr, tol, "Gaussian" );
  }

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::Meyer<double,void*> );

    justdoit( sherpa::fct_ptr( tstoptfct::MeyerInit<double> ),
	      fct, 3, pars, lo, hi, covarerr, tol, "Meyer" );
  }

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::GulfResearchDevelopment<double,void*> );

    justdoit( sherpa::fct_ptr( tstoptfct::GulfResearchDevelopmentInit<double> ),
	      fct, 3, pars, lo, hi, covarerr, tol, "GulfResearchDevelopment" );
  }

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::Box3d<double,void*> );

    justdoit( sherpa::fct_ptr( tstoptfct::Box3dInit<double> ),
	      fct, 3, pars, lo, hi, covarerr, tol, "Box3d" );
  }

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::PowellSingular<double,void*> );

    justdoit( sherpa::fct_ptr( tstoptfct::PowellSingularInit<double> ),
	      fct, 4*npars, pars, lo, hi, covarerr, tol, "PowellSingular" );
  }

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::Wood<double,void*> );

    justdoit( sherpa::fct_ptr( tstoptfct::WoodInit<double> ),
	      fct, 4, pars, lo, hi, covarerr, tol, "Wood" );
  }

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::KowalikOsborne<double,void*> );

    justdoit( sherpa::fct_ptr( tstoptfct::KowalikOsborneInit<double> ),
	      fct, 4, pars, lo, hi, covarerr, tol, "KowalikOsborne" );
  }

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::BrownDennis<double,void*> );

    justdoit( sherpa::fct_ptr( tstoptfct::BrownDennisInit<double> ),
	      fct, 4, pars, lo, hi, covarerr, tol, "BrownDennis" );
  }

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::Osborne1<double,void*> );

    justdoit( sherpa::fct_ptr( tstoptfct::Osborne1Init<double> ),
	      fct, 5, pars, lo, hi, covarerr, tol, "Osborne1" );
  }

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::Biggs<double,void*> );

    justdoit( sherpa::fct_ptr( tstoptfct::BiggsInit<double> ),
	      fct, 6, pars, lo, hi, covarerr, tol, "Biggs" );
  }

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::Osborne2<double,void*> );

    justdoit( sherpa::fct_ptr( tstoptfct::Osborne2Init<double> ),
	      fct, 11, pars, lo, hi, covarerr, tol, "Osborne2" );
  }

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::Watson<double,void*> );

    justdoit( sherpa::fct_ptr( tstoptfct::WatsonInit<double> ),
	      fct, 6, pars, lo, hi, covarerr, tol, "Watson" );
  }

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::PenaltyI<double,void*> );

    justdoit( sherpa::fct_ptr( tstoptfct::PenaltyIInit<double> ),
	      fct, 4, pars, lo, hi, covarerr, tol, "PenaltyI" );
  }

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::PenaltyII<double,void*> );

    justdoit( sherpa::fct_ptr( tstoptfct::PenaltyIIInit<double> ),
	      fct, 4, pars, lo, hi, covarerr, tol, "PenaltyII" );
  }

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::VariablyDimensioned<double,void*> );

    justdoit( sherpa::fct_ptr( tstoptfct::VariablyDimensionedInit<double> ),
	      fct, npars, pars, lo, hi, covarerr, tol, "VariablyDimensioned" );
  }

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::Trigonometric<double,void*> );

    justdoit( sherpa::fct_ptr( tstoptfct::TrigonometricInit<double> ),
	      fct, npars, pars, lo, hi, covarerr, tol, "Trigonometric" );
  }

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::BrownAlmostLinear<double,void*> );

    justdoit( sherpa::fct_ptr( tstoptfct::BrownAlmostLinearInit<double> ),
	      fct, npars, pars, lo, hi, covarerr, tol, "BrownAlmostLinear" );
  }

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::DiscreteBoundary<double,void*> );

    justdoit( sherpa::fct_ptr( tstoptfct::DiscreteBoundaryInit<double> ),
	      fct, npars, pars, lo, hi, covarerr, tol, "DiscreteBoundary" );
  }

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::DiscreteIntegral<double,void*> );

    justdoit( sherpa::fct_ptr( tstoptfct::DiscreteIntegralInit<double> ),
	      fct, npars, pars, lo, hi, covarerr, tol, "DiscreteIntegral" );
  }

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::BroydenTridiagonal<double,void*> );

    justdoit( sherpa::fct_ptr( tstoptfct::BroydenTridiagonalInit<double> ),
	      fct, npars, pars, lo, hi, covarerr, tol, "BroydenTridiagonal" );
  }

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::BroydenBanded<double,void*> );

    justdoit( sherpa::fct_ptr( tstoptfct::BroydenBandedInit<double> ),
	      fct, npars, pars, lo, hi, covarerr, tol, "BroydenBanded" );
  }

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::LinearFullRank<double,void*> );

    justdoit( sherpa::fct_ptr( tstoptfct::LinearFullRankInit<double> ),
	      fct, npars, pars, lo, hi, covarerr, tol, "LinearFullRank" );
  }

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::LinearFullRank1<double,void*> );

    justdoit( sherpa::fct_ptr( tstoptfct::LinearFullRank1Init<double> ),
	      fct, npars, pars, lo, hi, covarerr, tol, "LinearFullRank1" );
  }

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::LinearFullRank0cols0rows<double,void*> );

    justdoit( sherpa::fct_ptr( tstoptfct::LinearFullRank0cols0rowsInit<double> ),
	      fct, npars, pars, lo, hi, covarerr, tol, "LinearFullRank0cols0rows" );
  }

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::Chebyquad<double,void*> );

    justdoit( sherpa::fct_ptr( tstoptfct::ChebyquadInit<double> ),
	      fct, 9, pars, lo, hi, covarerr, tol, "Chebyquad" );
  }

  return 0;
  
}

/*
==21378== Memcheck, a memory error detector
==21378== Copyright (C) 2002-2009, and GNU GPL'd, by Julian Seward et al.
==21378== Using Valgrind-3.5.0 and LibVEX; rerun with -h for copyright info
==21378== Command: tstlm
==21378==
#
# sizeof(double) = 8
# A negative value for the nfev signifies that the optimization method did not converge
#
name	nfev	answer	fval	pars	err
S	N	N	N	N	N
lmdif_Rosenbrock	278	0	0	1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1	1,2.0026,1,2.0026,1,2.0026,1,2.0026,1,2.0026,1,2.0026,1,2.0026,1,2.0026
lmdif_FreudensteinRoth	-142	0	391.874	11.4123,-0.896842,11.4123,-0.896842,11.4123,-0.896842,11.4123,-0.896842,11.4123,-0.896842,11.4123,-0.896842,11.4123,-0.896842,11.4123,-0.896842	26456.6,1977.23,26456.6,1977.23,26456.6,1977.23,26456.6,1977.23,26456.6,1977.23,26456.6,1977.23,26456.6,1977.23,26456.6,1977.23
lmdif_PowellBadlyScaled	274	0	5.93125e-29	1.09816e-05,9.10615,1.09816e-05,9.10615,1.09816e-05,9.10615,1.09816e-05,9.10615,1.09816e-05,9.10615,1.09816e-05,9.10615,1.09816e-05,9.10615,1.09816e-05,9.10615	1.09816e-05,0,1.09816e-05,0,1.09816e-05,0,1.09816e-05,0,1.09816e-05,0,1.09816e-05,0,1.09816e-05,0,1.09816e-05,0
lmdif_BrownBadlyScaled	46	0	0	1e+06,2e-06	1,1e-06
lmdif_Beale	121	0	7.16896e-26	3,0.5,3,0.5,3,0.5,3,0.5,3,0.5,3,0.5,3,0.5,3,0.5	2.4998,0.653891,2.4998,0.653891,2.4998,0.653891,2.4998,0.653891,2.4998,0.653891,2.4998,0.653891,2.4998,0.653891,2.4998,0.653891
lmdif_JennrichSampson	213	994.896	994.897	0.25782,0.257829,0.25782,0.257829,0.25782,0.257829,0.25782,0.257829,0.25782,0.257829,0.25782,0.257829,0.25782,0.257829,0.25782,0.257829	474.176,474.133,474.176,474.133,474.176,474.133,474.176,474.133,474.176,474.133,474.176,474.133,474.176,474.133,474.176,474.133
lmdif_HelicalValley	492	0	3.92334e-45	1,2.49347e-32,-5.02346e-48,1,-4.21539e-44,-6.54827e-50,1,-5.02471e-42,7.72904e-45,1,-1.83647e-44,-1.05628e-47,1,1.26566e-43,3.41847e-48,1,-2.90901e-43,9.56088e-51,1,6.22264e-44,-5.06329e-47,1,-7.25041e-38,1.05794e-37,1,-2.07677e-43,4.97802e-50,1,1.80308e-44,3.94337e-50,1,6.94956e-44,4.5886e-47,1,-7.29113e-45,1.21534e-53,1,9.6301e-44,5.23862e-49,1,-7.28894e-44,8.1491e-89,1,-8.54792e-44,5.82166e-74,1,-3.93557e-24,3.87657e-63	0.1,0.631452,1,0.1,0.631452,1,0.1,0.631452,1,0.1,0.631452,1,0.1,0.631452,1,0.1,0.631457,1,0.1,0.631452,1,0.1,0.631452,1,0.1,0.631452,1,0.1,0.631452,1,0.1,0.631452,1,0.1,0.631269,1,0.1,0.0628319,1,0.1,0.0628319,1,0.1,0.0628319,1,0.1,0.0628319,1
lmdif_Bard	246	0.131438	0.131438	0.0824106,1.13304,2.34369,0.0824106,1.13304,2.34369,0.0824106,1.13304,2.34369,0.0824106,1.13304,2.34369,0.0824106,1.13304,2.34369,0.0824106,1.13304,2.34369,0.0824106,1.13304,2.34369,0.0824106,1.13304,2.34369,0.0824106,1.13304,2.34369,0.0824106,1.13304,2.34369,0.0824106,1.13304,2.34369,0.0824106,1.13304,2.34369,0.0824106,1.13304,2.34369,0.0824106,1.13304,2.34369,0.0824106,1.13304,2.34369,0.0824106,1.13304,2.34369	0.47294,11.7694,11.3256,0.47294,11.7694,11.3256,0.47294,11.7694,11.3256,0.47294,11.7694,11.3256,0.47294,11.7694,11.3256,0.47294,11.7694,11.3256,0.47294,11.7694,11.3256,0.47294,11.7694,11.3256,0.47294,11.7694,11.3256,0.47294,11.7694,11.3256,0.47294,11.7694,11.3256,0.47294,11.7694,11.3256,0.47294,11.7694,11.3256,0.47294,11.7694,11.3256,0.47294,11.7694,11.3256,0.47294,11.7694,11.3256
lmdif_Gaussian	13	1.12793e-08	1.12793e-08	0.398956,1.00002,-2.24701e-13	0.650508,3.76589,0
lmdif_Meyer	-387	87.9458	802.357	0.00732639,5960.66,337.728	1.88816e-06,0,0.00676482
lmdif_GulfResearchDevelopment	386	0	3.85848e-08	5.05075,36.9509,0.969567	35901.4,39786,1775.51
lmdif_Box3d	29	0	1.99834e-30	1,10,1	2.16802,37.3447,1.82013
lmdif_PowellSingular	4959	0	2.03141e-65	1.37045e-17,-1.37045e-18,6.53483e-18,6.53483e-18,1.44713e-17,-1.44713e-18,6.83409e-18,6.83409e-18,1.4092e-17,-1.4092e-18,6.84807e-18,6.84807e-18,1.46237e-17,-1.46237e-18,6.91695e-18,6.91695e-18,1.41974e-17,-1.41974e-18,6.64082e-18,6.64082e-18,1.40335e-17,-1.40335e-18,6.56943e-18,6.56943e-18,1.42634e-17,-1.42634e-18,6.75154e-18,6.75154e-18,1.43927e-17,-1.43927e-18,6.73499e-18,6.73499e-18,1.43962e-17,-1.43962e-18,6.8149e-18,6.8149e-18,1.43962e-17,-1.43962e-18,6.8149e-18,6.8149e-18,1.43962e-17,-1.43962e-18,6.8149e-18,6.8149e-18,1.43962e-17,-1.43962e-18,6.8149e-18,6.8149e-18,1.43962e-17,-1.43962e-18,6.8149e-18,6.8149e-18,1.43962e-17,-1.43962e-18,6.8149e-18,6.8149e-18,1.43962e-17,-1.43962e-18,6.8149e-18,6.8149e-18,1.43962e-17,-1.43962e-18,6.8149e-18,6.8149e-18	0,0.1,0,0.447214,0,0.1,0.447214,0,0,0.1,0.447214,0,0,0.1,0.447214,0,0,0.1,0.447214,0,0,0.1,0.447214,0,0,0.1,0.447214,0,0,0.1,0.447214,0,0,0.1,0.447214,0,0,0.1,0.447214,0,0,0.1,0.447214,0,0,0.1,0.447214,0,0,0.1,0.447214,0,0,0.1,0.447214,0,0,0.1,0.447214,0,0,0.1,0.447214,0
lmdif_Wood	326	0	3.6586e-28	1,1,1,1	0.534248,1.06653,0.534362,1.06653
lmdif_KowalikOsborne	82	0.000307505	0.000307506	0.192808,0.191257,0.123052,0.136051	1.72541,29.6237,12.1973,13.5845
lmdif_BrownDennis	-513	85822.2	101279	-7.85984,11.9153,-0.538758,0.228339	0.0708236,0.022317,0.403789,0.595972
lmdif_Osborne1	93	5.46489e-05	5.46489e-05	0.37541,1.93582,-1.46466,0.0128675,0.0221228	1.48313,157.678,158.709,0.321146,0.6403
lmdif_Biggs	7	0	0	1,10,1,5,4,3	0,204.531,64.0246,274.365,302.265,211.655
lmdif_Osborne2	148	0.0401377	0.0401683	1.30997,0.431458,0.633631,0.599303,0.753912,0.905583,1.36503,4.8248,2.39882,4.56887,5.67537	0.675507,0.859784,0.481429,1.30504,1.88633,5.02027,6.71767,17.2939,1.14824,1.05995,0.599553
lmdif_Watson	-36	0.00228767	0.00260576	-1.31773e-22,1.01364,-0.244317,1.37382,-1.68569,1.09812	1,0.759045,5.37423,14.8655,16.8653,6.70762
lmdif_PenaltyI	-130	9.37629e-06	2.24998e-05	0.250001,0.250017,0.25001,0.250001	273.85,273.867,273.864,273.866
lmdif_PenaltyII	514	9.37629e-06	9.37902e-06	0.199999,0.215021,0.467087,0.514754	1,1757.21,1771.18,2236.35
lmdif_VariablyDimensioned	154	0	1.14323e-25	1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1	0.999665,0.998661,0.996986,0.994636,0.991608,0.987895,0.983491,0.978384,0.972565,0.96602,0.958733,0.950688,0.941863,0.932237,0.921783,0.910471
lmdif_Trigonometric	595	0	0	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0	1.00005,0.999251,0.998452,0.997656,0.99686,0.996066,0.995273,0.994481,0.99369,0.992901,0.992113,0.991326,0.99054,0.989756,0.988973,0.988191
lmdif_BrownAlmostLinear	54	1	1	-0.0556601,-0.0556601,-0.0556601,-0.0556601,-0.0556601,-0.0556601,-0.0556601,-0.0556601,-0.0556601,-0.0556601,-0.0556601,-0.0556601,-0.0556601,-0.0556601,-0.0556601,17.8906	0.966227,0.966227,0.966227,0.966227,0.966227,0.966227,0.966227,0.966227,0.966227,0.966227,0.966227,0.966227,0.966227,0.966227,0.966227,0
lmdif_DiscreteBoundary	69	0	3.37602e-34	-0.0462988,-0.0908017,-0.0889284,-0.086703,-0.0840814,-0.0810162,-0.0774561,-0.0733461,-0.0686267,-0.0632337,-0.0570977,-0.0501438,-0.0422906,-0.0334499,-0.0235258,-0.0124139	1.86639,3.60664,3.48577,3.36158,3.2335,3.10082,2.96264,2.81783,2.66494,2.50202,2.3264,2.13427,1.9198,1.67333,1.3762,0.981057
lmdif_DiscreteIntegral	69	0	1.06592e-32	-0.0284861,-0.0550797,-0.0795978,-0.101833,-0.121548,-0.138475,-0.152302,-0.162673,-0.169172,-0.171318,-0.168542,-0.160174,-0.145417,-0.123314,-0.0927079,-0.0521848	0.995075,0.990626,0.986601,0.982949,0.979634,0.976625,0.973909,0.971491,0.969403,0.967718,0.966569,0.966179,0.966918,0.969379,0.974523,0.983895
lmdif_BroydenTridiagonal	86	0	2.09653e-26	-0.580265,-0.707105,-0.707103,-0.707097,-0.707083,-0.70705,-0.70697,-0.706777,-0.70631,-0.705184,-0.702468,-0.69593,-0.680249,-0.642988,-0.556421,-0.366025	0.206484,0.227554,0.227554,0.227556,0.227561,0.227571,0.227597,0.227659,0.227811,0.228179,0.229078,0.231288,0.23669,0.249267,0.273266,0.288681
lmdif_BroydenBanded	103	0	3.54099e-24	-0.428303,-0.476596,-0.519652,-0.558099,-0.592506,-0.624504,-0.623239,-0.62142,-0.619616,-0.618226,-0.617518,-0.617732,-0.617901,-0.617982,-0.618896,-0.586311	0.210528,0.185077,0.165399,0.150092,0.137958,0.127797,0.128203,0.128801,0.129398,0.129856,0.13009,0.130021,0.129966,0.129941,0.129613,0.140125
lmdif_LinearFullRank	35	0	1.97215e-31	-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1	1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1
lmdif_LinearFullRank1	35	3.63636	3.63636	-4883.01,-822.335,516.085,-518.777,53.0791,-21.7553,435.962,157.453,-22.3739,11.9037,267.798,-109.554,-80.2355,230.792,-143.096,63.6255	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.0016159
lmdif_LinearFullRank0cols0rows	35	5.13793	5.13793	1,-2412.19,-251.123,262.525,490.5,-327.764,375.448,374.867,-15.7958,-172.354,160.149,200.502,51.161,-173.107,-141.919,1	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.00209255,0
lmdif_Chebyquad	93	0	3.70419e-23	0.0442053,0.199491,0.235619,0.416047,0.5,0.583953,0.764381,0.800509,0.955795	0.370101,2.06122,1.55372,2.76901,2.92569,2.76925,1.55373,2.06184,0.369676
==21378==
==21378== HEAP SUMMARY:
==21378==     in use at exit: 0 bytes in 0 blocks
==21378==   total heap usage: 367 allocs, 367 frees, 302,876 bytes allocated
==21378==
==21378== All heap blocks were freed -- no leaks are possible
==21378==
==21378== For counts of detected and suppressed errors, rerun with: -v
==21378== ERROR SUMMARY: 0 errors from 0 contexts (suppressed: 4 from 4)
*/
#endif
