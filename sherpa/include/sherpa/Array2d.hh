#ifndef Array2d_hh
#define Array2d_hh

//_C++_INSERT_SAO_COPYRIGHT_HERE_(2007)_
//_C++_INSERT_GPL_LICENSE_HERE_

#include <iostream>
#include <stdexcept>
#include <vector>

//
// A simple 2d array class written for NelderMead/MultiDirSearch.
// Note Array2d as written, using vector of vector, does not guaranteed
// contiguous data. Check out one of the followings if contiguous 
// in memory is necessary:
//
//    http://math.nist.gov/tnt/index.html
//    http://www.oonumerics.org/blitz/
//    http://boost.org/libs/multi_array/doc/user.html
//
//
namespace sherpa {

  template < typename T >
  class Array2d {

    friend std::ostream& operator << ( std::ostream& os,
				       const Array2d< T >& a ) {
      return a.print( os );
    }

  public:

    virtual ~Array2d( ) { }

    Array2d( int r=0, int c=0 ) : nrow( r ), ncol( c ),
				  array( r, std::vector< T >( c ) ) { }

    std::vector< T >& operator[] ( int arg ) { 
      return array[ arg ];
    }
    const std::vector< T >& operator[] ( int arg ) const { 
      return array[ arg ];
    }

    int ncols( ) const { return ncol; }

    int nrows( ) const { return nrow; }

    static void copy_vector( int num, const std::vector< T >& from,
			     std::vector< T >& to ) {

      for ( int ii = 0; ii < num; ++ii )
	to[ ii ] = from[ ii ];

    }

    std::ostream& print( std::ostream& os ) const {

      for ( int ii = 0; ii < nrow; ++ii ) {
	os << array[ ii ][ 0 ];
	for ( int jj = 1; jj < ncol; ++jj )
	  os << ' ' << array[ ii ][ jj ];
	if ( nrow - 1 != ii )
	  os << '\n';
      }
      return os;
    }

    // resize the two dimentional array
    virtual void resize( int r, int c ) {
      array.resize( r );
      for( int ii = 0; ii < r; ++ii ) 
	array[ ii ].resize( c );
      nrow = r;
      ncol = c;
    }

  protected:

    T& get( int r, int c ) { return array[ r ][ c ]; }
    const T& get( int r, int c ) const { return array[ r ][ c ]; }

    void set( int r, int c, const T& val ) { array[ r ][ c ] = val; }

  private:

    int nrow, ncol;
    std::vector< std::vector< T > > array;

    Array2d& operator = (Array2d const&); // declare but, purposely, not define
    Array2d( Array2d const& );            // declare but, purposely, not define

  };

}

#endif

#ifdef testArray2d

#include <iostream>

#ifdef testArrayNd
#include "ArrayNd.hh"
#endif
#include "StopWatch.hh"

template < typename Type >
void initLaplacian( Type& uu, int r, int c ) {

  for ( int ii = 0; ii < r; ++ii ) {
    for ( int jj = 0; jj < c; ++jj )
      uu[ ii ][ jj ] = ii * ii + jj * jj;
  }

}

template < typename Type >
void timeme( Type& uu, Type& laplacian, int r, int c, const char* header ) {

  StopWatch stopwatch( header );

  initLaplacian( uu, r, c );

  for( int kk=0; kk < 10; kk++ ) {
    for( int ii = 1; ii < r - 1; ii++ ) {
      for( int jj = 1; jj < c - 1; jj++ ) {
        laplacian[ ii ][ jj ] = - uu[ ii - 1 ][ jj ] - uu[ ii ][ jj - 1 ] +
          4.0 * uu[ ii ][ jj ] - uu[ ii ][ jj + 1 ] - uu[ ii + 1 ][ jj ];
      }
    }
  }

  double sum = 0.0;
  for ( int ii = 1; ii < r - 1; ++ii )
    for ( int jj = 1; jj < c - 1; jj++ )
      sum += laplacian[ ii ][ jj ];

  fprintf( stdout, "%.12f\t", sum );

}

template < typename Type >
Type*** alloc3d( int nx, int ny, int nz ) {

  Type ***ptr = new int** [ nx ];

  for( int x = 0; x < nx; ++x ) {
    ptr[ x ] = new Type* [ ny ];

    for( int y = 0; y < ny; ++ y )
      ptr[ x ][ y ] = new Type[ nz ];
  }

  return ptr;

}

template < typename Type >
void del3d( Type*** ptr, int nx, int ny ) {
  for( int x = 0; x < nx; ++x ) {
    for( int y = 0; y < ny; ++y )
      delete [] ptr[ x ][ y ], ptr[ x ][ y ] = 0;
    delete [] ptr[ x ], ptr[ x ] = 0;
  }
  delete [] ptr;
  ptr = 0;
}

template < typename Type >
Type** alloc2d( int r, int c ) {
  Type** ptr = new Type*[ r ];
  for ( int ii = 0; ii < r; ++ii )
    ptr[ ii ] = new Type[ c ];
  return ptr;
}

template < typename Type >
void del2d( Type** ptr, int r ) {

  for ( int ii = 0; ii < r; ++ii )
    delete [] ptr[ ii ];
  delete [] ptr;

}

void timeclassic( int r, int c ) {

  double** uu = alloc2d< double >( r, c );
  double** laplacian = alloc2d< double >( r, c );

  timeme( uu, laplacian, r, c, "classic" );

  del2d( laplacian, r );
  del2d( uu, r );

}

void timebracket( int r, int c ) {

  sherpa::Array2d< double > uu( r, c ), laplacian( r, c );

  for ( int ii = 0; ii < r; ++ii ) {
    for ( int jj = 0; jj < c; ++jj )
      uu[ ii ][ jj ] = ii * ii + jj * jj;
  }

  timeme( uu, laplacian, r, c, "sherpa::Array2d[][]" );

}

#ifdef ARRAYND_HH
void timemulti_array( int r, int c ) {

  std::vector<size_t> dims( 2, r );
  multi_array< double, 2 > laplacian( dims ), uu( dims );

  timeme( uu, laplacian, r, c, "Multi_Array[][]" );

}

void timeArray( int r, int c ) {

  std::vector<size_t> dims( 2, r );
  Array< double, 2 > laplacian( dims ), uu( dims );

  timeme( uu, laplacian, r, c, "Array[][]" );

}

void timeBavestrelliArray( int r, int c ) {

  bavestrelli::Array<double, 2> laplacian( bavestrelli::ArraySizes(r)(c) ),
    uu(  bavestrelli::ArraySizes(r)(c) );

  timeme( uu, laplacian, r, c, "bavestrelli::Array[][]" );

}
#endif

template< typename A, typename B >
void cmpme( A& a, B& b, int r, int c ) {

  for ( int ii = 0; ii < r; ++ii )
    for ( int jj = 0; jj < c; ++jj )
      if ( a[ ii ][ jj ] != b[ ii ][ jj ] )
	printf( "a[%d][%d] = %g vs b[%d][%d] = %g\n",
		ii, jj, a[ ii ][ jj ], ii, jj, b[ ii ][ jj ] );

}


/*
void testme( int r, int c ) {

  sherpa::Array2d< int > foo( r, c );

  for ( int ii = 0; ii < r; ++ii )
    for ( int jj = 0; jj < c; ++jj )
      foo[ ii ][ jj ] = ii + jj;

  std::cout << foo << '\n';

  sherpa::Array2d< int >::Row bar = foo[ r / 2 ];
  int* ptr = &bar[ 0 ];
  std::cout << ptr[ 0 ];
  for ( int jj = 1; jj < c; ++jj )
    std::cout << '\t' << ptr[ jj ];
  std::cout << '\n';

}
*/

int main( int argc, char* argv[] ) {

  const int num = 4000;

  int row = num;
  int col = num;
  if ( argc == 3 ) {
    row = atoi( argv[1] );
    col = atoi( argv[2] );
  }


  timeclassic( row, col );
  timebracket( row, col );
#ifdef ARRAYND_HH
  timemulti_array( row, col );
  timeArray( row, col );
  timeBavestrelliArray( row, col );
#endif
  return 0;

}

#endif

//
//cp Array2d.hh tmp.cc; g++ -Wall -ansi -pedantic -O3 -DtestArray2d -DNDEBUG tmp.cc; a.out
//cp Array2d.hh tmp.cc; g++ -Wall -ansi -pedantic -O3 -DtestArrayNd -DtestArray2d -DNDEBUG tmp.cc; a.out
