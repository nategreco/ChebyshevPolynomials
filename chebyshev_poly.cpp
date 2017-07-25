/******************************************************************************************
  Date:    25.07.2017
  Author:  Nathan Greco (Nathan.Greco@gmail.com)

  Project:
      ChebyshevPolynomials: Library for fitting Chebyshev polynomials of the first type.
	  http://github.com/NateGreco/ChebyshevPolynomials.git

  License:
	  This software is licensed under GNU GPL v3.0
	  
******************************************************************************************/

//Standard libraries
#include <vector>
#include <stdexcept>

//Project libraries
#include "chebyshev_poly.h"


/*****************************************************************************************/
template <class T_Type> ChebyshevFit<T_Type>::ChebyshevFit(
	const std::vector<T_Type> x,
	const std::vector<T_Type> y,
	const int degrees):
	order{ degrees + 1 }
{
	//Check that x/y data points match 1:1
	if ( x.size() != y.size() )
	{
		throw std::invalid_argument("X and Y vectors must be same length.");
	}
	//Check for non-zero degree
	if ( degrees < 0 )
	{
		throw std::invalid_argument("Polynomial degree must be non-negative.");
	}
	//Check for minimum number of points
	if ( x.size() < 2 )
	{ //ToDo, check min number of pts
		throw std::invalid_argument("Not enough points.");
	}
	//Check enough points exist to fit degree
	if ( x.size() < degrees )
	{ //ToDo, not sure deg >= pts
		throw std::invalid_argument("Not enought points to fit polynomial degree.");
	}

	//Get Vandermore Matrix
	vandermatrix = VandermondeMatrix(x, degrees);
}

template <class T_Type> T_Type **ChebyshevFit<T_Type>::VandermondeMatrix(
	const std::vector<T_Type> x,
	const int degrees )
{
	//Check for non-zero degree
	if ( degrees < 0 )
	{
		throw std::invalid_argument("Polynomial degree must be non-negative.");
	}
	//Check enough points exist to fit degree
	if ( x.size() < degrees )
	{ //ToDo, not sure deg >= pts
		throw std::invalid_argument("Not enought points to fit polynomial degree.");
	}

	//Create a 2D matrix with C-arrays
	int rows = x.size();
	int cols = degrees + 1;
	T_Type** matrix = new T_Type*[rows];
	for( int i = 0; i < rows; ++i )
	{
		matrix[i] = new T_Type[cols];
	}

	//Generate values
	for( int i = 0; i < rows; ++i )
	{
		matrix[0][i] = static_cast<T_Type>( 1 );
	}
	
}

template <class T_Type> ChebyshevFit<T_Type>::~ChebyshevFit( )
{
	delete vandermatrix;
}