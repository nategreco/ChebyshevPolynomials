/******************************************************************************************
  Date:    25.07.2017
  Author:  Nathan Greco (Nathan.Greco@gmail.com)

  Project:
      ChebyshevPolynomials: Library for fitting Chebyshev polynomials of the first type.
	  http://github.com/NateGreco/ChebyshevPolynomials.git

  License:
	  This software is licensed under GNU GPL v3.0
	  
******************************************************************************************/

//Header guard
#ifndef CHEBYSHEV_POLY_H_INCLUDED
#define CHEBYSHEV_POLY_H_INCLUDED

//Standard libraries
#include <vector>
#include <stdexcept>

/*****************************************************************************************/
template <class T_Type> class ChebyshevFit
{
	public:
		ChebyshevFit( const std::vector<T_Type> x,
					  const std::vector<T_Type> y,
					  const int degrees );
		int order;
		std::vector<T_Type> coefs;
		~ChebyshevFit( );
	private:
		std::vector<std::vector<T_Type>>
			VandermondeMatrix( const std::vector<T_Type> x,
							   const int degrees );
		std::vector<std::vector<T_Type>>
			Transpose( const std::vector<std::vector<T_Type>> matrix );
};

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
	{
		throw std::invalid_argument("Not enough points.");
	}
	//Check enough points exist to fit degree
	if ( x.size() < degrees )
	{
		throw std::invalid_argument("Not enought points to fit polynomial degree.");
	}

	//Get Vandermonde Matrix
	std::vector<std::vector<T_Type>> vandermatrix{VandermondeMatrix(x, degrees) };

	//Get transpose of Vandermonde Matrix
	std::vector<std::vector<T_Type>> trans_matrix{Transpose(vandermatrix) };
	
}

template <class T_Type> std::vector<std::vector<T_Type>> 
	ChebyshevFit<T_Type>::VandermondeMatrix( const std::vector<T_Type> x,
											 const int degrees )
{
	//Check for non-zero degree
	if ( degrees < 0 )
	{
		throw std::invalid_argument("Polynomial degree must be non-negative.");
	}
	//Check enough points exist to fit degree
	if ( x.size() < degrees )
	{
		throw std::invalid_argument("Not enought points to fit polynomial degree.");
	}

	//Create a 2D matrix with vectors
	int rows = x.size();
	int cols = degrees + 1;
	std::vector<std::vector<T_Type>>  matrix;

	//Generate values
	for( int i = 0; i < x.size(); i++ )
	{
		std::vector<T_Type> row;
		row.push_back( static_cast<T_Type>( 1 ) );
		if (degrees == 0) continue;
		row.push_back(  x[i] );
		if (degrees == 1) continue;
		for( int j = 2; j < cols; j++ )
		{
			row.push_back(row[j - 1] * 2.0 * x[i] - row[j - 2] );
		}
		matrix.push_back(row);
	}
	
	return matrix;
	
}

template <class T_Type> std::vector<std::vector<T_Type>> 
	ChebyshevFit<T_Type>::Transpose( const std::vector<std::vector<T_Type>> matrix )
{
	//Define new matrix
	std::vector<std::vector<T_Type>> trans_matrix;

	//Copy values from old one
	for (int i = 0; i < matrix[0].size(); i++) {
		std::vector<T_Type> row;
		for (int j = 0; j < matrix.size(); j++) {
			row.push_back( matrix[j][i] );
		}
		trans_matrix.push_back( row );
	}
	return trans_matrix;
}

template <class T_Type> ChebyshevFit<T_Type>::~ChebyshevFit( )
{

}

#endif // CHEBYSHEV_POLY_H_INCLUDED
