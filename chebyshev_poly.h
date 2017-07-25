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


/*****************************************************************************************/
template <class T_Type> class ChebyshevFit
{
	public:
		ChebyshevFit( const std::vector<T_Type> x,
					  const std::vector<T_Type> y,
					  const int degrees );
		T_Type **VandermondeMatrix( const std::vector<T_Type> x,
									const int degrees );
		T_Type **vandermatrix;
		int order;
		std::vector<T_Type> coefs;
		~ChebyshevFit( );
};

#endif // CHEBYSHEV_POLY_H_INCLUDED
