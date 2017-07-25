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
#include <iostream>
#include <vector>

//Project libraries
#include "chebyshev_poly.h"


/*****************************************************************************************/
int main()
{
	std::vector<double> x{0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
	std::vector<double> y{1.0, 2.0, 0.0, 2.0, 3.0, 2.0, 8.0};
	ChebyshevFit<double> cheb(x, y, 4);
	for ( double &coef : cheb.coefs )
	{
		std::cout << coef << ", ";
	}

	return 0;
}