import numpy as np
from numpy.polynomial import chebyshev

# Create a few points
x = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
y = [1.0, 2.0, 0.0, 2.0, 3.0, 2.0, 8.0]

# Get chebyshev polynomial (first type) curve fit
cheb = chebyshev.chebfit(x, y, 4)
poly = chebyshev.cheb2poly(cheb)
cheb2 = chebyshev.poly2cheb(poly)

# Print
print(cheb)
print(poly)
print(cheb2)