# Reference documents

## Gauss.c / Guass.h

This is a (not so pretty) conversion of some c code that was used in a VME rack for flatness control system.

## chebyshev.py / ChebyshevTest.py

This is pulled directly from numpy's polynomial module, a great reference.

## Note

Both the old VME c code and the numpy produce the same results for the 7 point test data set, but the chebyshev series must be converted with numpy's cheb2poly() function.