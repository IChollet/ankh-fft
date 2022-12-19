## ANKH
The ankh-fft code is a c++ implementation of the ANKH-FFT method presented in https://arxiv.org/abs/2212.08284 . This library aims at efficiently perform energy computation arising in molecular dynamics.

# Dependencies
The ankh-fft code depends on both BLAS and FFTW3

# Install and test
In the root directory, change the BLAS and FFTW3 paths in ./tests/makefile
then run

$ cd ./tests

$ make analyze

$ ./analyze 2 4 15 6 5

Signification of the various parameters can be found in analyze.cpp (i.e. equispaced interpolation order / number of tree levels / number of far images in each direction / interpolation order for far iamges / Chebyshev interpolation order at leaves).

# Licence
The ankh library is under Lesser Gnu Public Licence 3 (LGPL).

# Author
Igor Chollet (Laboratoire Analyse Geometrie et Applications / Laboratoire de Chimie Theorique)
