# Package PPInterpolation
| Status | Coverage |
| :----: | :----: |
| [![Build Status](https://github.com/jherekhealy/PPInterpolation.jl/actions/workflows/julia-runtests.yml/badge.svg) | [![codecov.io](http://codecov.io/github/jherekhealy/PPInterpolation.jl/coverage.svg?branch=master)](http://codecov.io/github/jherekhealy/PPInterpolation.jl?branch=master) |


Piecewise polynomial interpolation in Julia following a straightforward implementation. This is mostly oriented towards various cubic spline interpolation:

* C2 cubic spline, eventually filtering out the first derivatives to ensure monotonicity.
* C1 Hermite spline with known first derivatives
* C1 Hermite spline with Bessel (parabolic) estimate of derivatives
* C1 Hermite spline using limiters such as VanLeer, VanAlbada, Huynh rational function, Fritsch-Butland (1980) harmonic mean, Fritch-Butland (1984) Brodlie formula.

with the following boundary conditions:

* first derivative at end point
* second derivative at end point
* not-a-knot
* first difference estimate at end point
 
Also supported are C3 Hermite quintic spline with  given first derivatives at all knots, plus second derivatives at end-points, cubic and quintic Lavery splines.


You may find other Julia interpolation packages relevant, such as Interpolations.jl


## Examples

