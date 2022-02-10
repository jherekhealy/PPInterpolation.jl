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
Interpolate f(x) at non-equistant knots with a cubic spline keeping continuity of class C2 (first and second derivatives continuous):

```julia
using PPInterpolation

x = [3.0, 5.0, 6.0, 8.0, 9.0, 11.0, 12.0, 14.0, 15.0]
f = [10.0, 10.0, 10.0, 10.0, 10.5, 15.0, 50.0, 60.0, 85.0]
spline = makeCubicPP(x, f, PPInterpolation.SECOND_DERIVATIVE, 0.0, PPInterpolation.SECOND_DERIVATIVE, 0.0, C2())
spline(4.0)
```

Various monotonicity preserving scheme are available by changing the last parameter `C2()` to

* `C2Hyman89()` : the first derivatives are adjusted to preserve monotonicity according to Dougherty and Hyman "Nonnegativity-, Monotonicity-, or Convexity-Preserving Cubic and Quintic Hermite Interpolation" (1989).
* `C2MP()` : adjust first derivatives according to the simplest monotonicity preserving scheme. See Huyn "Accurate Monotone Cubic Interpolation" (1991).
* `C2MP2()` : adjust first derivatives according to the extended monotonicity preserving scheme. See Huyn "Accurate Monotone Cubic Interpolation" (1991).
* `HuynRational()` : the first derivatives are estimated through a rational function which guarantees monotonicity. See Huyn "Accurate Monotone Cubic Interpolation" (1991).
* `VanAlbada()` : Van Albada type of approximation for the first derivatives. See Huyn "Accurate Monotone Cubic Interpolation" (1991).
* `VanLeer()` : Van Leer type of approximation for the first derivatives. See Huyn "Accurate Monotone Cubic Interpolation" (1991).
* `FritschButland()`  : the first derivatives are estimated through a weighted harmonic mean which guarantees monotonicity.
* `Brodlie()` : similar to FritchButland but takes into account the non-uniformity of the knots better, see Fritsch and Butland "A method for constructing local monotone piecewise cubic interpolants"  (1984).
