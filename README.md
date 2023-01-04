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
 
Also supported are C3 Hermite quintic spline with  given first derivatives at all knots, plus second derivatives at end-points. 
Cubic and quintic Lavery splines are available in the package [Lavery.jl](https://github.com/jherekhealy/LaverySpline.jl).


## Why another Julia interpolation package?
At the time this package was created (and still today), there are surprisingly few interpolation packages for Julia which support cubic splines. And the existing ones do not provide any control over the boundary conditions. The existing packages are not focused on cubic splines but tend to be more general.

* [Interpolations.jl](https://github.com/JuliaMath/Interpolations.jl) does not support C2 cubic splines with knots specified by a vector  (gridded in their terminology). 
* [DataInterpolations.jl](https://github.com/PumasAI/DataInterpolations.jl) does not support specifying endpoints conditions. 
* [Dierckx.jl](https://github.com/kbarbary/Dierckx.jl) is somewhat slow and supports only the standard C2 spline. It offers other features however, such as 2D splines.

This package support various boundary conditions (natural, not-a-knot, specific first or second derivatives), and not only C2 cubic splines but also C1 splines, monotonocity preserving splines. The flexible API allows to evaluate faster if we know the vector to evaluate is sorted.

Below is an example of performance using `x = sort(rand(500))` and `y=rand(500)` interpolating an array `z = collect(range(x[1],stop=x[end],length=500))`

| Package | Construction | Evaluation |
|:--------|--------------:|-----------------------:|
| Dierckx  Spline1D(x,y) | 48.692 μs (8 allocations: 97.11 KiB) |  21.405 μs (2 allocations: 4.09 KiB) |
| DataInterpolation CubicSpline(y,x) |  26.123 μs (524 allocations: 73.11 KiB) |21.591 μs (3 allocations: 4.14 KiB)|
| Interpolations interpolate(x,y,SteffenMonotonicInterpolation()) | 6.641 μs (5 allocations: 16.30 KiB) |   21.427 μs (3 allocations: 4.16 KiB)|
| PPInterpolation C2() | 21.610 μs (26 allocations: 89.41 KiB) | 15.503 μs (1 allocation: 48 bytes) |
| PPInterpolation VanAlbada() | 7.144 μs (13 allocations: 48.64 KiB) | 15.879 μs (1 allocation: 48 bytes) |
| PPInterpolation QuadraticLagrangePP(x,y) | 54.487 ns (1 allocation: 32 bytes) | 15.711 μs (1 allocation: 48 bytes) |
| PPInterpolation QuadraticLagrangePP(x,y) evaluateSorted | |  2.942 μs (0 allocations: 0 bytes)|



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


Quadratic Lagrange interpolation has fast evaluation API with mutable variable:

```julia
using PPInterpolation

x = [3.0, 5.0, 6.0, 8.0, 9.0, 11.0, 12.0, 14.0, 15.0]
f = [10.0, 10.0, 10.0, 10.0, 10.5, 15.0, 50.0, 60.0, 85.0]
spline = QuadraticLagrangePP(x, f)
spline(4.0)
y = Array{Float64}(undef, 100) #the output
z = collect(range(3.0,stop=15.0,length=100)) #this is sorted input
evaluateSorted!(spline,y,z)
```
