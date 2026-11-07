using PPInterpolation, Test, ForwardDiff, Random
const P = PPInterpolation

function allSplines()
	x = [0.0, 0.7, 1.5, 2.0, 3.2, 4.0, 5.5]
	y = [1.0, 2.2, 2.0, 0.8, 0.9, 1.6, 1.2]
	yp = [1.5, 0.4, -1.2, -0.6, 0.3, 0.2, -0.5]
	return [
		"C2" => makeCubicPP(x, y, C2()),
		"C2Hyman89" => makeCubicPP(x, y, C2Hyman89(), left = NotAKnot(), right = FirstDifference()),
		"VanLeer" => makeCubicPP(x, y, VanLeer(), left = FirstDerivative(1.0), right = SecondDerivative(0.5)),
		"Hermite" => makeCubicPP(x, y, Hermite(yp)),
		"Schumaker" => SchumakerQuadraticPP(x, y),
		"Lam" => SchumakerQuadraticPP(x, y, LamDerivative(0.5)),
		"Quintic" => makeHermiteQuinticPP(x, y, yp),
		"QuinticEnds" => makeHermiteQuinticPP(x, y, yp, 0.0, 1.0),
		"Linear" => makeLinearPP(x, y),
		"LagrangeMid" => convert(PP, QuadraticLagrangePP(x, y)),
		"LagrangeLeft" => convert(PP, QuadraticLagrangePP(x, y, knotStyle = P.LEFT_KNOT)),
	]
end

@testset "PPLayout" begin
	for (name, spl) in allSplines()
		n = length(spl.x)
		@test length(spl.a) == n
		@test length(spl.b) == n
		@test size(spl.c, 1) == n - 1
	end
	@test_throws DimensionMismatch PP(2, zeros(3), zeros(2), zeros(2, 1), [1.0, 2.0, 3.0])
	@test_throws DimensionMismatch PP(3, zeros(3), zeros(3), zeros(3, 2), [1.0, 2.0, 3.0])
end

@testset "DerivativeConsistency" begin
	#each derivative is the derivative of the lower order one, inside, at the knots and beyond, for every extrapolation
	fs = [evaluate, evaluateDerivative, evaluateSecondDerivative, evaluateThirdDerivative]
	for (name, spl) in allSplines(), ext in (LinearAutoExtrapolation(), ConstantAutoExtrapolation(), PieceAutoExtrapolation(), ZeroExtrapolation())
		zs = [spl.x[1] - 0.5; spl.x; spl.x[end] + 0.5; (spl.x[1:end-1] .+ spl.x[2:end]) ./ 2]
		for z in zs, k in 1:3
			#at a knot, the derivative from the right (or from the left at the last knot) is used, consistently with automatic differentiation
			expected = ForwardDiff.derivative(t -> fs[k](spl, t, extrapolation = ext), z)
			@test isapprox(fs[k+1](spl, z, extrapolation = ext), expected, atol = 1e-8, rtol = 1e-8)
		end
	end
end

@testset "StoredExtrapolation" begin
	x = [0.0, 1.0, 2.0, 3.0]
	y = [0.0, 1.0, 0.5, 2.0]
	zs = [-1.0, 0.5, 2.5, 4.0]
	linear = makeCubicPP(x, y, C2())
	for ext in (ConstantAutoExtrapolation(), PieceAutoExtrapolation(), ZeroExtrapolation())
		spl = makeCubicPP(x, y, C2(), extrapolation = ext)
		@test spl.(zs) == [evaluate(linear, z, extrapolation = ext) for z in zs]
		@test PP(linear, extrapolation = ext).(zs) == spl.(zs)
		@test evaluateSorted!(spl, zeros(4), zs) == spl.(zs)
	end
	@test makeCubicPP(x, y, C2(), extrapolation = ConstantAutoExtrapolation())(4.0) == 2.0
	@test makeCubicPP(x, y, C2(), extrapolation = ZeroExtrapolation())(-1.0) == 0.0
	@test SchumakerQuadraticPP(x, y, extrapolation = ConstantAutoExtrapolation())(4.0) == 2.0
	@test makeHermiteQuinticPP(x, y, zeros(4), extrapolation = ConstantAutoExtrapolation())(-1.0) == 0.0
	@test linear(4.0) ≈ 2.0 + linear.b[end]
end

@testset "BoundaryConditionTypes" begin
	x = [0.0, 1.0, 2.5, 3.0, 4.5]
	y = [1.0, 0.0, 2.0, 1.5, 3.0]
	pairs = [
		(Natural(), P.SECOND_DERIVATIVE, 0.0),
		(SecondDerivative(0.3), P.SECOND_DERIVATIVE, 0.3),
		(FirstDerivative(-1.0), P.FIRST_DERIVATIVE, -1.0),
		(FirstDifference(), P.FIRST_DIFFERENCE, 0.0),
		(NotAKnot(), P.NOT_A_KNOT, 0.0),
	]
	for kind in (C2(), C2MP(), Bessel(), FritschButland()), (bc, enum, value) in pairs
		spl = makeCubicPP(x, y, kind, left = bc, right = bc)
		ref = makeCubicPP(x, y, enum, value, enum, value, kind)
		@test spl.b == ref.b
		@test spl.c == ref.c
	end
	@test makeCubicPP(x, y).b == CubicSplineNatural(x, y).b
end

@testset "LagrangeConversion" begin
	rng = MersenneTwister(7)
	x = sort(rand(rng, 12))
	y = randn(rng, 12)
	z = x[1] .+ (x[end] - x[1]) .* rand(rng, 500)
	for style in (P.LEFT_KNOT, P.MID_KNOT)
		q = QuadraticLagrangePP(x, y, knotStyle = style)
		pp = convert(PP, q)
		@test isapprox(pp.(z), q.(z), atol = 1e-12)
		@test isapprox(evaluateDerivative.(pp, z), evaluateDerivative.(q, z), atol = 1e-10)
	end
end

@testset "SchumakerRightExtrapolation" begin
	x = [1.0, 2, 3, 4, 5]
	y = [1.0, 2, 3, 2.5, 4]
	s = SchumakerQuadraticPP(x, y)
	@test evaluateDerivative(s, 5.5) ≈ evaluateDerivative(s, 5.0)
	@test evaluateDerivative(s, 5.5) ≈ ForwardDiff.derivative(s, 5.5)
	@test evaluateThirdDerivative(s, 3.3) == 0.0
end

#data generators for the shape preservation properties
function shapeData(kind, rng)
	n = rand(rng, 4:10)
	x = cumsum(0.1 .+ rand(rng, n) .^ 2)
	y = if kind == :increasing
		cumsum(rand(rng, n) .^ 3)
	elseif kind == :decreasing
		-cumsum(rand(rng, n) .^ 3)
	elseif kind == :flat
		cumsum(rand(rng, n) .* (rand(rng, n) .< 0.5))
	else
		randn(rng, n)
	end
	return x, y
end

#true if the spline is monotone in the direction of the data on each interval, and constant on flat intervals
function isComonotone(spl, x, y; m = 40, tol = 1e-10)
	scale = maximum(abs, y) + 1
	for i in 1:length(x)-1
		v = spl.(range(x[i], x[i+1], length = m))
		dv = diff(v)
		S = y[i+1] - y[i]
		ok = S > 0 ? minimum(dv) >= -tol * scale : S < 0 ? maximum(dv) <= tol * scale : maximum(abs.(v .- y[i])) <= tol * scale
		ok || return false
	end
	return true
end

@testset "MonotonicityProperty" begin
	rng = MersenneTwister(2024)
	#monotone on monotone data
	for kind in (C2Hyman89(), C2MP(), C2MP2(), HuynRational(), VanLeer(), VanAlbada(), FritschButland(), Brodlie()), bc in (Natural(), FirstDifference())
		for dataKind in (:increasing, :decreasing, :flat), trial in 1:30
			x, y = shapeData(dataKind, rng)
			@test isComonotone(makeCubicPP(x, y, kind, left = bc, right = bc), x, y)
		end
	end
	#piecewise monotone on any data, when the derivative is zero at extrema
	for kind in (C2MP(), C2MP2(), HuynRational(), VanLeer(), FritschButland(), Brodlie()), trial in 1:30
		x, y = shapeData(:general, rng)
		@test isComonotone(makeCubicPP(x, y, kind, left = FirstDifference(), right = FirstDifference()), x, y)
	end
	for dataKind in (:increasing, :decreasing, :flat, :general), trial in 1:30
		x, y = shapeData(dataKind, rng)
		@test isComonotone(SchumakerQuadraticPP(x, y, LamDerivative(0.5)), x, y)
	end
end

@testset "ConvexityAndPositivityProperty" begin
	rng = MersenneTwister(99)
	for trial in 1:30
		x = cumsum(0.1 .+ rand(rng, 8))
		slopes = randn(rng) .+ cumsum(rand(rng, 7))
		y = [0.0; cumsum(slopes .* diff(x))]
		for spl in (SchumakerQuadraticPP(x, y), SchumakerQuadraticPP(x, y, LamDerivative(0.5)))
			z = range(x[1], x[end], length = 500)
			@test minimum(evaluateSecondDerivative.(spl, z)) >= -1e-9
		end
	end
	for trial in 1:30
		x = cumsum(0.1 .+ rand(rng, 10))
		y = rand(rng, 10) .^ 4
		spl = makeCubicPP(x, y, C2HymanNonNegative())
		@test minimum(spl.(range(x[1], x[end], length = 500))) >= -1e-12
	end
end
