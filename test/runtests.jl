using PPInterpolation, Test
using ForwardDiff

@testset "Akima" begin
	x = [3.0, 5.0, 6.0, 8.0, 9.0, 11.0, 12.0, 14.0, 15.0]
	f = [10.0, 10.0, 10.0, 10.0, 10.5, 15.0, 50.0, 60.0, 85.0]
	varieties = [C2Hyman89(), C2MP(), HuynRational(), VanLeer(), FritschButland(), Brodlie(), VanAlbada()]
	for v in varieties
		spline = makeCubicPP(x, f, PPInterpolation.SECOND_DERIVATIVE, 0.0, PPInterpolation.SECOND_DERIVATIVE, 0.0, v)
		nSample = 1000
		previousVal = spline(x[1])
		for i ∈ 1:nSample
			ti = x[1] + (x[end] - x[1]) * (i) / (nSample)
			zval = spline(ti)
			if i > 2
				@test zval >= previousVal
			end
			previousVal = zval
		end
	end


	# nSample = 10
	# spline = makeCubicPP(x, f, PPInterpolation.FIRST_DERIVATIVE,0.0, PPInterpolation.FIRST_DERIVATIVE,0.0, PPInterpolation.Fukasawa())
	# previousVal = spline(x[1])
	# for i = 1:nSample
	#     ti = x[1] + (x[end] - x[1]) * (i) / (nSample)
	#     zval = spline(ti)
	#     if i > 2
	#         @test zval >= previousVal
	#     end
	#     previousVal = zval
	# end

	#     varieties = [ C2Hyman89(), C2MP(), HuynRational(), VanLeer(), FritschButland(), Brodlie(), VanAlbada() ]^C
	#  p = plot()
	#  for v in varieties
	#        spline = makeCubicPP(x,f,PPInterpolation.SECOND_DERIVATIVE,0.0,PPInterpolation.SECOND_DERIVATIVE,0.0,v)
	#        plot!(p, t, spline.(t), label=string(v))
	#        end
	#  plot!(legend=:top)

end

@testset "FritschCarlson" begin
	x = [7.99, 8.09, 8.19, 8.7, 9.2, 10.0, 12.0, 15.0, 20.0]
	f = [0.0, 2.76429e-5, 4.37498e-2, 0.169183, 0.469428, 0.943740, 0.998636, 0.999919, 0.999994]

	varieties = [C2Hyman89(), C2MP(), HuynRational(), VanLeer(), FritschButland(), Brodlie(), VanAlbada()]
	for v in varieties
		spline = makeCubicPP(x, f, PPInterpolation.SECOND_DERIVATIVE, 0.0, PPInterpolation.SECOND_DERIVATIVE, 0.0, v)
		nSample = 1000
		previousVal = spline(x[1])
		for i ∈ 1:nSample
			ti = x[1] + (x[end] - x[1]) * (i) / (nSample)
			zval = spline(ti)
			if i > 2
				@test zval >= previousVal
			end
			previousVal = zval
		end
	end
end
@testset "HaganCurveAuto" begin
	#for quintic: first pass compute cubic first derivative, second pass, use cubic b, and optimize second derivative.
	t = [0.0, 0.1, 1.0, 2.0, 3.0, 4.0, 9.0, 20.0, 30.0]
	y = [0.06, 0.06, 0.06, 0.05, 0.05, 0.05, 0.065, 0.06, 0.06]
	t = [0.0, 0.1, 1.0, 4.0, 9.0, 20.0, 30.0]
	y = [0.081, 0.081, 0.07, 0.044, 0.07, 0.04, 0.03]
	t = [0, 0.1, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0, 3.5]
	y = [0.012, 0.008, 0.0104, 0.00995, 0.010633333333333333, 0.010849999999999999, 0.010179999999999998, 0.009816666666666665, 0.008612499999999999, 0.007789999999999998, 0.007741666666666665, 0.007635714285714285]

	# t = [0.0, 0.1, 1, 2, 3, 4, 5, 6, 7, 30.0]
	# f = [0.01, 0.01, 0.01, 0.01, 0.01, 0.03, 0.03, 0.03, 0.03, 0.02]
	# y = zeros(length(t))
	# y[1] = f[1]
	# y[2:end] = (y[1:end-1] .* t[1:end-1] + (t[2:end]-t[1:end-1]) .* f[2:end]) ./ t[2:end]
	# y[3:end-2] = f[3:end-2] ./ t[3:end-2]
	# println(y)
	yt = zeros(length(t))
	yt[1:end] = y[1:end] .* t[1:end]
	spline = makeCubicPP(t, yt, PPInterpolation.SECOND_DERIVATIVE, 0.0, PPInterpolation.SECOND_DERIVATIVE, 0.0, C2())
	nSample = 1000
	for i ∈ 1:nSample
		ti = (t[end]) * (i) / (nSample)
		val = evaluateDerivative(spline, ti)
		zval = spline(ti) / ti
		println(ti, " CubicPP ", val, " ", zval)
		autoval = ForwardDiff.derivative(spline, ti)
		@test isapprox(val, autoval, atol = 1e-15)
	end
	#AUTODIFF // where is autodiff useful (hagan wave read),  BigFloat
end

using Random
struct ReflectingAngle <: PPInterpolation.LimiterDerivative end

function PPInterpolation.fillDerivativeEstimate(limiter::ReflectingAngle, dx::AbstractArray{TX}, S::AbstractArray{T}, b::AbstractArray{T}) where {T, TX}
	# @. b[2:end-1] = (dx[1:end-1] * S[2:end] + dx[2:end] * S[1:end-1]) / (dx[1:end-1] + dx[2:end])
	n = length(S)
	for i ∈ 2:n
		tan2a = (S[i] + S[i-1]) / (1 - S[i] * S[i-1])
		sqrtDisc = sqrt(1 + tan2a^2)
		tana = (-1 - sqrtDisc) / (tan2a)
		tanb = (-1 + sqrtDisc) / (tan2a)
		b[i] = if tan2a == zero(T)
			zero(T)
		elseif sign(tana) == sign(S[i-1])
			tana
		else
			tanb
		end
	end
end

@testset "FiveHundredAngle" begin
	#  Random.seed!(1)
	#  x = sort(rand(500));
	#  y = rand(500);
	x = [3.0, 5.0, 6.0, 8.0, 9.0, 11.0, 12.0, 14.0, 15.0]
	y = [10.0, 10.0, 10.0, 10.0, 10.5, 15.0, 50.0, 60.0, 85.0]
	spline = makeCubicPP(x, y, PPInterpolation.FIRST_DERIVATIVE, 0.0, PPInterpolation.FIRST_DERIVATIVE, 0.0, Fukasawa())
	z = collect(range(x[1], stop = x[end], length = 500))
	yNew = Array{Float64}(undef, length(z))
	@time yNew[1:end] = spline(z)
	@test isapprox(y[5], spline(x[5]), atol = 1e-15)
	splineb = makeCubicPP(x, y, PPInterpolation.FIRST_DERIVATIVE, 0.0, PPInterpolation.FIRST_DERIVATIVE, 0.0, ReflectingAngle())
	ybNew = Array{Float64}(undef, length(z))
	@time ybNew[1:end] = splineb(z)
	for i ∈ 1:length(yNew)
		@test isapprox(yNew[i], ybNew[i], atol = 1e-8)
	end
end
using QuadGK
@testset "FiveHundred" begin
	Random.seed!(1)
	x = sort(rand(500))
	y = rand(500)
	spline = makeCubicPP(x, y, PPInterpolation.SECOND_DERIVATIVE, 0.0, PPInterpolation.SECOND_DERIVATIVE, 0.0, C2())
	z = collect(range(x[1], stop = x[end], length = 500))
	yNew = Array{Float64}(undef, length(x))
	@time yNew[1:end] = spline(z)
	@test isapprox(y[10], spline(x[10]), atol = 1e-15)
	yNewCopy = Array{Float64}(undef, length(x))
	@time for i ∈ 1:length(z)
		yNew[i] = spline(z[i])
	end
	@time evaluateSorted!(spline, yNewCopy, z)
	for i ∈ 1:length(yNew)
		@test isapprox(yNew[i], yNewCopy[i], atol = 0.0)
	end
	q = QuadraticLagrangePP(x, y)
	@time evaluateSorted!(q, yNew, z)
	@test isapprox(y[10], q(x[10]), atol = 1e-15)
	val = evaluateDerivative(q, 0.5)
	autoval = ForwardDiff.derivative(q, 0.5)
	@test isapprox(val, autoval, atol = 1e-12)
	# integration tests
	lower_bounds = rand(100) * (x[end] - x[1]) .+ x[1]
	upper_bounds = rand(100) * (x[end] - x[1]) .+ x[1]
	for (lower, upper) in zip(lower_bounds, upper_bounds)
		val = evaluateIntegral(spline, lower, upper)
		quadval = quadgk(spline, lower, upper)
		@test isapprox(val, quadval[1], atol = 5e-8)
	end

end

@testset "Quintic" begin
	Random.seed!(1)
	x = sort(rand(500))
	y = rand(500)
	yp = rand(500)
	spline = makeHermiteQuinticPP(x, y, yp)
	z = collect(range(x[1], stop = x[end], length = 500))
	yNew = Array{Float64}(undef, length(x))
	@time yNew[1:end] = spline(z)
	@test isapprox(y[10], spline(x[10]), atol = 1e-15)
	@test isapprox(yp[10], evaluateDerivative(spline, x[10]), atol = 1e-15)
end

@testset "SchumakerExample" begin
	t = Float64.([1, 2, 3, 4, 5])
	z = Float64.([1, 2, 3, 2, 1])
	ref =
		[1 1.08 1.16 1.24 1.32 1.4 1.48 1.56 1.6400000000000001 1.72 1.8 1.88 1.96 2.0408 2.1272 2.22 2.3192000000000004 2.4248000000000003 2.5368 2.6544 2.7600000000000002 2.8463999999999996 2.9135999999999997 2.9616 2.9904 3 2.9904 2.9616 2.9135999999999997 2.8464 2.7600000000000002 2.6544 2.5368 2.4248 2.3192 2.22 2.1272 2.0408 1.96 1.88 1.7999999999999998 1.7200000000000006 1.6400000000000006 1.5600000000000005 1.4800000000000004 1.4000000000000004 1.3200000000000003 1.2400000000000002 1.1600000000000001 1.08 1]

	interp = makeSchumakerQuadraticPP(t, z, SchumakerDerivative())
	for i ∈ 0:50
		ui = 1.0 + 4i / 50.0
		zi = interp(ui)
		println(ui, " ", zi, " ", evaluateDerivative(interp, ui))
		@test isapprox(ref[i+1], zi, atol = 1e-15)
	end
end

@testset "SchumakerRouah" begin
	strike = Float64.([300, 310, 320, 330, 340, 350, 360, 370, 380, 390, 400])
	put = [0.09, 0.136, 0.192, 0.288, 0.404, 0.53, 0.736, 1.232, 1.898, 3.104, 4.93]
	interp = makeSchumakerQuadraticPP(strike, put, SchumakerDerivative())
	interp2 = makeSchumakerQuadraticPP(strike, put, LamDerivative{Float64}())
	
    for i ∈ 0:50
		ui = strike[1] + (strike[end] - strike[1]) * i / 50.0
		zi = interp(ui)
		z2i = evaluateSecondDerivative(interp, ui)
		println(ui, " ", zi, " ", evaluateDerivative(interp, ui), " ", z2i)
		@test z2i >= 0
        z2i2 = evaluateSecondDerivative(interp2, ui)
		println(ui, " ", interp2(ui), " ", evaluateDerivative(interp2, ui), " ", z2i2)
		@test z2i2 >= 0
	
	end
end

include("conversion_tests.jl")
