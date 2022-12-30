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
        for i = 1:nSample
            ti = x[1] + (x[end] - x[1]) * (i) / (nSample)
            zval = spline(ti)
            if i > 2
                @test zval >= previousVal
            end
            previousVal = zval
        end
    end
    #     varieties = [ C2Hyman89(), C2MP(), HuynRational(), VanLeer(), FritschButland(), Brodlie(), VanAlbada() ]^C
    #  p = plot()
    #  for v in varieties
    #        spline = makeCubicPP(x,f,PPInterpolation.SECOND_DERIVATIVE,0.0,PPInterpolation.SECOND_DERIVATIVE,0.0,v)
    #        plot!(p, t, spline.(t), label=string(v))
    #        end
    #  plot!(legend=:top)

end

@testset "FrtschCarlson" begin
    x = [7.99, 8.09, 8.19, 8.7, 9.2, 10.0, 12.0, 15.0, 20.0]
    f = [0.0, 2.76429e-5, 4.37498e-2, 0.169183, 0.469428, 0.943740, 0.998636, 0.999919, 0.999994]

    varieties = [C2Hyman89(), C2MP(), HuynRational(), VanLeer(), FritschButland(), Brodlie(), VanAlbada()]
    for v in varieties
        spline = makeCubicPP(x, f, PPInterpolation.SECOND_DERIVATIVE, 0.0, PPInterpolation.SECOND_DERIVATIVE, 0.0, v)
        nSample = 1000
        previousVal = spline(x[1])
        for i = 1:nSample
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
    for i = 1:nSample
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
@testset "FourHundred" begin
    Random.seed!(1)
    x = sort(rand(500));
    y = rand(500);
    spline = makeCubicPP(x, y, PPInterpolation.SECOND_DERIVATIVE, 0.0, PPInterpolation.SECOND_DERIVATIVE, 0.0, C2())
    z = collect(range(start=x[1],stop=x[end],length=500));
    yNew = Array{Float64}(undef,length(x));
    @time yNew[1:end]= spline(z)
    @time for i=1:length(z) yNew[i] = spline(z[i]) end
    q = PPInterpolation.QuadraticLagrangeInterpolation(x,y)
    @time PPInterpolation.evaluateMid!(yNew, q, z)

end