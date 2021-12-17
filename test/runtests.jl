using PPInterpolation, Test
using ForwardDiff

@testset "HaganCurveAuto" begin
    #for quintic: first pass compute cubic first derivative, second pass, use cubic b, and optimize second derivative.
    t = [0.0, 0.1, 1.0, 2.0, 3.0, 4.0, 9.0, 20.0, 30.0]
    y = [0.06, 0.06, 0.06, 0.05, 0.05, 0.05, 0.065, 0.06, 0.06]
    t = [0.0, 0.1, 1.0, 4.0, 9.0, 20.0, 30.0]
    y = [0.081, 0.081, 0.07, 0.044, 0.07, 0.04, 0.03]
    t=[0, 0.1, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0, 3.5]
    y=[0.012, 0.008, 0.0104, 0.00995, 0.010633333333333333 ,0.010849999999999999, 0.010179999999999998, 0.009816666666666665, 0.008612499999999999, 0.007789999999999998, 0.007741666666666665, 0.007635714285714285]

    # t = [0.0, 0.1, 1, 2, 3, 4, 5, 6, 7, 30.0]
    # f = [0.01, 0.01, 0.01, 0.01, 0.01, 0.03, 0.03, 0.03, 0.03, 0.02]
    # y = zeros(length(t))
    # y[1] = f[1]
    # y[2:end] = (y[1:end-1] .* t[1:end-1] + (t[2:end]-t[1:end-1]) .* f[2:end]) ./ t[2:end]
    # y[3:end-2] = f[3:end-2] ./ t[3:end-2]
    # println(y)
    yt = zeros(length(t))
    yt[1:end] = y[1:end] .* t[1:end]
    spline = makeCubicPP(t,yt,PPInterpolation.second_derivative,0.0,PPInterpolation.second_derivative,0.0,C2())
    nSample = 1000
    for i=1:nSample
        ti = (t[end]) * (i) / (nSample)
        val = evaluateDerivative(spline, ti)
        zval = spline(ti)/ti
        println(ti," CubicPP ", val," ",zval)
        autoval = ForwardDiff.derivative(spline, ti)
        @test isapprox(val, autoval, atol=1e-15)
    end
    #AUTODIFF // where is autodiff useful (hagan wave read),  BigFloat
end

@testset "HaganCurveLavery" begin
    #for quintic: first pass compute cubic first derivative, second pass, use cubic b, and optimize second derivative.
    t = [0.0, 0.1, 1.0, 2.0, 3.0, 4.0, 9.0, 20.0, 30.0]
    y = [0.06, 0.06, 0.06, 0.05, 0.05, 0.05, 0.065, 0.06, 0.06]
    t = [0.0, 0.1, 1.0, 4.0, 9.0, 20.0, 30.0]
    y = [0.081, 0.081, 0.07, 0.044, 0.07, 0.04, 0.03]
    t=[0, 0.1, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0, 3.5]
    y=[0.012, 0.008, 0.0104, 0.00995, 0.010633333333333333 ,0.010849999999999999, 0.010179999999999998, 0.009816666666666665, 0.008612499999999999, 0.007789999999999998, 0.007741666666666665, 0.007635714285714285]

    # t = [0.0, 0.1, 1, 2, 3, 4, 5, 6, 7, 30.0]
    # f = [0.01, 0.01, 0.01, 0.01, 0.01, 0.03, 0.03, 0.03, 0.03, 0.02]
    # y = zeros(length(t))
    # y[1] = f[1]
    # y[2:end] = (y[1:end-1] .* t[1:end-1] + (t[2:end]-t[1:end-1]) .* f[2:end]) ./ t[2:end]
    # y[3:end-2] = f[3:end-2] ./ t[3:end-2]
    # println(y)
    yt = zeros(length(t))
    yt[1:end] = y[1:end] .* t[1:end]
    spline = newQuinticLaverySpline(t,yt,1e-3,1e-2,3)
    nSample = 1000
    for i=1:nSample
        ti = (t[end]) * (i) / (nSample)
        val = evaluateDerivative(spline, ti)
        zval = evaluate(spline, ti)/ti
        println(ti," QuinticLavery3 ", val," ",zval)
    end
end
