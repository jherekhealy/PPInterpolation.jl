using Test, BSplines, PPInterpolation

@testset "BSplinesConversion" begin
    Random.seed!(1234);
    x = sort(rand(10))
    a = rand(12)
    basis = BSplineBasis(4,x)
    t = BSplines.knots(basis)
    spl = BSplines.Spline(basis,a)
    pp = convert(PPInterpolation.PP{3,Float64,Float64}, spl)
    spl2 = convert(BSplines.Spline,pp)
    xt = range(start=x[1]*0.75, stop=x[end]*1.25, length=101)
    y = spl.(xt)
    y2 = spl2.(xt)
    @test isapprox(y,y2,atol=1e-15)        
end