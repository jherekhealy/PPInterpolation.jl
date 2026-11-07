# Generates the figures of README.md into docs/img.
# Usage: julia --project=examples examples/readme_plots.jl
using PPInterpolation, CairoMakie

const OUT = joinpath(@__DIR__, "..", "docs", "img")
mkpath(OUT)

const SERIES = [colorant"#2a78d6", colorant"#eb6834", colorant"#1baf7a"]
const INK = colorant"#0b0b0b"
const MUTED = colorant"#8a8985"

set_theme!(Theme(
    fontsize = 15,
    textcolor = INK,
    backgroundcolor = :white,
    Axis = (
        xgridcolor = (:black, 0.08), ygridcolor = (:black, 0.08),
        xminorgridvisible = false, yminorgridvisible = false,
        topspinevisible = false, rightspinevisible = false,
        leftspinecolor = MUTED, bottomspinecolor = MUTED,
        xtickcolor = MUTED, ytickcolor = MUTED,
        titlealign = :left, titlesize = 17,
    ),
    Legend = (framevisible = false, labelsize = 14),
    Lines = (linewidth = 2,),
))

#data with a sharp increase, where a C2 spline overshoots
const x = [3.0, 5.0, 6.0, 8.0, 9.0, 11.0, 12.0, 14.0, 15.0]
const f = [10.0, 10.0, 10.0, 10.0, 10.5, 15.0, 50.0, 60.0, 85.0]
const z = range(x[1], x[end], length = 1000)

function comparison(file, title, splines)
    fig = Figure(size = (760, 430))
    ax = Axis(fig[1, 1], title = title, xlabel = "x", ylabel = "f(x)")
    for (i, (label, spl)) in enumerate(splines)
        lines!(ax, z, spl.(z), color = SERIES[i], label = label)
    end
    scatter!(ax, x, f, color = INK, markersize = 9, label = "data",
        strokecolor = :white, strokewidth = 1.5)
    axislegend(ax, position = :lt)
    save(joinpath(OUT, file), fig, px_per_unit = 2)
end

cubic(kind) = makeCubicPP(x, f, PPInterpolation.SECOND_DERIVATIVE, 0.0, PPInterpolation.SECOND_DERIVATIVE, 0.0, kind)

comparison("cubic_c2.png", "C2 cubic spline and its monotone filtered variant", [
    "C2()" => cubic(C2()),
    "C2Hyman89()" => cubic(C2Hyman89()),
])

comparison("cubic_c1.png", "C1 Hermite splines with derivative limiters", [
    "Bessel()" => cubic(Bessel()),
    "VanAlbada()" => cubic(VanAlbada()),
    "FritschButland()" => cubic(FritschButland()),
])

comparison("quadratic.png", "Quadratic splines", [
    "QuadraticLagrangePP" => QuadraticLagrangePP(x, f),
    "SchumakerQuadraticPP" => SchumakerQuadraticPP(x, f),
    "SchumakerQuadraticPP, LamDerivative(0.5)" => SchumakerQuadraticPP(x, f, LamDerivative(0.5)),
])

#quintic versus cubic Hermite spline, with exact first derivatives, on g(x) = exp(-x/2) sin(2x)
g(t) = exp(-t / 2) * sin(2t)
g1(t) = exp(-t / 2) * (2cos(2t) - sin(2t) / 2)
g2(t) = exp(-t / 2) * (-15sin(2t) / 4 - 2cos(2t))
xq = [0.0, 0.4, 1.0, 1.5, 2.2, 2.8, 3.5, 4.0]
yq = g.(xq)
y1q = g1.(xq)
quintic = makeHermiteQuinticPP(xq, yq, y1q, g2(xq[1]), g2(xq[end]))
hermite = makeCubicPP(xq, yq, PPInterpolation.FIRST_DERIVATIVE, y1q[1], PPInterpolation.FIRST_DERIVATIVE, y1q[end], Hermite(y1q))
zq = range(xq[1], xq[end] - 1e-9, length = 1000)

fig = Figure(size = (760, 430))
ax = Axis(fig[1, 1], title = "Second derivative of Hermite splines with exact first derivatives",
    xlabel = "x", ylabel = "second derivative")
lines!(ax, zq, [evaluateSecondDerivative(hermite, t) for t in zq], color = SERIES[1], label = "cubic Hermite (C1)")
lines!(ax, zq, [evaluateSecondDerivative(quintic, t) for t in zq], color = SERIES[2], label = "makeHermiteQuinticPP (C3)")
lines!(ax, zq, g2.(zq), color = INK, linestyle = :dash, linewidth = 1.5, label = "exact")
vlines!(ax, xq, color = (:black, 0.12), linewidth = 1)
axislegend(ax, position = :rb)
save(joinpath(OUT, "quintic.png"), fig, px_per_unit = 2)

#boundary conditions of the C2 cubic spline
xb = [0.0, 1.0, 2.0, 3.5, 4.0, 5.0]
yb = [1.0, 2.2, 2.0, 0.8, 0.9, 1.6]
zb = range(xb[1], xb[end], length = 1000)
fig = Figure(size = (760, 430))
ax = Axis(fig[1, 1], title = "Boundary conditions of the C2 cubic spline", xlabel = "x", ylabel = "f(x)")
for (i, (label, bc, value)) in enumerate([
    ("SECOND_DERIVATIVE, 0 (natural)", PPInterpolation.SECOND_DERIVATIVE, 0.0),
    ("NOT_A_KNOT", PPInterpolation.NOT_A_KNOT, 0.0),
    ("FIRST_DERIVATIVE, 0 (clamped)", PPInterpolation.FIRST_DERIVATIVE, 0.0),
])
    spl = makeCubicPP(xb, yb, bc, value, bc, value, C2())
    lines!(ax, zb, spl.(zb), color = SERIES[i], label = label)
end
scatter!(ax, xb, yb, color = INK, markersize = 9, label = "data", strokecolor = :white, strokewidth = 1.5)
axislegend(ax, position = :rt)
save(joinpath(OUT, "boundary.png"), fig, px_per_unit = 2)

#extrapolation beyond the first and last knots
spl = makeCubicPP(xb, yb, PPInterpolation.NOT_A_KNOT, 0.0, PPInterpolation.NOT_A_KNOT, 0.0, C2())
ze = range(xb[1] - 1.5, xb[end] + 1.5, length = 1000)
fig = Figure(size = (760, 430))
ax = Axis(fig[1, 1], title = "Extrapolation of a cubic spline", xlabel = "x", ylabel = "f(x)")
vspan!(ax, [ze[1], xb[end]], [xb[1], ze[end]], color = (:black, 0.04))
for (i, (label, ext)) in enumerate([
    ("LinearAutoExtrapolation() (default)", LinearAutoExtrapolation()),
    ("ConstantAutoExtrapolation()", ConstantAutoExtrapolation()),
    ("PieceAutoExtrapolation()", PieceAutoExtrapolation()),
])
    lines!(ax, ze, [evaluate(spl, t, extrapolation = ext) for t in ze], color = SERIES[i], label = label)
end
scatter!(ax, xb, yb, color = INK, markersize = 9, label = "data", strokecolor = :white, strokewidth = 1.5)
axislegend(ax, position = :rb)
save(joinpath(OUT, "extrapolation.png"), fig, px_per_unit = 2)
