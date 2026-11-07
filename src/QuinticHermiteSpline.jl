using LinearAlgebra
#C^3 quintic spline with given first derivatives at all knots, plus second derivatives at end-points.
export makeHermiteQuinticPP

function computeDEF(
    c::AbstractMatrix{T},
    f1::AbstractVector{T},
    f2::AbstractVector{T},
    S::AbstractVector{T},
    dx::AbstractVector{T}
) where {T}
    for i = 1:length(f1)-1
        c[i, 1] = f2[i] / 2
        c[i, 2] = (f2[i+1] - 3 * f2[i]) / (2 * dx[i]) + 2 * (5 * S[i] - 3 * f1[i] - 2 * f1[i+1]) / dx[i]^2
        c[i, 3] = (3 * f2[i] - 2 * f2[i+1]) / (2 * dx[i]^2) + (8 * f1[i] + 7 * f1[i+1] - 15 * S[i]) / dx[i]^3
        c[i, 4] = (f2[i+1] - f2[i]) / (2 * dx[i]^3) + 3 * (2 * S[i] - f1[i+1] - f1[i]) / dx[i]^4
    end
end


#Tridiagonal system for the second derivatives at the interior knots, and the slopes. The boundary rows 1 and n are left for the caller.
function quinticSystem(x::AbstractVector{TX}, y::AbstractVector{T}, y1::AbstractVector{T}) where {T,TX}
    dx = @. x[2:end] - x[1:end-1]
    s = @. (y[2:end] - y[1:end-1]) / dx
    n = length(y)
    rhs = zeros(T, n)
    d = zeros(TX, n)
    dl = zeros(TX, n - 1)
    du = zeros(TX, n - 1)
    for i = 2:n-1
        rhs[i] = 2 * (5 * s[i] - 3 * y1[i] - 2 * y1[i+1]) / dx[i]^2 -
                 2 * (5 * s[i-1] - 3 * y1[i-1] - 2 * y1[i]) / dx[i-1]^2 -
                 4 * (8 * y1[i-1] + 7 * y1[i] - 15 * s[i-1]) / dx[i-1]^2 -
                 30 * (2 * s[i-1] - y1[i] - y1[i-1]) / dx[i-1]^2
        d[i] = 3 / (2 * dx[i]) + 3 / (2 * dx[i-1])
        dl[i-1] = -1 / (2 * dx[i-1])
        du[i] = -1 / (2 * dx[i])
    end
    return dx, s, rhs, dl, d, du
end

function makeQuinticPP(x::AbstractVector{TX}, y::AbstractVector{T}, y1::AbstractVector{T}, dx, s, rhs, dl, d, du, extrapolation::PPExtrapolation) where {T,TX}
    n = length(y)
    y2 = Tridiagonal(dl, d, du) \ rhs
    c = zeros(T, (n - 1, 4))
    #q''(xi+1)=f''i+1
    computeDEF(c, y1, y2, s, dx)
    return PP(5, copy(y), copy(y1), c, copy(x), extrapolation = extrapolation)
end

"""
    makeHermiteQuinticPP(x, y, y1[, y2left, y2right]; extrapolation = LinearAutoExtrapolation())

C3 quintic spline interpolating the values `y` and first derivatives `y1` at the knots `x`.
The second derivatives at the end-points are `y2left` and `y2right` when given, otherwise equal to the second derivatives at the next knots.
"""
function makeHermiteQuinticPP(x::AbstractVector{TX}, y::AbstractVector{T}, y1::AbstractVector{T}, y2left::T, y2right::T; extrapolation::PPExtrapolation = LinearAutoExtrapolation()) where {T,TX}
    dx, s, rhs, dl, d, du = quinticSystem(x, y, y1)
    n = length(y)
    rhs[1] = y2left
    d[1] = one(TX)
    d[n] = one(TX)
    rhs[n] = y2right
    return makeQuinticPP(x, y, y1, dx, s, rhs, dl, d, du, extrapolation)
end

function makeHermiteQuinticPP(x::AbstractVector{TX}, y::AbstractVector{T}, y1::AbstractVector{T}; extrapolation::PPExtrapolation = LinearAutoExtrapolation()) where {T,TX}
    dx, s, rhs, dl, d, du = quinticSystem(x, y, y1)
    n = length(y)
    d[1] = one(TX)
    du[1] = -one(TX) #y''(x0)= y''(x1)
    d[n] = one(TX)
    dl[n-1] = -one(TX)
    return makeQuinticPP(x, y, y1, dx, s, rhs, dl, d, du, extrapolation)
end
