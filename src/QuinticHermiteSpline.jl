using LinearAlgebra
#C^3 quintic spline with given first derivatives at all knots, plus second derivatives at end-points.
export makeHermiteQuinticPP

function computeDEF(
    c::AbstractMatrix{T},
    f1::Vector{T},
    f2::Vector{T},
    S::Vector{T},
    dx::Vector{T}
) where {T}
    for i = 1:length(f1)-1
        c[1, i] = f2[i] / 2
        c[2, i] = (f2[i+1] - 3 * f2[i]) / (2 * dx[i]) + 2 * (5 * S[i] - 3 * f1[i] - 2 * f1[i+1]) / dx[i]^2
        c[3, i] = (3 * f2[i] - 2 * f2[i+1]) / (2 * dx[i]^2) + (8 * f1[i] + 7 * f1[i+1] - 15 * S[i]) / dx[i]^3
        c[4, i] = (f2[i+1] - f2[i]) / (2 * dx[i]^3) + 3 * (2 * S[i] - f1[i+1] - f1[i]) / dx[i]^4
    end
end


function makeHermiteQuinticPP(x::Vector{TX}, y::Vector{T}, y1::Vector{T}, y2left::T, y2right::T) where {T,TX}
    dx = @. x[2:end] - x[1:end-1]
    s = @. (y[2:end] - y[1:end-1]) / dx
    n = length(y)
    rhs = zeros(T, n)
    dd = zeros(TX, n)
    dl = zeros(TX, n - 1)
    du = zeros(TX, n - 1)
    rhs[1] = y2left
    dd[1] = one(TX)
    dd[n] = one(TX)
    rhs[n] = y2right
    for i = 2:n-1
        rhs[i] = 2 * (5 * s[i] - 3 * y1[i] - 2 * y1[i+1]) / dx[i]^2 -
                 2 * (5 * s[i-1] - 3 * y1[i-1] - 2 * y1[i]) / dx[i-1]^2 -
                 4 * (8 * y1[i-1] + 7 * y1[i] - 15 * s[i-1]) / dx[i-1]^2 -
                 30 * (2 * s[i-1] - y1[i] - y1[i-1]) / dx[i-1]^2
        dd[i] = 3 / (2 * dx[i]) + 3 / (2 * dx[i-1])
        dl[i-1] = -1 / (2 * dx[i-1])
        du[i] = -1 / (2 * dx[i])
    end
    tri = Tridiagonal(dl, dd, du)
    y2 = tri \ rhs
    c = zeros(T, (4, n))
    #q''(xi+1)=f''i+1
    computeDEF(c, y1, y2, s, dx)
    return PP(5, copy(y), copy(y1), c, copy(x))
end

function makeHermiteQuinticPP(x::Vector{TX}, y::Vector{T}, y1::Vector{T}) where {T,TX}
    dx = @. x[2:end] - x[1:end-1]
    s = @. (y[2:end] - y[1:end-1]) / dx
    n = length(y)
    rhs = zeros(T, n)
    d = zeros(TX.n)
    dl = zeros(TX, n - 1)
    du = zeros(TX, n - 1)
    rhs[1] = zero(T)
    d[1] = one(TX)
    du[1] = -one(TX) #y''(x0)= y''(x1)
    d[n] = one(TX)
    dl[n-1] = -one(TX)
    rhs[n] = 0
    for i = 2:n-1
        rhs[i] = 2 * (5 * s[i] - 3 * y1[i] - 2 * y1[i+1]) / dx[i]^2 -
                 2 * (5 * s[i-1] - 3 * y1[i-1] - 2 * y1[i]) / dx[i-1]^2 -
                 4 * (8 * y1[i-1] + 7 * y1[i] - 15 * s[i-1]) / dx[i-1]^2 -
                 30 * (2 * s[i-1] - y1[i] - y1[i-1]) / dx[i-1]^2
        d[i] = 3 / (2 * dx[i]) + 3 / (2 * dx[i-1])
        dl[i-1] = -1 / (2 * dx[i-1])
        du[i] = -1 / (2 * dx[i])
    end
    tri = Tridiagonal(dl, d, du)
    y2 = tri \ rhs
    c = zeros(T, (4, n))
    #q''(xi+1)=f''i+1
    computeDEF(c, y1, y2, s, dx)
    return PP(5, copy(y), copy(y1), c, copy(x))
end

function evaluate(self::PP{5,T,TX}, z::TZ) where {T,TX,TZ}
    if z <= self.x[1]
        return self.b[1] * (z - self.x[1]) + self.a[1]
    elseif z >= self.x[end]
        rightSlope = self.b[end]
        return rightSlope * (z - self.x[end]) + self.a[end]
    end
    i = searchsortedfirst(self.x, z)  # x[i-1]<z<=x[i]
    if z != self.x[i] && i > 1
        i -= 1
    end
    h = z - self.x[i]
    return self.a[i] + h * (self.b[i] + h * (self.c[1, i] + h * (self.c[2, i] + h * (self.c[3, i] + h * self.c[4, i]))))
end

function evaluateDerivative(self::PP{5,T,TX}, z::TZ) where {T,TX,TZ}
    if z <= self.x[1]
        return self.b[1]
    elseif z >= self.x[end]
        rightSlope = self.b[end]
        return rightSlope
    end
    i = searchsortedfirst(self.x, z)  # x[i-1]<z<=x[i]
    if z != self.x[i] && i > 1
        i -= 1
    end
    h = z - self.x[i]
    return self.b[i] + h * (2 * self.c[1, i] + h * (3 * self.c[2, i] + h * (4 * self.c[3, i] + h * 5 * self.c[4, i])))
end
function evaluateSecondDerivative(self::PP{5,T,TX}, z::TZ) where {T,TX,TZ}
    if z <= self.x[1]
        return zero(TZ)
    elseif z >= self.x[end]
        return zero(TZ)
    end
    if z != self.x[i] && i > 1
        i -= 1
    end
    h = z - self.x[i]
    return 2 * self.c[1, i] + h * (3 * 2 * self.c[2, i] + h * (4 * 3 * self.c[3, i] + h * 5 * 4 * self.c[4, i]))
end
function evaluateThirdDerivative(self::PP{5,T,TX}, z::TZ) where {T,TX,TZ}
    if z <= self.x[1]
        return zero(TZ)
    elseif z >= self.x[end]
        return zero(TZ)
    end
    if z != self.x[i] && i > 1
        i -= 1
    end
    h = z - self.x[i]
    return ((3 * 2 * self.c[2, i] + h * (4 * 3 * 2 * self.c[3, i] + h * 5 * 4 * 3 * self.c[4, i])))
end
