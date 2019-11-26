using LinearAlgebra
#C^3 quintic spline with given first derivatives at all knots, plus second derivatives at end-points.

struct QuinticPP{T<:Real}
    a::Vector{T}
    b::Vector{T}
    c::Vector{T}
    d::Vector{T}
    e::Vector{T}
    f::Vector{T}
    x::Vector{T}
end
function computeDEF(
    d::Vector{T},
    e::Vector{T},
    f::Vector{T},
    f1::Vector{T},
    f2::Vector{T},
    S::Vector{T},
    dx::Vector{T}
) where {T}
    for i = 1:length(f1)-1
        d[i] = (f2[i+1] - 3 * f2[i]) / (2 * dx[i]) + 2 * (5 * S[i] - 3 * f1[i] - 2 * f1[i+1]) / dx[i]^2
        e[i] = (3 * f2[i] - 2 * f2[i+1]) / (2 * dx[i]^2) + (8 * f1[i] + 7 * f1[i+1] - 15 * S[i]) / dx[i]^3
        f[i] = (f2[i+1] - f2[i]) / (2 * dx[i]^3) + 3 * (2 * S[i] - f1[i+1] - f1[i]) / dx[i]^4
    end
end


function makeHermiteQuinticPP(x::Vector{T}, y::Vector{T}, y1::Vector{T}, y2left::T, y2right::T) where {T}
    dx = @. x[2:end] - x[1:end-1]
    s = @. (y[2:end] - y[1:end-1]) / dx
    n = length(y)
    rhs = zeros(n)
    dd = zeros(n)
    dl = zeros(n - 1)
    du = zeros(n - 1)
    rhs[1] = y2left
    dd[1] = 1
    dd[n] = 1
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
    d = zeros(n)
    e = zeros(n)
    f = zeros(n)
    #q''(xi+1)=f''i+1
    for i = 1:n-1
        d[i] = (y2[i+1] - 3 * y2[i]) / (2 * dx[i]) + 2 * (5 * s[i] - 3 * y1[i] - 2 * y1[i+1]) / dx[i]^2
        e[i] = (3 * y2[i] - 2 * y2[i+1]) / (2 * dx[i]^2) + (8 * y1[i] + 7 * y1[i+1] - 15 * s[i]) / dx[i]^3
        f[i] = (y2[i+1] - y2[i]) / (2 * dx[i]^3) + 3 * (2 * s[i] - y1[i+1] - y1[i]) / dx[i]^4
    end
    return QuinticPP(copy(y), copy(y1), y2 / 2, d, e, f, copy(x))
end

function makeHermiteQuinticPP(x::Vector{T}, y::Vector{T}, y1::Vector{T}) where {T}
    dx = @. x[2:end] - x[1:end-1]
    s = @. (y[2:end] - y[1:end-1]) / dx
    n = length(y)
    rhs = zeros(n)
    d = zeros(n)
    dl = zeros(n - 1)
    du = zeros(n - 1)
    rhs[1] = 0
    d[1] = 1
    du[1] = -1 #y''(x0)= y''(x1)
    d[n] = 1
    dl[n-1] = -1
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
    a = copy(y)
    d = zeros(length(a))
    e = zeros(length(a))
    f = zeros(length(a))
    computeDEF(d, e, f, y1, y2, s, dx)
    return QuinticPP(a, copy(y1), y2 / 2, d, e, f, copy(x))
end
function evaluate(self::QuinticPP, z::T) where {T}
    if z <= self.x[1]
        return self.b[1] * (z - self.x[1]) + self.a[1]
    elseif z >= self.x[end]
        rightSlope = self.b[end]
        return rightSlope * (z - self.x[end]) + self.a[end]
    end
    i = searchsortedfirst(self.x, z)  # x[i-1]<z<=x[i]
    if z==self.x[i]
        return self.a[i]
    end
    if i > 1
        i -= 1
    end
    h = z - self.x[i]
    return self.a[i] + h * (self.b[i] + h * (self.c[i] + h * (self.d[i] + h * (self.e[i] + h * self.f[i]))))
end

function evaluateDerivative(self::QuinticPP, z::T) where {T}
    if z <= self.x[1]
        return self.b[1]
    elseif z >= self.x[end]
        rightSlope = self.b[end]
        return rightSlope
    end
    i = searchsortedfirst(self.x, z)  # x[i-1]<z<=x[i]
    if i > 1
        i -= 1
    end
    h = z - self.x[i]
    return self.b[i] + h * (2 * self.c[i] + h * (3 * self.d[i] + h * (4 * self.e[i] + h * 5 * self.f[i])))
end
function evaluateSecondDerivative(self::QuinticPP, z::T) where {T}
    if z <= self.x[1]
        return self.b[1]
    elseif z >= self.x[end]
        rightSlope = self.b[end]
        return rightSlope
    end
    i = searchsortedfirst(self.x, z)  # x[i-1]<z<=x[i]
    if i > 1
        i -= 1
    end
    h = z - self.x[i]
    return 2 * self.c[i] + h * (3 * 2 * self.d[i] + h * (4 * 3 * self.e[i] + h * 5 * 4 * self.f[i]))
end
function evaluateThirdDerivative(self::QuinticPP, z::T) where {T}
    if z <= self.x[1]
        return self.b[1]
    elseif z >= self.x[end]
        rightSlope = self.b[end]
        return rightSlope
    end
    i = searchsortedfirst(self.x, z)  # x[i-1]<z<=x[i]
    if i > 1
        i -= 1
    end
    h = z - self.x[i]
    return ( (3 * 2 * self.d[i] + h * (4 * 3 * 2*self.e[i] + h * 5 * 4 *3* self.f[i])))
end
