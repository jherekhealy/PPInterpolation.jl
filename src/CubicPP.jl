using LinearAlgebra

export CubicPP, makeLinearCubicPP, makeC2CubicPP, makeHermiteCubicPP

struct CubicPP{T<:Real}
    a::Vector{T}
    b::Vector{T}
    c::Vector{T}
    d::Vector{T}
    x::Vector{T}
end

# abstract type PPBoundary <: Real end
# struct FirstDerivativeBoundary <: PPBoundary end
# struct SecondDerivativeBoundary <: PPBoundary end
@enum PPBoundary not_a_knot = 0 first_derivative = 1 second_derivative = 2

(spl::CubicPP)(x::Int) = CubicPP(zeros(x), zeros(x),zeros(x),zeros(x),zeros(x))

function makeLinearCubicPP(x::Vector{T}, y::Vector{T}) where {T}
    n = length(x)
    if n <= 1
        return CubicPP(copy(y), zeros(n), zeros(n), zeros(n), copy(x))
    elseif n == 2
        t = y[2] - y[1]
        b = zeros(n)
        if abs(x[2] - x[1]) > eps()
            b[1] = t / (x[2] - x[1])
        end
        b[2] = b[1]
        return CubicPP(copy(y), b, zeros(n), zeros(n), copy(x))
    else
		# on xi, xi+1, f(x)= yi (xi+1-x) + yi (x-xi) = A + B (x-xi) => B = (yi-yi+1)/(xi-xi+1)
        b = zeros(n)
        for i = 1:n-1
            b[i] = (y[i+1] - y[i]) / (x[i+1] - x[i])
        end
        return CubicPP(copy(y), b, zeros(n), zeros(n), copy(x))
    end
end

function makeC2CubicPP(
    x::Vector{T},
    y::Vector{T},
    leftBoundary::PPBoundary,
    leftValue::T,
    rightBoundary::PPBoundary,
    rightValue::T
) where {T}
  pp = CubicPP(length(y))
  computeC2CubicPP(pp, x, y, leftBoundary, leftValue, rightBoundary, rightValue)
  return pp
end

function computeC2CubicPP(pp:CubicPP{T},
    x::Vector{T},
    y::Vector{T},
    leftBoundary::PPBoundary,
    leftValue::T,
    rightBoundary::PPBoundary,
    rightValue::T
) where {T}
    n = length(y)
    if n <= 2
        return makeLinearCubicPP(x, y)
    end

    dx = x[2:end] - x[1:end-1]
    S = (y[2:end] - y[1:end-1]) ./ dx
    middle = zeros(n)
    alpha = zeros(n)
    lower = zeros(n - 1)
    upper = zeros(n - 1)
    for i = 2:n-1
        lower[i-1] = dx[i]
        upper[i] = dx[i-1]
        middle[i] = 2 * (dx[i] + dx[i-1])
        alpha[i] = 3.0 * (dx[i] * S[i-1] + dx[i-1] * S[i])
    end
    #middle[2:n-1] = 3 * (dx[2:n-1] + dx[1:n-2])
    #alpha[2:n-1] = 3 * (dx[2:n-1] .* S[1:n-2] + dx[1:n-2] .* S[2:n-1])
    if leftBoundary == not_a_knot
        middle[1] = dx[2] * (dx[2] + dx[1])
        upper[1] = (dx[2] + dx[1]) * (dx[2] + dx[1])
        alpha[1] = S[1] * dx[2] * (2.0 * dx[2] + 3.0 * dx[1]) + S[2] * dx[1]^2
    elseif leftBoundary == first_derivative
        middle[1] = 1.0
        upper[1] = 0.0
        alpha[1] = leftValue
    elseif leftBoundary == second_derivative
        middle[1] = 2.0
        upper[1] = 1.0
        alpha[1] = 3 * S[1] - leftValue * dx[1] / 2
    end
    if rightBoundary == not_a_knot
        lower[n-1] = -(dx[n-1] + dx[n-2]) * (dx[n-1] + dx[n-2])
        middle[n] = -dx[n-2] * (dx[n-2] + dx[n-1])
        alpha[n] = -S[n-2] * dx[n-1]^2 - S[n-1] * dx[n-2] * (3.0 * dx[n-1] + 2.0 * dx[n-2])
    elseif rightBoundary == first_derivative
        middle[n] = 1.0
        lower[n-1] = 0.0
        alpha[n] = rightValue
    elseif rightBoundary == second_derivative
        middle[n] = 2.0
        lower[n-1] = 1.0
        alpha[n] = 3 * S[n-1] - rightValue * dx[n-1] / 2
    end
    tri = LinearAlgebra.Tridiagonal(lower, middle, upper)
    fPrime = tri \ alpha
    c = (3 * S - fPrime[2:end] - 2 * fPrime[1:end-1]) ./ dx
    d = (fPrime[2:end] + fPrime[1:end-1] - 2 * S) ./ (dx.^2)

    pp.a[1:end] = y
	pp.b[1:end] = fPrime
	pp.c[1:end] = c
	pp.d[1:end] = d
	pp.x[1:end] = x
end

function makeHermiteCubicPP(
    x::Vector{T},
    y::Vector{T},
    y1::Vector{T},
    leftBoundary::PPBoundary,
    leftValue::T,
    rightBoundary::PPBoundary,
    rightValue::T
) where {T}
    n = length(y)
    if n <= 2
        return makeLinearCubicPP(x, y)
    end

    dx = x[2:end] - x[1:end-1]
    S = (y[2:end] - y[1:end-1]) ./ dx
    b = copy(y1)
    if leftBoundary == not_a_knot
        b[1] = S[2] * dx[1] / (dx[2] * (dx[2] + dx[1])) - S[1] * ((dx[2] / dx[1] + 2) / (dx[2] + dx[1]))
    elseif leftBoundary == first_derivative
        b[1] = leftValue
    elseif leftBoundary == second_derivative
        #c[1] = leftValue * 0.5
        b[1] = (-leftValue * 0.5 * dx[1] - b[2] + 3 * S[1]) * 0.5
    end
    if rightBoundary == not_a_knot
        b[n] = S[n-2] * dx[n-1] / (dx[n-2] * (dx[n-2] + dx[n-1])) -
               S[n-1] * ((dx[n-2] / dx[n-1] + 2) / (dx[n-2] + dx[n-1]))
    elseif rightBoundary == first_derivative
        b[n] = rightValue
    elseif rightBoundary == second_derivative
        b[n] = (rightValue * dx[n-1] + 6 * S[n-1] - 2 * b[n-1]) / 4
    end
    c = (3 * S - b[2:end] - 2 * b[1:end-1]) ./ dx
    d = (b[2:end] + b[1:end-1] - 2 * S) ./ (dx.^2)

    return CubicPP(copy(y), b, c, d, copy(x))
end

function evaluate(self::CubicPP, z::T) where {T}
    if z <= self.x[1]
        return self.b[1] * (z - self.x[1]) + self.a[1]
    elseif z >= self.x[end]
        rightSlope = self.b[end]
        return rightSlope * (z - self.x[end]) + self.a[end]
    end
    i = searchsortedfirst(self.x, z)  # x[i-1]<z<=x[i]
    if z == self.x[i]
        return self.a[i]
    end
    if i > 1
        i -= 1
    end
    h = z - self.x[i]
    return self.a[i] + h * (self.b[i] + h * (self.c[i] + h * (self.d[i])))
end

function evaluateDerivative(self::CubicPP, z::T) where {T}
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
    return self.b[i] + h * (2 * self.c[i] + h * (3 * self.d[i]))
end
function evaluateSecondDerivative(self::CubicPP, z::T) where {T}
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
    return 2 * self.c[i] + h * (3 * 2 * self.d[i])
end
