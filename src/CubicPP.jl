using LinearAlgebra
import SpecialFunctions: erfc
export PP, makeLinearCubicPP, makeCubicPP, C2, C2Hyman89, C2HymanNonNegative, C2MP, Bessel, HuynRational, VanAlbada, VanLeer, FritschButland, Brodlie, Hermite, Fukasawa
export evaluateDerivative, evaluateSecondDerivative, CubicSplineNatural, CubicSplineNotAKnot, evaluateSorted!, evaluateIntegral
export evaluatePiece, evaluateDerivativePiece
export LinearAutoExtrapolation, ConstantAutoExtrapolation, PieceAutoExtrapolation
abstract type DerivativeKind end
struct C2 <: DerivativeKind end
struct C2Hyman89 <: DerivativeKind end
struct C2HymanNonNegative <: DerivativeKind end
struct C2MP <: DerivativeKind end
struct C2MP2 <: DerivativeKind end
@enum PPBoundary NOT_A_KNOT = 0 FIRST_DERIVATIVE = 1 SECOND_DERIVATIVE = 2 FIRST_DIFFERENCE = 3
abstract type PPExtrapolation end
struct ZeroExtrapolation <: PPExtrapolation
end
"Extrapolation is constant beyond first and last knots, using the value at the endpoints"
struct ConstantAutoExtrapolation <: PPExtrapolation
end
"Extrapolation is linear beyond first and last knots, using the slope and value at the endpoints"
struct LinearAutoExtrapolation <: PPExtrapolation
end
struct ConstantExtrapolation <: PPExtrapolation
	value::Any
end
struct LinearExtrapolation <: PPExtrapolation
	slope::Any
end
"Extrapolation is using the first/last piecewise polynomial"
struct PieceAutoExtrapolation <: PPExtrapolation
end

abstract type LimiterDerivative <: DerivativeKind end
"The derivative at each knot is passed as parameter."
struct Hermite{T} <: LimiterDerivative
	b::Vector{T} #derivative values
end
"The derivative at each knot is such that the incidence angle and the reflecting angle are equal as in  Fukasawa et al (2011) 'Model-free implied volatility: from surface to index'."
struct Fukasawa <: LimiterDerivative end
"The derivative at each knot is the one of the 3-points parabola."
struct Bessel <: LimiterDerivative end
"The derivative at each knot is based on the rational limiter of Huyn. In particular it preserves monotonicity."
struct HuynRational <: LimiterDerivative end
struct VanLeer <: LimiterDerivative end
struct VanAlbada <: LimiterDerivative end
"Fritsch Butland 1980"
struct FritschButland <: LimiterDerivative end
"Fritch Butland 1984"
struct Brodlie <: LimiterDerivative end

struct PP{N, T, TX, VT <: AbstractArray{T}, MT <: AbstractMatrix{T}, VX <: AbstractArray{TX}}
	a::VT
	b::VT
	c::MT #c[:,i] = coeff of x^{i+1}
	x::VX
end

PP(N::Int, T, TX, n::Int) = PP{N, T, TX, Vector{T}, Matrix{T}, Vector{TX}}(zeros(T, n), zeros(T, n), zeros(T, (n - 1, N - 1)), zeros(TX, n))
PP(N::Int, a::AbstractArray{T}, b::AbstractArray{T}, c::AbstractMatrix{T}, x::AbstractArray{TX}) where {T <: Real, TX} =
	PP{N, T, TX, typeof(a), typeof(c), typeof(x)}(a, b, c, x)


Base.length(p::PP) = Base.length(p.x)
Base.size(p::PP) = Base.size(p.x)
Base.broadcastable(p::PP) = Ref(p)


function makeLinearPP(x::AbstractArray{TX}, y::AbstractArray{T}) where {T, TX}
	pp = PP(2, T, TX, length(y))
	computeLinearPP(pp, x, y)
	return pp
end

function computeLinearPP(pp::PP{N, T, TX}, x::AbstractArray{TX}, y::AbstractArray{T}) where {N, T, TX}
	n = length(x)
	if n <= 1
		pp.a[1:end] = y
		pp.b[1:end] = zeros(n)
		pp.c[1:end, 1:end] = zeros(n - 1, (N - 1))
		pp.x[1:end] = x
	elseif n == 2
		t = y[2] - y[1]
		b = zeros(n)
		if abs(x[2] - x[1]) > eps()
			b[1] = t / (x[2] - x[1])
		end
		b[2] = b[1]
		pp.a[1:end] = y
		pp.b[1:end] = b
		pp.c[1:end, 1:end] = zeros(n - 1, (N - 1))
		pp.x[1:end] = x
	else
		# on xi, xi+1, f(x)= yi (xi+1-x) + yi (x-xi) = A + B (x-xi) => B = (yi-yi+1)/(xi-xi+1)
		b = zeros(n)
		for i ∈ 1:n-1
			b[i] = (y[i+1] - y[i]) / (x[i+1] - x[i])
		end
		pp.a[1:end] = y
		pp.b[1:end] = b
		pp.c[1:end, 1:end] = zeros(n - 1, (N - 1))
		pp.x[1:end] = x
	end
end

CubicSplineNatural(x::AbstractArray{TX}, y::AbstractArray{T}) where {TX, T} = makeCubicPP(x, y, SECOND_DERIVATIVE, zero(T), SECOND_DERIVATIVE, zero(T), C2())

CubicSplineNotAKnot(x::AbstractArray{TX}, y::AbstractArray{T}) where {TX, T} = makeCubicPP(x, y, NOT_A_KNOT, zero(T), NOT_A_KNOT, zero(T), C2())

function makeCubicPP(
	x::AbstractArray{TX},
	y::AbstractArray{T},
	leftBoundary::PPBoundary,
	leftValue::T,
	rightBoundary::PPBoundary,
	rightValue::T,
	kind::DerivativeKind,
) where {T, TX}
	pp = PP(3, T, TX, length(y))
	computePP(pp, x, y, leftBoundary, leftValue, rightBoundary, rightValue, kind)
	return pp
end

function computePP(
	pp::PP{3, T, TX},
	x::AbstractArray{TX},
	y::AbstractArray{T},
	leftBoundary::PPBoundary,
	leftValue::T,
	rightBoundary::PPBoundary,
	rightValue::T,
	kind::Union{C2, C2Hyman89, C2HymanNonNegative, C2MP, C2MP2},
) where {T, TX}
	n = length(y)
	if n <= 2
		computeLinearPP(pp, x, y)
		return
	end

	dx = x[2:end] - x[1:end-1]
	S = (y[2:end] - y[1:end-1]) ./ dx
	middle = zeros(TX, n)
	alpha = zeros(T, n)
	lower = zeros(TX, n - 1)
	upper = zeros(TX, n - 1)
	for i ∈ 2:n-1
		lower[i-1] = dx[i]
		upper[i] = dx[i-1]
		middle[i] = 2 * (dx[i] + dx[i-1])
		alpha[i] = 3 * (dx[i] * S[i-1] + dx[i-1] * S[i])
	end
	#middle[2:n-1] = 3 * (dx[2:n-1] + dx[1:n-2])
	#alpha[2:n-1] = 3 * (dx[2:n-1] .* S[1:n-2] + dx[1:n-2] .* S[2:n-1])
	if leftBoundary == NOT_A_KNOT
		middle[1] = dx[2] * (dx[2] + dx[1])
		upper[1] = (dx[2] + dx[1]) * (dx[2] + dx[1])
		alpha[1] = S[1] * dx[2] * (2 * dx[2] + 3 * dx[1]) + S[2] * dx[1]^2
	elseif leftBoundary == FIRST_DERIVATIVE
		middle[1] = one(TX)
		upper[1] = zero(TX)
		alpha[1] = leftValue
	elseif leftBoundary == FIRST_DIFFERENCE
		middle[1] = one(TX)
		upper[1] = zero(TX)
		alpha[1] = S[1]
	elseif leftBoundary == SECOND_DERIVATIVE
		middle[1] = 2 * one(TX)
		upper[1] = one(TX)
		alpha[1] = 3 * S[1] - leftValue * dx[1] / 2
	end
	if rightBoundary == NOT_A_KNOT
		lower[n-1] = -(dx[n-1] + dx[n-2]) * (dx[n-1] + dx[n-2])
		middle[n] = -dx[n-2] * (dx[n-2] + dx[n-1])
		alpha[n] = -S[n-2] * dx[n-1]^2 - S[n-1] * dx[n-2] * (3 * dx[n-1] + 2 * dx[n-2])
	elseif rightBoundary == FIRST_DERIVATIVE
		middle[n] = one(TX)
		lower[n-1] = zero(TX)
		alpha[n] = rightValue
	elseif rightBoundary == FIRST_DIFFERENCE
		middle[n] = one(TX)
		lower[n-1] = zero(TX)
		alpha[n] = S[n-1]
	elseif rightBoundary == SECOND_DERIVATIVE
		middle[n] = 2 * one(TX)
		lower[n-1] = one(TX)
		alpha[n] = 3 * S[n-1] + rightValue * dx[n-1] / 2
	end
	tri = LinearAlgebra.Tridiagonal(lower, middle, upper)
	fPrime = tri \ alpha
	filterSlope(kind, y, fPrime, dx, S)
	c = pp.c
	for i ∈ 1:n-1
		c[i, 1] = (3 * S[i] - fPrime[i+1] - 2 * fPrime[i]) / dx[i]
		c[i, 2] = (fPrime[i+1] + fPrime[i] - 2 * S[i]) / (dx[i]^2)
	end
	pp.a[1:end] = y
	pp.b[1:end] = fPrime
	pp.x[1:end] = x
end


function computePP(
	pp::PP{3, T, TX},
	x::AbstractArray{TX},
	y::AbstractArray{T},
	leftBoundary::PPBoundary,
	leftValue::T,
	rightBoundary::PPBoundary,
	rightValue::T,
	limiter::LimiterDerivative,
) where {T, TX}
	n = length(y)
	if n <= 2
		computeLinearPP(pp, x, y)
		return
	end

	dx = x[2:end] - x[1:end-1]
	S = (y[2:end] - y[1:end-1]) ./ dx
	b = pp.b
	fillDerivativeEstimate(limiter, dx, S, b)
	if leftBoundary == NOT_A_KNOT
		b[1] = S[2] * dx[1] / (dx[2] * (dx[2] + dx[1])) - S[1] * ((dx[2] / dx[1] + 2) / (dx[2] + dx[1]))
	elseif leftBoundary == FIRST_DIFFERENCE
		b[1] = S[1]
	elseif leftBoundary == FIRST_DERIVATIVE
		b[1] = leftValue
	elseif leftBoundary == SECOND_DERIVATIVE
		b[1] = (-leftValue / 2 * dx[1] - b[2] + 3 * S[1]) / 2
	end
	if rightBoundary == NOT_A_KNOT
		b[n] =
			S[n-2] * dx[n-1] / (dx[n-2] * (dx[n-2] + dx[n-1])) -
			S[n-1] * ((dx[n-2] / dx[n-1] + 2) / (dx[n-2] + dx[n-1]))
	elseif rightBoundary == FIRST_DERIVATIVE
		b[n] = rightValue
	elseif rightBoundary == FIRST_DIFFERENCE
		b[n] = S[n-1]
	elseif rightBoundary == SECOND_DERIVATIVE
		b[n] = (rightValue * dx[n-1] + 6 * S[n-1] - 2 * b[n-1]) / 4
	end
	c = pp.c
	for i ∈ 1:n-1
		c[i, 1] = (3 * S[i] - b[i+1] - 2 * b[i]) / dx[i]
		c[i, 2] = (b[i+1] + b[i] - 2 * S[i]) / (dx[i]^2)
	end
	pp.a[1:end] = y
	pp.x[1:end] = x
end

"evaluate a sorted array `za` with piecewise cubic interpolation, putting the results in `v``. `za`` must be sorted in ascending order. For now it does not support extrapolation."
function evaluateSorted!(self::Union{PP{3, T, TX}, PP{2, T, TX}}, v::AbstractVector{T}, za::AbstractVector{TZ}) where {T, TX, TZ}
	ppIndex = 1
	for j ∈ 1:length(za)
		z = za[j]
		while (ppIndex < length(self.x) && (self.x[ppIndex] < z)) #Si[ppIndex]<=z<Si[ppIndex+1]  
			ppIndex += 1
		end
		ppIndex -= 1
		ppIndex = min(max(ppIndex, 1), length(self.x) - 1)
		v[j] = evaluatePiece(self, ppIndex, z)
	end
end

function evaluate(self::Union{PP{3, T, TX}, PP{2, T, TX}, PP{1, T, TX}}, z::TZ; extrapolation = LinearAutoExtrapolation()) where {T, TX, TZ}
	if z < self.x[1]
		return evaluateLeftExtrapolation(self, z, extrapolation)
	elseif z >= self.x[end]
		return evaluateRightExtrapolation(self, z, extrapolation)
	end
	i = searchsortedfirst(self.x, z)  # x[i-1]<z<=x[i]
	if z != self.x[i] && i > 1
		i -= 1
	end
	return evaluatePiece(self, i, z)
end


function evaluateLeftExtrapolation(self::Union{PP{3, T, TX}, PP{2, T, TX}, PP{1, T, TX}}, z::TZ, extrapolation::LinearAutoExtrapolation) where {T, TX, TZ}
	return self.b[1] * (z - self.x[1]) + self.a[1]
end

function evaluateRightExtrapolation(self::Union{PP{3, T, TX}, PP{2, T, TX}, PP{1, T, TX}}, z::TZ, extrapolation::LinearAutoExtrapolation) where {T, TX, TZ}
	return self.b[end] * (z - self.x[end]) + self.a[end]
end

function evaluateLeftExtrapolation(self::Union{PP{3, T, TX}, PP{2, T, TX}, PP{1, T, TX}}, z::TZ, extrapolation::ConstantAutoExtrapolation) where {T, TX, TZ}
	return self.a[1]
end

function evaluateRightExtrapolation(self::Union{PP{3, T, TX}, PP{2, T, TX}, PP{1, T, TX}}, z::TZ, extrapolation::ConstantAutoExtrapolation) where {T, TX, TZ}
	return self.a[end]
end

function evaluateLeftExtrapolation(self::Union{PP{3, T, TX}, PP{2, T, TX}, PP{1, T, TX}}, z::TZ, extrapolation::PieceAutoExtrapolation) where {T, TX, TZ}
	return evaluatePiece(self, 1, z)
end

function evaluateRightExtrapolation(self::Union{PP{3, T, TX}, PP{2, T, TX}, PP{1, T, TX}}, z::TZ, extrapolation::PieceAutoExtrapolation) where {T, TX, TZ}
	return evaluatePiece(self, length(self.x) - 2, z)
end

@inline function evaluatePiece(self::PP{3, T, TX}, i::Int, z::TZ) where {T, TX, TZ}
	h = z - self.x[i]
	return self.a[i] + h * (self.b[i] + h * (self.c[i, 1] + h * (self.c[i, 2])))
end

@inline function evaluatePiece(self::PP{2, T, TX}, i::Int, z::TZ) where {T, TX, TZ}
	h = z - self.x[i]
	return self.a[i] + h * (self.b[i] + h * (self.c[i, 1]))
end

@inline function evaluatePiece(self::PP{1, T, TX}, i::Int, z::TZ) where {T, TX, TZ}
	h = z - self.x[i]
	return self.a[i] + h * self.b[i]
end


function evaluateDerivativeLeftExtrapolation(self::Union{PP{3, T, TX}, PP{2, T, TX}, PP{1, T, TX}}, z::TZ, extrapolation::LinearAutoExtrapolation) where {T, TX, TZ}
	return self.b[1]
end

function evaluateDerivativeRightExtrapolation(self::Union{PP{3, T, TX}, PP{2, T, TX}, PP{1, T, TX}}, z::TZ, extrapolation::LinearAutoExtrapolation) where {T, TX, TZ}
	return self.b[end]
end

function evaluateDerivativeLeftExtrapolation(self::Union{PP{3, T, TX}, PP{2, T, TX}, PP{1, T, TX}}, z::TZ, extrapolation::ConstantAutoExtrapolation) where {T, TX, TZ}
	return zero(T)
end

function evaluateDerivativeRightExtrapolation(self::Union{PP{3, T, TX}, PP{2, T, TX}, PP{1, T, TX}}, z::TZ, extrapolation::ConstantAutoExtrapolation) where {T, TX, TZ}
	return zero(T)
end

function evaluateDerivativeLeftExtrapolation(self::Union{PP{3, T, TX}, PP{2, T, TX}, PP{1, T, TX}}, z::TZ, extrapolation::PieceAutoExtrapolation) where {T, TX, TZ}
	return evaluateDerivativePiece(self, 1, z)
end

function evaluateDerivativeRightExtrapolation(self::Union{PP{3, T, TX}, PP{2, T, TX}, PP{1, T, TX}}, z::TZ, extrapolation::PieceAutoExtrapolation) where {T, TX, TZ}
	return evaluateDerivativePiece(self, length(self.x) - 2, z)
end

function evaluateDerivative(self::Union{PP{3, T, TX}, PP{2, T, TX}, PP{1, T, TX}}, z::TZ; extrapolation = LinearAutoExtrapolation()) where {T, TX, TZ}
	if z < self.x[1]
		return evaluateDerivativeLeftExtrapolation(self, z, extrapolation)
	elseif z >= self.x[end]
		return evaluateDerivativeRightExtrapolation(self, z, extrapolation)
	end
	i = searchsortedfirst(self.x, z)  # x[i-1]<z<=x[i]
	if z != self.x[i] && i > 1
		i -= 1
	end
	return evaluateDerivativePiece(self, i, z)
end

@inline function evaluateDerivativePiece(self::PP{3, T, TX}, i::Int, z::TZ) where {T, TX, TZ}
	h = z - self.x[i]
	return self.b[i] + h * (2 * self.c[i, 1] + h * (3 * self.c[i, 2]))
end

@inline function evaluateDerivativePiece(self::PP{2, T, TX}, i::Int, z::TZ) where {T, TX, TZ}
	h = z - self.x[i]
	return self.b[i] + 2h * self.c[i, 1]
end

@inline function evaluateDerivativePiece(self::PP{1, T, TX}, i::Int, z::TZ) where {T, TX, TZ}
	return self.b[i]
end

function evaluateIntegral(self::PP{N, T, TX}, z1::TZ, z2::TZ) where {N, T, TX, TZ}
	# get final type of calculation
	U = promote_type(T, promote_type(TX, TZ))
	# check if z1 < z2
	if z1 < z2
		sign = one(U)
	elseif z2 < z1
		sign = -one(U)
		z1, z2 = z2, z1
	else
		return zero(U)
	end

	# for now throw error if z1 or z2 are out of bounds
	if z1 < self.x[1] || z2 > self.x[end]
		throw(
			AssertionError(
				"The integration bounds must be inside of the interpolation interval.",
			),
		)
	end

	# find edge polynomials
	i1 = searchsortedfirst(self.x, z1)
	i2 = searchsortedfirst(self.x, z2)
	if z1 != self.x[i1]
		# the edge case i1 = 1 can not happen, as this would mean z1 < x[1] and in this case an error would already been thrown
		i1 -= 1
	end

	# calculate integral
	integral::U = zero(U)
	# iterate over the involved polynomials
	for i in i1:i2-1
		# calculate integration bounds of current integral
		# this makes sure that the lower integration bound of the first polynomial, i.e. the case i = i1, is z1        
		lower_bound = max(z1, self.x[i])
		# smiliar this makes sure taht the upper bound of the last polynoimial, i.e. the case i = i2 - 1, is z2
		upper_bound = min(z2, self.x[i+1])
		h1 = upper_bound - self.x[i]
		h2 = lower_bound - self.x[i]
		# f(x) = a + b (x-xi) + c1 (x-xi)^2 + c2 (x-xi)^3 + ...
		# indefinite integral with integraion constant -a xi + ...
		# F(x) = a (x-xi) + b/2 (x-xi)^2 + c1/3 (x-xi)^3 + c2/4 (x-xi)^4
		integral += self.a[i] * (h1 - h2) + self.b[i] / 2 * (h1^2 - h2^2) + sum(self.c[i, j] / (j + 2) * (h1^(j + 2) - h2^(j + 2)) for j in 1:N-1)
	end
	# return the integral multiplied by the correct sign, if the integration bounds were swapped
	return integral * sign
end


function evaluateSecondDerivativeLeftExtrapolation(self::Union{PP{3, T, TX}, PP{2, T, TX}, PP{1, T, TX}}, z::TZ, extrapolation::LinearAutoExtrapolation) where {T, TX, TZ}
	return zero(T)
end

function evaluateSecondDerivativeRightExtrapolation(self::Union{PP{3, T, TX}, PP{2, T, TX}, PP{1, T, TX}}, z::TZ, extrapolation::LinearAutoExtrapolation) where {T, TX, TZ}
	return zero(T)
end

function evaluateSecondDerivativeLeftExtrapolation(self::Union{PP{3, T, TX}, PP{2, T, TX}, PP{1, T, TX}}, z::TZ, extrapolation::ConstantAutoExtrapolation) where {T, TX, TZ}
	return zero(T)
end

function evaluateSecondDerivativeRightExtrapolation(self::Union{PP{3, T, TX}, PP{2, T, TX}, PP{1, T, TX}}, z::TZ, extrapolation::ConstantAutoExtrapolation) where {T, TX, TZ}
	return zero(T)
end

function evaluateSecondDerivativeLeftExtrapolation(self::Union{PP{3, T, TX}, PP{2, T, TX}, PP{1, T, TX}}, z::TZ, extrapolation::PieceAutoExtrapolation) where {T, TX, TZ}
	return evaluateSecondDerivativePiece(self, 1, z)
end

function evaluateSecondDerivativeRightExtrapolation(self::Union{PP{3, T, TX}, PP{2, T, TX}, PP{1, T, TX}}, z::TZ, extrapolation::PieceAutoExtrapolation) where {T, TX, TZ}
	return evaluateSecondDerivativePiece(self, length(self.x) - 2, z)
end
function evaluateSecondDerivative(self::PP{N, T, TX}, z::TZ; extrapolation = LinearAutoExtrapolation()) where {N, T, TX, TZ}
	if z < self.x[1]
		return evaluateSecondDerivativeLeftExtrapolation(self, z, extrapolation)
	elseif z >= self.x[end]
		return evaluateSecondDerivativeRightExtrapolation(self, z, extrapolation)
	end
	i = searchsortedfirst(self.x, z) # x[i-1]<z<=x[i]
	if z != self.x[i] && i > 1
		i -= 1
	end
	return evaluateSecondDerivativePiece(self, i, z)
end

@inline function evaluateSecondDerivativePiece(self::PP{3, T, TX}, i::Int, z::TZ) where {T, TX, TZ}
	h = z - self.x[i]
	return 2 * self.c[i, 1] + h * (3 * 2 * self.c[i, 2])
end

@inline function evaluateSecondDerivativePiece(self::PP{2, T, TX}, i::Int, z::TZ) where {T, TX, TZ}
	return 2 * self.c[i, 1]
end

@inline function evaluateSecondDerivativePiece(self::PP{1, T, TX}, i::Int, z::TZ) where {T, TX, TZ}
	return zero(T)
end

function evaluateThirdDerivative(self::PP{N, T, TX}, z::TZ) where {N, T, TX, TZ}
	if z <= self.x[1]
		return zero(T)
	elseif z >= self.x[end]
		return zero(T)
	end
	i = searchsortedfirst(self.x, z)  # x[i-1]<z<=x[i]
	if z != self.x[i] && i > 1
		i -= 1
	end
	return evaluateThirdDerivativePiece(self, i, z)
end

@inline function evaluateThirdDerivativePiece(self::PP{3, T, TX}, i::Int, z::TZ) where {T, TX, TZ}
	return 3 * 2 * self.c[i, 2]
end

(spl::PP{N, T, TX})(x::TZ) where {N, T, TX, TZ} = evaluate(spl, x)
function (spl::PP{N, T, TX})(x::AbstractArray) where {N, T, TX}
	evaluate.(spl, x)
end


function filterSlope(kind::C2, y::AbstractArray{T}, b::AbstractArray{T}, dx::AbstractArray{TX}, S::AbstractArray{T}) where {T, TX}
	#do nothing
end


function filterSlope(kind::C2MP, y::AbstractArray{T}, b::AbstractArray{T}, dx::AbstractArray{TX}, S::AbstractArray{T}) where {T, TX}
	n = length(y)
	if S[1] > 0
		b[1] = min(max(0, b[1]), 3 * S[1])
	else
		b[1] = max(min(0, b[1]), 3 * S[1])
	end
	if S[n-1] > 0
		b[n] = min(max(0, b[n]), 3 * S[n-1])
	else
		b[n] = max(min(0, b[n]), 3 * S[n-1])
	end

	for i ∈ 2:n-1
		Sim = S[i-1]
		Sip = S[i]
		if Sim * Sip <= 0
			b[i] = 0
		elseif Sim > 0 && Sip > 0
			b[i] = min(max(0, b[i]), 3 * min(Sim, Sip))
		else
			b[i] = max(min(0, b[i]), 3 * max(Sim, Sip))
		end
	end
end

minmod(s::T, t::T) where {T} = s * t <= 0 ? zero(T) : sign(s) * min(abs(s), abs(t))


function filterSlope(kind::C2MP2, y::AbstractArray{T}, b::AbstractArray{T}, dx::AbstractArray{TX}, S::AbstractArray{T}) where {T, TX}
	#r = s/t and G = t*g(r)
	n = length(y)
	b[1] = minmod(b[1], 3 * S[1])
	b[n] = minmod(b[n], 3 * S[n-1])

	for i ∈ 2:n-1
		Sim = S[i-1]
		Sip = S[i]
		lowerBound = minmod(Sim, Sip)
		upperBound = (sign(Sim) + sign(Sip)) / 2 * min(max(abs(Sim), abs(Sip)), 3 * min(abs(Sim), abs(Sip)))
		b[i] = min(max(lowerBound, b[i]), upperBound)
	end
end


function filterSlope(kind::C2HymanNonNegative, y::AbstractArray{T}, b::AbstractArray{T}, dx::AbstractArray{TX}, S::AbstractArray{T}) where {T, TX}
	n = length(y)
	tau0 = sign(y[1])
	#Warning the paper (3.3) is wrong as it does not obey (3.1)
	b[1] = tau0 * max(-3 * tau0 * y[1] / dx[1], tau0 * b[1])
	for i ∈ 2:n-1
		taui = sign(y[i])
		b[i] = taui * min(3 * taui * y[i] / dx[i-1], max(-3 * taui * y[i] / dx[i], taui * b[i]))
	end
	taun = sign(y[n])
	b[n] = taun * min(3 * taun * y[n] / dx[n-1], taun * b[n])
end

function filterSlope(kind::C2Hyman89, y::AbstractArray{T}, b::AbstractArray{T}, dx::AbstractArray{TX}, S::AbstractArray{T}) where {T, TX}
	n = length(y)
	tmp = b

	local correction, pm, pu, pd, M
	if (tmp[1] * S[1]) > 0
		correction = sign(tmp[1]) * min(abs(tmp[1]), abs(3 * S[1]))
	else
		correction = zero(T)
	end
	if correction != tmp[1]
		tmp[1] = correction
	end
	for i ∈ 2:n-1
		pm = ((S[i-1] * dx[i]) + (S[i] * dx[i-1])) / (dx[i-1] + dx[i])
		M = 3 * min(min(abs(S[i]), abs(S[i-1])), abs(pm))
		if i > 2
			if ((S[i-1] - S[i-2]) * (S[i] - S[i-1])) > 0
				pd = ((S[i-1] * ((2 * dx[i-1]) + dx[i-2])) - (S[i-2] * dx[i-1])) / (dx[i-2] + dx[i-1])
				if (pm * pd) > 0 && (pm * (S[i-1] - S[i-2])) > 0
					M = max(M, 3 * min(abs(pm), abs(pd)) / 2)
				end
			end
		end
		if i < (n - 1)
			if ((S[i] - S[i-1]) * (S[i+1] - S[i])) > 0
				pu = ((S[i] * ((2 * dx[i]) + dx[i+1])) - (S[i+1] * dx[i])) / (dx[i] + dx[i+1])
				if ((pm * pu) > 0) && ((-pm * (S[i] - S[i-1])) > 0)
					M = max(M, 3 * min(abs(pm), abs(pu)) / 2)
				end
			end
		end
		if (tmp[i] * pm) > 0
			correction = sign(tmp[i]) * min(abs(tmp[i]), M)
		else
			correction = zero(T)
		end
		if correction != tmp[i]
			tmp[i] = correction
		end
	end
	if (tmp[n] * S[n-1]) > 0
		correction = sign(tmp[n]) * min(abs(tmp[n]), abs(3 * S[n-1]))
	else
		correction = zero(T)
	end
	if correction != tmp[n]
		tmp[n] = correction
	end
end


function estimateDerivativeParabolic(x::AbstractArray{TX}, y::AbstractArray{T}) where {T, TX}
	n = length(x)
	b = zeros(n)
	dx = x[2:end] - x[1:end-1]
	S = (y[2:end] - y[1:end-1]) ./ dx
	for i ∈ 2:n-1
		b[i] = (dx[i-1] * S[i] + dx[i] * S[i-1]) / (dx[i-1] + dx[i])
	end
	return b
end


function limit(::HuynRational, s::T, t::T) where {T}
	st = s * t
	if st <= 0
		return zero(T)
	end
	return st * 3 * (s + t) / (s^2 + 4 * st + t^2)
end

function limit(::VanAlbada, s::T, t::T) where {T}
	st = s * t
	if st == 0
		return s
	end
	return st * (s + t) / (s^2 + t^2)
end

function limit(::VanLeer, s::T, t::T) where {T}
	st = s * t
	if st <= 0
		return zero(T)
	end
	return 2 * st / (s + t) #warning, the product s*t can be zero even when s or t are not 0, this rewrite helps saving accuracy
end

function limit(::FritschButland, s::T, t::T) where {T}
	st = s * t
	if st <= 0
		return zero(T)
	end
	if abs(s) <= abs(t)
		return 3 * st / (2 * s + t)
	end
	return 3 * st / (s + 2 * t)   # (1+dxp / (dxp+dx)) * t + (2 - dxp / (dxp+dx)) * s
end

function fillDerivativeEstimate(limiter::LimiterDerivative, dx::AbstractArray{TX}, S::AbstractArray{T}, b::AbstractArray{T}) where {T, TX}
	n = length(S)
	for i ∈ 2:n
		b[i] = limit(limiter, S[i-1], S[i])
	end
end

function fillDerivativeEstimate(limiter::Fukasawa, dx::AbstractArray{TX}, S::AbstractArray{T}, b::AbstractArray{T}) where {T, TX}
	# @. b[2:end-1] = (dx[1:end-1] * S[2:end] + dx[2:end] * S[1:end-1]) / (dx[1:end-1] + dx[2:end])
	n = length(S)
	for i ∈ 2:n
		li = sqrt(dx[i]^2 + (S[i] * dx[i])^2)
		lim = sqrt(dx[i-1]^2 + (S[i-1] * dx[i-1])^2)
		denom = (S[i] * dx[i] / li - S[i-1] * dx[i-1] / lim)
		b[i] = denom == zero(T) ? S[i] : -(dx[i] / li - dx[i-1] / lim) / denom
	end
end


function fillDerivativeEstimate(limiter::Bessel, dx::AbstractArray{TX}, S::AbstractArray{T}, b::AbstractArray{T}) where {T, TX}
	# @. b[2:end-1] = (dx[1:end-1] * S[2:end] + dx[2:end] * S[1:end-1]) / (dx[1:end-1] + dx[2:end])
	n = length(S)
	for i ∈ 2:n
		b[i] = (dx[i-1] * S[i] + dx[i] * S[i-1]) / (dx[i-1] + dx[i])
	end
end

function fillDerivativeEstimate(limiter::Hermite{T}, dx::AbstractArray{TX}, S::AbstractArray{T}, b::AbstractArray{T}) where {T, TX}
	# @. b[2:end-1] = (dx[1:end-1] * S[2:end] + dx[2:end] * S[1:end-1]) / (dx[1:end-1] + dx[2:end])
	n = length(S)
	for i ∈ 2:n
		b[i] = limiter.b[i]
	end
end

function fillDerivativeEstimate(limiter::Brodlie, dx::AbstractArray{TX}, S::AbstractArray{T}, b::AbstractArray{T}) where {T, TX}
	n = length(S)
	for i ∈ 2:n
		s, t = S[i-1], S[i]
		st = s * t
		if st == 0
			b[i] = s
		else
			α = (dx[i-1] + 2 * dx[i]) / (3 * (dx[i-1] + dx[i]))
			b[i] = (st) / (α * t + (1 - α) * s)
		end
	end
end


@inline function evaluate(self::PPInterpolation.PP{3, T, TX}, i::Int, z::TZ) where {T, TX, TZ}
	h = z - self.x[i]
	return self.a[i] + h * (self.b[i] + h * (self.c[i, 1] + h * (self.c[i, 2])))
end

normpdf(x) = exp(-x^2 / 2) / sqrt(2π)
normcdf(x) = erfc(-x / sqrt(2)) / 2

#Integral from z to Infty. i is index of z in pp representation.
#Assumes linear extrapolation for now.
function evaluateHermiteIntegralBounded(pp::PPInterpolation.PP{2, T, U}, i::Int, z::TZ)::T where {T, U, TZ}
	n = length(pp.x)
	if i > length(pp.x)
		i -= 1
	end
	integral = zero(z)
	l = z
	pdfL = normpdf(l)
	cdfL = normcdf(l)
	if z <= pp.x[1]
		#should be up to x[1]
		r = pp.x[1]
		pdfR = normpdf(r)
		cdfR = normcdf(r)
		slope = pp.b[1]
		moment = (pp.a[1] - pp.x[1] * slope) * (cdfR - cdfL) - slope * (pdfR - pdfL)
		integral += moment
		pdfL = pdfR
		cdfL = cdfR
		i = 0
	elseif z < pp.x[n]
		r = pp.x[i+1]
		pdfR = normpdf(r)
		cdfR = normcdf(r)
		moment = (pp.a[i] - pp.x[i] * (pp.b[i]) + pp.c[i, 1] * (pp.x[i]^2 + 1)) * (cdfR - cdfL)
		moment += pdfL * (pp.b[i] - pp.c[i, 1] * (2 * pp.x[i] - l)) - pdfR * (pp.b[i] - pp.c[i, 1] * (2 * pp.x[i] - r))
		integral += moment
		pdfL = pdfR
		cdfL = cdfR
	end
	for j ∈ i+1:n-1
		l = pp.x[j]
		r = pp.x[j+1]
		pdfR = normpdf(r)
		cdfR = normcdf(r)
		moment = (pp.a[j] - pp.x[j] * (pp.b[j]) + pp.c[j, 1] * (pp.x[j]^2 + 1)) * (cdfR - cdfL)
		moment += pdfL * (pp.b[j] - pp.c[j, 1] * (2 * pp.x[j] - l)) - pdfR * (pp.b[j] - pp.c[j, 1] * (2 * pp.x[j] - r))
		integral += moment
		pdfL = pdfR
		cdfL = cdfR
	end
	l = max(pp.x[n], z)
	slope = PPInterpolation.evaluateDerivativePiece(pp, n - 1, pp.x[n])
	yn = PPInterpolation.evaluatePiece(pp, n - 1, pp.x[n])
	moment = (yn - slope * pp.x[n]) * (1 - cdfL) + slope * pdfL
	integral += moment
	if isnan(integral) || isinf(integral)
		println(z, " ", ckIndex, " inf integral ", pp)
		throw(DomainError("infinite integral"))
	end
	return integral
end

function hermiteIntegralRightExtrapolation(pp::PPInterpolation.PP{N, T, U}, pdfL, cdfL, rightExtrapolation::ZeroExtrapolation) where {N, T, U}
	return zero(cdfL)
end
function hermiteIntegralLeftExtrapolation(pp::PPInterpolation.PP{N, T, U}, pdfL, cdfL, leftExtrapolation::ZeroExtrapolation) where {N, T, U}
	return zero(cdfL)
end

function hermiteIntegralLeftExtrapolation(pp::PPInterpolation.PP{N, T, U}, pdfL, cdfL, leftExtrapolation::LinearAutoExtrapolation) where {N, T, U}
	l = pp.x[1]
	slope = pp.b[1]
	(pp.a[1] - pp.x[1] * slope) * cdfL - slope * pdfL
end

function hermiteIntegralRightExtrapolation(pp::PPInterpolation.PP{N, T, U}, pdfL, cdfL, rightExtrapolation::LinearAutoExtrapolation) where {N, T, U}
	n = length(pp.x)
	l = pp.x[n]
	slope = evaluateDerivativePiece(pp, n - 1, pp.x[n])
	yn = evaluatePiece(pp, n - 1, pp.x[n])
	(yn - slope * pp.x[n]) * (1 - cdfL) + slope * pdfL
end

function hermiteIntegralLeftExtrapolation(pp::PPInterpolation.PP{N, T, U}, pdfL, cdfL, leftExtrapolation::ConstantAutoExtrapolation) where {N, T, U}
	pp.a[1] * cdfL
end

function hermiteIntegralRightExtrapolation(pp::PPInterpolation.PP{N, T, U}, pdfL, cdfL, rightExtrapolation::ConstantAutoExtrapolation) where {N, T, U}
	n = length(pp.x)
	l = pp.x[n]
	yn = evaluatePiece(pp, n - 1, pp.x[n])
	yn * (1 - cdfL)
end

"Integral of pp(x)*exp(-x^2 / 2) from -Infty to +Infty."
function evaluateHermiteIntegral(pp::PPInterpolation.PP{3, T, U}; leftExtrapolation = ConstantAutoExtrapolation(), rightExtrapolation = ConstantAutoExtrapolation()) where {T, U}
	n = length(pp.x)
	l = pp.x[1]
	pdfL = normpdf(l)
	cdfL = normcdf(l)
	integral = hermiteIntegralLeftExtrapolation(pp, pdfL, cdfL, leftExtrapolation)
	for i ∈ 1:n-1
		l = pp.x[i]
		r = pp.x[i+1]
		pdfR = normpdf(r)
		cdfR = normcdf(r)
		moment = (pp.a[i] - pp.x[i] * (pp.b[i] + pp.c[i, 2] * (pp.x[i]^2 + 3)) + pp.c[i, 1] * (pp.x[i]^2 + 1)) * (cdfR - cdfL)
		moment += pdfL * (pp.b[i] - pp.c[i, 1] * (2 * pp.x[i] - l) + pp.c[i, 2] * (3 * pp.x[i]^2 - 3 * l * pp.x[i] + l^2 + 2)) - pdfR * (pp.b[i] - pp.c[i, 1] * (2 * pp.x[i] - r) + pp.c[i, 2] * (3 * pp.x[i]^2 - 3 * pp.x[i] * r + r^2 + 2))
		integral += moment
		pdfL = pdfR
		cdfL = cdfR
	end
	integral += hermiteIntegralRightExtrapolation(pp, pdfL, cdfL, rightExtrapolation)
	if isnan(integral) || isinf(integral)
		println(z, " ", ckIndex, " inf integral ", pp)
		throw(DomainError("infinite integral"))
	end
	return integral
end
