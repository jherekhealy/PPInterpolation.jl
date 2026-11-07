export PP, evaluate, evaluateDerivative, evaluateSecondDerivative, evaluateThirdDerivative, evaluateSorted!, evaluateIntegral
export evaluatePiece, evaluateDerivativePiece, evaluateSecondDerivativePiece, evaluateThirdDerivativePiece
export LinearAutoExtrapolation, ConstantAutoExtrapolation, PieceAutoExtrapolation, ZeroExtrapolation

abstract type PPExtrapolation end
"Extrapolation is zero beyond first and last knots"
struct ZeroExtrapolation <: PPExtrapolation end
"Extrapolation is constant beyond first and last knots, using the value at the endpoints"
struct ConstantAutoExtrapolation <: PPExtrapolation end
"Extrapolation is linear beyond first and last knots, using the slope and value at the endpoints"
struct LinearAutoExtrapolation <: PPExtrapolation end
"Extrapolation is using the first/last piecewise polynomial"
struct PieceAutoExtrapolation <: PPExtrapolation end

"""
	PP{N,T,TX}

Piecewise polynomial of degree `N` on the knots `x`. On `[x[i], x[i+1]]`, it equals
`a[i] + b[i] h + c[i,1] h^2 + ... + c[i,N-1] h^N` with `h = z - x[i]`.

`a` and `b` have one entry per knot: the value and first derivative at the knot (at the last knot, those of the last piece).
`c` has one row per piece. Beyond the first and last knots, the polynomial is evaluated according to `extrapolation`.
"""
struct PP{N, T, TX, VT <: AbstractArray{T}, MT <: AbstractMatrix{T}, VX <: AbstractArray{TX}, E <: PPExtrapolation}
	a::VT
	b::VT
	c::MT #c[i,j] = coeff of (x-x[i])^{j+1}
	x::VX
	extrapolation::E
	function PP{N, T, TX, VT, MT, VX, E}(a, b, c, x, extrapolation) where {N, T, TX, VT, MT, VX, E}
		n = length(x)
		if length(a) != n || length(b) != n
			throw(DimensionMismatch("a and b must have one entry per knot, got $(length(a)) and $(length(b)) for $n knots"))
		end
		if size(c) != (max(n - 1, 0), N - 1)
			throw(DimensionMismatch("c must be of size $((max(n - 1, 0), N - 1)), got $(size(c))"))
		end
		new{N, T, TX, VT, MT, VX, E}(a, b, c, x, extrapolation)
	end
end

PP(N::Int, T, TX, n::Int; extrapolation::PPExtrapolation = LinearAutoExtrapolation()) =
	PP{N, T, TX, Vector{T}, Matrix{T}, Vector{TX}, typeof(extrapolation)}(zeros(T, n), zeros(T, n), zeros(T, (max(n - 1, 0), N - 1)), zeros(TX, n), extrapolation)
PP(N::Int, a::AbstractArray{T}, b::AbstractArray{T}, c::AbstractMatrix{T}, x::AbstractArray{TX}; extrapolation::PPExtrapolation = LinearAutoExtrapolation()) where {T <: Real, TX} =
	PP{N, T, TX, typeof(a), typeof(c), typeof(x), typeof(extrapolation)}(a, b, c, x, extrapolation)
"A copy of the piecewise polynomial `pp` sharing its coefficients, with a different `extrapolation`."
PP(pp::PP{N}; extrapolation::PPExtrapolation) where {N} = PP(N, pp.a, pp.b, pp.c, pp.x; extrapolation = extrapolation)

Base.length(p::PP) = Base.length(p.x)
Base.size(p::PP) = Base.size(p.x)
Base.broadcastable(p::PP) = Ref(p)

(spl::PP)(z) = evaluate(spl, z)
(spl::PP)(z::AbstractArray) = evaluate.(spl, z)

"index `i` of the piece such that `x[i] <= z < x[i+1]`, clamped to the first and last pieces."
@inline function findPiece(self::PP, z)
	x = self.x
	i = searchsortedfirst(x, z)  # x[i-1] < z <= x[i]
	if i > length(x) || z != x[i]
		i -= 1
	end
	return clamp(i, 1, max(length(x) - 1, 1))
end

#Horner evaluation through evalpoly, unrolled since the number of coefficients is known at compile time
@inline function evaluatePiece(self::PP{N}, i::Int, z) where {N}
	h = z - self.x[i]
	N == 1 && return self.a[i] + h * self.b[i]
	c = self.c
	return self.a[i] + h * (self.b[i] + h * evalpoly(h, ntuple(j -> c[i, j], Val(N - 1))))
end

@inline function evaluateDerivativePiece(self::PP{N}, i::Int, z) where {N}
	N == 1 && return self.b[i]
	h = z - self.x[i]
	c = self.c
	return self.b[i] + h * evalpoly(h, ntuple(j -> (j + 1) * c[i, j], Val(N - 1)))
end

@inline function evaluateSecondDerivativePiece(self::PP{N}, i::Int, z) where {N}
	N == 1 && return zero(self.a[i])
	h = z - self.x[i]
	c = self.c
	return evalpoly(h, ntuple(j -> (j + 1) * j * c[i, j], Val(N - 1)))
end

@inline function evaluateThirdDerivativePiece(self::PP{N}, i::Int, z) where {N}
	N <= 2 && return zero(self.a[i])
	h = z - self.x[i]
	c = self.c
	return evalpoly(h, ntuple(k -> (k + 2) * (k + 1) * k * c[i, k+1], Val(N - 2)))
end

#explicit forms of the most common quadratic and cubic cases, faster than the generic loops above
@inline function evaluatePiece(self::PP{3}, i::Int, z)
	h = z - self.x[i]
	return self.a[i] + h * (self.b[i] + h * (self.c[i, 1] + h * self.c[i, 2]))
end
@inline function evaluateDerivativePiece(self::PP{3}, i::Int, z)
	h = z - self.x[i]
	return self.b[i] + h * (2 * self.c[i, 1] + h * (3 * self.c[i, 2]))
end
@inline evaluateSecondDerivativePiece(self::PP{3}, i::Int, z) = 2 * self.c[i, 1] + (z - self.x[i]) * (6 * self.c[i, 2])
@inline function evaluatePiece(self::PP{2}, i::Int, z)
	h = z - self.x[i]
	return self.a[i] + h * (self.b[i] + h * self.c[i, 1])
end
@inline evaluateDerivativePiece(self::PP{2}, i::Int, z) = self.b[i] + 2 * (z - self.x[i]) * self.c[i, 1]

#derivative of order 0 to 3 of the piece i
@inline derivativePiece(self::PP, i::Int, z, ::Val{0}) = evaluatePiece(self, i, z)
@inline derivativePiece(self::PP, i::Int, z, ::Val{1}) = evaluateDerivativePiece(self, i, z)
@inline derivativePiece(self::PP, i::Int, z, ::Val{2}) = evaluateSecondDerivativePiece(self, i, z)
@inline derivativePiece(self::PP, i::Int, z, ::Val{3}) = evaluateThirdDerivativePiece(self, i, z)

#derivative of order 0 to 3 of the extrapolation beyond the end knot k, the adjacent piece being p
extrapolate(self::PP, z, ::ZeroExtrapolation, k::Int, p::Int, ::Val) = zero(self.a[k])
extrapolate(self::PP, z, ::ConstantAutoExtrapolation, k::Int, p::Int, ::Val{0}) = self.a[k]
extrapolate(self::PP, z, ::ConstantAutoExtrapolation, k::Int, p::Int, ::Val) = zero(self.a[k])
extrapolate(self::PP, z, ::LinearAutoExtrapolation, k::Int, p::Int, ::Val{0}) = self.a[k] + self.b[k] * (z - self.x[k])
extrapolate(self::PP, z, ::LinearAutoExtrapolation, k::Int, p::Int, ::Val{1}) = self.b[k]
extrapolate(self::PP, z, ::LinearAutoExtrapolation, k::Int, p::Int, ::Val) = zero(self.a[k])
extrapolate(self::PP, z, ::PieceAutoExtrapolation, k::Int, p::Int, order::Val) = derivativePiece(self, p, z, order)

#value and first derivative at the last knot are the stored ones, higher derivatives are the ones of the last piece.
#the (zero) term in z - x[end] keeps the derivative for automatic differentiation.
atLastKnot(self::PP, z, ::Val{0}) = self.a[end] + self.b[end] * (z - self.x[end])
atLastKnot(self::PP, z, ::Val{1}) = self.b[end] + evaluateSecondDerivativePiece(self, length(self.x) - 1, z) * (z - self.x[end])
atLastKnot(self::PP, z, order::Val) = derivativePiece(self, length(self.x) - 1, z, order)

@inline function evaluateOrder(self::PP, z, extrapolation::PPExtrapolation, order::Val)
	x = self.x
	n = length(x)
	if z < x[1]
		return extrapolate(self, z, extrapolation, 1, 1, order)
	elseif z > x[n]
		return extrapolate(self, z, extrapolation, n, max(n - 1, 1), order)
	elseif z == x[n]
		return atLastKnot(self, z, order)
	end
	return derivativePiece(self, findPiece(self, z), z, order)
end

"Value of the piecewise polynomial at `z`. The `extrapolation` defaults to the one of the polynomial."
evaluate(self::PP, z; extrapolation::PPExtrapolation = self.extrapolation) = evaluateOrder(self, z, extrapolation, Val(0))
"First derivative of the piecewise polynomial at `z`."
evaluateDerivative(self::PP, z; extrapolation::PPExtrapolation = self.extrapolation) = evaluateOrder(self, z, extrapolation, Val(1))
"Second derivative of the piecewise polynomial at `z`."
evaluateSecondDerivative(self::PP, z; extrapolation::PPExtrapolation = self.extrapolation) = evaluateOrder(self, z, extrapolation, Val(2))
"Third derivative of the piecewise polynomial at `z`."
evaluateThirdDerivative(self::PP, z; extrapolation::PPExtrapolation = self.extrapolation) = evaluateOrder(self, z, extrapolation, Val(3))

"evaluate a sorted array `za` with the piecewise polynomial, putting the results in `v`. `za` must be sorted in ascending order."
function evaluateSorted!(self::PP, v::AbstractVector, za::AbstractVector; extrapolation::PPExtrapolation = self.extrapolation)
	x = self.x
	x1, xn = x[1], x[end]
	i = 1
	for j ∈ eachindex(za)
		z = za[j]
		if x1 <= z < xn
			#stops at the last piece at the latest, since z < x[end]
			while x[i+1] <= z
				i += 1
			end
			v[j] = evaluatePiece(self, i, z)
		else
			v[j] = evaluate(self, z, extrapolation = extrapolation)
		end
	end
	return v
end

"Integral of the piecewise polynomial between `z1` and `z2`, which must be within the knots."
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
	i1 = findPiece(self, z1)
	i2 = searchsortedfirst(self.x, z2)

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
		integral += self.a[i] * (h1 - h2) + self.b[i] / 2 * (h1^2 - h2^2) + sum((self.c[i, j] / (j + 2) * (h1^(j + 2) - h2^(j + 2)) for j in 1:N-1); init = zero(U))
	end
	# return the integral multiplied by the correct sign, if the integration bounds were swapped
	return integral * sign
end
