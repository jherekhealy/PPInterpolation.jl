struct QuadraticPP{T<:Number,U<:Real}
    a::Vector{T}
    b::Vector{T}
    c::Vector{T}
    x::Vector{U}
    QuadraticPP(T, U, n::Int) = new{T,U}(zeros(T, n), zeros(T, n), zeros(T, n - 1), zeros(U, n))
    QuadraticPP(a::Vector{T}, b::Vector{T}, c::Vector{T}, x::Vector{U}) where {T<:Real,U<:Real} = new{T,U}(a, b, c, x)
end

(p::QuadraticPP)(x::T) where {T} = evaluate(p, x)
Base.length(pp::QuadraticPP) = length(pp.a)

function evaluate(self::QuadraticPP, z::T) where {T}
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
    return self.a[i] + h * (self.b[i] + h * (self.c[i]))
end


function evaluateDerivative(self::QuadraticPP, z::T) where {T}
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
    return self.b[i] + 2h * self.c[i]
end