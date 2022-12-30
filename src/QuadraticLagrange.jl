#Quadratic Lagrange Interpolation. Warning, x and y are not copied. Assumes that x is sorted in ascending order.
export QuadraticLagrangePP

struct QuadraticLagrangePP{TX, T}
    x::Vector{TX}
    y::Vector{T}
end

Base.length(p::QuadraticLagrangePP) = Base.length(p.x)
Base.size(p::QuadraticLagrangePP) = Base.size(p.x)
Base.broadcastable(p::QuadraticLagrangePP) = Ref(p)


function evaluateMid!(v::AbstractArray{T}, pp::QuadraticLagrangePP{TX,T}, za::AbstractArray{T}) where {TX,T}
    evaluateMid!(v, pp.x, pp.y, za)
end

function evaluateMid(pp::QuadraticLagrangePP{TX,T}, z::AbstractArray{T}) where {TX,T}
    v = Array{T}(undef, length(z))
    evaluateMid!(v, pp, pp.x, pp.y, z)
end

function evaluate(pp::QuadraticLagrangePP{TX,T}, z::T) where {TX,T}
    ppIndex = findIndex(pp, z)
    evaluate(pp,ppIndex,z)
end

@inline function findIndex(pp::QuadraticLagrangePP{TX,T}, z::T) where {TX,T}
    ppIndex = searchsortedlast(pp.x,z)
    ppIndex = min(max(ppIndex, 2), length(pp.x) - 1)
    ppIndex
end

@inline function evaluate(pp::QuadraticLagrangePP{TX,T}, ppIndex::Int, z::T) where {TX,T}
    return pp.y[ppIndex] * (pp.x[ppIndex-1] - z) * (pp.x[ppIndex+1] - z) / ((pp.x[ppIndex-1] - pp.x[ppIndex]) * (pp.x[ppIndex+1] - pp.x[ppIndex])) + pp.y[ppIndex-1] * (pp.x[ppIndex] - z) * (pp.x[ppIndex+1] - z) / ((pp.x[ppIndex] - pp.x[ppIndex-1]) * (pp.x[ppIndex+1] - pp.x[ppIndex-1])) + pp.y[ppIndex+1] * (pp.x[ppIndex-1] - z) * (pp.x[ppIndex] - z) / ((pp.x[ppIndex-1] - pp.x[ppIndex+1]) * (pp.x[ppIndex] - pp.x[ppIndex+1]))
end


function evaluateMid!(v::AbstractArray{T}, x::AbstractArray{T},y::AbstractArray{T}, za::AbstractArray{T}) where {T}
    ppIndex = 2
    for j = 1:length(za)
        z = za[j]
        while (ppIndex < length(x) && (x[ppIndex] + x[ppIndex-1] < 2z)) #Si[ppIndex]<=z<Si[ppIndex+1]  
            ppIndex += 1
        end
        ppIndex -= 1
        ppIndex = min(max(ppIndex, 2), length(x) - 1)
        v[j] =  y[ppIndex] * (x[ppIndex-1] - z) * (x[ppIndex+1] - z) / ((x[ppIndex-1] - x[ppIndex]) * (x[ppIndex+1] - x[ppIndex])) + y[ppIndex-1] * (x[ppIndex] - z) * (x[ppIndex+1] - z) / ((x[ppIndex] - x[ppIndex-1]) * (x[ppIndex+1] - x[ppIndex-1])) + y[ppIndex+1] * (x[ppIndex-1] - z) * (x[ppIndex] - z) / ((x[ppIndex-1] - x[ppIndex+1]) * (x[ppIndex] - x[ppIndex+1]))

    end
end
