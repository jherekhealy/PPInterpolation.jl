#Quadratic Lagrange Interpolation. Warning, x and y are not copied. Assumes that x is sorted in ascending order.
export QuadraticLagrangePP, evaluate, evaluate!, evaluateAtLeft!, evaluateAtMid!, findIndex

@enum LagrangeKnot LEFT_KNOT = 0 MID_KNOT = 1 

struct QuadraticLagrangePP{TX, T}
    x::Vector{TX}
    y::Vector{T}
    knotStyle::LagrangeKnot
    QuadraticLagrangePP(x::AbstractArray{TX},y::AbstractArray{T};knotStyle=MID_KNOT) where {TX,T}= new{TX,T}(x,y,knotStyle)
end

#TODO add method to convert to PP{2,T,TX}

Base.length(p::QuadraticLagrangePP) = Base.length(p.x)
Base.size(p::QuadraticLagrangePP) = Base.size(p.x)
Base.broadcastable(p::QuadraticLagrangePP) = Ref(p)

(spl::QuadraticLagrangePP{TX,T})(x::TZ) where {TX,T,TZ} = evaluate(spl, x)
function (spl::QuadraticLagrangePP{TX,T})(x::AbstractArray) where {TX,T}
    evaluate.(spl, x)
end

function evaluate!(pp::QuadraticLagrangePP{TX,T}, v::AbstractArray{T}, za::AbstractArray{TZ}) where {TX,T,TZ}  
    if pp.knotStyle == LEFT_KNOT
        evaluateAtLeft!(v,pp.x, pp.y, za)
    else
        evaluateAtMid!(v,pp.x, pp.y, za)
    end
end

function evaluate(pp::QuadraticLagrangePP{TX,T}, z::AbstractArray{TZ}) where {TX,T,TZ}
    v = Array{T}(undef, length(z))
    evaluate!(v, pp.x, pp.y, za)   
end

function evaluate(pp::QuadraticLagrangePP{TX,T}, z::TZ) where {TX,T,TZ}
    ppIndex = findIndex(pp, z)
    evaluate(pp,ppIndex,z)
end



function evaluateDerivative(pp::QuadraticLagrangePP{TX,T}, z::TZ) where {TX,T,TZ}
    ppIndex = findIndex(pp, z)
    evaluateDerivative(pp,ppIndex,z)
end

function evaluateSecondDerivative(pp::QuadraticLagrangePP{TX,T}, z::TZ) where {TX,T,TZ}
    ppIndex = findIndex(pp, z)
    evaluateSecondDerivative(pp,ppIndex,z)
end

@inline function findIndex(pp::QuadraticLagrangePP{TX,T}, z::TZ) where {TX,T,TZ}
    if pp.knotStyle == LEFT_KNOT
        ppIndex = searchsortedlast(pp.x,z) #   x[i]<=z<x[i+1]
    else
        # we want  x[i]+x[i-1]<=2z<x[i]+x[i+1]  .i-1   .i    .i+1
        ppIndex = searchsortedlast(pp.x,z) 
        if ppIndex > 0 && ppIndex < length(pp.x) 
            if (pp.x[ppIndex] + pp.x[ppIndex+1] < 2z)
                ppIndex += 1
            end
        end
    end
    ppIndex = min(max(ppIndex, 2), length(pp.x) - 1)
    ppIndex
end

evaluate(pp::QuadraticLagrangePP{TX,T}, ppIndex::Int, z::TZ) where {TX,T,TZ} = evaluate(ppIndex,pp.x,pp.y,z)
evaluateDerivative(pp::QuadraticLagrangePP{TX,T}, ppIndex::Int, z::TZ) where {TX,T,TZ} = evaluateDerivative(ppIndex,pp.x,pp.y,z)
evaluateSecondDerivative(pp::QuadraticLagrangePP{TX,T}, ppIndex::Int, z::TZ) where {TX,T,TZ} = evaluateSecondDerivative(ppIndex,pp.x,pp.y,z)

@inline function evaluate(ppIndex::Int, x::AbstractArray{TX},y::AbstractArray{T},z::TZ) where {TX,T,TZ} 
   return y[ppIndex] * (x[ppIndex-1] - z) * (x[ppIndex+1] - z) / ((x[ppIndex-1] - x[ppIndex]) * (x[ppIndex+1] - x[ppIndex])) + y[ppIndex-1] * (x[ppIndex] - z) * (x[ppIndex+1] - z) / ((x[ppIndex] - x[ppIndex-1]) * (x[ppIndex+1] - x[ppIndex-1])) + y[ppIndex+1] * (x[ppIndex-1] - z) * (x[ppIndex] - z) / ((x[ppIndex-1] - x[ppIndex+1]) * (x[ppIndex] - x[ppIndex+1]))
end

@inline function evaluateDerivative(ppIndex::Int, x::AbstractArray{TX},y::AbstractArray{T},z::TZ) where {TX,T,TZ} 
    return -y[ppIndex] * ((x[ppIndex-1] - z) + (x[ppIndex+1] - z)) / ((x[ppIndex-1] - x[ppIndex]) * (x[ppIndex+1] - x[ppIndex])) - y[ppIndex-1] * ((x[ppIndex] - z) + (x[ppIndex+1] - z)) / ((x[ppIndex] - x[ppIndex-1]) * (x[ppIndex+1] - x[ppIndex-1])) - y[ppIndex+1] * ((x[ppIndex-1] - z) + (x[ppIndex] - z)) / ((x[ppIndex-1] - x[ppIndex+1]) * (x[ppIndex] - x[ppIndex+1]))
 end

 @inline function evaluateSecondDerivative(ppIndex::Int, x::AbstractArray{TX},y::AbstractArray{T},z::TZ) where {TX,T,TZ} 
    return 2y[ppIndex] / ((x[ppIndex-1] - x[ppIndex]) * (x[ppIndex+1] - x[ppIndex])) +2y[ppIndex-1]  / ((x[ppIndex] - x[ppIndex-1]) * (x[ppIndex+1] - x[ppIndex-1])) + 2y[ppIndex+1] / ((x[ppIndex-1] - x[ppIndex+1]) * (x[ppIndex] - x[ppIndex+1]))
 end
 
function evaluateAtMid!(v::AbstractArray{T}, x::AbstractArray{TX},y::AbstractArray{T}, za::AbstractArray{TZ}) where {TX,T,TZ}
    ppIndex = 2
    for j = 1:length(za)
        z = za[j]
        while (ppIndex < length(x) && (x[ppIndex] + x[ppIndex-1] < 2z)) #Si[ppIndex]<=z<Si[ppIndex+1]  
            ppIndex += 1
        end
        ppIndex -= 1
        ppIndex = min(max(ppIndex, 2), length(x) - 1)
        v[j] = evaluate(ppIndex,x,y,z)

    end
end

function evaluateAtLeft!(v::AbstractArray{T}, x::AbstractArray{TX},y::AbstractArray{T}, za::AbstractArray{TZ}) where {TX,T,TZ}
    ppIndex = 1
    for j = 1:length(za)
        z = za[j]
        while (ppIndex < length(x) && (x[ppIndex] < z)) #Si[ppIndex]<=z<Si[ppIndex+1]  
            ppIndex += 1
        end
        ppIndex -= 1
        ppIndex = min(max(ppIndex, 2), length(x) - 1)
        v[j] = evaluate(ppIndex,x,y,z)
    end
end
