#Quadratic Lagrange Interpolation. Warning, x and y are not copied. Assumes that x is sorted in ascending order.
export QuadraticLagrangePP, evaluate, evaluateSortedAtLeft!, evaluateSortedAtMid!, findIndex, evaluateSorted, evaluateSorted!, evaluatePiece, evaluateDerivativePiece, evaluateSecondDerivativePiece

"""
Controls which piece is used for a given interpolation point `z`. The index `i` of the piece is used to interpolate using x[i-1],x[i],x[i+1]
LEFT_KNOT, we chose the index such that  x[i] <= z < x[i+1]
MID_KNOT, we chose the index such that  (x[i]+x[i-1])/2 <= z < (x[i]+x[i+1])/2
When z <= x[2], the index is 2, when z >= x[end-1], the index is length(x)-1
"""
@enum LagrangeKnot LEFT_KNOT = 0 MID_KNOT = 1 

"""A piecewise quadratic lagrange interpolation based on the abscissae `x` and ordinates `y`. """
struct QuadraticLagrangePP{TX, T}
    x::Vector{TX}
    y::Vector{T}
    knotStyle::LagrangeKnot
    QuadraticLagrangePP(x::AbstractArray{TX},y::AbstractArray{T};knotStyle=MID_KNOT) where {TX,T}= new{TX,T}(x,y,knotStyle)
end

"""
    convert(PP{2,T,TX}, q::QuadraticLagrangePP)

Exact piecewise quadratic representation of the quadratic Lagrange interpolation `q`. Its breakpoints are the boundaries of the Lagrange pieces,
which depend on the `knotStyle`. The interpolation is not continuous at those breakpoints.
"""
function Base.convert(::Type{PP{2, T, TX}}, q::QuadraticLagrangePP) where {T, TX}
    x, y = q.x, q.y
    n = length(x)
    n >= 3 || throw(ArgumentError("quadratic Lagrange interpolation requires at least 3 knots"))
    # breakpoints and Lagrange index of each piece, consistent with findIndex
    if q.knotStyle == LEFT_KNOT
        breaks = TX[x[1]; x[3:end]]
        lagrangeIndex = j -> j == 1 ? 2 : j + 1
    else
        breaks = TX[x[1]; [(x[i] + x[i+1]) / 2 for i in 2:n-2]; x[end]]
        lagrangeIndex = j -> j + 1
    end
    m = length(breaks)
    a = zeros(T, m)
    b = zeros(T, m)
    c = zeros(T, (m - 1, 1))
    for j in 1:m-1
        k = lagrangeIndex(j)
        a[j] = evaluate(k, x, y, breaks[j])
        b[j] = evaluateDerivative(k, x, y, breaks[j])
        c[j, 1] = evaluateSecondDerivative(k, x, y, breaks[j]) / 2
    end
    k = lagrangeIndex(m - 1)
    a[m] = evaluate(k, x, y, breaks[m])
    b[m] = evaluateDerivative(k, x, y, breaks[m])
    return PP(2, a, b, c, breaks)
end

Base.convert(::Type{PP}, q::QuadraticLagrangePP{TX, T}) where {TX, T} = convert(PP{2, T, TX}, q)

Base.length(p::QuadraticLagrangePP) = Base.length(p.x)
Base.size(p::QuadraticLagrangePP) = Base.size(p.x)
Base.broadcastable(p::QuadraticLagrangePP) = Ref(p)

(spl::QuadraticLagrangePP{TX,T})(x::TZ) where {TX,T,TZ} = evaluate(spl, x)
function (spl::QuadraticLagrangePP{TX,T})(x::AbstractArray) where {TX,T}
    evaluate.(spl, x)
end

"evaluate a sorted array `za` with quadratic lagrange interpolation, putting the results in `v``. `za`` must be sorted in ascending order"
function evaluateSorted!(pp::QuadraticLagrangePP{TX,T}, v::AbstractArray{T}, za::AbstractArray{TZ}) where {TX,T,TZ}  
    if pp.knotStyle == LEFT_KNOT
        evaluateSortedAtLeft!(v,pp.x, pp.y, za)
    else
        evaluateSortedAtMid!(v,pp.x, pp.y, za)
    end
end

"evaluate a sorted array `za` with quadratic lagrange interpolation, returning the result of the interpolation as a new array. `za` must be sorted in ascending order"
function evaluateSorted(pp::QuadraticLagrangePP{TX,T}, z::AbstractArray{TZ}) where {TX,T,TZ}
    v = Array{T}(undef, length(z))
    evaluateSorted!(pp, v, z)
    return v
end

"evaluate the quadratic lagrange interpolation at a given point `z`"
function evaluate(pp::QuadraticLagrangePP{TX,T}, z::TZ) where {TX,T,TZ}
    ppIndex = findIndex(pp, z)
    evaluatePiece(pp,ppIndex,z)
end


"evaluate the quadratic lagrange interpolation first derivative at a given point `z`"
function evaluateDerivative(pp::QuadraticLagrangePP{TX,T}, z::TZ) where {TX,T,TZ}
    ppIndex = findIndex(pp, z)
    evaluateDerivativePiece(pp,ppIndex,z)
end

function evaluateSecondDerivative(pp::QuadraticLagrangePP{TX,T}, z::TZ) where {TX,T,TZ}
    ppIndex = findIndex(pp, z)
    evaluateSecondDerivativePiece(pp,ppIndex,z)
end

"find the index of the piece of the quadratic lagrange interpolation for the given point `z`. The `knotStyle` will control which piece is returned."
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

"evaluate the quadratic lagrange interpolation piece corresponding to index `ppIndex` at a given point `z`"
evaluatePiece(pp::QuadraticLagrangePP{TX,T}, ppIndex::Int, z::TZ) where {TX,T,TZ} = evaluate(ppIndex,pp.x,pp.y,z)
evaluateDerivativePiece(pp::QuadraticLagrangePP{TX,T}, ppIndex::Int, z::TZ) where {TX,T,TZ} = evaluateDerivative(ppIndex,pp.x,pp.y,z)
evaluateSecondDerivativePiece(pp::QuadraticLagrangePP{TX,T}, ppIndex::Int, z::TZ) where {TX,T,TZ} = evaluateSecondDerivative(ppIndex,pp.x,pp.y,z)

@inline function evaluate(ppIndex::Int, x::AbstractArray{TX},y::AbstractArray{T},z::TZ) where {TX,T,TZ} 
   return y[ppIndex] * (x[ppIndex-1] - z) * (x[ppIndex+1] - z) / ((x[ppIndex-1] - x[ppIndex]) * (x[ppIndex+1] - x[ppIndex])) + y[ppIndex-1] * (x[ppIndex] - z) * (x[ppIndex+1] - z) / ((x[ppIndex] - x[ppIndex-1]) * (x[ppIndex+1] - x[ppIndex-1])) + y[ppIndex+1] * (x[ppIndex-1] - z) * (x[ppIndex] - z) / ((x[ppIndex-1] - x[ppIndex+1]) * (x[ppIndex] - x[ppIndex+1]))
end

@inline function evaluateDerivative(ppIndex::Int, x::AbstractArray{TX},y::AbstractArray{T},z::TZ) where {TX,T,TZ} 
    return -y[ppIndex] * ((x[ppIndex-1] - z) + (x[ppIndex+1] - z)) / ((x[ppIndex-1] - x[ppIndex]) * (x[ppIndex+1] - x[ppIndex])) - y[ppIndex-1] * ((x[ppIndex] - z) + (x[ppIndex+1] - z)) / ((x[ppIndex] - x[ppIndex-1]) * (x[ppIndex+1] - x[ppIndex-1])) - y[ppIndex+1] * ((x[ppIndex-1] - z) + (x[ppIndex] - z)) / ((x[ppIndex-1] - x[ppIndex+1]) * (x[ppIndex] - x[ppIndex+1]))
 end

@inline function evaluateSecondDerivative(ppIndex::Int, x::AbstractArray{TX},y::AbstractArray{T},z::TZ) where {TX,T,TZ} 
    return 2y[ppIndex] / ((x[ppIndex-1] - x[ppIndex]) * (x[ppIndex+1] - x[ppIndex])) +2y[ppIndex-1]  / ((x[ppIndex] - x[ppIndex-1]) * (x[ppIndex+1] - x[ppIndex-1])) + 2y[ppIndex+1] / ((x[ppIndex-1] - x[ppIndex+1]) * (x[ppIndex] - x[ppIndex+1]))
end
 
function evaluateSortedAtMid!(v::AbstractArray{T}, x::AbstractArray{TX},y::AbstractArray{T}, za::AbstractArray{TZ}) where {TX,T,TZ}
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

function evaluateSortedAtLeft!(v::AbstractArray{T}, x::AbstractArray{TX},y::AbstractArray{T}, za::AbstractArray{TZ}) where {TX,T,TZ}
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
