using BSplines


function Base.convert(::Type{PP{2,T,TX}}, spl::BSplines.Spline) where {T,TX}
    t = BSplines.knots(spl.basis)
    n = length(spl.basis) - 1
    a = zeros(T, n)
    b = zeros(T, n)
    c = zeros(T, (n - 1, 1))
    x = Vector(spl.basis.breakpoints)
    α = spl.coeffs
    for i = 1:n
        a[i] = (t[i+2] - t[i+1]) / (t[i+3] - t[i+1]) * (α[i+1] - α[i]) + α[i]
        b[i] = 2 * (α[i+1] - α[i]) / (t[i+3] - t[i+1])
    end
    for i = 1:n-1
        c[i, 1] = ((α[i+2] - α[i+1]) / (t[i+4] - t[i+2]) + (α[i] - α[i+1]) / (t[i+3] - t[i+1])) / (t[i+3] - t[i+2])
    end
    return PP(2, a, b, c, x)
end

function convertG(::Type{PP{2,T,TX}}, spl::BSplines.Spline) where {T,TX}
    t = BSplines.knots(spl.basis)
    α = spl.coeffs
    n = length(α)
    breakA, coef, l = bspp(t, α, n, 3)
    x = breakA
    a = @view(coef[1, :])
    b = @view(coef[2, :])
    c = coef[3:3, 1:l]'
    dx = (x[l+1] - x[l])
    a[l+1] = a[l] + dx * (b[l] + dx * (c[l, 1]))
    b[l+1] = b[l] + dx * (2 * c[l, 1])
    return PP(2, a, b, c, x)
end
function Base.convert(::Type{BSplines.Spline}, pp::PP{2,T,TX}) where {T,TX}
    basis = BSplines.BSplineBasis(3, pp.x)
    t = BSplines.knots(basis)
    n = length(pp.x)
    α = zeros(T, n + 1)
    for i = 1:n
        α[i] = pp.a[i] - pp.b[i] / 2 * (t[i+2] - t[i+1])
    end
    α[n+1] = pp.a[n]

    return BSplines.Spline(basis, α)
end

function Base.convert(::Type{BSplines.Spline}, pp::PP{3,T,TX}) where {T,TX}
    basis = BSplines.BSplineBasis(4, pp.x)
    t = BSplines.knots(basis)
    n = length(pp.x)
    α = zeros(T, n + 2)
    α[1] = pp.a[1]
    #spl.(t[3:end-2]) + spl.(t[3:end-2],Derivative(1))/3 .* (t[4:end-1]-t[3:end-2]-(t[3:end-2]-t[2:end-3])) - spl.(t[3:end-2],Derivative(2))/6 .* (t[4:end-1]-t[3:end-2]).*(t[3:end-2]-t[2:end-3])
    for i = 1:n-1
        dt1 = t[i+3] - t[i+2]
        dt2 = t[i+4] - t[i+3]
        α[i+1] = pp.a[i] + pp.b[i] / 3 * (dt2 - dt1) - pp.c[i, 1] / 3 * dt2 * dt1
    end
    α[n+1] = pp.a[n] + pp.b[n] / 3 * (t[n+4] - 2 * t[n+3] + t[n+2])
    α[n+2] = pp.a[n]
    return BSplines.Spline(basis, α)
end


function Base.convert(::Type{PPInterpolation.PP{3,T,TX}}, spl::BSplines.Spline) where {T,TX}
    t = BSplines.knots(spl.basis)
    α = spl.coeffs
    n = length(α)
    breakA, coef, l = bspp(t, α, n, 4)
    x = breakA
    a = @view(coef[1, :])
    b = @view(coef[2, :])
    c = coef[3:4, 1:l]'
    dx = (x[l+1] - x[l])
    a[l+1] = a[l] + dx * (b[l] + dx * (c[l, 1] + dx * (c[l, 2])))
    b[l+1] = b[l] + dx * (2 * c[l, 1] + dx * 3 * c[l, 2])
    return PPInterpolation.PP(3, a, b, c, x)
end


function bspp(t::AbstractArray{T}, bcoef::AbstractArray{B}, n::Int, k::Int) where {T,B}
    #=  !! BSPP converts from b-spline representation to pp representation
      !
      !     input
      !       t     knot sequence of length n+k
      !       bcoef b-spline coefficient sequence of length n
      !       n     length of bcoef
      !       k     order of the b-splines
      !
      !     output
      !       break breakpoint sequence, of length l+1, containing
      !             (in increasing order) the distinct points of the
      !             sequence t(k),...,t(n+1).
      !       coef  kxl matrix where coef(i,j) = (i-1)st right derivative
      !             of the pp at break(j) divided by factorial(i-1).
      !       l     number of polynomials which form the pp
      !
      !     work area
      !       wk    2-dimensional array of dimension (k,k+1)
      !     ------------------
        real t(*),bcoef(n),break(*),coef(k,*),wk(k,*)
      !     ------------------=#
    l = 0
    breakA = zeros(T, n + 1)
    coef = zeros(B, (k, n))
    wk = zeros(B, (k, k + 1))
    breakA[1] = t[k]
    if k != 1
        km1 = k - 1
        kp1 = k + 1
        #          general k-th order case

        for left = k:n
            if t[left] != t[left+1]
                l = l + 1
                breakA[l+1] = t[left+1]
                for j = 1:k
                    jj = left - k + j
                    wk[j, 1] = bcoef[jj]
                end
                for j = 1:km1
                    jp1 = j + 1
                    kmj = k - j
                    for i = 1:kmj
                        il = i + left
                        ilkj = il - kmj
                        diff = t[il] - t[ilkj]
                        wk[i, jp1] = (wk[i+1, j] - wk[i, j]) / diff
                    end
                end
                wk[1, kp1] = one(B)
                x = t[left]
                coef[k, l] = wk[1, k]
                a = one(B)
                for j = 1:km1
                    jp1 = j + 1
                    s = zero(B)
                    for i = 1:j
                        il = i + left
                        ilj = il - j
                        term = wk[i, kp1] / (t[il] - t[ilj])
                        wk[i, kp1] = s + (t[il] - x) * term
                        s = (x - t[ilj]) * term
                    end
                    wk[jp1, kp1] = s
                    s = zero(B)
                    kmj = k - j
                    for i = 1:jp1
                        s = s + wk[i, kmj] * wk[i, kp1]
                    end
                    a = (a * kmj) / j
                    coef[kmj, l] = a * s

                end
            end
        end
        return @view(breakA[1:l+1]), @view(coef[:, 1:l+1]), l
    end
    #!          piecewise constant case
    for left = k:n
        if (t[left] != t[left+1])
            l += 1
            breakA[l+1] = t[left+1]
            coef[1, l] = bcoef[left]
        end
    end
    return @view(breakA[1:l+1]), @view(coef[:, 1:l+1]), l
end
