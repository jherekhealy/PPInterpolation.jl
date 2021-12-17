using JuMP
using Cbc
#using GLPK
export LaverySpline, QuinticLaverySpline, newQuinticLaverySpline, evaluate, evaluateDerivative


mutable struct QuinticLaverySpline{T <: Real}
    a::Vector{T}
    b::Vector{T}
    c::Vector{T}
    d::Vector{T}
    e::Vector{T}
    f::Vector{T}
    x::Vector{T}
end


function evaluate(self::QuinticLaverySpline, z::T) where {T}
    if z <= self.x[1]
        return self.b[1] * (z - self.x[1]) + self.a[1]
    elseif z >= self.x[end]
        rightSlope = self.b[end]
        return rightSlope * (z - self.x[end]) + self.a[end]
    end
    i = searchsortedfirst(self.x, z)  # x[i-1]<z<=x[i]
    if i > 1
        i -= 1
    end
    h = z - self.x[i]
    return self.a[i] + h * (self.b[i] + h * (self.c[i] + h*(self.d[i]+h*(self.e[i]+h*self.f[i]))))
end

(spl::QuinticLaverySpline)(x::Real) = evaluate(spl, x)

function evaluateDerivative(self::QuinticLaverySpline, z::T) where {T}
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
    return self.b[i] + h * (2*self.c[i] + h*(3*self.d[i]+h*(4*self.e[i]+h*5*self.f[i])))
end

function newLaverySpline(x::Vector{T}, y::Vector{T}) where {T}
    pp = CubicPP(T, length(y))
    computeLaverySpline(pp, x, y)
    return pp
end

function computeLaverySpline(pp::CubicPP{T}, x::Vector{T}, z::Vector{T}) where {T}
    h = x[2:end]-x[1:end-1] #h1 = x2-x1
    dz = (z[2:end]-z[1:end-1]) ./ h
    size = length(z)
    nIntervals = 10
    firstSize = (size-1)*nIntervals #first sum from 1 to size-1 (included) and interior sum from 1 to nInterval
    #1,1 1,2, ..., 1,nIntervals ...,size-1,1, size-1,2, ..., size-1,nIntervals
    eps = 1e-4
    m = Model(with_optimizer(Cbc.Optimizer))
    # m = Model(with_optimizer(GLPK.Optimizer, tm_lim = 60.0, msg_lev = GLPK.OFF))
    @variable(m, bv[1:size])
    @variable(m, s0[1:firstSize] >= 0) #first slack
    # @objective(m, Min, sum(s0))

    @variable(m, s1[1:size-2] >= 0) #second slack
    @objective(m, Min, sum(s0)+eps*sum(s1)) #s0 and s1 don't have the same size! => add s0 to s1 elements! or expand s0 and s1
    @constraint(m, s1 .>= bv[2:size-1])
    @constraint(m, s1  .>= -bv[2:size-1])
deltatk = 1.0/(nIntervals)

    for k=1:nIntervals
    for i=1:size-1
        dzi = dz[i]
                   tk = -0.5 + (k-0.5)/nIntervals
            @constraint(m,s0[k+(i-1)*nIntervals] >= ((-1+6*tk)*bv[i]+(1+6*tk)*bv[i+1]-12*dzi*tk)*deltatk)
            @constraint(m,s0[k+(i-1)*nIntervals] >= -((-1+6*tk)*bv[i]+(1+6*tk)*bv[i+1]-12*dzi*tk)*deltatk)
        end
    end
    @constraint(m, bv[1] == dz[1]-0.5*(bv[2]-dz[1]))
    @constraint(m, bv[size] == dz[size-1]-0.5*(bv[size-1]-dz[size-1]))
   optimize!(m)

    println("Objective: ", objective_value(m))

    b = zeros(size)

    for i=1:size
        println("b $i: ", value(bv[i]))
        b[i] = value(bv[i])
    end
    #b[1] = dz[1]-0.5*(b[2]-dz[1]) #natural BC
    #b[end] = dz[end]-0.5*(b[end-1]-dz[end]) #natural BC

    c = (-(2*b[1:end-1]+b[2:end])+3*dz) ./ h
    d = (b[1:end-1]+b[2:end]-2*dz) ./ (h.^2)
    pp.a[1:end] = z
    pp.b[1:end] = b
    pp.c[1:end] = c
    pp.d[1:end] = d
    pp.x[1:end]  = x
end


function newQuinticLaverySpline(x::Vector{T}, z::Vector{T}, eps1=1e-3, eps2=1e-3, order=3) where {T}
    h = x[2:end]-x[1:end-1] #h1 = x2-x1
    dz = (z[2:end]-z[1:end-1]) ./ h
    size = length(z)
    nIntervals = 10
    firstSize = (size-1)*nIntervals #first sum from 1 to size-1 (included) and interior sum from 1 to nInterval
    #1,1 1,2, ..., 1,nIntervals ...,size-1,1, size-1,2, ..., size-1,nIntervals

    #first: compute slopes for a cubic.
    m = Model(with_optimizer(Cbc.Optimizer))
    # m = Model(with_optimizer(GLPK.Optimizer, tm_lim = 60.0, msg_lev = GLPK.OFF))
    @variable(m, bv[1:size])
    @variable(m, s0[1:firstSize] >= 0) #first slack
    # @objective(m, Min, sum(s0))

    @variable(m, s1[1:size-2] >= 0) #second slack
    @objective(m, Min, sum(s0)+eps1*sum(s1)) #s0 and s1 don't have the same size! => add s0 to s1 elements! or expand s0 and s1
    @constraint(m, s1 .>= bv[2:size-1])
    @constraint(m, s1  .>= -bv[2:size-1])
deltatk = 1.0/(nIntervals)

    for k=1:nIntervals
    for i=1:size-1
        dzi = dz[i]
                   tk = -0.5 + (k-0.5)/nIntervals
            @constraint(m,s0[k+(i-1)*nIntervals] >= ((-1+6*tk)*bv[i]+(1+6*tk)*bv[i+1]-12*dzi*tk)*deltatk)
            @constraint(m,s0[k+(i-1)*nIntervals] >= -((-1+6*tk)*bv[i]+(1+6*tk)*bv[i+1]-12*dzi*tk)*deltatk)
        end
    end
    @constraint(m, bv[1] == dz[1]-0.5*(bv[2]-dz[1]))
    @constraint(m, bv[size] == dz[size-1]-0.5*(bv[size-1]-dz[size-1]))
   optimize!(m)

    println("Objective: ", objective_value(m))

    b = zeros(size)

    for i=1:size
        println("b $i: ", value(bv[i]))
        b[i] = value(bv[i])
    end

    #second: compute second derivative for quintic

    m = Model(with_optimizer(Cbc.Optimizer))
    # m = Model(with_optimizer(GLPK.Optimizer, tm_lim = 60.0, msg_lev = GLPK.OFF))
    @variable(m, cv[1:size])
    @variable(m, s0[1:firstSize] >= 0) #first slack
    @variable(m, s1[1:size-2] >= 0) #second slack
    @objective(m, Min, sum(s0)+eps2*sum(s1)) #s0 and s1 don't have the same size! => add s0 to s1 elements! or expand s0 and s1
    @constraint(m, s1 .>= cv[2:size-1])
    @constraint(m, s1  .>= -cv[2:size-1])

if order == 1
    for k=1:nIntervals
        for i=1:size-1
            deltatk = (x[i+1]-x[i])/(nIntervals)
    hk  = deltatk*(k-0.5)
    @constraint(m, s0[k+(i-1)*nIntervals] >= (b[i] + hk*cv[i] + 3*hk^2*((cv[i+1]-3*cv[i])/(2*h[i]) + 2*(5*dz[i]-3*b[i]-2*b[i+1])/(h[i]^2))+4*hk^3*((3*cv[i]-2*cv[i+1])/(2*h[i]^2)+(8*b[i]+7*b[i+1]-15*dz[i])/(h[i]^3))+5*(hk^4)*((cv[i+1]-cv[i])/(2*h[i]^3)+3*(2*dz[i]-b[i+1]-b[i])/(h[i]^4)))*deltatk)
    @constraint(m, -s0[k+(i-1)*nIntervals] <= (b[i] + hk*cv[i] + 3*hk^2*((cv[i+1]-3*cv[i])/(2*h[i]) + 2*(5*dz[i]-3*b[i]-2*b[i+1])/(h[i]^2))+4*hk^3*((3*cv[i]-2*cv[i+1])/(2*h[i]^2)+(8*b[i]+7*b[i+1]-15*dz[i])/(h[i]^3))+5*(hk^4)*((cv[i+1]-cv[i])/(2*h[i]^3)+3*(2*dz[i]-b[i+1]-b[i])/(h[i]^4)))*deltatk)
            end
        end
elseif order == 2
#second der constraint
for k=1:nIntervals
    for i=1:size-1
        deltatk = (x[i+1]-x[i])/(nIntervals)
        hk  = deltatk*(k-0.5)
    @constraint(m, s0[k+(i-1)*nIntervals] >= (cv[i] + 6*hk*((cv[i+1]-3*cv[i])/(2*h[i]) + 2*(5*dz[i]-3*b[i]-2*b[i+1])/(h[i]^2))+12*hk^2*((3*cv[i]-2*cv[i+1])/(2*h[i]^2)+(8*b[i]+7*b[i+1]-15*dz[i])/(h[i]^3))+20*(hk^3)*((cv[i+1]-cv[i])/(2*h[i]^3)+3*(2*dz[i]-b[i+1]-b[i])/(h[i]^4)))*deltatk)
@constraint(m, -s0[k+(i-1)*nIntervals] <= (cv[i] + 6*hk*((cv[i+1]-3*cv[i])/(2*h[i]) + 2*(5*dz[i]-3*b[i]-2*b[i+1])/(h[i]^2))+12*hk^2*((3*cv[i]-2*cv[i+1])/(2*h[i]^2)+(8*b[i]+7*b[i+1]-15*dz[i])/(h[i]^3))+20*(hk^3)*((cv[i+1]-cv[i])/(2*h[i]^3)+3*(2*dz[i]-b[i+1]-b[i])/(h[i]^4)))*deltatk)
        end
    end
#third der
elseif order == 3
for k=1:nIntervals
    for i=1:size-1
        deltatk = (x[i+1]-x[i])/(nIntervals)
        hk  = deltatk*(k-0.5)
    @constraint(m, s0[k+(i-1)*nIntervals] >= (6*((cv[i+1]-3*cv[i])/(2*h[i]) + 2*(5*dz[i]-3*b[i]-2*b[i+1])/(h[i]^2))+24*hk*((3*cv[i]-2*cv[i+1])/(2*h[i]^2)+(8*b[i]+7*b[i+1]-15*dz[i])/(h[i]^3))+60*(hk^2)*((cv[i+1]-cv[i])/(2*h[i]^3)+3*(2*dz[i]-b[i+1]-b[i])/(h[i]^4)))*deltatk)
@constraint(m, -s0[k+(i-1)*nIntervals] <= (6*((cv[i+1]-3*cv[i])/(2*h[i]) + 2*(5*dz[i]-3*b[i]-2*b[i+1])/(h[i]^2))+24*hk*((3*cv[i]-2*cv[i+1])/(2*h[i]^2)+(8*b[i]+7*b[i+1]-15*dz[i])/(h[i]^3))+60*(hk^2)*((cv[i+1]-cv[i])/(2*h[i]^3)+3*(2*dz[i]-b[i+1]-b[i])/(h[i]^4)))*deltatk)
        end
    end
elseif order == 4
    for k=1:nIntervals
        for i=1:size-1
            deltatk = (x[i+1]-x[i])/(nIntervals)
            hk  = deltatk*(k-0.5)
        @constraint(m, s0[k+(i-1)*nIntervals] >= (24*((3*cv[i]-2*cv[i+1])/(2*h[i]^2)+(8*b[i]+7*b[i+1]-15*dz[i])/(h[i]^3))+120*(hk)*((cv[i+1]-cv[i])/(2*h[i]^3)+3*(2*dz[i]-b[i+1]-b[i])/(h[i]^4)))*deltatk)
    @constraint(m, -s0[k+(i-1)*nIntervals] <= (24*((3*cv[i]-2*cv[i+1])/(2*h[i]^2)+(8*b[i]+7*b[i+1]-15*dz[i])/(h[i]^3))+120*(hk)*((cv[i+1]-cv[i])/(2*h[i]^3)+3*(2*dz[i]-b[i+1]-b[i])/(h[i]^4)))*deltatk)
            end
        end
end
    @constraint(m, cv[1] == 0)
    @constraint(m, cv[size] == 0)
   optimize!(m)

    println("Objective: ", objective_value(m))

    c = zeros(size)

    for i=1:size
        println("c $i: ", value(cv[i]))
        c[i] = value(cv[i])/2
    end

    #b[1] = dz[1]-0.5*(b[2]-dz[1]) #natural BC
    #b[end] = dz[end]-0.5*(b[end-1]-dz[end]) #natural BC
    d = (c[2:end]-3*c[1:end-1]) ./ h + 2*(5*dz - 3*b[1:end-1]-2*b[2:end]) ./ (h.^2)
    e = (3*c[1:end-1]-2*c[2:end]) ./ (h.^2) + (8*b[1:end-1]+7*b[2:end]-15*dz) ./ (h.^3)
    f = (c[2:end]-c[1:end-1]) ./ (h.^3) + 3*(2*dz-b[2:end]-b[1:end-1]) ./ (h.^4)
    return QuinticLaverySpline(copy(z), b, c, d, e,f,copy(x))
end
