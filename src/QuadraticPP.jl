export SchumakerDerivative, LamDerivative, makeSchumakerQuadraticPP

abstract type FirstDerivativeVariant end
struct SchumakerDerivative <: FirstDerivativeVariant end
struct LamDerivative{T} <: FirstDerivativeVariant
	ξ::T
end
LamDerivative{T}() where {T} = LamDerivative(one(T) / 2)



function makeSchumakerQuadraticPP(
	t::AbstractArray{TX},
	z::AbstractArray{T},
	variant::FirstDerivativeVariant;
 epsilon = sqrt(eps(T))
) where {T, TX}
	dt = t[2:end] - t[1:end-1]
	dz = z[2:end] - z[1:end-1]
	d = dz ./ dt
	s = computeFirstDerivatives(variant, t, z, epsilon=epsilon)
	j = 0
	x = Vector{TX}(undef, 0)
	a = Vector{T}(undef, 0)
	b = Vector{T}(undef, 0)
	c = Vector{T}(undef, 0)
    n = length(t)
	for i ∈ 1:n-1
		if abs(s[i] + s[i+1] - 2 * d[i]) <= epsilon * abs(d[i])
			j += 1
			x = push!(x, t[i])
			a = push!(a, z[i])
			b = push!(b, s[i])
			c = push!(c, (s[i+1] - s[i]) / (2 * (t[i+1] - t[i])))
		else
			aa = s[i] - d[i]
			bb = s[i+1] - d[i]
			ei = 	if aa * bb >= zero(T)
				(t[i+1] + t[i]) / 2
			else
				if abs(aa) > abs(bb)
					t[i] + bb * (t[i+1] - t[i]) / (s[i+1] - s[i])
				else
					t[i+1] + aa * (t[i+1] - t[i]) / (s[i+1] - s[i])
				end
			end
			if ei == t[i]
				#can happen when bb close to zero and machine epsilon makes ei=ti
				ei = t[i] + (t[i+1] - t[i]) / 2
			end
			sbi = 2 * d[i] - ((ei - t[i]) * s[i] + (t[i+1] - ei) * s[i+1]) / (t[i+1] - t[i])
			etai = (sbi - s[i]) / (ei - t[i])

			if abs(ei - t[i]) > zero(T)
				j += 1
				x = push!(x, t[i])
				a = push!(a, z[i])
				b = push!(b, s[i])
				c = push!(c, etai / 2)
			end
			if abs(ei - t[i+1]) > zero(T)
				j += 1
				x = push!(x, ei)
				a = push!(a, z[i] + s[i] * (ei - t[i]) + etai * (ei - t[i])^2 / 2)
				b = push!(b, sbi)
				c = push!(c, (s[i+1] - sbi) / (2 * (t[i+1] - ei)))
			end
		end
	end
	x = push!(x, t[end])
	a = push!(a, z[end])
	pp = PP(2, a, b, reshape(c,length(c),1), x)
	return pp
end

function computeFirstDerivatives(variant::LamDerivative{TZ}, t::AbstractArray{TX}, z::AbstractArray{TZ}; epsilon = sqrt(eps(TZ))) where {TX, TZ}
	ξ = variant.ξ
	η = one(TX) - ξ

	n = length(t)
	dt = t[2:end] - t[1:end-1]
	dz = z[2:end] - z[1:end-1]
	d = dz ./ dt
	s = zeros(TZ, n)
	for i ∈ 2:n-1
		if d[i] * d[i-1] > zero(TZ)
			if (abs(d[i-1]) - abs(d[i])) * (ξ - one(TX) / 2) >= zero(TZ)
				s[i] = d[i-1] * d[i] / (ξ * d[i-1] + η * d[i])
			else
				s[i] = d[i-1] * d[i] / (ξ * d[i] + η * d[i-1])
			end
		end
	end
	if d[1] * (2 * d[1] - s[2]) > zero(TZ)
		s[1] = (2 * d[1] - s[2])
	else
		s[1] = s[1] * epsilon
	end
	if d[n-1] * (2 * d[n-1] - s[n-1]) > zero(TZ)
		s[n] = (2 * d[n-1] - s[n-1])
	else
		s[n] = s[n-1] * epsilon

	end
	return s
end

function computeFirstDerivatives(variant::SchumakerDerivative, t::AbstractArray{TX}, z::AbstractArray{TZ}; epsilon = zero(TZ)) where {TX, TZ}
	n = length(t)
	dt = t[2:end] - t[1:end-1]
	dz = z[2:end] - z[1:end-1]
	d = dz ./ dt
	l = @. sqrt(dt^2 + dz^2)
	#derivatives not provided, compute sb.
	s = zeros(TZ, n)
	for i ∈ 2:n-1
		if d[i] * d[i-1] > zero(TZ)
			s[i] = (l[i-1] * d[i-1] + l[i] * d[i]) / (l[i-1] + l[i])
		end
	end
	s[1] = (3 * d[1] - s[2]) / 2 #neg sign?
	s[n] = (3 * d[n-1] - s[n-1]) / 2
	if s[1] * s[2] < zero(TZ)
		alpha = one(TZ) / 4
		for i ∈ 1:32
			if s[1] * s[2] >= zero(TZ)
				break
			end
			s[1] = (1 + alpha) * d[1] - alpha * s[2]
			alpha /= 2
		end

	end
	if s[n] * s[n-1] < zero(TZ)
		alpha = one(TZ) / 4
		for i ∈ 1:32
			if s[n] * s[n-1] >= zero(TZ)
				break
			end

			s[n] = (1 + alpha) * d[n-1] - alpha * s[n-1]
			alpha /= 2
		end
	end
	return s
end
