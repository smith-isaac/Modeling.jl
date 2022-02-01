module Modeling

export StateSpace, DiscreteSS

using LinearAlgebra

struct StateSpace
    A::Matrix{Float64}
    B::Vector{Float64}
    C::Matrix{Float64}
end

mutable struct DiscreteSS
    Ad::Matrix{Float64}
    Bd::Vector{Float64}
    C::Matrix{Float64}

    function DiscreteSS(ss::StateSpace, dt::Float64)
	this = new()
	this.Ad = exp(ss.A * dt)
	if det(ss.A) == 0
	    this.Bd = dt / 2 * (this.Ad + I) * ss.B
	else
	    this.Bd = inv(ss.A) * (this.Ad - I) * ss.B
	end
	this.C = ss.C
	return this
    end
end

end # module
