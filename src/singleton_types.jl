#Old way of organizing package, should be deprecated in the future
# Dielectrics
struct Sic <: OptProp end
struct Cbn <: OptProp end

# Conductors
struct Al <: OptProp end
struct Au <: OptProp
    mfp :: Float64
    a   :: Float64
end
Au(mfp) = Au(mfp,1.0)
Au() = Au(1.0,0.0)

struct Au_latella{T <: Real} <: OptProp
    mfp :: T
    a   :: T
end
Au_latella(mfp) = Au_latella(mfp,1.0)
Au_latella() = Au_latella(1.0,0.0)


function permittivity(material::Cbn,w) :: Complex{Float64}
    eps_fin = 4.46 + 0.0*im
    w_lo    = 2.451e14 # rad/s
    w_to    = 1.985e14 # rad/s
    gamma   = 9.934e11 # rad/s
    return eps_fin*(w^2-w_lo^2 + im*gamma*w)/(w^2-w_to^2 + im*gamma*w)
end

function permittivity(material::Sic,w) :: Complex{Float64}
    eps_fin = 6.7 + 0.0*im
    w_lo    = 1.827e14 # rad/s
    w_to    = 1.495e14 # rad/s
    gamma   = 8.971e11 # rad/s
    return eps_fin*(w^2-w_lo^2 + im*gamma*w)/(w^2-w_to^2 + im*gamma*w)
end


function permittivity(material::Al,w) :: Complex{Float64}
    alum = Model(1.0,2.24e16,0.0,1.22e14)
    return permittivity(alum,w)
end

function permittivity(material::Au,w) :: Complex{Float64}
    gamma1 = material.a*vf_au/material.mfp
    gold   = Model(9.4,13584.25e12,0.0,109.96e12,gamma1)
    return permittivity(gold,w)
end

function permittivity(material::Au_latella,w) :: Complex{Float64}
    gamma1 = material.a*vf_au/material.mfp
    gold   = Model(1.0,1.37e16,0.0,5.32e13,gamma1)
    return permittivity(gold,w)
end
