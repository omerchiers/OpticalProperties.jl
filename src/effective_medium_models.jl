
# Types for effective medium models
struct Bruggeman{T,U,V} <: PropertyEffective
    phase1::T
    phase2::U
    volfrac::V
end


"""
    permittivity(model::Bruggeman, w)

Compute the effective dielectric permittivity of a 2 phase mixture using the Bruggeman model.
Here, the model takes as an iput value, the model or data for some material.

# Arguments
* `model :: Bruggeman` : is a type for the material you want to pass
* `w :: Real` : pulsation in rad s^(-1)
"""
function permittivity(model::Bruggeman, w)
    (; phase1, phase2, volfrac) = model
    eps1 = permittivity(phase1, w)
    eps2 = permittivity(phase2, w)
    p = sqrt(eps1) / sqrt(eps2)
    b = ((3.0 * volfrac - 1.0) * (1.0 / p - p) + p) / 4.0
    z = b + sqrt(b * b + 0.5)
    return z * sqrt(eps2) * sqrt(eps1)
end


struct MaxwellGarnett{T,U,V} <: PropertyEffective
    particle::T
    matrix::U
    volfrac::V
end

"""
    permittivity(model::MaxwellGarnett, w)

the effective dielectric permittivity for spherical particles much smaller than the incident wavelength 
dispersed in a matrix using the MaxwellGarnett model.
Here, the model takes as an iput value, the model or data for some material.

# Arguments
* `model :: MaxwellGarnett` : is a type for the material you want to pass
* `w :: Real` : pulsation in rad s^(-1)
"""
function permittivity(model::MaxwellGarnett, w)
    (; particle, matrix, volfrac) = model
    epsm = permittivity(matrix, w)
    epsp = permittivity(particle, w)
    f = volfrac
    num = 3.0 * epsm + (1.0 + 2.0 * f) * (epsp - epsm)
    den = 3.0 * epsm + (1.0 - f) * (epsp - epsm)
    return epsm * num / den
end



struct ParallelMixing{T,U,V} <: PropertyEffective
    phase1::T
    phase2::U
    volfrac::V
end

"""
    permittivity(model::ParallelMixing, w)

Compute the effective dielectric permittivity of two layers in serie. 
The electric field is parallel to the interface between both layers.
Here, the model takes as an iput value, the model or data for some material.

# Arguments
* `model :: ParallelMixing` : is a type for the material you want to pass
* `w :: Real` : pulsation in rad s^(-1)
"""
function permittivity(model::ParallelMixing, w)
    (; phase1, phase2, volfrac) = model
    eps1 = permittivity(phase1, w)
    eps2 = permittivity(phase2, w)
    f = volfrac
    return f * eps1 + (1 - f) * eps2
end



struct PerpendicularMixing{T,U,V} <: PropertyEffective
    phase1::T
    phase2::U
    volfrac::V
end

"""
    permittivity(model::PerpendicularMixing, w)

Compute the effective dielectric permittivity of two layers in in parallel. 
The electric field is perpendicular to the interface between both layers.

# Arguments
* `model :: PerpendicularMixing` : is a type for the material you want to pass
* `w :: Real` : pulsation in rad s^(-1)
"""
function permittivity(model::PerpendicularMixing, w)
    (; phase1, phase2, volfrac) = model
    eps1 = permittivity(phase1, w)
    eps2 = permittivity(phase2, w)
    f = volfrac
    return 1 / (f / eps1 + (1 - f) / eps2)
end
