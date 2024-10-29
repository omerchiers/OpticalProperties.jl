

struct DrudeLorentz{T,U,V} <: PropertyModel
    eps0::T
    wp::U
    w0::U
    gamma0::V
    gamma1::V
end

DrudeLorentz(eps0, wp, w0, gamma0) = DrudeLorentz(eps0, wp, w0, gamma0, 0.0)
gamma1(mfp=1.0, a=0.0; vf=0.0) = a * vf / mfp


"""
    permittivity(model::DrudeLorentz, w)

Compute the dielectric permittivity for when given a DrudeLorentz model

# Arguments
* `model :: DrudeLorentz` : is a type for the material you want to pass
* `w :: Real` : pulsation in rad s^(-1)
"""
function permittivity(model::DrudeLorentz, w)
    (; eps0, wp, w0, gamma0, gamma1) = model
    return eps0 + wp^2 / (w0 * w0 - w * w - im * w * (gamma0 + gamma1))
end


# Gold
# Fit to Johnson-Christy data 
Au(mfp=1.0, a=0.0) = DrudeLorentz(9.4, 13584.25e12, 0.0, 109.96e12, gamma1(mfp, a; vf=vf_au))
# From Latella
Au_latella(mfp, a) = DrudeLorentz(1.0, 1.37e16, 0.0, 5.32e13, gamma1(mfp, a; vf=vf_au))

#Aluminium
Al(mfp=1.0, a=0.0) = DrudeLorentz(1.0, 2.24e16, 0.0, 1.22e14, gamma1(mfp, a; vf=vf_ag))

"Drude model of Tungsten: fitted from data obtained from Ordal et al."
const W = DrudeLorentz(6.0, 6.2e15, 0.0, 2.6e14, 0.0)

"Drude model of Titanium: fitted from data obtained from Ordal et al."
const Ti = DrudeLorentz(1.0, 3.82e15, 0.0, 719.6e11, 0.0)

"Drude model of Titanium Tungsten: obtained from Varpula et al."
const TiW = DrudeLorentz(1.0, sqrt(7.75e5 * 1e15 / 8.85e-12), 0.0, 1.0e15, 0.0)
const TiW_v2 = DrudeLorentz(1.0, sqrt(5.0e5 * 1e15 / 8.85e-12), 0.0, 1.0e15, 0.0)

"Drude model of Titanium Nitride obtained from Varpula et al."
const TiN = DrudeLorentz(1.0, sqrt(7.57e5 / 1.5e-15 / 8.85e-12), 0.0, 1.0 / (1.5e-15), 0.0)

"""
    gamma0_TiN_cethil(thickness)

    Drude model for TiN based on reflectance measurements performed at CETHIL 04/2021.
# Arguments
* `thickness`: layer thickness in meter
"""
function gamma0_TiN_cethil(thickness)
    th = thickness * 1e9
    if th <= 20
        gamma = 0.0087e15 * th + 0.8939e15
    else
        gamma = 1.13e15
    end
end


"""
    TiN_cethil(thickness)

    Returns the Model object for a Drude model with the parameters obtained from reflectance measurements performed at CETHIL 04/2021.
# Arguments
* `thickness`: layer thickness in meter
"""
TiN_cethil(thickness) = DrudeLorentz(1.0, 6.93e15, 0.0, gamma0_TiN_cethil(thickness), 0.0)

" Special case of TiN_cethil(thickness) where gamma does not depend on thickness"
const TiN_cethil_v2 = DrudeLorentz(1.0, 6.93e15, 0.0, 1.13e15, 0.0)




struct LorentzTerm{T} <: PropertyModel
    f::T
    wp::T
    w0::T
    gamma0::T
end


function susceptibility(model::LorentzTerm, w)
    (; f, wp, w0, gamma0) = model
    return f * wp^2 / (w^2 - w0^2 + im * w * gamma0)
end


struct DrudeTerm{T,U} <: PropertyModel
    wp::T
    gamma0::U
    gamma1::U
end

DrudeTerm(wp, gamma0) = DrudeTerm(wp, gamma0, 0.0)

function susceptibility(model::DrudeTerm, w)
    (; wp, gamma0, gamma1) = model
    return wp^2 / (w * (w + im * (gamma0 + gamma1)))
end


struct DrudeLorentzSerie{T,U} <: PropertyModel
    eps0::T
    drude::DrudeTerm{U}
    lorentz::Vector{LorentzTerm{U}}
end



function permittivity(model::DrudeLorentzSerie, w)
    (; eps0, drude, lorentz) = model
    ϵ = eps0 - susceptibility(drude, w)
    for term in lorentz
        ϵ -= susceptibility(term, w)
    end
    return ϵ
end



# Create interpolations
# Cu
knots = (cu_data[:, 1],)
const itp_cu_re = interpolate(knots, cu_data[:, 2], Gridded(Linear()))
const itp_cu_im = interpolate(knots, cu_data[:, 3], Gridded(Linear()))

const Cu = PropertyInterpolated(:refr_ind, :wavelength, itp_cu_re, itp_cu_im)
