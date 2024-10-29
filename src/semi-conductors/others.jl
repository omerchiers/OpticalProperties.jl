struct Polariton{T,U,V} <: PropertyModel
    eps_fin::T
    w_lo::U # rad/s
    w_to::U # rad/s
    gamma::V # 1/s
end

const SiC = Polariton(6.7, 1.827e14, 1.495e14, 8.971e11)
const cBN = Polariton(4.46, 2.451e14, 1.985e14, 9.934e11)


"""
    permittivity(model::Polariton, w)

Compute the dielectric permittivity for when given a Polariton model

# Arguments
* `model :: Polariton` : is a type for the material you want to pass
* `w :: Real` : pulsation in rad s^(-1)
"""
function permittivity(model::Polariton, w)::Complex{Float64}
    (; eps_fin, w_lo, w_to, gamma) = model
    return eps_fin * (w^2 - w_lo^2 + im * gamma * w) / (w^2 - w_to^2 + im * gamma * w)
end




struct Sellmeier{T} <: PropertyData
    values::Array{Tuple{T,T},1}
end

const Si = Sellmeier([(10.6684293, 0.301516485), (0.0030434748, 1.13475115), (1.54133408, 1104.0)])
const SiN = Sellmeier([(3.0249, 0.1353406), (40314.0, 1239.842)])

"""
    permittivity(model::Sellmeier, w)

for a Sellmeier model.
This model is used for transpararent dielectrics.
w is still given in rad/s but the data requires the wavelength in microns. The permittivity method
does the conversion to microns.

"""
function permittivity(model::Sellmeier, w)
    (; values) = model
    lamb = 2.0 * pi * c0 / w * 1e6
    eps = 1.0 + im * 0.0
    for val in values
        eps += val[1] / (1 - (val[2] / lamb)^2)
    end
    return eps
end