
# From file containing measurements
struct PropertyInterpolated{T,U,V} <: PropertyData
    realp::V
    imagp::V
end


# Outer constructor to dispatch on options
PropertyInterpolated(T, U, re::V, im::V) where {V} = PropertyInterpolated{T,U,V}(re, im)


"""
    permittivity(model::PropertyInterpolated{:refr_ind,:wavelength}, w)

Compute the dielectric permittivity using data interpolation when data is given for the
complex refractive index and wavelength in microns. These options are specified using type parameters given 
as the symbols :refr_ind and :wavelength respectively.  
w is given in rad/s but the data requires the wavelength in microns. The permittivity method
does the conversion to microns. The frequancy range that is available depends on each material.

# Arguments
* `model :: PropertyInterpolated{:refr_ind,:wavelength}`
* `w :: Real` : pulsation in rad s^(-1)
"""
function permittivity(data::PropertyInterpolated{:refr_ind,:wavelength}, w)
    (; realp, imagp) = data
    lamb = 2.0 * pi * c0 / w * 1e6 # wavelength in microns
    return (realp(lamb)^2 - imagp(lamb)^2) + im * 2 * imagp(lamb) * realp(lamb)
end


"""
    permittivity(model::PropertyInterpolated{:permittivity,:frequency}, w)

When data interpolation is given for the
complex permittivity and frequency in rad/s. These options are specified using type parameters given 
as the symbols :permittivity and :frequancy respectively. w is given in rad/s.
The frequancy range that is available depends on each material.

# Arguments
* `model :: PropertyInterpolated{:permittivity,:frequency}`
* `w :: Real` : pulsation in rad s^(-1)
"""
function permittivity(data::PropertyInterpolated{:permittivity,:frequency}, w)
    (; realp, imagp) = data
    return realp(w) + im * imagp(w)
end



function Base.show(io::IO, mf::PropertyInterpolated{T,U}) where {T,U}
    println(io, "property: $T")
    println(io, "parameter: $U")
end