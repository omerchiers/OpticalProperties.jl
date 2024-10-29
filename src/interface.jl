abstract type AbstractProperty end

abstract type OpticalProperty <: AbstractProperty end
abstract type PropertyModel <: OpticalProperty end
abstract type PropertyEffective <: OpticalProperty end
abstract type PropertyData <: OpticalProperty end

abstract type Interband <: OpticalProperty end
abstract type Intraband <: OpticalProperty end


abstract type ElectricalProperty <: AbstractProperty end
abstract type MobilityModel <: ElectricalProperty end
abstract type ResistivityData <: ElectricalProperty end


"""
    permittivity(material::OpticalProperty, w)

Compute the dielectric permittivity for a material. w is the pulsation in rad/s.
To use the whole interface, it is required to implement this method of a given material.

# Arguments
* `material :: OpticalProperty` : object containing the data in order to compute the permittivity
* `w :: Real` : pulsation in rad s^(-1)

"""
function permittivity(material::OpticalProperty, w)
end



"""
    refractive_index(permittivity)

give the refractive index when given the complex permittivity of a material

"""
refractive_index(permittivity) = sqrt(permittivity)


"""
    refractive_index(permittivity::OpticalProperty, w)

when the permittivity of a material

# Arguments
* `material :: OpticalProperty` : object containing the data in order to compute the permittivity
* `w :: Real` : pulsation in rad s^(-1)
"""
refractive_index(material::OpticalProperty, w) = refractive_index(permittivity(material, w))



"""
    skin_depth(permittivity, w)

give the skin depth when given the complex permittivity of a material and the frequency

# Arguments
* `permittivity :: Complex` : value of permittivity
* `w :: Real` : pulsation in rad s^(-1)
"""
skin_depth(permittivity, w) = c0 / imag(refractive_index(permittivity)) / w


"""
    skin_depth(material::OpticalProperty, w)

give the skin depth when given the complex permittivity of a material and the frequency

# Arguments
* `material :: OpticalProperty` : object containing the data in order to compute the skin depth
* `w :: Real` : pulsation in rad s^(-1)
"""
skin_depth(material::OpticalProperty, w) = c0 / imag(refractive_index(material, w)) / w
