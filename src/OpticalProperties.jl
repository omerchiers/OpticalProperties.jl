module OpticalProperties

using MyPhysicalConstants, Interpolations

# Functions
export convert_prop, skin_depth,
       permittivity, refractive_index,eps_si,eps_sin


# Types
export OptProp,
       Model, Bruggeman, MaxwellGarnett,SimpleMixing
       MaterialFile, Sellmeier,
       Cbn, Sic, Si_cst, Al,
       Au,Au_latella,Cst

# Constants
export Cu,SiO2,Si,SiN

abstract type OptProp end

# Generic models
struct Model <: OptProp
    eps0   :: Float64
    wp     :: Float64
    w0     :: Float64
    gamma0 :: Float64
    gamma1 :: Float64
end
Model(eps0,wp,w0,gamma0) = Model(eps0,wp,w0,gamma0,0.0)

struct Sellmeier{T} <: OptProp
    values :: Array{Tuple{T,T},1}
end

const Si  = Sellmeier([(10.6684293 , 0.301516485), (0.0030434748 , 1.13475115), (1.54133408 , 1104.0)])
const SiN = Sellmeier([(3.0249 , 0.1353406), (40314.0 , 1239.842)])


# Effective medium models
struct Bruggeman{T, U , V} <: OptProp
   phase1  :: T
   phase2  :: U
   volfrac :: V
end

struct MaxwellGarnett{T, U , V} <: OptProp
   particle :: T
   matrix   :: U
   volfrac  :: V
end

struct SimpleMixing{T,U,V} <: OptProp
   phase1  :: T
   phase2  :: U
   volfrac :: V
end

# From file containing measurements
struct MaterialFile{T,U} <: OptProp
    re :: T
    im :: U
end



# Dielectrics
struct Sic <: OptProp end
struct Cbn <: OptProp end
struct Si_cst <: OptProp end

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

# Constant permittivity
struct Cst <: OptProp
    val :: Complex128
end
Cst() = Cst(1.0+im*0.0)

# Refractive index
refractive_index(material :: OptProp, w) = sqrt(permittivity(material,w))

# Skin depth
skin_depth(material :: OptProp , w ) = c0/imag(refractive_index(material,w))/w


# Create interpolations
include("interpolation_objects.jl")

"""
    permittivity(material,w)

Compute the dielectric permittivity for a material.

# Arguments
* `material :: OptProp` : is a type for the material you want to pass
* `w :: Real` : pulsation in rad s^(-1)

# Example
```julia
julia> permittivity(Sic(),1e13)
6.805820438080644 + 0.002847538251166107im
```
"""
function permittivity() end


function permittivity(material::Model,w) :: Complex128
     eps0   = material.eps0
     wp     = material.wp
     w0     = material.w0
     gamma0 = material.gamma0
     gamma1 = material.gamma1
    return  eps0 + wp^2/(w0*w0 - w*w - im*w*(gamma0+gamma1))
end



function permittivity(material::Cbn,w) :: Complex128
    eps_fin = 4.46 + 0.0*im
    w_lo    = 2.451e14 # rad/s
    w_to    = 1.985e14 # rad/s
    gamma   = 9.934e11 # rad/s
    return eps_fin*(w^2-w_lo^2 + im*gamma*w)/(w^2-w_to^2 + im*gamma*w)
end

function permittivity(material::Sic,w) :: Complex128
    eps_fin = 6.7 + 0.0*im
    w_lo    = 1.827e14 # rad/s
    w_to    = 1.495e14 # rad/s
    gamma   = 8.971e11 # rad/s
    return eps_fin*(w^2-w_lo^2 + im*gamma*w)/(w^2-w_to^2 + im*gamma*w)
end


function permittivity(material::Al,w) :: Complex128
    alum = Model(1.0,2.24e16,0.0,1.22e14)
    return permittivity(alum,w)
end

function permittivity(material::Au,w) :: Complex128
    gamma1 = material.a*vf_au/material.mfp
    gold   = Model(9.4,13584.25e12,0.0,109.96e12,gamma1)
    return permittivity(gold,w)
end

function permittivity(material::Au_latella,w) :: Complex128
    gamma1 = material.a*vf_au/material.mfp
    gold   = Model(1.0,1.37e16,0.0,5.32e13,gamma1)
    return permittivity(gold,w)
end

function permittivity(material::Si_cst,w) :: Complex128
    return 11.7+im*0.0
end

function permittivity(material::Cst,w) :: Complex128
    return material.val
end

function permittivity(material::Bruggeman,w) :: Complex128
    eps1 = permittivity(material.phase1,w)
    eps2 = permittivity(material.phase2,w)
    p    = sqrt(eps1)/sqrt(eps2)
    b    = ((3.0*material.volfrac-1.0)*(1.0/p - p)+p)/4.0
    z    = b+sqrt(b*b + 0.5)
    return transpose(z*sqrt(eps2)*sqrt(eps1))
end

function permittivity(material::MaxwellGarnett,w) :: Complex128
    epsm = permittivity(material.matrix,w)
    epsp = permittivity(material.particle,w)
    f    = material.volfrac
    num  = 3.0*epsm + (1.0 + 2.0*f)*(epsp - epsm)
    den  = 3.0*epsm + (1.0 - f)*(epsp - epsm)
    return epsm*num/den
end

function permittivity(material::SimpleMixing,w) 
    eps1 = permittivity(material.phase1,w)
    eps2 = permittivity(material.phase2,w)
    f    = material.volfrac
    return f*eps1 + (1-f)*eps2
end

function permittivity(material :: MaterialFile, w)
    lamb = 2.0*pi*c0/w*1e6 # wavelength in microns
    return (material.re[lamb]^2 - material.im[lamb]^2) + im*2*material.im[lamb]*material.re[lamb]
end

function permittivity(material :: Sellmeier, w)
    lamb   = 2.0*pi*c0/w*1e6
    eps    = 1.0 + im*0.0
    for val in material.values
        eps += val[1]/(1-(val[2]/lamb)^2)
    end
    return eps
end
eps_si(x)  = 1+10.6684293./(1-(0.301516485./x).^2)+0.0030434748./(1-(1.13475115./x).^2)+1.54133408./(1-(1104./x).^2)
eps_sin(x) = 1+3.0249./(1-(0.1353406./x).^2)+40314./(1-(1239.842./x).^2)
end # module
