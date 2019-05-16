module OpticalProperties

using MyPhysicalConstants, Interpolations, Parameters, QuadGK

# Functions
export convert_prop, skin_depth,
       permittivity, refractive_index,eps_si,eps_sin,
       resistivity


# Types
export OptProp, ElectricalProperties,
       Model, Polariton,Cst,
       Bruggeman, MaxwellGarnett,
       SimpleMixingPar,SimpleMixingPerp,
       MaterialFile, Sellmeier,
       ResistivityFile, MobilityModel,
       Cbn, Sic, Al,
       Au,Au_latella,
       Au2,Au_latella2,
       Si_n_doped,Si_p_doped

# Constants
export Cu,SiO2,Si,SiN,Al2,Vacuum,
       pSi_masetti,nSi_masetti,
       pSi_sze,nSi_sze,
       SiC,cBN,Si_cst

abstract type AbstractMaterial end
abstract type OptProp <: AbstractMaterial end
abstract type ElectricalProperties <: AbstractMaterial end


# Generic models
struct Model <: OptProp
    eps0   :: Float64
    wp     :: Float64
    w0     :: Float64
    gamma0 :: Float64
    gamma1 :: Float64
end
Model(eps0,wp,w0,gamma0) = Model(eps0,wp,w0,gamma0,0.0)


# Gold
gamma1(mfp=1.0,a=0.0) = a*vf_au/mfp
Au2(mfp=1.0,a=0) = Model(9.4,13584.25e12,0.0,109.96e12,gamma1(mfp,a))
Au_latella2(mfp,a) = Model(1.0,1.37e16,0.0,5.32e13,gamma1(mfp,a))

#Aluminium
const Al2 = Model(1.0,2.24e16,0.0,1.22e14)

struct Polariton{T} <: OptProp
    eps_fin :: T
    w_lo    :: T # rad/s
    w_to    :: T # rad/s
    gamma   :: T # rad/s
end
SiC = Polariton(6.7,1.827e14,1.495e14,8.971e11)
cBN = Polariton(4.46,2.451e14,1.985e14,9.934e11)

struct Sellmeier{T} <: OptProp
    values :: Array{Tuple{T,T},1}
end

const Si  = Sellmeier([(10.6684293 , 0.301516485), (0.0030434748 , 1.13475115), (1.54133408 , 1104.0)])
const SiN = Sellmeier([(3.0249 , 0.1353406), (40314.0 , 1239.842)])


# Constant permittivity
struct Cst <: OptProp
    val :: Complex{Float64}
end

Cst()  = Cst(1.0+im*0.0)
const Vacuum = Cst()
const Si_cst = Cst(11.7 + im*0.0)



# Types for effective medium models
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

struct SimpleMixingPar{T,U,V} <: OptProp
   phase1  :: T
   phase2  :: U
   volfrac :: V
end

struct SimpleMixingPerp{T,U,V} <: OptProp
   phase1  :: T
   phase2  :: U
   volfrac :: V
end

# From file containing measurements
struct MaterialFile{T,U} <: OptProp
    re :: T
    im :: U
end

struct ResistivityFile{T} <: ElectricalProperties
    rho :: T
end

struct MobilityModel{T} <: ElectricalProperties
    μ1   :: T
    μ2   :: T
    μmax :: T
    cr   :: T
    cs   :: T
    α    :: T
    β    :: T
    pc   :: T
    charge :: T
end

const nSi_masetti = MobilityModel(68.5,56.1,1414.0,9.2e16,3.41e20,0.711,1.98,0.0,0.0)
const pSi_masetti = MobilityModel(44.9,29.0,470.5,2.23e17,6.10e20,0.719,2.0,9.23e16,1.0)


# Create interpolations
include("interpolation_objects.jl")
include("SiN.jl")

#To be deprecated
include("singleton_types.jl")

"""
    mobility(material,N)

Compute the electric mobility for a material in cm^2 V^{-1} s^{-1}.

# Arguments
* `material :: MobilityModel` : the material for which you need the resisitivity
* `N :: Real` : electron or hole concentration in cm^{-3}

# Example
```julia
julia> mobility(nSi,1e15)
4.593579841190333
```
"""

function mobility(material :: ElectricalProperties, N)
    @unpack μ1, μ2, μmax, cr, cs, α, β, pc, charge = material
    return μ1*exp(-charge*pc/N) + (μmax - (1.0 - charge)*μ1)/(1 + (N/cr)^α) - μ2/(1 + (cs/N)^β)
end



"""
    resistivity(material,N)

Compute the electric resistivity for a material in Ohm cm.

# Arguments
* `material :: ResistivityFile` : the material for which you need the resisitivity
* `N :: Real` : impurity concentration in cm^{-3}

# Example
```julia
julia> resistivity(nSi,1e15)
4.593579841190333
```
"""
resistivity(material :: ResistivityFile, N) = material.rho(N)
resistivity(material :: ElectricalProperties, N) = 1/(N*electron*mobility(material,N))



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


function permittivity(material::Model,w) :: Complex{Float64}
     eps0   = material.eps0
     wp     = material.wp
     w0     = material.w0
     gamma0 = material.gamma0
     gamma1 = material.gamma1
    return  eps0 + wp^2/(w0*w0 - w*w - im*w*(gamma0+gamma1))
end

function permittivity(material::Polariton,w) :: Complex{Float64}
    eps_fin = material.eps_fin
    w_lo    = material.w_lo
    w_to    = material.w_to
    gamma   = material.gamma
    return  eps_fin*(w^2-w_lo^2 + im*gamma*w)/(w^2-w_to^2 + im*gamma*w)
end

function permittivity(material::Cst,w) :: Complex{Float64}
    return material.val
end


function permittivity(material::Bruggeman,w) :: Complex{Float64}
    eps1 = permittivity(material.phase1,w)
    eps2 = permittivity(material.phase2,w)
    p    = sqrt(eps1)/sqrt(eps2)
    b    = ((3.0*material.volfrac-1.0)*(1.0/p - p)+p)/4.0
    z    = b+sqrt(b*b + 0.5)
    return transpose(z*sqrt(eps2)*sqrt(eps1))
end

function permittivity(material::MaxwellGarnett,w) :: Complex{Float64}
    epsm = permittivity(material.matrix,w)
    epsp = permittivity(material.particle,w)
    f    = material.volfrac
    num  = 3.0*epsm + (1.0 + 2.0*f)*(epsp - epsm)
    den  = 3.0*epsm + (1.0 - f)*(epsp - epsm)
    return epsm*num/den
end

function permittivity(material::SimpleMixingPar,w)
    eps1 = permittivity(material.phase1,w)
    eps2 = permittivity(material.phase2,w)
    f    = material.volfrac
    return f*eps1 + (1-f)*eps2
end

function permittivity(material::SimpleMixingPerp,w)
    eps1 = permittivity(material.phase1,w)
    eps2 = permittivity(material.phase2,w)
    f    = material.volfrac
    return 1/(f/eps1 + (1-f)/eps2)
end


function permittivity(material :: MaterialFile, w)
    lamb = 2.0*pi*c0/w*1e6 # wavelength in microns
    return (material.re(lamb)^2 - material.im(lamb)^2) + im*2*material.im(lamb)*material.re(lamb)
end

function permittivity(material :: Sellmeier, w)
    lamb   = 2.0*pi*c0/w*1e6
    eps    = 1.0 + im*0.0
    for val in material.values
        eps += val[1]/(1-(val[2]/lamb)^2)
    end
    return eps
end

# Refractive index
refractive_index(material :: OptProp, w) = sqrt(permittivity(material,w))

# Skin depth
skin_depth(material :: OptProp , w ) = c0/imag(refractive_index(material,w))/w


end # module
