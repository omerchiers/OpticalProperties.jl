
# Definitions of resisitivity and mobility should be moved to separate package ElectricalProperties.jl
# the content of this file should be moved to ElectricalProperties.jl

@enum DopingType begin
    n_doping
    p_doping
end

struct ResistivityInterpolated{T,U} <: ResistivityData
    rho::T
    effective_mass::U
    doping_level::U
end

# Resitivity n-Si using Sze data
knots = (nSi_data[:, 1],)
const itp_nSi = interpolate(knots, nSi_data[:, 2], Gridded(Linear()))
nSi_sze(Nn) = ResistivityInterpolated(itp_nSi, 0.27 * m0, Nn)

# Resitivity p-Si using Sze data
knots = (pSi_data[:, 1],)
const itp_pSi = interpolate(knots, pSi_data[:, 2], Gridded(Linear()))
pSi_sze(Np) = ResistivityInterpolated(itp_pSi, 0.34 * m0, Np)


struct SiDopedMasetti{T} <: MobilityModel
    μ1::T
    μ2::T
    μmax::T
    cr::T
    cs::T
    α::T
    β::T
    pc::T
    charge::T
    effective_mass::T
    doping_level::T
end

nSi_masetti(Nn) = SiDopedMasetti(68.5, 56.1, 1414.0, 9.2e16, 3.41e20, 0.711, 1.98, 0.0, 0.0, 0.27 * m0, Nn)
pSi_masetti(Np) = SiDopedMasetti(44.9, 29.0, 470.5, 2.23e17, 6.10e20, 0.719, 2.0, 9.23e16, 1.0, 0.34 * m0, Np)


doping_level(model) = model.doping_level
effective_mass(model) = model.effective_mass


"""
    mobility(resistivity, doping_level)

Compute the electric mobility for a material in cm^2 V^{-1} s^{-1} using the formula
 μ = 1/(n ⋅ e ⋅ ρ) where ρ is the resistivity (Ohm cm),  n the doping level (cm^{-3}), e the elementary charge (C) 
 and μ the carrier mobility (cm^2 V^{-1} s^{-1}).
"""
mobility(resistivity, doping_level) = 1 / (doping_level * electron * resistivity)


"""
     mobility(data)

When given the interpolated data.

# Example
```julia
julia> mobility(Si_sze_n(1e15)) #1e15 is the doping level in cm^{-3}
1358.7461946416095
```
"""

mobility(data) = mobility(resistivity(data), doping_level(data))


"""
     mobility(model::SiDopedMasetti)

When given the Masetti model for Silicon.

# Example
```julia
julia> mobility(nSi_masetti(1e15)) #1e15 is the doping level in cm^{-3}
1362.0563795023859
```
"""

function mobility(model::SiDopedMasetti)
    (; μ1, μ2, μmax, cr, cs, α, β, pc, charge, doping_level) = model
    N = doping_level
    return μ1 * exp(-charge * pc / N) + (μmax - (1.0 - charge) * μ1) / (1 + (N / cr)^α) - μ2 / (1 + (cs / N)^β)
end



"""
    resistivity(mobility, doping_level)

 Compute the electric resistivity for a material in Ohm cm using the formula
 ρ = 1/(n ⋅ e ⋅ μ) where ρ is the resistivity (Ohm cm),  n the doping level (cm^{-3}), e the elementary charge (C) 
 and μ the carrier mobility  (cm^2 V^{-1} s^{-1}).
"""
resistivity(mobility, doping_level) = 1 / (doping_level * electron * mobility)


"""
    resistivity(data::ResistivityInterpolated)

When given a interpolated resisitivity data.

# Example
```julia
julia> resistivity(nSi_sze(1e15))
4.593579841190333
```
"""
resistivity(data::ResistivityData) = data.rho(doping_level(data))


"""
    resistivity(model::MobilityModel)

When given a MobilityModel.

# Example
```julia
julia> resistivity(nSi_masetti(1e15))
4.582416134110431
```
"""
resistivity(model::MobilityModel) = resistivity(mobility(model), doping_level(model))


"""
    DrudeLorentz(material::ElectricalProperty)

Return DrudeLorentz object for the optical properties when given a material
with mobility, doping level and effective mass. 

The plasma frequency as obtained using the relation ωp = sqrt(n⋅e²⋅10⁶/me/ϵ0) where n is the electron density in
cm^{-3}, e is the elementary charge in C, me is the effective mass in kg and ϵ0 the vacuum permittivity.

The damping constant γ0 = n⋅e²⋅10⁴⋅ρ / me where \rho is the resistivity in Ω cm
 """
function DrudeLorentz(material::ElectricalProperty)
    mu = mobility(material)
    N = carrier_concentration(material)
    me = effective_mass(material)
    return DrudeLorentz(11.7,
        sqrt(N * 1e6 * electron^2 / me / epsilon0),
        0.0,
        N * 1e4 * electron^2 * resistivity(mu, N) / me,
        0.0)
end

nSi_masetti_drude(Nn) = DrudeLorentz(nSi_masetti(Nn))
pSi_masetti_drude(Np) = DrudeLorentz(pSi_masetti(Np))

nSi_sze_drude(Nn) = DrudeLorentz(nSi_sze(Nn))
pSi_sze_drude(Np) = DrudeLorentz(pSi_sze(Np))



