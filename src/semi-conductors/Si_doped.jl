

struct ResistivityInterpolated{T} <: ResistivityData
    rho::T
end

# Resitivity n-Si
knots = (nSi_data[:, 1],)
const itp_nSi = interpolate(knots, nSi_data[:, 2], Gridded(Linear()))
const nSi_sze = ResistivityInterpolated(itp_nSi)

# Resitivity p-Si
knots = (pSi_data[:, 1],)
const itp_pSi = interpolate(knots, pSi_data[:, 2], Gridded(Linear()))
const pSi_sze = ResistivityInterpolated(itp_pSi)




struct Masetti{T} <: MobilityModel
    μ1::T
    μ2::T
    μmax::T
    cr::T
    cs::T
    α::T
    β::T
    pc::T
    charge::T
end

const Si_ndoped = Masetti(68.5, 56.1, 1414.0, 9.2e16, 3.41e20, 0.711, 1.98, 0.0, 0.0)
const Si_pdoped = Masetti(44.9, 29.0, 470.5, 2.23e17, 6.10e20, 0.719, 2.0, 9.23e16, 1.0)




"""
    mobility(model, doping_level)

Compute the electric mobility for a material in cm^2 V^{-1} s^{-1}.
"""

function mobility end


"""
     mobility(model::Masetti, doping_level)

when given the Masetti model and a doping level.
# Arguments
* `model :: Masetti` : model using the Masetti data
* `doping_level :: Real` : electron or hole concentration in cm^{-3}

# Example
```julia
julia> mobility(Si_ndoped, 1e15)
4.593579841190333
```
"""

function mobility(model::Masetti, doping_level)
    (; μ1, μ2, μmax, cr, cs, α, β, pc, charge) = model
    N = doping_level
    return μ1 * exp(-charge * pc / N) + (μmax - (1.0 - charge) * μ1) / (1 + (N / cr)^α) - μ2 / (1 + (cs / N)^β)
end



"""
    resistivity(mobility, doping_level)

Compute the electric resistivity for a material in Ohm cm using the formula
 ρ = 1/(n ⋅ e ⋅ μ) where ρ is the resistivity (Ohm cm),  n the doping level (cm^{-3}), e the elementary charge (C) 
 and μ the carrier mobility in cm^2 V^{-1} s^{-1}.

# Arguments
* `mobility` : the material for which you need the resisitivity
* `doping_level` : impurity concentration in cm^{-3}
"""
resistivity(mobility, doping_level) = 1 / (doping_level * electron * mobility)


"""
    resistivity(model::MobilityModel, doping_level)

when given a MobilityModel

# Arguments
* `model` : the model describing the mobility
* `doping_level` : impurity concentration in cm^{-3}

# Example
```julia
julia> resistivity(nSi,1e15)
4.593579841190333
```
"""
resistivity(model::MobilityModel, doping_level) = resistivity(mobility(model, doping_level), doping_level)


"""
    resistivity(data::ResistivityInterpolated, doping_level)

when given a interpolated resisitivity data

# Arguments
* `data` : data describing the resisitivity
* `doping_level` : impurity concentration in cm^{-3}

# Example
```julia
julia> resistivity(nSi,1e15)
4.593579841190333
```
"""
resistivity(data::ResistivityInterpolated, doping_level) = data.rho(doping_level)



"""
    Si_n_doped(mobility, doping_level)

Return DrudeLorentz object for n doped Silicon doped Silicon for the optical properties. 
doping_level varies between 3e19 and 5e20 cm^{-3}.

# Arguments
* `mobility` : passed as a value in cm^2 V^{-1} s^{-1} or as MobilityModel object.
* `doping_level` : impurity concentration in cm^{-3}
"""

Si_n_doped(mobility, doping_level) = DrudeLorentz(11.7, sqrt(doping_level * 1e6 * electron^2 / (0.27 * m0) / epsilon0), 0.0, doping_level * 1e4 * electron^2 * resistivity(mobility, doping_level) / (0.27 * m0), 0.0)

"""
    Si_p_doped(mobility, doping_level)

Return DrudeLorentz object for p doped Silicon doped Silicon for the optical properties. 
doping_level varies between 3e19 and 5e20 cm^{-3}.

# Arguments
* `mobility` : passed as a value in cm^2 V^{-1} s^{-1} or as MobilityModel object.
* `doping_level` : impurity concentration in cm^{-3}
"""
Si_p_doped(mobility, doping_level) = DrudeLorentz(11.7, sqrt(doping_level * 1e6 * electron^2 / (0.34 * m0) / epsilon0), 0.0, doping_level * 1e4 * electron^2 * resistivity(mobility, doping_level) / (0.34 * m0), 0.0)



