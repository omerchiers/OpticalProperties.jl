

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

Si_n_masetti_drude(Nn) = DrudeLorentz(Si_n_masetti(Nn))
Si_p_masetti_drude(Np) = DrudeLorentz(Si_p_masetti(Np))

Si_n_sze_drude(Nn) = DrudeLorentz(Si_n_sze(Nn))
Si_p_sze_drude(Np) = DrudeLorentz(Si_p_sze(Np))

Si_n_regianni_drude(Nd, Na, T) = DrudeLorentz(Si_n_regianni(Nd, Na, T))
Si_p_regianni_drude(Nd, Na, T) = DrudeLorentz(Si_p_regianni(Nd, Na, T))

