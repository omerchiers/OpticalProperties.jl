
# Constant permittivity
struct Cst{T} <: PropertyModel
    permittivity::Complex{T}
end

const Vacuum = Cst(1.0 + im * 0.0)
const Si_cst = Cst(11.7 + im * 0.0)

"""
    permittivity(model::Cst, w)

Give the dielectric permittivity independent of frequency.
To obtain the permittivity, the frequency is not needed. However, it is still required to pass the value
to the function to satisfy the interface.
"""
permittivity(model::Cst, w) = model.permittivity


# SiO2
knots = (sio2_data[:,1],)
const itp_sio2_re = interpolate(knots, sio2_data[:,2], Gridded(Linear()))
const itp_sio2_im = interpolate(knots, sio2_data[:,3], Gridded(Linear()))
const SiO2 = PropertyInterpolated(:refr_ind, :wavelength, itp_sio2_re,itp_sio2_im)



# Water
knots = (h2o_data[:,1],)
const itp_h2o_re = interpolate(knots, h2o_data[:,2], Gridded(Linear()))
const itp_h2o_im = interpolate(knots, h2o_data[:,3], Gridded(Linear()))
const H2O = PropertyInterpolated(:refr_ind, :wavelength, itp_h2o_re,itp_h2o_im)

