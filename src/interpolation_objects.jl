
# Material files
include("./data/Cu.jl")
include("./data/SiO2.jl")
include("./data/SiO2_T.jl")
include("./data/resistivity_n-Si.jl")
include("./data/resistivity_p-Si.jl")

# Create interpolations
# Cu
knots = (cu_data[:,1],)
const itp_cu_re = interpolate(knots, cu_data[:,2], Gridded(Linear()))
const itp_cu_im = interpolate(knots, cu_data[:,3], Gridded(Linear()))

const Cu = MaterialFile(:refr_ind, :wavelength, itp_cu_re,itp_cu_im)

# SiO2
knots = (sio2_data[:,1],)
const itp_sio2_re = interpolate(knots, sio2_data[:,2], Gridded(Linear()))
const itp_sio2_im = interpolate(knots, sio2_data[:,3], Gridded(Linear()))

const SiO2 = MaterialFile(:refr_ind, :wavelength, itp_sio2_re,itp_sio2_im)

## SiO2 for T = 416, 511, 626, 691, 789 K
knots = (sio2_T_data[:,1],)
const itp_sio2_416_re = interpolate(knots, sio2_T_data[:,2], Gridded(Linear()))
const itp_sio2_416_im = interpolate(knots, sio2_T_data[:,3], Gridded(Linear()))
const SiO2_416 = MaterialFile(:permittivity, :frequency, itp_sio2_416_re,itp_sio2_416_im)

const itp_sio2_511_re = interpolate(knots, sio2_T_data[:,4], Gridded(Linear()))
const itp_sio2_511_im = interpolate(knots, sio2_T_data[:,5], Gridded(Linear()))
const SiO2_511 = MaterialFile(:permittivity, :frequency, itp_sio2_511_re,itp_sio2_511_im)

const itp_sio2_626_re = interpolate(knots, sio2_T_data[:,6], Gridded(Linear()))
const itp_sio2_626_im = interpolate(knots, sio2_T_data[:,7], Gridded(Linear()))
const SiO2_626 = MaterialFile(:permittivity, :frequency, itp_sio2_626_re,itp_sio2_626_im)

const itp_sio2_691_re = interpolate(knots, sio2_T_data[:,8], Gridded(Linear()))
const itp_sio2_691_im = interpolate(knots, sio2_T_data[:,9], Gridded(Linear()))
const SiO2_691 = MaterialFile(:permittivity, :frequency, itp_sio2_691_re,itp_sio2_691_im)

const itp_sio2_789_re = interpolate(knots, sio2_T_data[:,10], Gridded(Linear()))
const itp_sio2_789_im = interpolate(knots, sio2_T_data[:,11], Gridded(Linear()))
const SiO2_789 = MaterialFile(:permittivity, :frequency, itp_sio2_789_re,itp_sio2_789_im)



# Resitivity n-Si
knots = (nSi_data[:,1],)
const itp_nSi = interpolate(knots, nSi_data[:,2], Gridded(Linear()))
const nSi_sze = ResistivityFile(itp_nSi)

# Resitivity p-Si
knots = (pSi_data[:,1],)
const itp_pSi = interpolate(knots, pSi_data[:,2], Gridded(Linear()))
const pSi_sze = ResistivityFile(itp_pSi)
