
# Material files
include("./data/Cu.jl")
include("./data/SiO2.jl")

# Create interpolations
# Cu
knots = (cu_data[:,1],)
const itp_cu_re = interpolate(knots, cu_data[:,2], Gridded(Linear()))
const itp_cu_im = interpolate(knots, cu_data[:,3], Gridded(Linear()))

const Cu = MaterialFile(itp_cu_re,itp_cu_im)

# SiO2
knots = (sio2_data[:,1],)
const itp_sio2_re = interpolate(knots, sio2_data[:,2], Gridded(Linear()))
const itp_sio2_im = interpolate(knots, sio2_data[:,3], Gridded(Linear()))

const SiO2 = MaterialFile(itp_sio2_re,itp_sio2_im)
