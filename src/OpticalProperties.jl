module OpticalProperties


using MyPhysicalConstants, Interpolations, QuadGK, StaticArrays

# Interface functions
export permittivity,
    refractive_index,
    skin_depth

# Abstract types
export OpticalProperty, ElectricalProperty

# Structs
export DrudeLorentz,
    Polariton,
    Cst,
    Bruggeman,
    MaxwellGarnett,
    ParallelMixing,
    PerpendicularMixing,
    Sellmeier,
    ResistivityFile, MobilityModel,
    cBN, SiC,
    Si_n_doped, Si_p_doped

# Constants
export Cu, SiO2, Si, SiN, Vacuum,
    Au, Au_latella, Al,
    W, Ti, TiW, TiW_v2,
    TiN, TiN_cethil,
    pSi_masetti, nSi_masetti,
    pSi_sze, nSi_sze,
    SiC, cBN, Si_cst,
    AlGaAs


#data necessary for interpolations
include("./data/SiO2.jl")
include("./data/Cu.jl")
include("./data/resistivity_n-Si.jl")
include("./data/resistivity_p-Si.jl")


include("interface.jl") #defines the functions and methods of the interface
include("interpolated.jl") #defines the methods for interpolating data stored in ./data folder
include("metals.jl")
include("dielectrics.jl")
include("semi-conductors/AlGaAs.jl")
include("semi-conductors/Si_doped.jl")
include("semi-conductors/SiN.jl")
include("semi-conductors/others.jl")



end # module
