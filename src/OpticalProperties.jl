module OpticalProperties

using MyPhysicalConstants, Interpolations, QuadGK, StaticArrays, ElectricalProperties

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
    ResistivityFile,
    MobilityModel

# Instances
export Au, Au_latella, Al, AlGaAs,
    cBN, Cu,
    H2O,
    Si, Si_cst,
    pSi_masetti_drude, nSi_masetti_drude,
    pSi_sze_drude, nSi_sze_drude,
    SiC, SiN, SiN20_300, SiO2,
    Ti, TiN, TiN_cethil, TiW, TiW_v2,
    W,
    Vacuum

#data necessary for interpolations
include("./data/SiO2.jl")
include("./data/Cu.jl")
include("./data/resistivity_n-Si.jl")
include("./data/resistivity_p-Si.jl")
include("./data/H2O.jl")


include("interface.jl") #defines the functions and methods of the interface
include("interpolated.jl") #defines the methods for interpolating data stored in ./data folder
include("metals.jl")
include("dielectrics.jl")
include("semi-conductors/AlGaAs.jl")
include("semi-conductors/Si_doped.jl")
include("semi-conductors/SiN.jl")
include("semi-conductors/others.jl")
include("effective_medium_models.jl")
include("list_materials.jl")



end # module
