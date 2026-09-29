module OpticalProperties

using MyPhysicalConstants, Interpolations, QuadGK, StaticArrays
import ElectricalProperties

# Interface functions
export permittivity,
    refractive_index,
    skin_depth

# Abstract types
export OpticalProperty

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
    SiC, SiN, SiN20_300, SiO2,
    Ti, TiN, TiN_cethil, TiW, TiW_v2,
    W,
    Vacuum

# the istances for those models require prefixing with OpticalProperties.
# This is to keep the same convetion as for the same functions from 
# ElectricalProperties which have exactly the same naming and arguments but give back instead 
# a mobility model. 
public Si_p_masetti, Si_n_masetti,
    Si_p_sze, Si_n_sze,
    Si_p_regianni, Si_n_regianni


#data necessary for interpolations
include("./data/SiO2.jl")
include("./data/Cu.jl")
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
