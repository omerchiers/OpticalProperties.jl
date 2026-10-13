# OpticalProperties
This julia package computes permittivities and refractive indices for usual materials.
This is very much work in progress.

## Warning
This package is still in development and the interface might break regularly.

## Installation
From version 0.4.0 you will need Julia version 1.11 or higher. 

To install from the github repository, type in the Julia repl
```julia
julia> using Pkg
julia> Pkg.add("https://github.com/omerchiers/OpticalProperties.jl")
```

To install from cethil gitlab, add first CETHIL_REGISTRY, then type in the Julia repl
```julia
julia> using Pkg
julia> Pkg.add("OpticalProperties")
```


## Available materials
The following materials have been implemented:
Au, Al, AlGaAs, cBN, Cu, H2O, Si, doped Si, SiC, SiN, SiO2, Ti, TiN, TiN, TiW, TiW, W, Vacuum

## Quickstart
This is an example of the usage for gold. The functions presented take as arguments the material, here gold and the frequency in rad/s.
```julia
julia> permittivity(Au(),1e15)
-172.92729201753394 + 20.04870903024803im
julia> refractive_index(Au(),1e15)
0.7610244413149028 + 13.17218471696372im
julia> skin_depth(Au(),1e15)
2.2759509105115573e-8
```