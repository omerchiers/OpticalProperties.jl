# OpticalProperties

This julia package computes permittivities and refractive indices for usual materials.
This is very much work in progress.

## To do
- Introduce struct for Au, Sic, Al
- Making the names consistent with their chemical formula

## Warning
The API is not very consistent yet. Especially the way materials are exported.
Some are exported as singleton types (Al(),Au(),Sic()), others as constants (Cu).

## Installation
You need Julia version 0.6.4 or higher.

To install on julia 0.7 or higher, enter the Pkg mode by typing ]
and then when installing from gitlab

```julia
(v1.0)> https://git-cethil.insa-lyon.fr/omerchiers/myphysicalconstants.jl
(v1.0)> https://git-cethil.insa-lyon.fr/omerchiers/opticalproperties.jl
```
or if you use the github repository

```julia
(v1.0)> add https://github.com/omerchiers/MyPhysicalProperties.jl
(v1.0)> add https://github.com/omerchiers/OpticalProperties.jl
```


## Examples
This is an example of the usage for gold. The functions presented take as arguments the material, here gold and the frequency in rad/s.
```julia
julia> permittivity(Au(),1e15)
-172.92729201753394 + 20.04870903024803im
julia> refractive_index(Au(),1e15)
0.7610244413149028 + 13.17218471696372im
julia> skin_depth(Au(),1e15)
2.2759509105115573e-8
```
