# OpticalProperties

This julia package computes permittivities and refractive indices for usual materials.
This is very much work in progress.

## To do
- Introduce struct for Au, Sic, Al
- Making the names consistent with their chemical formula

## Warning
The API is not very consistent yet. Especially the way materials are exported.
Some are exported as singleton types (Al(),Au(),Sic()), others as constants (Cu).


## Examples

```julia
julia> permittivity(Au(),1e15)
-172.92729201753394 + 20.04870903024803im
```
