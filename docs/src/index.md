# PiMTs.jl

This is the documentation of the awesome `PiMTs.jl` Julia package.

!!! note

    This package is still under development!

## Installation
`PiMTs.jl` is **not an officially registered Julia package** but it's available on
the **[ECAP Julia registry](https://git.ecap.work/common/julia-registry)**. To add
the ECAP Julia registry to your local Julia registry list, follow the
instructions in its
[README](https://git.ecap.work/common/julia-registry#adding-the-registry) or simply do

    git clone https://git.ecap.work/common/julia-registry ~/.julia/registries/ECAP
    
After that, you can install `PiMTs.jl` just like any other Julia package:

    julia> import Pkg; Pkg.add("PiMTs")
    
## Quickstart

``` julia-repl
julia> using PiMTs
julia> λ, q₀, σ₀, w, c₀, μ, σ, kmax = 0.8, 1.0, 0.2, 0.3, 1.0, 6.0, 2.0, 10
julia> cs = ChargeSpectrum(λ, q₀, σ₀, w, c₀, μ, σ, kmax)
```
