# PMT Charge Spectrumdocs
The photomultiplier charge spectrum implemented in this package is based on this [paper](https://iopscience.iop.org/article/10.1088/1748-0221/15/02/P02001) by Milind Diwan. 
It is assumed that the spectrum can be described as a convolution of three processes. Electronics noise, additional noise (probably mainly from dark pulses, altough I am not sure here) and the amplification process of photoelectrons released due to an outside optical signal.
The electronics noise is modeled as a Gaussian. Additional noise, if present is assumed to be exponentially distributed (not yet understood the deeper assumption here) and the multiplication is approximated by a Gaussian.

The following show how we can initialize a PMT charge spectrum for some given parameters
```@example usage
using PiMTs

λ, q₀, σ₀, w, c₀, μ, σ, kmax = 3.0, 1.0, 0.2, 0.3, 10.0, 5.0, 2.0, 10
cs = ChargeSpectrum(λ, q₀, σ₀, w, c₀, μ, σ, kmax)
```

`ChargeSpectrum` is extending the `Distributions.jl` ecosystem and has full functionality that is expected from a distribution.
```@example usage
rand(cs)
```
By loading one of the Makie plotting backend we can also easily visualize the spectrum out of the box.
```@example usage
using CairoMakie
lines(cs)
```
