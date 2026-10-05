# Fitting
The following shows how we can estimate parameters of the charge spectrum, by fitting the spectrum to a histogram.
Here we use the standard tool of high energy physicists when performing fits: Minuit (`NativeMinuit.jl`).
## Fitting on simulation data
### Sampling data
First we need to loading packages and define a charge spectrum where we sample our data from.

```@example usage
using PMTools, NativeMinuit, StatsBase

λ, q₀, σ₀, w, c₀, μ, σ, kmax = 0.8, 1.0, 0.2, 0.3, 1.0, 6.0, 2.0, 10
cs = ChargeSpectrum(λ, q₀, σ₀, w, c₀, μ, σ, kmax)
```

Now we sample from this distribution
```@example usage
Qs = rand(cs, 10_000)
```

In an experiment we usually can infer the average noise and standard deviation from regions in the time series where not pulse was present. We simulate this by sampling from a charge spectrum with $\lambda$ set to 0.0.
```@example usage
cs = ChargeSpectrum(0.0, q₀, σ₀, w, c₀, μ, σ, kmax)
Qs_noise = rand(cs, 10_000)
```
The sample mean and standard deviation can directly be used as starting parameters. 
Based on thes a rough guess on the average number of photoelectrons is possible.
The expectation value of the charge spectrum is given by the expectation value of the electronics noise, noise from spurious pulses and by the expected number of PE multiplied with the gain, i.e. `X = q₀ + w/c₀ + λ·μ`.
If we like to set a starting value for `μ = (X - q₀ + w/c₀)/λ`, we need some assumptions on `w, c₀` and `λ`. `λ` depends mainly on the intensity of the illumination source (and the detection efficiency, but lets not worry about this "second order" effect).
Frequently the intensity of the light source is fixed and thus a rough starting value on `λ` can be given.
By additionally setting w=0 in the expectation value, we get a starting parameter for the gain.
Once we have `μ`, the standard deviation of the gain can often be assumed to be 1/3 of the average.
`w, c₀` could in part be estimated by dark rate measurements without illumination, but here often typical values just work out of the box.
We set `w` to 0.5 and `c₀` to `1/μ` as one often expects dark pulses similar to single PE pulses.
```@example usage
q₀_init, σ₀_init = mean(Qs_noise), std(Qs_noise)
λ_init = 1.0
μ_init = (mean(Qs) -  q₀_init)/λ_init
σ_init = 0.3 * μ_init
w_init = 0.5
c_init = 1/μ_init

p0 = [λ_init, q₀_init, σ₀_init, w_init, 1/μ_init, μ_init, σ_init]
```

### Binning and likelihood definition
Now we need to bin our data and construct the likelihood function.
```@example usage
N_bins = 200
bin_edges = range(extrema(Qs)..., length=N_bins+1)
h = fit(Histogram, Qs, bin_edges)
counts = h.weights

bll = BinnedNLL(counts, bin_edges, cs)
```

### Minuit fit
What follows is standard `NativeMinuit.jl` code
```@example usage
m = Minuit(bll, p0,
        names = ["λ", "q₀", "σ₀", "w", "c₀", "μ", "σ"],
        limits = [(0.0, 10.0), (nothing, nothing), (0.0, nothing),
                    (0.0, 1.0), (0.0, nothing), (0.0, nothing),
                    (0.0, 10.0)],
        )

migrad!(m)
```

We can also calculate the goodness of fit
```@example usage
goodness_of_fit(bll, m)
```

### Visualization
We can visualize the data and the fit.
```@example usage
using CairoMakie
fig = Figure()
ax = Axis(fig[1,1], xlabel="Charge / arb")
stephist!(ax, Qs, bins=N_bins, color=:black, normalization=:pdf)
p_fit = m.values
lines!(ax, ChargeSpectrum(p_fit..., kmax), color=:red)
fig
```

## Fitting on real data
In the example above data was created based on samples taken from the model that was assumed in the fit.
To prove that this also works on real data, we have a helper function which provides real data from a Hamamatsu 10 inch PMT.
```@example usage
using PMTools, NativeMinuit, StatsBase, CairoMakie
Qs, Qs_noise = PMTools.get_test_data()

q₀_init, σ₀_init = mean(Qs_noise), std(Qs_noise)
λ_init = 1.0
μ_init = (mean(Qs) -  q₀_init)/λ_init
σ_init = 0.3 * μ_init
w_init = 0.5
c_init = 1/μ_init



N_bins = 200
bin_edges = range(extrema(Qs)..., length=N_bins+1)
h = fit(Histogram, Qs, bin_edges)
ys = h.weights

kmax = 10
p0 = [λ_init, q₀_init, σ₀_init, w_init, c_init, μ_init, σ_init]
cs = ChargeSpectrum(p0..., kmax)

bll = BinnedNLL(ys, bin_edges, cs)


m = Minuit(bll, p0,
        names = ["λ", "q₀", "σ₀", "w", "c₀", "μ", "σ"],
        limits = [(0.0, 10.0), (nothing, nothing), (0.0, nothing),
                    (0.0, 1.0), (0.0, nothing), (0.0, nothing),
                    (0.0, 10.0)],
        )

migrad!(m)
```

And we can again visulize the result
```@example usage
using CairoMakie
fig = Figure()
ax = Axis(fig[1,1], xlabel="Charge / arb")
stephist!(ax, Qs, bins=N_bins, color=:black, normalization=:pdf)
p_fit = m.values
lines!(ax, ChargeSpectrum(p_fit..., kmax), color=:red)
fig
```
