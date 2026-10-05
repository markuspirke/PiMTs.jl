"""
    ExGaussian(μ::T, σ::T, c::T) where {T <: Real}

The *exponentially modified Gaussian* (ex-Gaussian) distribution: the distribution of

```math
X = Z + E, \\qquad Z \\sim \\mathrm{Normal}(\\mu, \\sigma), \\quad E \\sim \\mathrm{Exponential}(\\text{rate} = c),
```

with ``Z`` and ``E`` independent. It is the convolution of a normal and an exponential
density.

The density is

```math
f(x; \\mu, \\sigma, c) = \\frac{c}{2}
    \\exp\\!\\left(\\frac{c}{2}\\left(2\\mu + c\\sigma^2 - 2x\\right)\\right)
    \\operatorname{erfc}\\!\\left(\\frac{\\mu + c\\sigma^2 - x}{\\sqrt{2}\\,\\sigma}\\right).
```
# Parameters
- `μ`: mean of the Gaussian component (location; any real number).
- `σ`: standard deviation of the Gaussian component; must satisfy `σ ≥ 0`, though the
  density is only defined for `σ > 0`.
- `c`: **rate** of the exponential component, `c > 0`. The exponential component has mean
  `1/c`, so larger `c` means a shorter right tail. As `c → ∞` the distribution tends to
  `Normal(μ, σ)`.
"""
struct ExGaussian{T <: Real} <: ContinuousUnivariateDistribution
    μ::T
    σ::T
    c::T
    function ExGaussian(μ::T, σ::T, c::T) where {T <: Real}
        @argcheck σ ≥ 0
        @argcheck c > 0
        return new{T}(μ, σ, c)
    end
end
"""
    params(d::ExGaussian)

Return the parameters of `d` as the tuple `(μ, σ, c)`, where `μ` and `σ` are the mean and
standard deviation of the Gaussian component and `c` is the **rate** of the exponential
component.

# Examples
```julia
julia> params(ExGaussian(0.0, 1.0, 2.0))
(0.0, 1.0, 2.0)
```
"""
params(d::ExGaussian) = (d.μ, d.σ, d.c)

"""
    rand(rng::AbstractRNG, d::ExGaussian)

Draw a single variate from `d` by exploiting its definition as a sum: sample
`Normal(μ, σ)` and `Exponential(1/c)` independently and add them.

The array and count forms (`rand(rng, d, n)`, `rand(d, dims...)`) come for free from the
`Distributions.jl` fallbacks for `UnivariateDistribution`.

# Examples
```julia
julia> using Random

julia> rand(ExGaussian(0.0, 1.0, 2.0));

julia> rand(ExGaussian(0.0, 1.0, 2.0), 5);
```
"""
function Random.rand(rng::AbstractRNG, d::ExGaussian)
    μ, σ, c = params(d)
    return rand(rng, Normal(μ, σ)) + rand(rng, Exponential(1/c))
end

"""
    mean(d::ExGaussian)

Mean of `d`, `μ + 1/c`: the mean of the Gaussian component. The exponential part only shifts mass to the right, so the mean
always exceeds the mode and the location parameter `μ`.

# Examples
```julia
julia> mean(ExGaussian(0.0, 1.0, 2.0))
0.5
```
"""
Statistics.mean(d::ExGaussian) = d.μ + 1/d.c # prob wrong 1/τ -> τ

"""
    var(d::ExGaussian)

Variance of `d`, `σ² + 1/c²`. Because the Gaussian and exponential components are
independent, their variances add; the exponential contributes `1/c²`.

# Examples
```julia
julia> var(ExGaussian(0.0, 1.0, 2.0))
1.25
```
"""
Statistics.var(d::ExGaussian) = d.σ^2 + 1/d.c^2
"""
    std(d::ExGaussian)

Standard deviation of `d`, `sqrt(var(d)) = sqrt(σ² + 1/c²)`.
"""
Statistics.std(d::ExGaussian) = sqrt(var(d))


function Distributions.logpdf(d::ExGaussian, x::Real)
    μ, σ, c = params(d)
    z = (μ + c*σ^2 - x) / (sqrt(2) * σ)
    return log(c) + c*(μ - x) + c^2*σ^2/2 + logcdf(Normal(), (x - μ)/σ - c*σ)
end

"""
    pdf(d::ExGaussian, x::Real)

Probability density of `d` at `x`, evaluated in closed form:

```math
f(x) = \\frac{c}{2}
    \\exp\\!\\left(\\frac{c}{2}\\left(2\\mu + c\\sigma^2 - 2x\\right)\\right)
    \\operatorname{erfc}\\!\\left(\\frac{\\mu + c\\sigma^2 - x}{\\sqrt{2}\\,\\sigma}\\right).
```
"""
Distributions.pdf(d::ExGaussian, x::Real) = exp(logpdf(d, x))

"""
    cdf(d::ExGaussian, x::Real)

Cumulative distribution function of `d` at `x`, i.e. ``P(X \\le x)``, computed from the
identity

```math
F(x) = \\Phi\\!\\left(\\frac{x - \\mu}{\\sigma}\\right) - \\frac{f(x)}{c},
```
"""
function Distributions.cdf(d::ExGaussian, x::Real)
    μ, σ, c = params(d)

    return cdf(Normal(μ, σ), x) - pdf(d, x) / c
end


"""
    function ChargeSpectrum(λ::T, q₀::T, σ₀::T, w::T, c₀::T, μ::T, σ::T, kmax::Int) where {T <: Real}

The distribution of a PMT Charge spectrum according to this paper [[https://arxiv.org/abs/1909.05373]].
Nearly full Distributions.jl functionality is added.
The charge spectrum comes from a sum of random variables.
On the one hand due to background processes B, which can be further divided in noise
from the electronics chain Z and due to dark pulses D.
On top of these background processes we have the amplification
 of photoelectrons emitted at the kathode Q.

```math
Y = Z + D + Q
```

For the amplification process, this distribution assumes a normally distributed gain at the
dynodes.
The electronics noise `Z` is Gaussian with mean `q₀` and width `σ₀`.  The dark-pulse
contribution `D` is zero with probability `1 - w` and exponentially distributed with rate
`c₀` with probability `w`, so `Z + D` is a two-component mixture of a Normal and an
exponentially modified Gaussian.  The amplification term `Q` is a compound Poisson sum:
`k ~ Poisson(λ)` photoelectrons each contribute a Gaussian charge with mean `μ` and width
`σ`, giving `μₖ = q₀ + k·μ` and `σₖ = sqrt(σ₀² + k·σ²)` for the `k`-photoelectron peak.
The infinite sum over `k` is truncated at `kmax`.

# Arguments
- `λ`: mean number of photoelectrons emitted at the photocathode per trigger (Poisson
  intensity of the light source together with the quantum efficiency).
- `q₀`: pedestal position, i.e. the mean of the electronics noise in charge units.
- `σ₀`: pedestal width, i.e. the standard deviation of the electronics noise. Must be `≥ 0`.
- `w`: probability that a dark pulse occurs within the integration gate, so the weight of
  the exponential background component. Must lie in `[0, 1]`.
- `c₀`: rate (inverse decay constant) of the exponential dark-pulse charge distribution;
  its mean contribution is `w/c₀`. Must be `> 0`.
- `μ`: mean charge of a single photoelectron, i.e. the gain of the dynode chain. Must be `> 0`.
- `σ`: standard deviation of the single-photoelectron charge. Must be `≥ 0`.
- `kmax`: truncation order of the Poisson sum. Should be chosen well above `λ` so that the
  neglected terms are negligible; `kmax ≥ 1` is enforced.
"""
struct ChargeSpectrum{T<:Real} <: ContinuousUnivariateDistribution
    λ::T
    q₀::T
    σ₀::T
    w::T
    c₀::T
    μ::T
    σ::T
    kmax::Int

    function ChargeSpectrum(λ::T, q₀::T, σ₀::T, w::T, c₀::T, μ::T, σ::T, kmax::Int) where {T <: Real}
        @argcheck σ₀ ≥ 0
        @argcheck c₀ > 0
        @argcheck 1 ≥ w ≥ 0
        @argcheck σ ≥ 0
        @argcheck kmax ≥ 1
        return new{T}(λ, q₀, σ₀, w, c₀, μ, σ, kmax)
    end
end

"""
    params(d::ChargeSpectrum)

Return the parameters of `d` as the tuple `(λ, q₀, σ₀, w, c₀, μ, σ, kmax)`, in the same
order as expected by the constructor.
"""
params(d::ChargeSpectrum) = (d.λ, d.q₀, d.σ₀, d.w, d.c₀, d.μ, d.σ, d.kmax)

"""
    mean(d::ChargeSpectrum)

Mean of the charge spectrum, `q₀ + w/c₀ + λ·μ`: the pedestal position plus the mean
dark-pulse charge plus the mean amplified charge.  This is exact, i.e. independent of the
truncation `kmax`.
"""
Statistics.mean(d::ChargeSpectrum) = d.q₀ + d.w/d.c₀ + d.λ * d.μ

"""
    var(d::ChargeSpectrum)

Variance of the charge spectrum.  The three summands `Z`, `D` and `Q` are independent, so
the variances add: `σ₀²` from the electronics noise, `(w/c₀)²` from the dark pulses and
`λ·(σ² + μ²)` from the compound Poisson amplification.
"""
Statistics.var(d::ChargeSpectrum) = d.σ₀^2 + (d.w/d.c₀)^2 + d.λ * (d.σ^2 + d.μ^2)

"""
    std(d::ChargeSpectrum)

Standard deviation of the charge spectrum, `sqrt(var(d))`.
"""
Statistics.std(d::ChargeSpectrum) = sqrt(var(d))

"""
    minimum(d::ChargeSpectrum)

Theoretical minimal value of the charge spectrum is returned.
"""
Base.minimum(d::ChargeSpectrum) = -Inf

"""
    maximum(d::ChargeSpectrum)

Theoretical maximal value of the charge spectrum is returned.
"""
Base.maximum(d::ChargeSpectrum) = Inf

"""
    insupport(d::ChargeSpectrum, x::Real)

Checks whether the distribution is defined for a given number.
"""
Distributions.insupport(d::ChargeSpectrum, x::Real) = minimum(d) <= x <= maximum(d)

"""
    quantile(d::ChargeSpectrum, x::Real)

Returns the x-quantile of the charge spectrum.
As there is no analyitcal inverse of the cdf, this is estimated
numerically based on `Distributions.quantile_bisect`.
"""
function Distributions.quantile(d::ChargeSpectrum, x::Real)
    qs = rand(d, 1_000_000) # needs improvements in future

    return quantile(qs, x)
end

"""
    pdf(d::ChargeSpectrum, x::Real)

Probability density of the charge spectrum at `x`, evaluated as the Poisson-weighted sum
over the photoelectron peaks

```math
f(x) = \\sum_{k=0}^{k_\\mathrm{max}} \\frac{\\lambda^k e^{-\\lambda}}{k!}
       \\left[(1-w)\\, \\mathcal{N}(x; \\mu_k, \\sigma_k)
            + w\\, \\mathrm{EMG}(x; \\mu_k, \\sigma_k, c_0)\\right]
```

with `μₖ = q₀ + k·μ` and `σₖ = sqrt(σ₀² + k·σ²)`.  Terms whose Poisson weight falls below
`1e-12` are skipped.  Because the sum is truncated at `kmax`, the density does not
integrate to exactly one; the missing mass is the Poisson tail above `kmax`.
"""
function Distributions.pdf(d::ChargeSpectrum, x::Real)
    λ, q₀, σ₀, w, c₀, μ, σ, kmax = params(d)

    y = 0.0

    pois = Poisson(λ)
    for k in 0:kmax
        pk = pdf(pois, k)
        pk < 1e-12 && continue
        μk = q₀ + k * μ
        σk = sqrt(σ₀^2 + k * σ^2)
        Nk = Normal(μk, σk)
        Ek = ExGaussian(μk, σk, c₀)

        y += pk * ((1 - w) * pdf(Nk, x) + w * pdf(Ek, x))
    end

    return y
end

"""
    cdf(d::ChargeSpectrum, x::Real)

Cumulative distribution function of the charge spectrum at `x`.  Same Poisson-weighted
mixture as [`pdf`](@ref), with the component densities replaced by their CDFs, and subject
to the same truncation at `kmax`.
"""
function Distributions.cdf(d::ChargeSpectrum, x::Real)
    λ, q₀, σ₀, w, c₀, μ, σ, kmax = params(d)

    y = 0.0

    pois = Poisson(λ)
    for k in 0:kmax
        pk = pdf(pois, k)
        pk < 1e-12 && continue
        μk = q₀ + k * μ
        σk = sqrt(σ₀^2 + k * σ^2)
        Nk = Normal(μk, σk)
        Ek = ExGaussian(μk, σk, c₀)

        y += pk * ((1 - w) * cdf(Nk, x) + w * cdf(Ek, x))
    end

    return y
end

"""
    rand(rng::AbstractRNG, d::ChargeSpectrum)

Draw a single sample from the charge spectrum.  A photoelectron count `k ~ Poisson(λ)` is
sampled first, which fixes `μₖ = q₀ + k·μ` and `σₖ = sqrt(σ₀² + k·σ²)`.  With probability
`w` a dark pulse is added and the sample is drawn from `ExGaussian(μₖ, σₖ, c₀)`, otherwise
from `Normal(μₖ, σₖ)`.  Unlike [`pdf`](@ref) and [`cdf`](@ref) this is not truncated at
`kmax`.

# Example

```julia
julia> λ, q₀, σ₀, w, c₀, μ, σ, kmax = 3.0, 1.0, 0.2, 0.3, 10.0, 5.0, 2.0, 10
julia> charge_spectrum = ChargeSpectrum(λ, q₀, σ₀, w, c₀, μ, σ, kmax)

julia> rand(charge_spectrum)
julia> rand(charge_spectrum, 100)
```
"""
function Random.rand(rng::AbstractRNG, d::ChargeSpectrum)
    λ, q₀, σ₀, w, c₀, μ, σ, kmax = params(d)
    k = rand(rng, Poisson(λ))
    if rand(Bernoulli(w))
        μk = q₀ + k*μ
        σk = sqrt(σ₀^2 + k*σ^2)
        return rand(rng, ExGaussian(μk, σk, c₀))
    else
        μk = q₀ + k*μ
        σk = sqrt(σ₀^2 + k*σ^2)
        return rand(rng, Normal(μk, σk))
    end

end

"""
    function peak2valley(cs::ChargeSpectrum)

Calculates the peak to valley ratio for a given charge spectrum.
This is accomplished by searching for a minimum between the pedestal peak
and the single photoelectron peak.
The search is based on finding the minimum with Minuit.
"""
function peak2valley(cs::ChargeSpectrum)
    xmin = cs.q₀
    xmax = cs.μ
    x0 = (xmin + xmax)/2

    m = Minuit(x -> pdf(cs, x[1]), [x0]; limits=[(xmin, xmax)])
    migrad!(m)

    return pdf(cs, cs.μ)/pdf(cs, m.values[1])
end
