"""
    BinnedNLL(counts, bin_edges, d::ChargeSpectrum)

Using `BinnedNLL` together with a `ChargeSpectrum` in the background a multinomial
likelihood ratio is constructed based on the observed counts, choosen bin size and
the CDF of the charge spectrum.

```math
\\text{BinnedNLL} = - \\log \\lambda = \\sum_i n_i \\log{n_i/y_i}
```

``y_i`` are the individual counts in each bin.

# Arguments
- `counts` binned integrated charges
- `bin_edges` as the variable name suggests bin edges
- `d` a given ChargeSpectrum distribution with selected `k_max`
"""
function NativeMinuit.BinnedNLL(counts, bin_edges, d::ChargeSpectrum)
    return BinnedNLL(counts, bin_edges, _cdf_for_minuit(d))
end

function _cdf_for_minuit(d::ChargeSpectrum)
    λ, q₀, σ₀, w, c₀, μ, σ, kmax = params(d)
    p0 = λ, q₀, σ₀, w, c₀, μ, σ

    f(x, p, kmax) = cdf(ChargeSpectrum(p..., kmax), x)

    return (x, p) -> f(x, p, kmax)
end


"""
    function goodness_of_fit(bll::BinnedNLL, m::Minuit)

If `BinnedNLL` is used as a cost function for Minuit the likelihood is constructed based on
a multinomial distribution ``L_m``.

```math
L_m(y; n) = N! \\prod_i p_i^{n^i} / n_i! = N! N^N \\prod_i y_i^{n_i} / n_i!
```

In this setting Minuit does not directly minimize the negative log likelihood, but a log likelihood ratio ``-2\\log\\lambda`` with

```math
\\lambda = L_m(y; n) / L_m(n; n)
```

This is equivalent to maximizing the likelihood, as ``-2\\log\\lambda`` splits into a sum of two parts,
where the second part is not dependent on ``y``.

```math
\\chi^2_{\\lambda} = -2\\log\\lambda = -2\\log L_m(y; n) - 2\\log L_m(n; n)
```

As the variable definition ``\\chi^2_{\\lambda}`` hints, it will be, in the asymptotic limit,
 distributed according to a chi-square distribution, which can also be used as a test statistic.

This is what is used here. The binned likelihood function (i.e. the ratio) is evaluated at the
maximum likelihood estimate. Then this is divided by the number of free parameters.

More details can be found in `Clarification of the use of CHI-square and likelihood functions in fits to histograms` by Baker and Cousins.
"""
function goodness_of_fit(bll::BinnedNLL, m::Minuit)
    χ² = 2 * bll(m.values)
    ndf = length(bll.n) - n_free(m.params) - 1

    return χ² / ndf
end

"""
    save_fit_results(fname::String, bll::BinnedNLL, m::Minuit)

Saves parameters names, values and errors into a csv file.
Additionally the goodness of fit value is add as a last line with
an error set to 0.0
"""
function save_fit_results(fname::String, bll::BinnedNLL, m::Minuit, kmax::Int)
    gof_name = "gof"
    gof = goodness_of_fit(bll, m)
    gof_err = 0.0

    cs = ChargeSpectrum(m.values..., kmax)
    p2v_name = "p2v"
    p2v = peak2valley(cs)
    p2v_err = 0.0



    open(fname, "w") do io
        join(io, ["parameter", "value", "err"], ",")
        println(io)

        for (name, value, err) in zip(m.parameters, m.values, m.errors)
            err = round(err, sigdigits=2)
            err_str = "$(err)"
            n_digits = length.(split(err_str, "."))[end]
            value = round(value, digits=n_digits)
            join(io, [name, value, err], ",")
            println(io)
        end

        join(io, [gof_name, round(gof, digits=2), gof_err], ",")
        println(io)
        join(io, [p2v_name, round(p2v, digits=2), p2v_err], ",")
        println(io)
    end
end
