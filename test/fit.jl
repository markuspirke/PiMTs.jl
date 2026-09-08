using PMTools
using NativeMinuit
using StatsBase
using Random
using Test



@testset "BLL Fit" begin
    # simulated data
    λ, q₀, σ₀, w, c₀, μ, σ, kmax = 0.8, 1.0, 0.2, 0.3, 10.0, 5.0, 2.0, 10
    cs = ChargeSpectrum(λ, q₀, σ₀, w, c₀, μ, σ, kmax)
    cs_0 = ChargeSpectrum(0.0, q₀, σ₀, w, c₀, μ, σ, kmax) # this we have in our data


    Random.seed!(123)
    Qs = rand(cs, 10_000)
    Qs_0 = rand(cs_0, 10_000)
    q₀_init, σ₀_init = mean(Qs_0), std(Qs_0)
    λ_init = 1.0
    μ_init = (mean(Qs) -  q₀_init)/λ_init
    σ_init = 0.3 * μ_init

    N_bins = 200
    bin_edges = range(extrema(Qs)..., length=N_bins+1)
    h = fit(Histogram, Qs, bin_edges)
    counts = h.weights

    bll = BinnedNLL(counts, bin_edges, cs)

    p0 = [λ_init, q₀_init, σ₀_init, w, c₀, μ_init, σ_init]

    m = Minuit(bll, p0,
            names = ["λ", "q₀", "σ₀", "w", "c₀", "μ", "σ"],
            limits = [(0.0, 10.0), (nothing, nothing), (0.0, nothing),
                        (0.0, 1.0), (0.0, nothing), (0.0, nothing),
                        (0.0, 10.0)],
            )

    migrad!(m)
    p_fit = m.values

    @test isapprox(0.799171968302097, p_fit[1], rtol=0.01)
    @test isapprox(0.9281309744858782, p_fit[2], rtol=0.1)
    @test 7 == n_free(m.params)

    @test isapprox(0.9702282241081267, goodness_of_fit(bll, m), rtol=0.1)

    # real data

    Qs, Qs_noise = PMTools.get_test_data()

    q₀_init, σ₀_init = mean(Qs_noise), std(Qs_noise)
    λ_init = 1.0
    μ_init = (mean(Qs) -  q₀_init)/λ_init
    σ_init = 0.3 * μ_init

    N_bins = 200
    bin_edges = range(extrema(Qs)..., length=N_bins+1)
    h = fit(Histogram, Qs, bin_edges)
    ys = h.weights

    kmax = 10
    p0 = [λ_init, q₀_init, σ₀_init, 0.5, μ_init, μ_init, σ_init]
    cs = ChargeSpectrum(p0..., kmax)

    bll = BinnedNLL(ys, bin_edges, cs)


    m = Minuit(bll, p0,
            names = ["λ", "q₀", "σ₀", "w", "c₀", "μ", "σ"],
            limits = [(0.0, 10.0), (nothing, nothing), (0.0, nothing),
                        (0.0, 1.0), (0.0, nothing), (0.0, nothing),
                        (0.0, 10.0)],
            )

    migrad!(m)
    p_fit = m.values

    @test isapprox(0.7240891558747176, p_fit[1], rtol=0.01)
    @test isapprox(-0.9820798224945543, p_fit[2], rtol=0.1)
    @test 7 == n_free(m.params)

    @test isapprox(0.727896333582715, goodness_of_fit(bll, m), rtol=0.1)



    PMTools.save_fit_results("fit_results.csv", bll, m)
    vals, errs = Float64[], Float64[]
    for line in readlines("fit_results.csv")[2:end]
        param, val, err = split(line, ",")
        val, err = parse(Float64, val), parse(Float64, err)
        push!(vals, val)
        push!(errs, err)
    end

    @test vals[1] == 0.724
    @test errs[1] == 0.049
    rm("fit_results.csv")
end
