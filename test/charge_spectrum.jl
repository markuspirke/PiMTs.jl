using PiMTs
using Test


@testset "ExGaussian" begin
    μ, σ, c = 0.0, 1.0, 10.0

    emg = PiMTs.ExGaussian(μ, σ, c)
    xs = rand(emg, 100_000)
    @test (μ, σ, c) == params(emg)
    @test 0.1 ≈ mean(emg)
    @test 1.01 ≈ var(emg)
    @test sqrt(1.01) ≈ std(emg)
    @test isreal(rand(emg))
    @test isapprox(mean(emg), mean(xs), rtol=0.1)
    @test isapprox(var(emg), var(xs), rtol=0.1)
    @test isapprox(std(emg), std(xs), rtol=0.1)
    @test 0.3950669410138599 ≈ pdf(emg, 0.0) # values from scipy.stats.exponnorm.pdf
    @test 0.26565308308461183 ≈ pdf(emg, 1.0)
    @test 1.0 ≈ cdf(emg, 1000.0)
    @test 0.8147794377600818 ≈ cdf(emg, 1.0)
    @test 0.9980291814654276 ≈ cdf(emg, 3.0)

end

@testset "ChargeSpectrum" begin
    λ, q₀, σ₀, w, c₀, μ, σ, kmax = 3.0, 1.0, 0.2, 0.3, 10.0, 5.0, 2.0, 10
    cs = ChargeSpectrum(λ, q₀, σ₀, w, c₀, μ, σ, kmax)
    Qs = rand(cs, 100_000)
    A_pedestal = exp(-cs.λ) * 1/sqrt(2π)/cs.σ₀

    @test (λ, q₀, σ₀, w, c₀, μ, σ, kmax) == params(cs)
    @test 1.0 + 3.0*5.0 + 0.3/10.0 ≈ mean(cs)
    @test isapprox(mean(cs), mean(Qs), rtol=0.1)
    @test isapprox(var(cs), var(Qs), rtol=0.1)
    @test isapprox(A_pedestal, pdf(cs, q₀ + w/c₀), rtol=0.1)
    @test isapprox(cdf(cs, 1e6), 1.0, rtol=1e-3)
    @test isapprox(cdf(cs, -1e2), 0.0, atol=1e-10)
    @test isreal(rand(cs))
    @test -Inf == minimum(cs)
    @test Inf == maximum(cs)
    @test insupport(cs, 10.0)
    @test insupport(cs, -10.0)
    @test isapprox(quantile(Qs, 0.5), quantile(cs, 0.5), rtol=0.1)

end
