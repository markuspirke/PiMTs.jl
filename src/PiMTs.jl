module PiMTs

import Base: minimum, maximum
import Distributions: pdf, cdf, insupport
using Distributions
import Statistics: mean, var, std, quantile
using Statistics
using Random
using ArgCheck
import NativeMinuit: BinnedNLL
using NativeMinuit
using StatsBase

export ExGaussian, params, ChargeSpectrum, mean, var, std, pdf, cdf, insupport, quantile, peak2valley
export BinnedNLL, goodness_of_fit

include("utils.jl")
include("charge_spectrum.jl")
include("fit.jl")

end
