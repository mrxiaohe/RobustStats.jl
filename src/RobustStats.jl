__precompile__()

module RobustStats
using Compat

using DataFrames
using Rmath
using StatsBase
using PyPlot
import Base.show
import DataFrames.complete_cases
import DataFrames.deleterows!


export
    outOutput,
    tmean,
    winval,
    winmean,
    winvar,
    wincov,
    trimse,
    trimci,
    idealf,
    pbvar,
    bivar,
    tauloc,
    tauvar,
    outbox,
    msmedse,
    binomci,
    acbinomci,
    sint,
    sintv2,
    seq,
    hpsi,
    onestep,
    bootstrapci,
    mom,
    momci,
    contam_randn,
    trimpb,
    pcorb,

    trimcibt,
    bootse,
    yuend,

    bisquareWM,
    huberWM,
    trimean,
    iqrn,
    shorthrange,
    shorthrange!,
    scaleS,
    scaleS!,
    scaleQ,
    scaleQ!


include("types.jl")
include("utils.jl")
include("show.jl")
include("functions.jl")
include("location_estimators.jl")
include("dispersion_estimators.jl")
end
