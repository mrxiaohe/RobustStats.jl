__precompile__()

module RobustStats

using Compat
using DataFrames
using Rmath
using StatsBase
import Base.show

export
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
    hpsi,
    onestep,
    bootstrapci,
    bootstrapse,
    mom,
    momci,
    contam_randn,
    trimpb,
    pcorb,
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
