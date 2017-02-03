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

    akerd,
    rdplot,

    sint,
    sintv2,
    seq,
    cnorm,
    hpsi,
    onestep,
    bootstrapci,
    mom,
    momci,
    trimpb,
    trimcibt,
    indt,
    pcorb,
    bootse,
    indirectTest,
    yuend,
    t1way,
    bisquareWM,
    huberWM,
    trimean,
    iqrn,
    shorthrange,
    shorthrange!,
    scaleS,
    scaleS!,
    scaleQ,
    scaleQ!,

# Below here in unmaintained.jl
    stein1,
    stein2,
    stein1_tr,
    stein2_tr


include("types.jl")
include("utils.jl")
include("show.jl")
include("functions.jl")
include("location_estimators.jl")
include("dispersion_estimators.jl")
include("unmaintained.jl")
end
