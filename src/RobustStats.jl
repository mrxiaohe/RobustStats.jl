module RobustStats
using Compat

using DataFrames
using Rmath
using StatsBase
import Base.show
import DataFrames.complete_cases
import DataFrames.deleterows!


export
    outOutput,
    tmean,
    winval,
    winmean,
    winvar,
    trimse,
    trimci,
    stein1,
    stein2,
    idealf,
    pbvar,
    bivar,
    tauloc,
    tauvar,
    outbox,
    akerd,
    akerd_C,
    rdplot,
    sint,
    binomci,
    acbinomci,
    near,
    msmedse,
    stein1_tr,
    stein2_tr,
    sintv2,
    seq,
    cnorm,
    hpsi,
    onestep,
    onesampb,
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
include("data.jl")
end
