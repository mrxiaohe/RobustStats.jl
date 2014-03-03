using DataFrames
using Rmath
using Distributions
using GLM
using StatsBase
using Winston
using DataArrays

module RobustStats
using DataFrames
using Rmath
using Winston
using Distributions
using StatsBase
using GLM
using DataArrays
import Base.show
import DataFrames.complete_cases
import DataFrames.deleterows!


export 
    outOutput,
    idealfOutput,
    testOutput,
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
    t1way
    
include("types.jl")
include("utilis.jl")
include("show.jl")
include("functions.jl")
include("data.jl")
end


