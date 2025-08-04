# JPEC.jl
module JPEC

include(joinpath(@__DIR__, "..", "deps", "build.jl"))

include("Splines/Splines.jl")
import .SplinesMod as Spl
export SplinesMod, Spl

include("Equilibrium/Equilibrium.jl")
export Equilibrium

end # module JPEC