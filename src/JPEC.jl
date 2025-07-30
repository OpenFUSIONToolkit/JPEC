# JPEC.jl
module JPEC

include(joinpath(@__DIR__, "..", "deps", "build.jl"))

include("Splines/Splines.jl")
import .SplinesMod as Spl
export SplinesMod, Spl

include("Equilibrium/Equilibrium.jl")
import .Equilibrium as Equilibrium
export Equilibrium

include("DCON/dcon_mod.jl")
import .DconMod as DconMod
export DconMod

end # module JPEC