# JPEC.jl
module JPEC

include("Splines/Splines.jl")
import .SplinesMod as Spl
export SplinesMod, Spl

include("Equilibrium/Equilibrium.jl")
import .Equilibrium as Equilibrium
export Equilibrium

include("DCON/dcon_mod.jl")
import .DconMod as DconMod
export DconMod

include("Vacuum/Vacuum.jl")
import .VacuumMod as VacuumMod
export VacuumMod

end # module JPEC