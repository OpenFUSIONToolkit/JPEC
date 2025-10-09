# JPEC.jl
module JPEC

include("Splines/Splines.jl")
import .SplinesMod as Spl
export SplinesMod, Spl

include("Equilibrium/Equilibrium.jl")
import .Equilibrium as Equilibrium
export Equilibrium

include("Vacuum/Vacuum.jl")
import .VacuumMod as VacuumMod
export VacuumMod

include("DCON/DCON.jl")
import .DCON as DCON
export DCON

end # module JPEC