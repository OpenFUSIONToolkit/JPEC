module DconMod

# All imports and includes for the DCON module
using LinearAlgebra
using LinearAlgebra.LAPACK #required for banded matrix operations
using BandedMatrices #for banded matrix operations
using TOML
using FFTW
using DifferentialEquations
import ..Equilibrium
import ..Spl
import ..SplinesMod
using Printf


# Include all necessary files
include("dcon_structs.jl")
include("dcon.jl")
include("Ode.jl")
include("sing.jl")
include("fourfit_wip.jl")
include("ode_output.jl")
include("utils.jl")

const Version = "JULIA-PORT-1.0"


end