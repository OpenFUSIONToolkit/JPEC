# Vacuum_data.jl
# This file defines the data structures for the Vacuum module in Julia.



# Vacuum data structures. 
# These are variables that are used globally in the 
# fortran code. 

struct VacuumInputType
    r::Vector{Float64}  # Plasma boundary r-coordinate as a function of poloidal angle
    z::Vector{Float64}  # Plasma boundary z-coordinate as a function of poloidal angle
    delta::Vector{Float64}  # -dphi/qa; 0 for coordinate systems using machine angle, e.g. PEST
    mlow::Int # Lower poloidal mode number
    mhigh::Int # Upper poloidal mode number
    n::Int # Toroidal mode number
    qa::Float64 # Safety factor at the plasma boundary
    mtheta_in::Int # Number of poloidal angles in the input
end

struct VacuumGlobals

end

# Vacuum settings struct

# Each namelist group is represented as a separate struct.
# The main VacuumSettingsType struct contains them as fields.
# This corresponds to the Fortran input namelist structure.

struct Modes
    mth::Int # Even integer, Number of poloidal grid points used for the calculation. The values of the needed quantities on these points are interpolated from those on the plasma boundary from DCON.
    xiin::NTuple{9, Int}
    lsymz::Bool # Whether to symmetrize the vacuum energy matrix
    leqarcw::Int # Turns on equal arcs distribution of the nodes on the shell. Best results unless the wall is very close to the plasma. See ishape=6 option.
end

struct Vacdat
    ishape::Int # 
    aw::Float64
    bw::Float64
    cw::Float64
    dw::Float64
    tw::Float64
    nsing::Int
    epsq::Float64
    noutv::Int
    idgt::Int
    idot::Int
    idsk::Int
    delg::Float64
    delfac::Float64
    cn0::Int
end

struct Shape
    ipshp::Int
    xpl::Float64
    apl::Float64
    a::Float64
    b::Float64
    bpl::Float64
    dpl::Float64
    r::Float64
    abulg::Float64
    bbulg::Float64
    tbulg::Float64
    qain::Float64
end

struct Diagns
    lkdis::Bool
    ieig::Int
    iloop::Int
    lpsub::Int
    nloop::Int
    nloopr::Int
    nphil::Int
    nphse::Int
    xofsl::Float64
    ntloop::Int
    aloop::Float64
    bloop::Float64
    dloop::Float64
    rloop::Float64
    deloop::Float64
    mx::Int
    mz::Int
    nph::Int
    nxlpin::Int
    nzlpin::Int
    epslp::Float64
    xlpmin::Float64
    xlpmax::Float64
    zlpmin::Float64
    zlpmax::Float64
    linterior::Int
end

struct Sprk
    nminus::Int
    nplus::Int
    mphi::Int
    lwrt11::Int
    civ::Float64
    sp2sgn1::Int
    sp2sgn2::Int
    sp2sgn3::Int
    sp2sgn4::Int
    sp2sgn5::Int
    sp3sgn1::Int
    sp3sgn2::Int
    sp3sgn3::Int
    sp3sgn4::Int
    sp3sgn5::Int
    lff::Int
    ff::Float64
    fv::Vector{Float64}
end

struct VacuumSettingsType
    modes::Modes
    vacdat::Vacdat
    shape::Shape
    diagns::Diagns
    sprk::Sprk
end

# Example instantiation of the DIIID ideal example input namelist
# VacuumSettings = VacuumSettingsType(
#     Modes(480, (0,0,0,0,0,0,0,1,0), true, 1, 0, 0, 0),
#     Debugs(false, false, false, false, false, false, 0, false),
#     Vacdat(6, 0.05, 1.5, 0.0, 0.5, 0.05, 500, 1e-5, 37, 6, 0, 0, 15.01, 0.001, 1),
#     Shape(0, 100.0, 1.0, 20.0, 170.0, 1.0, 0.0, 1.0, 0.932, 17.0, 0.02, 2.5),
#     Diagns(false, 0, 0, 1, 128, 0, 3, 1, 0.0, 32, 0.01, 1.6, 0.5, 1.0, 0.001, 21, 21, 0, 6, 11, 0.02, 0.7, 2.7, -1.5, 1.5, 2),
#     Sprk(0, 0, 16, 0, 0.0, 1, 1, 1, 1, 1, -1, -1, -1, 1, 1, 0, 1.6,
#         [1.6, 1.6, 1.6, 1.6, 1.6, 1.0, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6])
# )


