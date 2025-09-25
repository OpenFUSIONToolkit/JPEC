# TODO: add descriptions of what all variables are and/or explanation of defaults
@kwdef mutable struct ResistType
    e::Float64 = 0.0
    f::Float64 = 0.0
    h::Float64 = 0.0
    m::Float64 = 0.0
    g::Float64 = 0.0
    k::Float64 = 0.0
    eta::Float64 = 0.0
    rho::Float64 = 0.0
    taua::Float64 = 0.0
    taur::Float64 = 0.0
end

# TODO: ideally, everything is allocated at construction, but mpert is determined after
# since these are allocated in sing_find. What's the best way to handle this?
# For now, leaving as nothings
# Something simple could be not creating sing_types within sing_find, but instead saving the data,
# then creating them all at once after mpert is determined in dcon.jl
# This wouldn't be as clean, but would allow preallocation. Does this greatly impact performance?
@kwdef mutable struct SingType
    m::Int = 0
    r1::Vector{Int} = [0]
    r2::Vector{Int} = [0, 0]
    n1::Vector{Int} = Int[]
    n2::Vector{Int} = Int[]
    psifac::Float64 = 0.0
    rho::Float64 = 0.0
    q::Float64 = 0.0
    q1::Float64 = 0.0
    di::Float64 = 0.0
    alpha::ComplexF64 = 0.0 + 0.0im
    power::Union{Nothing, Vector{ComplexF64}} = nothing
    vmat::Union{Nothing, Array{ComplexF64,4}} = nothing
    mmat::Union{Nothing, Array{ComplexF64,4}} = nothing
    m0mat::Union{Nothing, Matrix{ComplexF64}} = zeros(ComplexF64, 2, 2)
    restype::ResistType = ResistType()
end
# @kwdef mutable struct SingType
#     mpert::Int
#     order::Int
#     m::Int = 0
#     psifac::Float64 = 0.0
#     rho::Float64 = 0.0
#     q::Float64 = 0.0
#     q1::Float64 = 0.0
#     di::Float64 = 0.0
#     alpha::ComplexF64 = 0.0 + 0.0im
#     r1::Vector{Int} = [0]
#     r2::Vector{Int} = [0, 0]
#     n1::Vector{Int} = Vector{Int}(undef, mpert - 1)
#     n2::Vector{Int} = Vector{Int}(undef, 2 * mpert - 2)
#     power::Vector{ComplexF64} = Vector{ComplexF64}(undef, 2 * mpert)
#     vmat::Array{ComplexF64,4} = Array{ComplexF64}(undef, mpert, 2 * mpert, 2, order + 1)
#     mmat::Array{ComplexF64,4} = Array{ComplexF64}(undef, mpert, 2 * mpert, 2, order + 3)
#     m0mat::Union{Nothing, Matrix{ComplexF64}} = zeros(ComplexF64, 2, 2)
#     restype::ResistType = ResistType()
# end
# # Constructor to allocate matrices
# SingType(mpert::Int, order::Int; kwargs...) = SingType(; mpert, order, kwargs...)

#TODO: I am assuming this will change - needs discussion of what output files should be
# we should have a working base case first
# @kwdef mutable struct DconFileNames
#     out_bal1_unit::String = "bal1.out"
#     out_bal2_unit::String = "bal2.out"
#     bin_bal1_unit::String = "bal1.bin"
#     bin_bal2_unit::String = "bal2.bin"
#     fourfit_out_unit::String = "imats.out"
#     fourfit_bin_unit::String = "imats.bin"
#     evals_out_unit::String = "feval.out"
#     evals_bin_unit::String = "feval.bin"
#     crit_out_unit::String = "crit.out"
#     crit_bin_unit::String = "crit.bin"
#     euler_bin_unit::String = "euler.bin"
#     init_out_unit::String = "init.out"
#     reinit_out_unit::String = "reinit.out"
#     dcon_unit::String = "dcon.out"
#     unorm_unit::String = "unorm.bin"
#     ca_unit::String = "ca.out"
#     err_unit::String = "error.bin"
# end

@kwdef mutable struct DconInternal
    dir_path::String = ""
    mlow::Int = 0 #Copy of delta_mlow, but with limits enforced
    mhigh::Int = 0 #Copy of delta_mhigh, but with limits enforced
    mpert::Int = 0 #mpert = mhigh-mlow+1
    mband::Int = 0 #mband = mpert-1-delta_mband
    vac_memory::Bool = false
    keq_out::Bool = false
    theta_out::Bool = false
    xlmda_out::Bool = false
    fkg_kmats_flag::Bool = false
    sol_base::Int = 50
    locstab::Union{Nothing, Spl.CubicSpline{Float64}} = nothing
    msing::Int = 0
    kmsing::Int = 0
    sing::Union{Nothing, Vector{SingType}} = SingType[]
    kinsing::Union{Nothing, Vector{SingType}} = SingType[]
    psilim::Float64 = 0.0
    qlim::Float64 = 0.0
    q1lim::Float64 = 0.0
    size_edge::Int = 0
    pre_edge::Int = 1
    i_edge::Int = 1
    q_edge::Union{Nothing, Vector{Float64}} = nothing
    psi_edge::Union{Nothing, Vector{Float64}} = nothing
    dw_edge::Union{Nothing, Vector{ComplexF64}} = nothing
end

@kwdef mutable struct DconControl
    verbose::Bool = true
    bal_flag::Bool = false
    mat_flag::Bool = false
    ode_flag::Bool = false
    vac_flag::Bool = false
    mer_flag::Bool = false
    res_flag::Bool = false
    fft_flag::Bool = false
    node_flag::Bool = false
    mthvac::Int = 480
    sing_start::Float64 = 0.0
    nn::Int = 0
    delta_mlow::Int = 0
    delta_mhigh::Int = 0
    delta_mband::Int = 0
    thmax0::Float64 = 1.0
    nstep::Int = typemax(Int)
    ksing::Int = -1
    tol_nr::Float64 = 1e-5
    tol_r::Float64 = 1e-5
    crossover::Float64 = 1e-2
    ucrit::Float64 = 1e4
    singfac_min::Float64 = 0.0
    singfac_max::Float64 = 0.0
    cyl_flag::Bool = false
    dmlim::Float64 = 0.5
    lim_flag::Bool = false
    sas_flag::Bool = false
    sing_order::Int = 2
    sort_type::String = "absm"
    termbycross_flag::Bool = false
    qhigh::Float64 = 1e3
    kin_flag::Bool = false
    con_flag::Bool = false
    kinfac1::Float64 = 1.0
    kinfac2::Float64 = 1.0
    kingridtype::Int = 0
    ktanh_flag::Bool = false
    passing_flag::Bool = false
    trapped_flag::Bool = true
    ion_flag::Bool = true
    electron_flag::Bool = false
    ktc::Float64 = 0.1
    ktw::Float64 = 50.0
    qlow::Float64 = 0.0
    use_classic_splines::Bool = false
    reform_eq_with_psilim::Bool = false
    psiedge::Float64 = 1.0
    nperq_edge::Int = 20
    wv_farwall_flag::Bool = false
    dcon_kin_threads::Int = 1
    parallel_threads::Int = 1
    diagnose::Bool = false
    diagnose_ca::Bool = false
end

# @kwdef mutable struct DconOutput
#     interp::Bool = false
#     crit_break::Bool = true
#     out_bal1::Bool = false
#     bin_bal1::Bool = false
#     out_bal2::Bool = false
#     bin_bal2::Bool = false
#     out_metric::Bool = false
#     bin_metric::Bool = false
#     feval_flag::Bool = false
#     out_fmat::Bool = false
#     bin_fmat::Bool = false
#     out_gmat::Bool = false
#     bin_gmat::Bool = false
#     out_kmat::Bool = false
#     bin_kmat::Bool = false
#     out_sol::Bool = false
#     out_sol_min::Int = 0
#     out_sol_max::Int = 0
#     bin_sol::Bool = false
#     bin_sol_min::Int = 0
#     bin_sol_max::Int = 0
#     out_fl::Bool = false
#     bin_fl::Bool = false
#     out_evals::Bool = false
#     bin_evals::Bool = false
#     bin_euler::Bool = false
#     euler_stride::Int = 1
#     bin_vac::Bool = false
#     ahb_flag::Bool = false
#     mthsurf0::Float64 = 1.0
#     msol_ahb::Int = 0
#     netcdf_out::Bool = true
#     out_fund::Bool = false
#     out_ahg2msc::Bool = true
# end

# TODO: how can we initialize the splines to not be nothings?
@kwdef mutable struct FourFitVars
    mpert::Int

    # Spline matrices
    amats::Union{Nothing,Spl.CubicSpline{ComplexF64}} = nothing
    bmats::Union{Nothing,Spl.CubicSpline{ComplexF64}} = nothing
    cmats::Union{Nothing,Spl.CubicSpline{ComplexF64}} = nothing
    dmats::Union{Nothing,Spl.CubicSpline{ComplexF64}} = nothing
    emats::Union{Nothing,Spl.CubicSpline{ComplexF64}} = nothing
    hmats::Union{Nothing,Spl.CubicSpline{ComplexF64}} = nothing
    dbats::Union{Nothing,Spl.CubicSpline{ComplexF64}} = nothing
    ebats::Union{Nothing,Spl.CubicSpline{ComplexF64}} = nothing
    fbats::Union{Nothing,Spl.CubicSpline{ComplexF64}} = nothing
    fmats::Union{Nothing,Spl.CubicSpline{ComplexF64}} = nothing
    kmats::Union{Nothing,Spl.CubicSpline{ComplexF64}} = nothing
    gmats::Union{Nothing,Spl.CubicSpline{ComplexF64}} = nothing
    # kaats::Union{Nothing,Spl.CubicSpline{ComplexF64}} = nothing
    # gaats::Union{Nothing,Spl.CubicSpline{ComplexF64}} = nothing
    # f0mats::Union{Nothing,Spl.CubicSpline{ComplexF64}} = nothing
    # pmats::Union{Nothing,Spl.CubicSpline{ComplexF64}} = nothing
    # paats::Union{Nothing,Spl.CubicSpline{ComplexF64}} = nothing
    # kkmats::Union{Nothing,Spl.CubicSpline{ComplexF64}} = nothing
    # kkaats::Union{Nothing,Spl.CubicSpline{ComplexF64}} = nothing
    # r1mats::Union{Nothing,Spl.CubicSpline{ComplexF64}} = nothing
    # r2mats::Union{Nothing,Spl.CubicSpline{ComplexF64}} = nothing
    # r3mats::Union{Nothing,Spl.CubicSpline{ComplexF64}} = nothing
    # akmats::Union{Nothing,Spl.CubicSpline{ComplexF64}} = nothing
    # bkmats::Union{Nothing,Spl.CubicSpline{ComplexF64}} = nothing
    # ckmats::Union{Nothing,Spl.CubicSpline{ComplexF64}} = nothing

    # k0s::JPEC.Spl.SplineType

    # ipiva::Union{Nothing,Vector{Int}} = nothing
    asmat::Matrix{ComplexF64} = Matrix{ComplexF64}(undef, mpert, mpert)
    bsmat::Matrix{ComplexF64} = Matrix{ComplexF64}(undef, mpert, mpert)
    csmat::Matrix{ComplexF64} = Matrix{ComplexF64}(undef, mpert, mpert)
    # jmat::Union{Nothing,Vector{ComplexF64}} = nothing

    parallel_threads::Int = 0
    dcon_kin_threads::Int = 0

    # kinetic ABCDEH mats for sing_mod
    # kwmats::Vector{JPEC.Spl.CubicSpline{ComplexF64}} = [JPEC.Spl.CubicSpline{ComplexF64}() for _ in 1:6]
    # ktmats::Vector{JPEC.Spl.CubicSpline{ComplexF64}} = [JPEC.Spl.CubicSpline{ComplexF64}() for _ in 1:6]
end

FourFitVars(mpert::Int) = FourFitVars(; mpert)
