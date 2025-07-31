#TODO: add descriptions of what all variables are and/or explanation of defaults
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

@kwdef mutable struct SingType
    m::Int = 0
    order::Int = 0
    r1::NTuple{1,Int} = (0,)
    r2::NTuple{2,Int} = (0,0)
    n1::Vector{Int} = Int[]
    n2::Vector{Int} = Int[]
    psifac::Float64 = 0.0
    rho::Float64 = 0.0
    q::Float64 = 0.0
    q1::Float64 = 0.0
    di::Float64 = 0.0
    alpha::Complex64 = 0.0 + 0.0im
    power::Union{Nothing, Vector{Complex64}} = nothing
    vmat::Union{Nothing, Array{Complex64,4}} = nothing
    mmat::Union{Nothing, Array{Complex64,4}} = nothing
    restype::ResistType = ResistType()
end

#TODO: I am assuming this will change - needs discussion of what output files should be
# we should have a working base case first
@kwdef mutable struct DconFileNames
    out_bal1_unit::String = "bal1.out"
    out_bal2_unit::String = "bal2.out"
    bin_bal1_unit::String = "bal1.bin"
    bin_bal2_unit::String = "bal2.bin"
    fourfit_out_unit::String = "imats.out"
    fourfit_bin_unit::String = "imats.bin"
    evals_out_unit::String = "feval.out"
    evals_bin_unit::String = "feval.bin"
    crit_out_unit::String = "crit.out"
    crit_bin_unit::String = "crit.bin"
    euler_bin_unit::String = "euler.bin"
    init_out_unit::String = "init.out"
    reinit_out_unit::String = "reinit.out"
    dcon_unit::String = "dcon.out"
    unorm_unit::String = "unorm.bin"
    ca_unit::String = "ca.out"
    err_unit::String = "error.bin"
end

@kwdef mutable struct DconInternal
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
    locstab::Any = nothing # Replace Any with actual type if available
    msing::Int = 0
    kmsing::Int = 0
    sing::Union{Nothing, Vector{SingType}} = SingType[]
    kinsing::Union{Nothing, Vector{SingType}} = nothing
    psilim::Float64 = 0.0
    qlim::Float64 = 0.0
    q1lim::Float64 = 0.0
    size_edge::Int = 0
    pre_edge::Int = 1
    i_edge::Int = 1
    q_edge::Union{Nothing, Vector{Float64}} = nothing
    psi_edge::Union{Nothing, Vector{Float64}} = nothing
    dw_edge::Union{Nothing, Vector{Complex64}} = nothing
    ud::Union{Nothing, Array{Float64,3} } = nothing #storage for u-derivative and xss # JMH - is there overlap with OdeState here? 
    # JMH - changed previous line from Matrix to Array due to error - might have to modify depending on what ud looks like
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
    sing_order::Int = 0
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

@kwdef mutable struct DconOutput
    interp::Bool = false
    crit_break::Bool = true
    out_bal1::Bool = false
    bin_bal1::Bool = false
    out_bal2::Bool = false
    bin_bal2::Bool = false
    out_metric::Bool = false
    bin_metric::Bool = false
    feval_flag::Bool = false
    out_fmat::Bool = false
    bin_fmat::Bool = false
    out_gmat::Bool = false
    bin_gmat::Bool = false
    out_kmat::Bool = false
    bin_kmat::Bool = false
    out_sol::Bool = false
    out_sol_min::Int = 0
    out_sol_max::Int = 0
    bin_sol::Bool = false
    bin_sol_min::Int = 0
    bin_sol_max::Int = 0
    out_fl::Bool = false
    bin_fl::Bool = false
    out_evals::Bool = false
    bin_evals::Bool = false
    bin_euler::Bool = false
    euler_stride::Int = 1
    bin_vac::Bool = false
    ahb_flag::Bool = false
    mthsurf0::Float64 = 1.0
    msol_ahb::Int = 0
    netcdf_out::Bool = true
    out_fund::Bool = false
    out_ahg2msc::Bool = true
end

# JMH - what is a CsplineType? Can't find it anywhere. Commenting out so it will compile
# @kwdef struct FourFitVars

#     # Spline matrices
#     amats::CsplineType
#     bmats::CsplineType
#     cmats::CsplineType
#     dmats::CsplineType
#     emats::CsplineType
#     hmats::CsplineType
#     dbats::CsplineType
#     ebats::CsplineType
#     fbats::CsplineType
#     fmats::CsplineType
#     kmats::CsplineType
#     gmats::CsplineType
#     kaats::CsplineType
#     gaats::CsplineType
#     f0mats::CsplineType
#     pmats::CsplineType
#     paats::CsplineType
#     kkmats::CsplineType
#     kkaats::CsplineType
#     r1mats::CsplineType
#     r2mats::CsplineType
#     r3mats::CsplineType
#     akmats::CsplineType
#     bkmats::CsplineType
#     ckmats::CsplineType

#     k0s::SplineType

#     ipiva::Vector{Int}
#     asmat::Array{ComplexF64,2}
#     bsmat::Array{ComplexF64,2}
#     csmat::Array{ComplexF64,2}
#     jmat::Vector{ComplexF64}

#     parallel_threads::Int = 0
#     dcon_kin_threads::Int = 0

#     # kinetic ABCDEH mats for sing_mod
#     kwmats::Vector{CsplineType} = [CsplineType() for _ in 1:6]
#     ktmats::Vector{CsplineType} = [CsplineType() for _ in 1:6]
# end

@kwdef mutable struct SingVars

      msol::Int = 0
      sing_order::Int = 2
      det_max::Float64 = 0
      r1::Vector{Int} = Int[]
      r2::Vector{Int} = Int[]
      n1::Vector{Int} = Int[]
      n2::Vector{Int} = Int[]

    #   m0mat::Union{Nothing, Matrix{ComplexF64,2}} = nothing
    #   sing_detf::Union{Nothing, Matrix{ComplexF64,2}} = nothing
      m0mat::Union{Nothing, Array{ComplexF64,2}} = nothing # JMH - same as DconInternal - changing to array so it will compile, might be wrong
      sing_detf::Union{Nothing, Array{ComplexF64,2}} = nothing
end