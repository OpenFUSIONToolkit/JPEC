using JPEC

struct ResistType
    e::Float64
    f::Float64
    h::Float64
    m::Float64
    g::Float64
    k::Float64
    eta::Float64
    rho::Float64
    taua::Float64
    taur::Float64
end

struct SingType
    m::Int
    order::Int
    r1::NTuple{1,Int}
    r2::NTuple{2,Int}
    n1::Vector{Int}
    n2::Vector{Int}
    psifac::Float64
    rho::Float64
    q::Float64
    q1::Float64
    di::Float64
    alpha::Complex{Float64}
    power::Union{Nothing, Vector{Complex{Float64}}}
    vmat::Union{Nothing, Array{Complex{Float64},4}}
    mmat::Union{Nothing, Array{Complex{Float64},4}}
    restype::ResistType
end

struct DconFileNames
    # Output file names
    out_bal1_unit::Str
    out_bal2_unit::Str
    bin_bal1_unit::Str
    bin_bal2_unit::Str
    fourfit_out_unit::Str
    fourfit_bin_unit::Str
    evals_out_unit::Str
    evals_bin_unit::Str
    crit_out_unit::Str
    crit_bin_unit::Str
    euler_bin_unit::Str
    init_out_unit::Str
    dcon_unit::Str
    unorm_unit::Str
    ca_unit::Str
    err_unit::Str
end

struct DconInternal
    # Parameter
    sol_base::Int

    # Spline type (replace Any with actual type if available)
    locstab::CubicSplineType

    # Singularity arrays and counters
    msing::Int
    kmsing::Int
    sing::Union{Nothing, Vector{SingType}}
    kinsing::Union{Nothing, Vector{SingType}}

    # Real variables not in namelists
    psilim::Float64
    qlim::Float64
    q1lim::Float64

    # Boundary diagnostics
    size_edge::Int
    pre_edge::Int
    i_edge::Int
    q_edge::Union{Nothing, Vector{Float64}}
    psi_edge::Union{Nothing, Vector{Float64}}
    dw_edge::Union{Nothing, Vector{Complex{Float64}}}
end


struct DconControl
    bal_flag::Bool
    mat_flag::Bool
    ode_flag::Bool
    vac_flag::Bool
    mer_flag::Bool
    res_flag::Bool
    fft_flag::Bool
    node_flag::Bool
    mthvac::Int
    sing_start::Float64
    nn::Int
    delta_mlow::Int
    delta_mhigh::Int
    delta_mband::Int
    thmax0::Float64
    nstep::Int
    ksing::Int
    tol_nr::Float64
    tol_r::Float64
    crossover::Float64
    ucrit::Float64
    singfac_min::Float64
    singfac_max::Float64
    cyl_flag::Bool
    dmlim::Float64
    lim_flag::Bool
    sas_flag::Bool
    sing_order::Int
    sort_type::Int
    termbycross_flag::Bool
    qhigh::Float64
    kin_flag::Bool
    con_flag::Bool
    kinfac1::Float64
    kinfac2::Float64
    kingridtype::Int
    ktanh_flag::Bool
    passing_flag::Bool
    trapped_flag::Bool
    ion_flag::Bool
    electron_flag::Bool
    ktc::Float64
    ktw::Float64
    qlow::Float64
    use_classic_splines::Bool
    reform_eq_with_psilim::Bool
    psiedge::Float64
    nperq_edge::Int
    wv_farwall_flag::Bool
    dcon_kin_threads::Int
    parallel_threads::Int
end

struct DconOutput
    interp::Bool
    crit_break::Bool
    out_bal1::Bool
    bin_bal1::Bool
    out_bal2::Bool
    bin_bal2::Bool
    out_metric::Bool
    bin_metric::Bool
    out_fmat::Bool
    bin_fmat::Bool
    out_gmat::Bool
    bin_gmat::Bool
    out_kmat::Bool
    bin_kmat::Bool
    out_sol::Bool
    out_sol_min::Int
    out_sol_max::Int
    bin_sol::Bool
    bin_sol_min::Int
    bin_sol_max::Int
    out_fl::Bool
    bin_fl::Bool
    out_evals::Bool
    bin_evals::Bool
    bin_euler::Bool
    euler_stride::Int
    bin_vac::Bool
    ahb_flag::Bool
    mthsurf0::Float64
    msol_ahb::Int
    netcdf_out::Bool
    out_fund::Bool
    out_ahg2msc::Bool
end

# Initialization of defaults
# --- Control Namelist ---
dcon_control_defaults = DconControl(
    bal_flag=false,
    mat_flag=false,
    ode_flag=false,
    vac_flag=false,
    mer_flag=false,
    res_flag=false,
    fft_flag=false,
    node_flag=false,
    mthvac=480,
    sing_start=0.0,
    nn=0,
    delta_mlow=0,
    delta_mhigh=0,
    delta_mband=0,
    thmax0=1.0,
    nstep=typemax(Int),
    ksing=-1,
    tol_nr=1e-5,
    tol_r=1e-5,
    crossover=1e-2,
    ucrit=1e4,
    singfac_min=0.0,
    singfac_max=0.0,
    cyl_flag=false,
    dmlim=0.5,
    lim_flag=false,
    sas_flag=false,
    sing_order=0,
    sort_type=0,
    termbycross_flag=false,
    qhigh=1e3,
    kin_flag=false,
    con_flag=false,
    kinfac1=1.0,
    kinfac2=1.0,
    kingridtype=0,
    ktanh_flag=false,
    passing_flag=false,
    trapped_flag=true,
    ion_flag=true,
    electron_flag=false,
    ktc=0.1,
    ktw=50.0,
    qlow=0.0,
    use_classic_splines=false,
    reform_eq_with_psilim=false,
    psiedge=1.0,
    nperq_edge=20,
    wv_farwall_flag=false,
    dcon_kin_threads=0,
    parallel_threads=0
)

# --- Output Namelist ---
dcon_output_defaults = DconOutput(
    interp=false,
    crit_break=true,
    out_bal1=false,
    bin_bal1=false,
    out_bal2=false,
    bin_bal2=false,
    out_metric=false,
    bin_metric=false,
    out_fmat=false,
    bin_fmat=false,
    out_gmat=false,
    bin_gmat=false,
    out_kmat=false,
    bin_kmat=false,
    out_sol=false,
    out_sol_min=0,
    out_sol_max=0,
    bin_sol=false,
    bin_sol_min=0,
    bin_sol_max=0,
    out_fl=false,
    bin_fl=false,
    out_evals=false,
    bin_evals=false,
    bin_euler=false,
    euler_stride=1,
    bin_vac=false,
    ahb_flag=false,
    mthsurf0=1.0,
    msol_ahb=0,
    netcdf_out=true,
    out_fund=false,
    out_ahg2msc=true
)

# --- File Names ---
dcon_file_names = DconFileNames(
    out_bal1_unit="bal1.out",
    out_bal2_unit="bal2.out",
    bin_bal1_unit="bal1.bin",
    bin_bal2_unit="bal2.bin",
    fourfit_out_unit="imats.out",
    fourfit_bin_unit="imats.bin",
    evals_out_unit="feval.out",
    evals_bin_unit="feval.bin",
    crit_out_unit="crit.out",
    crit_bin_unit="crit.bin",
    euler_bin_unit="euler.bin",
    init_out_unit="init.out",
    reinit_out_unit="reinit.out",
    dcon_unit="dcon.out",
    unorm_unit="unorm.bin",
    ca_unit="ca.out",
    err_unit="error.bin"
)
    
# --- Internal Variables ---
dcon_internal_defaults = DconInternal(
    sol_base=50,
    locstab=nothing, # Replace with actual spline type if available
    msing=0,
    kmsing=0,
    sing=nothing,
    kinsing=nothing,
    psilim=0.0,
    qlim=0.0,
    q1lim=0.0,
    size_edge=0,
    pre_edge=1,
    i_edge=1,
    q_edge=nothing,
    psi_edge=nothing,
    dw_edge=nothing
)
