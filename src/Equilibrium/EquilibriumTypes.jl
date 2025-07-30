#=
This file is the one stop shop for all the fundemental structures used in 
    creating equilibrium descriptions for the DCON ODE to use.
=#

using Base: @kwdef

# --- Helper function --- 


function symbolize_keys(dict::Dict{String, Any})
    return Dict(Symbol(k) => v for (k, v) in dict)
end


# --- Main Structures for the Equilibrium Code ---


@kwdef mutable struct EquilControl
    eq_type::String = "efit"
    eq_filename::String = "mypath"

    jac_type::String = "hamada"
    power_bp::Int = 0
    power_b::Int = 0
    power_r::Int = 0

    grid_type::String = "ldp"
    psilow::Float64 = 1e-2
    psihigh::Float64 = 0.994
    mpsi::Int = 128
    mtheta::Int = 256

    newq0::Int = 0
    etol::Float64 = 1e-7
    use_classic_splines::Bool = false

    input_only::Bool = false
    use_galgrid::Bool = true

    """
    Modified internal constructor that enforces self consistency within the inputs
    """
    function EquilControl(eq_type, eq_filename, jac_type, power_bp, power_b, power_r,
        grid_type, psilow, psihigh, mpsi, mtheta, newq0, etol, use_classic_splines,
        input_only,use_galgrid)
        if jac_type == "hamada"
            @info "Forcing hamada coordinate jacobian exponents: power_*"
            power_b=0
            power_bp=0
            power_r=0
        elseif jac_type == "pest"
            @info "Forcing pest coordinate jacobian exponents: power_*"
            power_b=0
            power_bp=0
            power_r=2
        elseif jac_type == "equal_arc"
            @info "Forcing equal_arc coordinate jacobian exponents: power_*"
            power_b=0
            power_bp=1
            power_r=0
        elseif jac_type == "boozer"
            @info "Forcing boozer coordinate jacobian exponents: power_*"
            power_b=2
            power_bp=0
            power_r=0
        elseif jac_type == "park"
            @info "Forcing park coordinate jacobian exponents: power_*"
            power_b=1
            power_bp=0
            power_r=0
        elseif jac_type == "other"
            @info "Using manual jacobian exponents: power b, bp, r = $(power_b), $(power_bp), $(power_r)"
        elseif jac_type != "other"
            error("Cannot recognize jac_type = $(jac_type)")
        end
        return new(eq_type, eq_filename, jac_type, power_bp, power_b, power_r,
        grid_type, psilow, psihigh, mpsi, mtheta, newq0, etol, use_classic_splines,
        input_only,use_galgrid)
    end
end

@kwdef mutable struct EquilOutput
    gse_flag::Bool = false
    out_eq_1d::Bool = false
    bin_eq_1d::Bool = false
    out_eq_2d::Bool = false
    bin_eq_2d::Bool = true
    out_2d::Bool = false
    bin_2d::Bool = false
    dump_flag::Bool = false
end


"""
    EquilConfig(...)

A container struct that bundles all necessary configuration settings originally specified in the equil
    fortran namelsits.
"""
@kwdef mutable struct EquilConfig
    control::EquilControl = EquilControl()
    output::EquilOutput = EquilOutput()
end
#

"""
Constructor that allows users to form a EquilConfig struct from dictionaries
    for convinience when most of the defaults are fine.
"""
function EquilConfig(control::Dict, output::Dict)
    construct = EquilControl(;control...)
    outstruct = EquilOutput(;output...)
    return EquilConfig(control=construct, output=outstruct)
end

"""
Outer constructor for EquilConfig that enables a toml file 
    interface for specifying the configuration settings
"""
# if this also have default, then conflicts with @kwdef mutable struct EquilConfig.
function EquilConfig(path::String)
    raw = TOML.parsefile(path)

    # Extract EQUIL_CONTROL with default fallback
    control_data = get(raw, "EQUIL_CONTROL", Dict())
    output_data  = get(raw, "EQUIL_OUTPUT", Dict())

    # Check for required fields in control_data
    required_keys = ("eq_filename", "eq_type")
    missingkeys = filter(k -> !haskey(control_data, k), required_keys)

    if !isempty(missingkeys)
        error("Missing required key(s) in [EQUIL_CONTROL]: $(join(missing, ", "))")
    end

    # Construct validated structs
    control = EquilControl(; symbolize_keys(control_data)...)
    if !isabspath(control.eq_filename)
        control.eq_filename = normpath(joinpath(dirname(path), control.eq_filename))
    end
    output  = EquilOutput(; symbolize_keys(output_data)...)

    return EquilConfig(control=control, output=output)
end



"""
    LargeAspectRatioConfig(...)  

A mutable struct holding parameters for the Large Aspect Ratio (LAR) plasma equilibrium model.

## Fields:

- `lar_r0`: The major radius of the plasma [m].
- `lar_a`: The minor radius of the plasma [m].
- `beta0`: The beta value on axis (normalized pressure).
- `q0`: The safety factor on axis.
- `p_pres`: The exponent for the pressure profile, defined as `p00 * (1 - (r / a)^2)^p_pres`.
- `p_sig`: The exponent that determines the shape of the current-related function profile.
- `sigma_type`: The type of sigma profile, can be "default" or "wesson". If "wesson", the sigma profile is defined as `sigma0 * (1 - (r / a)^2)^p_sig`.
- `mtau`: The number of grid points in the poloidal direction.
- `ma`: The number of grid points in the radial direction.
- `zeroth`: If set to true, it neglects the Shafranov shift
"""
@kwdef mutable struct LargeAspectRatioConfig
    lar_r0::Float64 = 10.0    # Major radius of the plasma
    lar_a::Float64 = 1.0      # Minor radius of the plasma

    beta0::Float64 = 1e-3     # beta on axis
    q0::Float64 = 1.5         # q (safety factor) on axis

    p_pres::Float64 = 2.0     # p00 * (1-(r/a)**2)**p_pres
    p_sig::Float64 = 1.0      # The exponent that determines the shape of the current-related function profile

    sigma_type::String = "default" # can be 'default' or 'wesson'. If 'wesson', switch sigma profile to sigma0*(1-(r/a)**2)**p_sig

    mtau::Int = 128       # the number of grid points in the poloidal direction
    ma::Int = 128         # the number of grid points in the radial direction

    zeroth ::Bool = false     #  If set to true, it neglects the Shafranov shift, creating an ideal concentric circular cross-section.
end

"""
Outer constructor for LargeAspectRatioConfig that enables a toml file 
    interface for specifying the configuration settings
"""
function LargeAspectRatioConfig(path::String)
    raw = TOML.parsefile(path)
    input_data = get(raw, "LAR_INPUT", Dict())
    return LargeAspectRatioConfig(; symbolize_keys(input_data)...)
end



"""
    SolevevConfig(...)  

A mutable struct holding parameters for the Solev'ev (SOL) plasma equilibrium model.

## Fields:

- `mr`: number of radial grid zones
- `mz`: number of axial grid zones
- `ma`: number of flux grid zones
- `e`:  elongation
- `a`: minor radius
- `r0`: major radius
- `q0`: safety factor at the o-point
- `p0fac`: scale on-axis pressure (P-> P+P0*p0fac. beta changes. Phi,q constant)
- `b0fac`: scale toroidal field at constant beta (s*Phi,s*f,s^2*P. bt changes. Shape,beta constant)
- `f0fac`: scale toroidal field at constant pressure (s*f. beta,q changes. Phi,p,bp constant)
"""

@kwdef mutable struct SolevevConfig
    mr::Int = 128      # number of radial grid zones
    mz::Int = 128      # number of axial grid zones
    ma::Int = 128      # number of flux grid zones
    e::Float64 = 1.6       # elongation
    a::Float64 = 0.33      # minor radius
    r0::Float64 = 1.0      # major radius
    q0::Float64 = 1.9      # safety factor at the o-point
    p0fac::Float64 = 1       # scale on-axis pressure (P-> P+P0*p0fac. beta changes. Phi,q constant)
    b0fac::Float64 = 1       # scale toroidal field at constant beta (s*Phi,s*f,s^2*P. bt changes. Shape,beta constant)
    f0fac::Float64 = 1       # scale toroidal field at constant pressure (s*f. beta,q changes. Phi,p,bp constant)
end

"""
Outer constructor for LarConfig that enables a toml file 
    interface for specifying the configuration settings
"""
function SolevevConfig(path::String) # if we use @kwdef, it generates SolevevConfig() so it conflicts with this line.
    raw = TOML.parsefile(path)
    input_data = get(raw, "SOL_INPUT", Dict())
    return SolevevConfig(; symbolize_keys(input_data)...)
end


"""
    DirectRunInput(...)

A container struct that bundles all necessary inputs for the `direct_run` function.
It is created by the `_read_efit` function after parsing the raw equilibrium file
and preparing the initial splines.

## Fields:
- `equil_input`: The original `EquilInput` object.
- `sq_in`
        # x value: psin
        # Quantity 1: F = R*Bt  [m T]
        # Quantity 2: mu0 * Pressure (non-negative) [nt^2 / m^2 * mu0 = T^2]
        # Quantity 3: q-profile 
        # Quantity 4: sqrt(psi_norm)
- `psi_in`:
        # x, y value: R, Z [m]
        # z value : poloidal flux adjusted to be zero at the boundary [Weber/radian]
            # 1. ψ(R,Z) = ψ_boundary - ψ(R,Z)
            # 2. if ψ = ψ * sign(ψ(centerR,centerZ))
- `rmin`: Minimum R-coordinate of the computational grid [m].
- `rmax`: Maximum R-coordinate of the computational grid [m].
- `zmin`: Minimum Z-coordinate of the computational grid [m].
- `zmax`: Maximum Z-coordinate of the computational grid [m].
- `psio`: The total flux difference `abs(ψ_axis - ψ_boundary)` [Weber / radian].
"""
mutable struct DirectRunInput
    config::EquilConfig
    sq_in::Any       # 1D profile spline (CubicSplineType)
    psi_in::Any      # 2D flux spline (BicubicSplineType)
    rmin::Float64
    rmax::Float64
    zmin::Float64
    zmax::Float64
    psio::Float64
end

"""
    InverseRunInput(...)

A container struct for inputs to the `inverse_run` function.

## Fields:
- `equil_input`: The original `EquilInput` object.
"""
mutable struct InverseRunInput
    config::EquilConfig

    sq_in::Any           # 1D spline input profile (e.g. F*Bt, Pressure, q)
    rz_in::Any           # 2D bicubic spline input for (R,Z) geometry

    ro::Float64          # R axis location
    zo::Float64          # Z axis location
    psio::Float64        # Total flux difference |psi_axis - psi_boundary|
end



Base.@kwdef mutable struct EquilibriumParameters
    ro::Union{Nothing, Float64} = nothing
    zo::Union{Nothing, Float64} = nothing
    psio::Union{Nothing, Float64} = nothing
    rsep::Union{Nothing, Vector{Float64}} = nothing
    zsep::Union{Nothing, Vector{Float64}} = nothing
    rext::Union{Nothing, Vector{Float64}} = nothing
    zext::Union{Nothing, Vector{Float64}} = nothing
    psi0::Union{Nothing, Float64} = nothing
    b0::Union{Nothing, Float64} = nothing
    q0::Union{Nothing, Float64} = nothing
    qmin::Union{Nothing, Float64} = nothing
    qmax::Union{Nothing, Float64} = nothing
    qa::Union{Nothing, Float64} = nothing
    q95::Union{Nothing, Float64} = nothing
    qextrema_psi::Union{Nothing, Vector{Float64}} = nothing
    qextrema_q::Union{Nothing, Vector{Float64}} = nothing
    psi_norm::Union{Nothing, Float64} = nothing
    b_norm::Union{Nothing, Float64} = nothing
    psi_axis::Union{Nothing, Float64} = nothing
    psi_boundary::Union{Nothing, Float64} = nothing
    psi_boundary_norm::Union{Nothing, Float64} = nothing
    psi_axis_norm::Union{Nothing, Float64} = nothing
    psi_boundary_offset::Union{Nothing, Float64} = nothing
    psi_axis_offset::Union{Nothing, Float64} = nothing
    psi_boundary_sign::Union{Nothing, Int} = nothing
    psi_axis_sign::Union{Nothing, Int} = nothing
    psi_boundary_zero::Union{Nothing, Bool} = nothing
    rmean::Union{Nothing, Float64} = nothing
    amean::Union{Nothing, Float64} = nothing
    aratio::Union{Nothing, Float64} = nothing
    kappa::Union{Nothing, Float64} = nothing
    delta1::Union{Nothing, Float64} = nothing
    delta2::Union{Nothing, Float64} = nothing
    bt0::Union{Nothing, Float64} = nothing
    crnt::Union{Nothing, Float64} = nothing
    bwall::Union{Nothing, Float64} = nothing
end





"""
    PlasmaEquilibrium(...)

The final, self-contained result of the equilibrium reconstruction. This object
provides a complete representation of the processed plasma equilibrium in flux coordinates.

## Fields:
- `equil_input`: The original `EquilInput` object used for the reconstruction.
- `sq`: The final 1D profile spline (`RealSplineType`).
        # x value: normalized psi
        # Quantity 1: Toroidal Field Function * 2π, `F * 2π` (where `F = R * B_toroidal`)
        # Quantity 2: Pressure * μ₀, `P * μ₀`.
        # Quantity 3: dVdpsi
        # Quantity 4: q
- `rzphi`: The final 2D flux-coordinate mapping spline (`BicubicSplineType`).
        # x value: normlized psi
        # y value: SFL poloidal angle [0,1]
        # Quantity 1: r_coord² = (R - ro)² + (Z - zo)²
        # Quantity 2: Offset between the geometric poloidal angle (η) and the new angle (θ_new)
                     `η / (2π) - θ_new
        # Quantity 3: ν in ϕ=2πζ+ν(ψ,θ)
        # Quantity 4: Jacobian.
-   `eqfun`: A 2D spline storing local physics and geometric quantities that vary across the flux surfaces.
        # These are pre-calculated for efficient use in subsequent stability and transport codes.
        # x value: Normalized poloidal flux, ψ_norm ∈ [0, 1].
        # y value: SFL poloidal angle, θ_new ∈ [0, 1].
        # Quantity 1: Total magnetic field strength, B [T]
        # Quantity 2: (e₁⋅e₂ + q⋅e₃⋅e₁) / (J⋅B²).
        # Quantity 3: (e₂⋅e₃ + q⋅e₃⋅e₃) / (J⋅B²).
- `ro`: R-coordinate of the magnetic axis [m].
- `zo`: Z-coordinate of the magnetic axis [m].
- `psio`: Total flux difference `|Ψ_axis - Ψ_boundary|` [Weber / radian].
"""
mutable struct PlasmaEquilibrium
    config::EquilConfig
    params::EquilibriumParameters  # Parameters for the equilibrium
    sq::Spl.CubicSplineType                     # Final 1D profile spline
    rzphi::Spl.BicubicSplineType                # Final 2D coordinate mapping spline
    eqfun::Spl.BicubicSplineType
    ro::Float64
    zo::Float64
    psio::Float64
end
