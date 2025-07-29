# src/Equilibrium/EquilibriumTypes.jl

"""
- `EquilInput`:          User-facing input parameters.
- `DirectRunInput`:      Internal data structure for the direct solver.
- `InverseRunInput`:     Internal data structure for the inverse solver.
- `PlasmaEquilibrium`:   The final, user-facing output object.
"""

"""
    EquilInput(...)

A mutable struct holding the input parameters that control the equilibrium
reconstruction process. These parameters are typically set by the user.

## Arguments:
- `eq_filename`: A string with the path to the equilibrium file (e.g., g-file).
- `eq_type`: A string specifying the equilibrium file format (e.g., "efit").
- `jac_type`: A string specifying the Jacobian type for coordinate transformation.
- `psilow`: The lower bound of the normalized psi grid for output (typically 0.0).
- `psihigh`: The upper bound of the normalized psi grid for output (typically 1.0).
- `mpsi`: The number of grid points in the `psi` direction for the output.
- `mtheta`: The number of grid points in the `theta` direction for the output.

## Keyword Arguments:
- `newq0`: Target q-value on axis. If non-zero, triggers q-profile revision (Default: 0.0).
"""
# it is mutable because power_bp _r _b is defined and then modified by jac_type in input.
# Maybe if we organize the order a bit more, we can change it to a struct.
mutable struct EquilInput
    eq_filename::String
    eq_type::String
    jac_type::String
    power_bp::Int
    power_r::Int
    power_b::Int
    grid_type::String
    psilow::Float64
    psihigh::Float64
    mpsi::Int
    mtheta::Int
    newq0::Float64

    function EquilInput(
        eq_filename::String, eq_type::String, jac_type::String,
        psilow::Float64, psihigh::Float64, mpsi::Int, mtheta::Int;
        newq0::Float64=0.0
    )
        new(eq_filename, eq_type, jac_type, 0, 0, 0, "ldp", psilow, psihigh, mpsi, mtheta, newq0)
    end
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
    equil_input::EquilInput
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
    equil_input::EquilInput

    sq_in::Any           # 1D spline input profile (e.g. F*Bt, Pressure, q)
    rz_in::Any           # 2D bicubic spline input for (R,Z) geometry

    ro::Float64          # R axis location
    zo::Float64          # Z axis location
    psio::Float64        # Total flux difference |psi_axis - psi_boundary|
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
    equil_input::EquilInput
    sq::Any             # Final 1D profile spline
    rzphi::Any          # Final 2D coordinate mapping spline
    eqfun::Any
    ro::Float64
    zo::Float64
    psio::Float64
end


"""
    LarInput(...)  

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

mutable struct LarInput
    lar_r0::Float64     # Major radius of the plasma
    lar_a::Float64      # Minor radius of the plasma

    beta0::Float64      # beta on axis
    q0::Float64         # q (safety factor) on axis

    p_pres::Float64     # p00 * (1-(r/a)**2)**p_pres
    p_sig::Float64      # The exponent that determines the shape of the current-related function profile

    sigma_type::String  # can be 'default' or 'wesson'. If 'wesson', switch sigma profile to sigma0*(1-(r/a)**2)**p_sig

    mtau::Int      # the number of grid points in the poloidal direction
    ma::Int        # the number of grid points in the radial direction

    zeroth ::Bool      #  If set to true, it neglects the Shafranov shift, creating an ideal concentric circular cross-section.
end
