"""
    DirectBField

Internal mutable struct to hold B-field components and their derivatives at a point.
It is used as a temporary workspace to avoid allocations in tight loops.
"""
@kwdef mutable struct DirectBField
    psi::Float64 = 0.0
    psir::Float64 = 0.0   # d(psi)/dr
    psiz::Float64 = 0.0   # d(psi)/dz
    psirz::Float64 = 0.0  # d2(psi)/drdz
    psirr::Float64 = 0.0  # d2(psi)/drdr
    psizz::Float64 = 0.0  # d2(psi)/dzdz
    f::Float64 = 0.0      # F = R*Bt
    f1::Float64 = 0.0     # dF/d(psi_norm)
    p::Float64 = 0.0      # mu0*Pressure
    p1::Float64 = 0.0     # dP/d(psi_norm)
    br::Float64 = 0.0     # Br
    bz::Float64 = 0.0     # Bz
    brr::Float64 = 0.0    # d(Br)/dr
    brz::Float64 = 0.0    # d(Br)/dz
    bzr::Float64 = 0.0    # d(Bz)/dr
    bzz::Float64 = 0.0    # d(Bz)/dz
end

"""
    FieldLineDerivParams

A struct to hold constant parameters for the ODE integration, making them
easily accessible within the derivative function `direct_fl_der!`.
"""
struct FieldLineDerivParams
    ro::Float64
    zo::Float64
    psi_in::Spl.BicubicSpline
    sq_in::Spl.CubicSpline
    psio::Float64
    power_bp::Int
    power_b::Int
    power_r::Int
    bfield::DirectBField
end

"""
    direct_get_bfield!(bf_out, r, z, psi_in, sq_in, psio; derivs=0)

Calculates the magnetic field and its derivatives at a given (R,Z) point.
The results are stored in-place in the `bf_out` object. This is equivalent
to the `direct_get_bfield` subroutine in the Fortran code, adapted for
the Julia spline implementation.

## Arguments:

  - `bf_out`: A mutable `DirectBField` struct to store the results
  - `r`: R-coordinate to evaluate at
  - `z`: Z-coordinate to evaluate at
  - `psi_in`: 2D bicubic spline for poloidal flux `ψ(R,Z)`
  - `sq_in`: 1D cubic spline for profiles `F(ψ_norm)` and `P(ψ_norm)`
  - `psio`: total toroidal flux
  - `derivs`: An integer specifying number of derivatives to compute (0, 1, or 2)
"""
function direct_get_bfield!(
    bf_out::DirectBField,
    r::Float64,
    z::Float64,
    psi_in::Spl.BicubicSpline,
    sq_in::Spl.CubicSpline,
    psio::Float64;
    derivs::Int=0
)
    # Evaluate 2D spline for psi(r,z) and its derivatives
    if derivs == 0
        f_psi = Spl.bicube_eval!(psi_in, r, z)
        bf_out.psi = f_psi[1]
    elseif derivs == 1
        f_psi, fx_psi, fy_psi = Spl.bicube_deriv1!(psi_in, r, z)
        bf_out.psi = f_psi[1]
        bf_out.psir = fx_psi[1]
        bf_out.psiz = fy_psi[1]
    else # derivs >= 2
        f_psi, fx_psi, fy_psi, fxx_psi, fxy_psi, fyy_psi = Spl.bicube_deriv2!(psi_in, r, z)
        bf_out.psi = f_psi[1]
        bf_out.psir = fx_psi[1]
        bf_out.psiz = fy_psi[1]
        bf_out.psirr = fxx_psi[1]
        bf_out.psirz = fxy_psi[1]
        bf_out.psizz = fyy_psi[1]
    end

    # Evaluate magnetic fields from equilibrium profiles
    psi_norm = (psio > 1e-12) ? (1.0 - bf_out.psi / psio) : 0.0
    psi_norm = clamp(psi_norm, 0.0, 1.0)

    f_sq, f1_sq = Spl.spline_deriv1!(sq_in, psi_norm)
    bf_out.f = f_sq[1]  # F = R*Bt
    bf_out.f1 = f1_sq[1] # dF/dψ
    bf_out.p = f_sq[2]  # μ0*Pressure
    bf_out.p1 = f1_sq[2] # dP/dψ

    (derivs == 0) && return

    # Evaluate B-field derivative components
    bf_out.br = bf_out.psiz / r # Br = (1/R) * ∂ψ/∂Z
    bf_out.bz = -bf_out.psir / r # Bz = -(1/R) * ∂ψ/∂R

    (derivs == 1) && return

    # Evaluate more derivatives of B-field components
    bf_out.brr = (bf_out.psirz - bf_out.br) / r
    bf_out.brz = bf_out.psizz / r
    bf_out.bzr = -(bf_out.psirr + bf_out.bz) / r
    bf_out.bzz = -bf_out.psirz / r
end

"""
    direct_position(raw_profile)

Finds the key geometric locations of the equilibrium: the magnetic axis (O-point)
and the inboard/outboard separatrix crossings on the midplane. This performs the same
overall function as the Fortran `direct_position` subroutine with better iteration
control and error handling. We have also added a helper function for separatrix finding.

## Arguments:

  - `raw_profile`: A `DirectRunInput` object containing splines and parameters.

## Returns:

  - `ro`: R-coordinate of the magnetic axis [m].
  - `zo`: Z-coordinate of the magnetic axis [m].
  - `rs1`: R-coordinate of the inboard separatrix crossing [m].
  - `rs2`: R-coordinate of the outboard separatrix crossing [m].
  - `psi_in_new` : returns psi_in renormalized by * psio/psi(ro,zo)
"""
function direct_position(raw_profile::DirectRunInput)

    bfield = DirectBField()
    max_iterations = 200

    # For an axis initial guess, we find zero crossing of Bz on the midplane
    # Note the Fortran had this wrapped in a if ro == 0 block, which we omit here
    # because I don't think it's ever used. If needed, it can be re-added.
    r = (raw_profile.rmax + raw_profile.rmin) / 2.0
    z = (raw_profile.zmax + raw_profile.zmin) / 2.0
    dr = (raw_profile.rmax - raw_profile.rmin) / 20.0

    for _ in 1:max_iterations
        direct_get_bfield!(bfield, r, z, raw_profile.psi_in, raw_profile.sq_in, raw_profile.psio; derivs=1)
        if bfield.bz >= 0.0
            break
        end
        r += dr
    end

    # If we never exited early, the loop failed to find bz = 0
    !(bfield.bz >= 0) && error("Took too many iterations to get bz=0.")

    # Now, use Newton iteration to find the O-point (magnetic axis) where Br=0 and Bz=0
    dr, dz = 0.0, 0.0
    for _ in 1:max_iterations
        direct_get_bfield!(bfield, r, z, raw_profile.psi_in, raw_profile.sq_in, raw_profile.psio; derivs=2)
        det = bfield.brr * bfield.bzz - bfield.brz * bfield.bzr
        if abs(det) < 1e-20
            error("Jacobian matrix is singular near ($r, $z).")
        end
        # Δx = -J⁻¹ F
        dr = (bfield.brz * bfield.bz - bfield.bzz * bfield.br) / det
        dz = (bfield.bzr * bfield.br - bfield.brr * bfield.bz) / det
        r += dr
        z += dz
        if abs(dr) <= 1e-12 * abs(r) && abs(dz) <= 1e-12 * abs(r)
            @printf("   Magnetic axis found at R = %.5f, Z = %.5f\n", r, z)
            break
        end
    end

    if !(abs(dr) <= 1e-12 * abs(r) && abs(dz) <= 1e-12 * abs(r))
        error("Failed to find magnetic axis after $max_iterations iterations.")
    else
        ro, zo = r, z
    end

    # Renormalize psi based on the value at the magnetic axis
    direct_get_bfield!(bfield, ro, zo, raw_profile.psi_in, raw_profile.sq_in, raw_profile.psio; derivs=0)
    fac = raw_profile.psio / bfield.psi
    new_psi_fs = raw_profile.psi_in.fs .* fac
    x_coords = Vector(raw_profile.psi_in.xs)
    y_coords = Vector(raw_profile.psi_in.ys)
    # Because DirectRunInput is a mutable struct, we can update the spline here
    raw_profile.psi_in = Spl.BicubicSpline(x_coords, y_coords, new_psi_fs; bctypex=3, bctypey=3)

    # Helper function for robust Newton-Raphson search with restarts
    function find_separatrix_crossing(start_r, end_r, label)
        local r_sol::Float64
        found = false
        for ird in 0:5 # 6 restart attempts
            r_sep = (start_r * (3.0 - 0.5 * ird) + end_r) / (4.0 - 0.5 * ird)
            for _ in 1:max_iterations

                direct_get_bfield!(bfield, r_sep, zo, raw_profile.psi_in, raw_profile.sq_in, raw_profile.psio; derivs=1)
                if abs(bfield.psir) < 1e-14
                    @warn "d(psi)/dr is near zero."
                    break
                end
                dr = -bfield.psi / bfield.psir
                r_sep += dr
                if abs(dr) <= 1e-12 * abs(r_sep)
                    r_sol = r_sep
                    found = true
                    break
                end
            end
            if found
                break
            end
        end
        !found && error("Could not find $label separatrix after all attempts.")
        println("   $label separatrix found at R = $(r_sol).")
        return r_sol
    end

    # Find inboard (rs1) and outboard (rs2) separatrix positions
    rs1 = find_separatrix_crossing(ro, raw_profile.rmin, "Inboard")
    rs2 = find_separatrix_crossing(ro, raw_profile.rmax, "Outboard")

    return ro, zo, rs1, rs2
end

"""
    direct_fieldline_int(psifac, raw_profile, ro, zo, rs2)

Performs the field-line integration for a single flux surface. This is a Julia adaptation
of the Fortran `direct_fl_int` subroutine. Note that the array `y_out` is now indexed
from 1:5 rather than 0:4 as in Fortran.

## Arguments:

  - `psifac`: normalized psi value for the surface (ψ_norm).
  - `raw_profile`: `DirectRunInput` object containing splines and parameters.
  - `ro`, `zo`: Coordinates of the magnetic axis [m].
  - `rs2`: R-coordinate of the outboard separatrix [m].

## Returns:

    - `y_out`: A matrix containing the integrated quantities vs. the geometric angle `η`.
        - `y_out[:, 1]`: η (geometric poloidal angle)
        - `y_out[:, 2]`: ∫(dl/Bp)
        - `y_out[:, 3]`: rfac (radial distance from magnetic axis)
        - `y_out[:, 4]`: ∫(dl/(R²Bp))
        - `y_out[:, 5]`: ∫(jac*dl/Bp)

  - `bfield`: A `DirectBField` object with values at the integration start point.
"""
function direct_fieldline_int(psifac::Float64, raw_profile::DirectRunInput, ro::Float64, zo::Float64, rs2::Float64)

    # Find the starting point on the flux surface (outboard midplane)
    psi0 = raw_profile.psio * (1.0 - psifac)
    r = ro + sqrt(psifac) * (rs2 - ro)
    z = zo
    bfield = DirectBField()

    # Refine starting R using Newton's method
    dr = 0.0
    for _ in 1:10
        direct_get_bfield!(bfield, r, z, raw_profile.psi_in, raw_profile.sq_in, raw_profile.psio; derivs=1)
        dr = (psi0 - bfield.psi) / bfield.psir
        r += dr
        if abs(dr) <= 1e-12 * r
            break
        end
    end
    !(abs(dr) <= 1e-12 * r) && error("Failed to refine starting R on flux surface.")

    direct_get_bfield!(bfield, r, z, raw_profile.psi_in, raw_profile.sq_in, raw_profile.psio; derivs=2)
    psi0 = bfield.psi

    # Set up and solve the ODE for fieldline following
    # Initial condition
    u0 = zeros(Float64, 4)
    u0[2] = sqrt((r - ro)^2 + (z - zo)^2)

    bfield = DirectBField()
    equil_input = raw_profile.config.control
    params = FieldLineDerivParams(ro, zo, raw_profile.psi_in, raw_profile.sq_in, raw_profile.psio,
        equil_input.power_bp, equil_input.power_b, equil_input.power_r, bfield)

    # Use a callback to refine the solution at each step to stay on the flux surface
    function refine_affect!(integrator)
        integrator.u[2] = direct_refine(integrator.u[2], integrator.t, psi0, params)
    end

    # We save the solution at each step before refinement (before=true, after=false) to match Fortran
    callback = DiscreteCallback((u, t, i) -> true, refine_affect!; save_positions=(true, false))

    prob = ODEProblem{true}(direct_fieldline_der!, u0, (0.0, 2π), params)
    sol = solve(prob, Tsit5(); callback=callback, reltol=1e-6, abstol=1e-8, dt=2π / 200, adaptive=true)

    if sol.retcode != :Success && sol.retcode != :Terminated
        error("ODE integration failed for psi = $psifac with code: $(sol.retcode)")
    end

    return hcat(sol.t, hcat(sol.u...)'), bfield
end

"""
    direct_fieldline_der!(dy, y, params, eta)

The derivative function for the field-line integration ODE. This is passed to
the `DifferentialEquations.jl` solver. This is a Julia adaptation of the Fortran
`direct_fl_der` subroutine, with an added safeguard against division by zero.

## Arguments:

  - `dy`: The derivative vector (output, modified in-place).
  - `y`: The state vector `[∫(dl/Bp), rfac, ∫(dl/(R²Bp)), ∫(jac*dl/Bp)]`.
  - `params`: A `FieldLineDerivParams` struct with all necessary parameters.
  - `eta`: The independent variable (geometric angle `η`).
"""
function direct_fieldline_der!(dy, y, params::FieldLineDerivParams, eta)

    cos_eta, sin_eta = cos(eta), sin(eta)
    r = params.ro + y[2] * cos_eta
    z = params.zo + y[2] * sin_eta
    direct_get_bfield!(params.bfield, r, z, params.psi_in, params.sq_in, params.psio; derivs=1)

    bp = sqrt(params.bfield.br^2 + params.bfield.bz^2)
    bt = params.bfield.f / r
    b = sqrt(bp^2 + bt^2)
    jac = (bp^params.power_bp) * (b^params.power_b) / (r^params.power_r)

    # Denominator for d(l_pol)/d(eta) = rfac |B_pol|/denominator
    denominator = params.bfield.bz * cos_eta - params.bfield.br * sin_eta
    if abs(denominator) < 1e-14
        fill!(dy, 1e20) # Return large derivatives for solver to handle stiffness
        @warn "Denominator in direct_fieldline_der! near zero at eta=$eta."
        return
    end

    # Compute derivatives
    # d/dη [∫(dl/Bp)] = 1/|B_P| dl/d(eta) = rfac/denominator
    dy[1] = y[2] / denominator
    # d(rfac)/d(eta) = rfac/denom *( Br cos(eta) + Bz sin(eta) )
    dy[2] = dy[1] * (params.bfield.br * cos_eta + params.bfield.bz * sin_eta)
    # d/dη [∫(dl/(R²Bp))]
    dy[3] = dy[1] / (r^2)
    # d/dη [∫(jac*dl/Bp)]
    dy[4] = dy[1] * jac
end

"""
    direct_refine(rfac, eta, psi0, params)

Refines the radial distance `rfac` at a given angle `eta` to ensure the
point lies exactly on the target flux surface `psi0`. This performs the
same function as the Fortran `direct_refine` subroutine, with more clear
iteration control and error handling.

## Arguments:

  - `rfac`: The current guess for the radial distance from the magnetic axis.
  - `eta`: The geometric poloidal angle.
  - `psi0`: The target `ψ` value for the flux surface.
  - `params`: A `FieldLineDerivParams` struct.

## Returns:

  - The refined `rfac` value.
"""
function direct_refine(rfac::Float64, eta::Float64, psi0::Float64, params::FieldLineDerivParams; max_iter::Int=50)::Float64

    cos_eta, sin_eta = cos(eta), sin(eta)
    r = params.ro + rfac * cos_eta
    z = params.zo + rfac * sin_eta
    direct_get_bfield!(params.bfield, r, z, params.psi_in, params.sq_in, params.psio; derivs=1)
    dpsi = params.bfield.psi - psi0

    for _ in 1:max_iter
        # Newton's method derivative: d(psi)/d(rfac)
        dpsi_drfac = params.bfield.psir * cos_eta + params.bfield.psiz * sin_eta
        if abs(dpsi_drfac) < 1e-14
            @warn "Refinement failed at eta=$eta: d(psi)/d(rfac) is zero."
            return rfac # Return current best guess
        end
        drfac = -dpsi / dpsi_drfac
        rfac += drfac
        r = params.ro + rfac * cos_eta
        z = params.zo + rfac * sin_eta
        direct_get_bfield!(params.bfield, r, z, params.psi_in, params.sq_in, params.psio; derivs=1)
        dpsi = params.bfield.psi - psi0

        if abs(dpsi) <= 1e-12 * psi0 || abs(drfac) <= 1e-12 * abs(rfac)
            return rfac
        end
    end

    error("direct_refine did not converge after $max_iter iterations at eta=$eta.")
end

"""
    equilibrium_solver(raw_profile)

The main driver for the direct equilibrium reconstruction. It orchestrates the entire
process from finding the magnetic axis to integrating along field lines and
constructing the final coordinate and physics quantity splines. This performs the same
overall function as the Fortran `direct_run` subroutine, with better checks for numerical
robustness.

## Arguments:

  - `raw_profile`: A `DirectRunInput` object containing the initial splines (`psi_in`, `sq_in`)
    and run parameters (`equil_input`).

## Returns:

  - A `PlasmaEquilibrium` object containing the final, processed equilibrium data,
    including the profile spline (`sq`), the coordinate mapping spline (`rzphi`), and
    the physics quantity spline (`eqfun`).
"""
function equilibrium_solver(raw_profile::DirectRunInput)

    # Shorthand
    equil_params = raw_profile.config.control
    sq_in = raw_profile.sq_in
    psi_in = raw_profile.psi_in
    psio = raw_profile.psio
    mtheta = equil_params.mtheta
    mpsi = equil_params.mpsi
    psilow = equil_params.psilow
    psihigh = equil_params.psihigh

    # Warn if psihigh is too close to 1.0
    if psihigh >= 1 - 1e-6
        @warn "Warning: direct equilibrium with psihigh = $(psihigh) could hang on separatrix."
    end

    # TODO: there's some fortran logic for grid_type = original that should be added when needed.

    # Set up radial and poloidal grid
    psi_nodes = Array{Float64}(undef, mpsi + 1)
    if equil_params.grid_type == "ldp"
        psi_nodes .= [psilow + (psihigh - psilow) * sin((ipsi / mpsi) * (π / 2))^2 for ipsi in 0:mpsi]
    else
        # TODO: add additional grid types
        error("Unsupported grid_type: $(equil_params.grid_type)")
    end
    theta_nodes = range(0.0, 1.0; length=mtheta + 1)

    # Find radial position of magnetic axis and separatrix
    ro, zo, rs1, rs2 = direct_position(raw_profile)

    # Loop over flux surfaces from outermost to innermost, integrating over field lines
    sq_fs_nodes = zeros(Float64, mpsi + 1, 4)
    rzphi_fs_nodes = zeros(Float64, mpsi + 1, mtheta + 1, 4)
    for ipsi in (mpsi+1):-1:1
        # Integrate along the field line for this surface
        y_out, bfield = direct_fieldline_int(psi_nodes[ipsi], raw_profile, ro, zo, rs2)

        # Fit data into temporary straight fieldline poloidal angle splines
        ff_x_nodes = y_out[:, 5] ./ y_out[end, 5]
        ff_fs_nodes = hcat(
            y_out[:, 3] .^ 2,
            y_out[:, 1] / (2π) .- ff_x_nodes,
            bfield.f * (y_out[:, 4] .- ff_x_nodes .* y_out[end, 4]),
            y_out[:, 2] ./ y_out[end, 2] .- ff_x_nodes
        )
        ff = Spl.CubicSpline(ff_x_nodes, ff_fs_nodes; bctype="periodic")

        # Interpolate `ff` onto the uniform `theta` grid for `rzphi`
        for itheta in 1:(mtheta+1)
            f, f1 = Spl.spline_deriv1!(ff, theta_nodes[itheta])
            @views rzphi_fs_nodes[ipsi, itheta, 1:3] = f[1:3]
            jac_term = (1.0 + f1[4]) * y_out[end, 2] * 2π * psio
            rzphi_fs_nodes[ipsi, itheta, 4] = jac_term
        end

        # Store surface-averaged quantities for the `sq` spline
        sq_fs_nodes[ipsi, 1] = bfield.f * 2π
        sq_fs_nodes[ipsi, 2] = bfield.p
        sq_fs_nodes[ipsi, 3] = y_out[end, 2] * 2π * psio
        sq_fs_nodes[ipsi, 4] = y_out[end, 4] * bfield.f / (2π)
    end

    # Fit 1D profile spline `sq` and perform q-profile revision if needed
    sq = Spl.CubicSpline(psi_nodes, sq_fs_nodes; bctype="extrap")
    q0 = sq.fs[1, 4] - sq.fs1[1, 4] * sq.xs[1]
    if equil_params.newq0 == -1
        equil_params.newq0 = -q0
    end
    if equil_params.newq0 != 0.0
        println("Revising q-profile for newq0 = $(equil_params.newq0)...")
        f0 = sq.fs[1, 1] - sq.fs1[1, 1] * sq.xs[1]
        f0fac = f0^2 * ((equil_params.newq0 / q0)^2 - 1.0)
        for i in 1:(mpsi+1)
            ffac = sqrt(1.0 + f0fac / sq.fs[i, 1]^2) * sign(equil_params.newq0)
            sq_fs_nodes[i, 1] *= ffac
            sq_fs_nodes[i, 4] *= ffac
            rzphi_fs_nodes[i, :, 3] .*= ffac
        end
        # Re-create the spline with the revised data
        sq = Spl.CubicSpline(psi_nodes, sq_fs_nodes; bctype="extrap")
    end

    # Fit the 2D geometric spline `rzphi`. Periodic in theta (y-dimension)
    rzphi = Spl.BicubicSpline(psi_nodes, collect(theta_nodes), rzphi_fs_nodes; bctypex="extrap", bctypey="periodic")

    # Calculate physics quantities (B-field, metric components, etc.) in 2D spline `eqfun`
    # for use in stability and transport codes
    eqfun_fs_nodes = zeros(Float64, mpsi + 1, mtheta + 1, 3)
    v = @MMatrix zeros(Float64, 2, 3)
    for ipsi in 1:(mpsi+1)
        psi_norm = psi_nodes[ipsi]
        fsq = Spl.spline_eval!(sq, psi_norm)
        q = fsq[4]
        f_val = fsq[1]
        for itheta in 1:(mtheta+1)
            theta_norm = theta_nodes[itheta]
            f, fx, fy = Spl.bicube_deriv1!(rzphi, psi_norm, theta_norm)
            rfac = sqrt(max(0.0, f[1])) # add in protection just in case of small negative due to numerical error
            eta = 2π * (theta_norm + f[2])
            r = ro + rfac * cos(eta)
            jacfac = f[4]

            v[1, 1] = (rfac > 0) ? fx[1] / (2.0 * rfac) : 0.0       # 1/(2rfac) * d(rfac)/d(psi_norm)
            v[1, 2] = fx[2] * 2π * rfac                             # 2π*rfac * d(eta)/d(psi_norm)
            v[1, 3] = fx[3] * r                                     # r * d(phi_s)/d(psi_norm)
            v[2, 1] = (rfac > 0) ? fy[1] / (2.0 * rfac) : 0.0       # 1/(2rfac) d(rfac)/d(theta_new)
            v[2, 2] = (1.0 + fy[2]) * 2π * rfac                     # 2π*rfac * d(eta)/d(theta_new)
            v[2, 3] = fy[3] * r                                     # r * d(phi_s)/d(theta_new)
            v33 = 2π * r
            w11 = (jacfac != 0) ? (1.0 + fy[2]) * (2π)^2 * rfac * r / jacfac : 0.0
            w12 = (jacfac * rfac != 0) ? -fy[1] * π * r / (rfac * jacfac) : 0.0
            delpsi_norm = sqrt(w11^2 + w12^2)
            modB = sqrt(((2π * psio * delpsi_norm)^2 + f_val^2) / (2π * r)^2)

            # Fill in eqfun nodes
            eqfun_fs_nodes[ipsi, itheta, 1] = modB
            denom = jacfac * modB^2
            if abs(denom) > 1e-20
                # Gyrokinetic coefficient C1
                numerator_2 = dot(v[1, :], v[2, :]) + q * v33 * v[1, 3]
                eqfun_fs_nodes[ipsi, itheta, 2] = numerator_2 / denom
                # Gyrokinetic coefficient C2
                numerator_3 = v[2, 3] * v33 + q * v33^2
                eqfun_fs_nodes[ipsi, itheta, 3] = numerator_3 / denom
            else
                eqfun_fs_nodes[ipsi, itheta, 2] = 0.0
                eqfun_fs_nodes[ipsi, itheta, 3] = 0.0
            end
        end
    end
    eqfun = Spl.BicubicSpline(psi_nodes, collect(theta_nodes), eqfun_fs_nodes; bctypex="extrap", bctypey="periodic")
    return PlasmaEquilibrium(raw_profile.config, EquilibriumParameters(), sq, rzphi, eqfun, ro, zo, psio)
end