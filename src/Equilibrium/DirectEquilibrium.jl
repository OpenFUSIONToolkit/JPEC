#=
This file contains the logic for the "direct" equilibrium reconstruction method.
It takes parsed data and splines from the IO module and calculates the final flux-coordinate
representation of the plasma equilibrium.
=#


# --- Internal Helper Structs ---
"""
    DirectBField

Internal mutable struct to hold B-field components and their derivatives at a point.
It is used as a temporary workspace to avoid allocations in tight loops.
"""
mutable struct DirectBField
    psi::Float64
    psir::Float64   # d(psi)/dr
    psiz::Float64   # d(psi)/dz
    psirz::Float64  # d2(psi)/drdz
    psirr::Float64  # d2(psi)/drdr
    psizz::Float64  # d2(psi)/dzdz
    f::Float64      # F = R*Bt
    f1::Float64     # dF/d(psi_norm)
    p::Float64      # mu0*Pressure
    p1::Float64     # dP/d(psi_norm)
    br::Float64     # Br
    bz::Float64     # Bz
    brr::Float64    # d(Br)/dr
    brz::Float64    # d(Br)/dz
    bzr::Float64    # d(Bz)/dr
    bzz::Float64    # d(Bz)/dz

    DirectBField() = new(zeros(Float64, 16)...)
end

"""
    FieldLineDerivParams

A struct to hold constant parameters for the ODE integration, making them
easily accessible within the derivative function `direct_fl_der!`.
"""
struct FieldLineDerivParams
    ro::Float64
    zo::Float64
    psi_in::Spl.BicubicSplineType
    sq_in::Spl.CubicSpline
    psio::Float64
    power_bp::Int
    power_b::Int
    power_r::Int
end


"""
    direct_get_bfield!(bf_out, r, z, psi_in, sq_in, psio; derivs=0)

Calculates the magnetic field and its derivatives at a given (R,Z) point.
The results are stored in-place in the `bf_out` object.

## Arguments:
- `bf_out`: A mutable `DirectBField` struct to store the results.
- `r`: The R-coordinate [m].
- `z`: The Z-coordinate [m].
- `psi_in`: The 2D bicubic spline for poloidal flux `ψ(R,Z)`.
- `sq_in`: The 1D cubic spline for profiles `F(ψ_norm)` and `P(ψ_norm)`.
- `psio`: The total flux difference `|ψ_axis - ψ_boundary|`.

## Keyword Arguments:
- `derivs`: An integer specifying the derivative level to compute.
    - `0`: Calculates `ψ`, `F`, and `P`.
    - `1`: Also calculates 1st derivatives of `ψ` and `B` field components.
    - `2`: Also calculates 2nd derivatives of `ψ` and 1st derivatives of `B`.

## Returns:
- `nothing`. The `bf_out` object is modified in-place.
"""
function direct_get_bfield!(
    bf_out::DirectBField,
    r::Float64,
    z::Float64,
    psi_in::Spl.BicubicSplineType,
    sq_in::Spl.CubicSpline,
    psio::Float64;
    derivs::Int=0
)
    # 1. Evaluate 2D spline for psi(r,z) and its derivatives
    if derivs == 0
        f_psi = Spl.bicube_eval(psi_in, r, z, 0)
        bf_out.psi = f_psi[1]
    elseif derivs == 1
        f_psi, fx_psi, fy_psi = Spl.bicube_eval(psi_in, r, z, 1)
        bf_out.psi  = f_psi[1]
        bf_out.psir = fx_psi[1]
        bf_out.psiz = fy_psi[1]
    else # derivs >= 2
        f_psi, fx_psi, fy_psi, fxx_psi, fxy_psi, fyy_psi = Spl.bicube_eval(psi_in, r, z, 2)
        bf_out.psi   = f_psi[1]
        bf_out.psir  = fx_psi[1]
        bf_out.psiz  = fy_psi[1]
        bf_out.psirr = fxx_psi[1]
        bf_out.psirz = fxy_psi[1]
        bf_out.psizz = fyy_psi[1]
    end

    # 2. Evaluate 1D splines for F(psi_norm) and P(psi_norm)
    # psi_norm = 0 at axis, 1 at boundary.
    psi_norm = (psio > 1e-12) ? (1.0 - bf_out.psi / psio) : 0.0
    psi_norm = clamp(psi_norm, 0.0, 1.0)

    f_sq, f1_sq = Spl.spline_eval(sq_in, psi_norm, 1)
    bf_out.f  = f_sq[1]  # F = R*Bt
    bf_out.f1 = f1_sq[1] # dF/d(psi_norm)
    bf_out.p  = f_sq[2]  # mu0*Pressure
    bf_out.p1 = f1_sq[2] # dP/d(psi_norm)

    (derivs < 1) && return

    # 3. Evaluate magnetic field components (derivs >= 1)
    bf_out.br =  bf_out.psiz / r # Br = (1/R) * ∂ψ/∂Z
    bf_out.bz = -bf_out.psir / r # Bz = -(1/R) * ∂ψ/∂R

    (derivs < 2) && return

    # 4. Evaluate derivatives of B-field components (derivs >= 2)
    bf_out.brr = (bf_out.psirz - bf_out.br) / r #d(Br)/dr
    bf_out.brz =  bf_out.psizz / r #d(Br)/dz
    bf_out.bzr = -(bf_out.psirr + bf_out.bz) / r #d(Bz)/dr
    bf_out.bzz = -bf_out.psirz / r # d(Bz)/dz

    return
end


"""
    direct_position(psi_in, sq_in, psio, ro_guess, zo_guess, rmin, rmax)

Finds the key geometric locations of the equilibrium: the magnetic axis (O-point)
and the inboard/outboard separatrix crossings on the midplane.

## Arguments:
- `psi_in`: The 2D `ψ(R,Z)` spline.
- `sq_in`: The 1D profile spline.
- `psio`: Initial flux difference.
- `ro_guess`, `zo_guess`: Initial guess for the magnetic axis location [m].
- `rmin`, `rmax`: Radial bounds of the computational domain [m].

## Returns:
- `ro`: R-coordinate of the magnetic axis [m].
- `zo`: Z-coordinate of the magnetic axis [m].
- `rs1`: R-coordinate of the inboard separatrix crossing [m].
- `rs2`: R-coordinate of the outboard separatrix crossing [m].
- `psi_in_new` : returns psi_in renormalized by * psio/psi(ro,zo)
"""
function direct_position(
    psi_in::Spl.BicubicSplineType,
    sq_in::Spl.CubicSpline,
    psio::Float64,
    ro_guess::Float64,
    zo_guess::Float64,
    rmin::Float64,
    rmax::Float64
)
    bf_temp = DirectBField()
    r, z = ro_guess, zo_guess
    max_newton_iter = 200

    # 1. Find the magnetic axis (O-point) using Newton-Raphson
    println("Finding magnetic axis...")
    # (Initial coarse search for r if guess is 0.0)
    # find o point coarse guess by finding point psi extrema by Bz.
    if r ≈ 0.0
        r_mid, z_mid = (rmax + rmin) / 2.0, zo_guess
        dr_scan = (rmax - rmin) / 20.0
        r = r_mid
        for _ in 1:20
            direct_get_bfield!(bf_temp, r, z_mid, psi_in, sq_in, psio, derivs=1)
            if bf_temp.bz >= 0.0; break; end
            r += dr_scan
        end
        z = z_mid
    end

    for iter in 1:20
        direct_get_bfield!(bf_temp, r, z, psi_in, sq_in, psio, derivs=2)
        det = bf_temp.brr * bf_temp.bzz - bf_temp.brz * bf_temp.bzr
        if abs(det) < 1e-20; error("Jacobian matrix is singular near ($r, $z)."); end
        # Δx = -J^-1 F
        dr = (bf_temp.brz * bf_temp.bz - bf_temp.bzz * bf_temp.br) / det
        dz = (bf_temp.bzr * bf_temp.br - bf_temp.brr * bf_temp.bz) / det
        r += dr; z += dz
        @printf "  Iter %2d: R = %.6f, Z = %.6f, |ΔR|=%.2e, |ΔZ|=%.2e\n" iter r z abs(dr) abs(dz)
        if abs(dr) <= 1e-12 * abs(r) && abs(dz) <= 1e-12 * abs(r)
            println("Magnetic axis found at R=$(r), Z=$(z).")
            break
        end
        (iter == 20) && @warn "O-point search did not converge."
    end
    ro, zo = r, z

    # 2. Psi Renormalization
    # renormalize psi_in by psi_in * psio/psi(ro,zo)
    # Redefined to avoid direct modification of the psi_rz structure.
    # However, if it is determined that redefining psi_rz will not cause major problems,
    # modification may be permitted.
    direct_get_bfield!(bf_temp, ro, zo, psi_in, sq_in, psio, derivs=0)
    psi_at_axis = bf_temp.psi
    fac = psio/psi_at_axis
    new_psi_fs = psi_in.fs * fac
    x_coords = Vector(psi_in.xs)
    y_coords = Vector(psi_in.ys)
    psi_in_new = Spl.bicube_setup(x_coords, y_coords, new_psi_fs, bctypex=4, bctypey=4)

    if abs(psi_at_axis - psio) / psio > 1e-3
        @warn "Psi at located axis (O-point) differs from expected psio. " *
              "Psi(axis)=$(psi_at_axis), psio=$(psio). Proceeding without re-fitting."
    end

    # 3. Find inboard (rs1) and outboard (rs2) separatrix positions
    # Helper function for robust Newton-Raphson search with restarts
    function find_separatrix_crossing(start_r, end_r, label)
        println("Finding $label separatrix crossing...")
        local r_sol::Float64
        found = false
        for ird in 0:5 # 6 restart attempts
            r_sep = (start_r * (3.0 - 0.5 * ird) + end_r) / (4.0 - 0.5 * ird)
            @printf "  Restart attempt %d/6 with initial R = %.6f\n" (ird + 1) r_sep
            for _ in 1:max_newton_iter
                
                direct_get_bfield!(bf_temp, r_sep, zo, psi_in_new, sq_in, psio, derivs=1)
                if abs(bf_temp.psir) < 1e-14; @warn "d(psi)/dr is near zero."; break; end
                dr = -bf_temp.psi / bf_temp.psir
                r_sep += dr
                if abs(dr) <= 1e-12 * abs(r_sep); r_sol = r_sep; found = true; break; end
            end
            if found; break; end
        end
        !found && error("Could not find $label separatrix after all attempts.")
        println("$label separatrix found at R=$(r_sol).")
        return r_sol
    end

    rs1 = find_separatrix_crossing(ro, rmin, "inboard")
    rs2 = find_separatrix_crossing(ro, rmax, "outboard")

    # 4. Return the new spline
    return ro, zo, rs1, rs2, psi_in_new
end


"""
    direct_fl_int(ipsi, psifac, raw_profile, ro, zo, rs2)

Performs the field-line integration for a single flux surface.

## Arguments:
- `ipsi`: The index of the current flux surface.
- `psifac`: The normalized psi value for the surface (ψ_norm).
- `raw_profile`: The `DirectRunInput` object containing splines and parameters.
- `ro`, `zo`: Coordinates of the magnetic axis [m].
- `rs2`: R-coordinate of the outboard separatrix [m].

## Returns:
- `y_out`: A matrix containing the integrated quantities vs. the geometric angle `η`.
- `bf_start`: A `DirectBField` object with values at the integration start point.
"""
function direct_fl_int(
    ipsi::Int,
    psifac::Float64,
    raw_profile::DirectRunInput,
    ro::Float64,
    zo::Float64,
    rs2::Float64
)
    # 1. Find the starting point on the flux surface (outboard midplane)
    psi0_target = raw_profile.psio * (1.0 - psifac)
    r_start = ro + sqrt(psifac) * (rs2 - ro)
    z_start = zo

    bf_start = DirectBField()
    for _ in 1:10 # Refine starting R using Newton's method
        direct_get_bfield!(bf_start, r_start, z_start, raw_profile.psi_in, raw_profile.sq_in, raw_profile.psio, derivs=1)
        dr = (psi0_target - bf_start.psi) / bf_start.psir
        r_start += dr
        if abs(dr) <= 1e-12 * r_start; break; end
    end
    direct_get_bfield!(bf_start, r_start, z_start, raw_profile.psi_in, raw_profile.sq_in, raw_profile.psio, derivs=2)
    psi0_actual = bf_start.psi

    # 2. Set up and solve the ODE for field line following
    u0 = zeros(Float64, 4)
    #[1]:∫(dl/Bp)
    #[2]: rfac (radial distance from magnetic axis)
    #[3]: ∫(dl/(R²Bp)) 
    #[4]: ∫(jac*dl/Bp) 
    u0[2] = sqrt((r_start - ro)^2 + (z_start - zo)^2) # Initial rfac
    tspan = (0.0, 2.0 * pi)

    equil_input = raw_profile.config.control

    params = FieldLineDerivParams(ro, zo, raw_profile.psi_in, raw_profile.sq_in, raw_profile.psio,
    equil_input.power_bp, equil_input.power_b, equil_input.power_r)

    # Use a callback to refine the solution at each step to stay on the flux surface
    saved_values = Vector{Vector{Float64}}()
    function refine_and_save_affect!(integrator)
        rfac_refined = direct_refine(integrator.u[2], integrator.t, psi0_actual, params)
        integrator.u[2] = rfac_refined
        push!(saved_values, [integrator.t; integrator.u])
    end
    callback = DiscreteCallback((u,t,i) -> true, refine_and_save_affect!, save_positions=(false,false))

    # Add the initial state
    push!(saved_values, [0.0; u0])

    prob = ODEProblem(direct_fl_der!, u0, tspan, params)
    sol = solve(prob, Tsit5(), callback=callback, reltol=1e-6, abstol=1e-8, dt=2*pi/200, adaptive=true)

    if sol.retcode != :Success && sol.retcode != :Terminated
        error("ODE integration failed for ipsi=$ipsi with code: $(sol.retcode)")
    end

    # 3. Finalize output
    y_out = hcat(saved_values...)' # Convert vector of vectors to a matrix
    return y_out, bf_start
end


"""
    direct_fl_der!(dy, y, params, eta)

The derivative function for the field-line integration ODE. This is passed to
the `DifferentialEquations.jl` solver.

## Arguments:
- `dy`: The derivative vector (output, modified in-place).
- `y`: The state vector `[∫(dl/Bp), rfac, ∫(dl/(R²Bp)), ∫(jac*dl/Bp)]`.
- `params`: A `FieldLineDerivParams` struct with all necessary parameters.
- `eta`: The independent variable (geometric angle `η`).
"""
function direct_fl_der!(dy, y, params::FieldLineDerivParams, eta)
    cos_eta, sin_eta = cos(eta), sin(eta)
    r = params.ro + y[2] * cos_eta
    z = params.zo + y[2] * sin_eta

    bf_temp = DirectBField()
    direct_get_bfield!(bf_temp, r, z, params.psi_in, params.sq_in, params.psio, derivs=1)

    bp_sq = bf_temp.br^2 + bf_temp.bz^2
    bp = bp_sq > 1e-28 ? sqrt(bp_sq) : 1e-14
    bt = bf_temp.f / r
    b_sq = bp_sq + bt^2
    b = sqrt(b_sq)
    jacfac = (bp^params.power_bp) * (b^params.power_b) / (r^params.power_r)

    # Denominator for d(l_pol)/d(eta) = rfac |B_pol|/denominator
    denominator = bf_temp.bz * cos_eta - bf_temp.br * sin_eta
    if abs(denominator) < 1e-14
        fill!(dy, 1e20) # Return large derivatives for solver to handle stiffness
        @warn "Denominator in direct_fl_der! near zero at eta=$eta."
        return
    end

    # dy/d(eta)
    dy[1] = y[2] / denominator  # d/dη [∫(dl/Bp)] = 1/|B_P| dl/d(eta) = rfac/denominator
    dy[2] = dy[1] * (bf_temp.br * cos_eta + bf_temp.bz * sin_eta)
    # d(rfac)/d(eta) = rfac/denom *( Br cos(eta) + Bz sin(eta) ) 
    dy[3] = dy[1] / (r^2)       # d/dη [∫(dl/(R²Bp))]
    dy[4] = dy[1] * jacfac      # d/dη [∫(jac*dl/Bp)]
    return
end


"""
    direct_refine(rfac, eta, psi0, params)

Refines the radial distance `rfac` at a given angle `eta` to ensure the
point lies exactly on the target flux surface `psi0`.

## Arguments:
- `rfac`: The current guess for the radial distance from the magnetic axis.
- `eta`: The geometric poloidal angle.
- `psi0`: The target `ψ` value for the flux surface.
- `params`: A `FieldLineDerivParams` struct.

## Returns:
- The refined `rfac` value.
"""
function direct_refine(rfac::Float64, eta::Float64, psi0::Float64, params::FieldLineDerivParams; max_iter::Int=10)::Float64
    cos_eta, sin_eta = cos(eta), sin(eta)
    bf_temp = DirectBField()

    for _ in 1:max_iter
        r_current = params.ro + rfac * cos_eta
        z_current = params.zo + rfac * sin_eta
        direct_get_bfield!(bf_temp, r_current, z_current, params.psi_in, params.sq_in, params.psio, derivs=1)
        dpsi = bf_temp.psi - psi0

        # Newton's method derivative: d(psi)/d(rfac)
        dpsi_drfac = bf_temp.psir * cos_eta + bf_temp.psiz * sin_eta
        if abs(dpsi_drfac) < 1e-14
            @warn "Refinement failed at eta=$eta: d(psi)/d(rfac) is zero."
            return rfac # Return current best guess
        end

        drfac = -dpsi / dpsi_drfac
        rfac += drfac

        if abs(dpsi) <= 1e-12 || abs(drfac) <= 1e-12 * abs(rfac)
            return rfac # Converged
        end
    end

    @warn "direct_refine did not converge after $max_iter iterations at eta=$eta."
    return rfac
end


"""
    equilibrium_solver(raw_profile)

The main driver for the direct equilibrium reconstruction. It orchestrates the entire
process from finding the magnetic axis to integrating along field lines and
constructing the final coordinate and physics quantity splines.

## Arguments:
- `raw_profile`: A `DirectRunInput` object containing the initial splines (`psi_in`, `sq_in`)
  and run parameters (`equil_input`).

## Returns:
- A `PlasmaEquilibrium` object containing the final, processed equilibrium data,
  including the profile spline (`sq`), the coordinate mapping spline (`rzphi`), and
  the physics quantity spline (`eqfun`).
"""
function equilibrium_solver(raw_profile::DirectRunInput)
    println("--- Starting Direct Equilibrium Processing ---")

    # 1. Unpack initial data
    equil_params = raw_profile.config.control
    sq_in = raw_profile.sq_in
    psi_in = raw_profile.psi_in
    psio = raw_profile.psio
    mtheta = equil_params.mtheta

    # 2. Setup the new output radial grid (`ψ_norm`)
    mpsi = equil_params.mpsi
    psilow = equil_params.psilow
    psihigh = equil_params.psihigh

    sq_x_nodes = zeros(Float64, mpsi + 1)
    if equil_params.grid_type == "ldp"
        for i in 0:mpsi
            x_param = Float64(i) / mpsi
            sq_x_nodes[i+1] = psilow + (psihigh - psilow) * sin(x_param * pi / 2.0)^2
        end
    else
        error("Unsupported grid_type: $(equil_params.grid_type)")
    end
    # need to add more grid type

    #2pi*F, mu0P, dvdpsi, q
    sq_fs_nodes = zeros(Float64, mpsi + 1, 4)

    ro_guess = (raw_profile.rmin + raw_profile.rmax) / 2.0
    zo_guess = (raw_profile.zmin + raw_profile.zmax) / 2.0



    # 3. Find key geometric positions and perform normalization
    ro, zo, rs1, rs2, psi_in_norm = direct_position(psi_in, sq_in, psio,
                                        ro_guess, zo_guess, # initial guess
                                       raw_profile.rmin, raw_profile.rmax)



    normalized_profile = DirectRunInput(
        raw_profile.config,
        raw_profile.sq_in,
        psi_in_norm, 
        raw_profile.rmin,
        raw_profile.rmax,
        raw_profile.zmin,
        raw_profile.zmax,
        raw_profile.psio
    )

    # 4. Main integration loop over flux surfaces
    local rzphi::Spl.BicubicSplineType
    local eqfun::Spl.BicubicSplineType
    local rzphi_fs_nodes 

    println("Starting loop over flux surfaces...")
    # Loop from edge inwards (index mpsi+1 down to 1)
    for ipsi in (mpsi+1):-1:1
        psi_norm_surf = sq_x_nodes[ipsi]
        @printf "--> Processing surface ipsi = %d / %d (ψ_norm = %.4f)\n" (ipsi-1) mpsi psi_norm_surf

        # a. Integrate along the field line for this surface
        y_out, bf_start = direct_fl_int(ipsi, psi_norm_surf, normalized_profile, ro, zo, rs2)
        #y_out
        #[1]: eta
        #[2]:∫(dl/Bp)
        #[3]: rfac (radial distance from magnetic axis)
        #[4]: ∫(dl/(R²Bp)) 
        #[5]: ∫(jac*dl/Bp) 

        # b. Process integration results into a temporary periodic spline `ff(θ_new)`
        theta_new_nodes = y_out[:, 5] ./ y_out[end, 5] #SFL angle θ_new
        ff_fs_nodes = hcat(
            y_out[:, 3].^2,                                            # 1: rfac²
            y_out[:, 1] / (2*pi) .- theta_new_nodes,                   # 2: η/(2π) - θ_new
            bf_start.f * (y_out[:, 4] .- theta_new_nodes .* y_out[end, 4]), # 3: Toroidal stream function term (term for calculate q)
            y_out[:, 2] ./ y_out[end, 2] .- theta_new_nodes            # 4: Jacobian-related term
        )
        ff = Spl.spline_setup(theta_new_nodes, ff_fs_nodes, bctype="periodic")

        # c. On first iteration, allocate the main output data array
        if ipsi == (mpsi+1)
            rzphi_fs_nodes = zeros(Float64, mpsi + 1, mtheta + 1, 4)
        end

        # d. Interpolate `ff` onto the uniform `theta` grid for `rzphi`
        rzphi_y_nodes = range(0.0, 1.0, length=mtheta + 1)
        for itheta in 1:(mtheta + 1)
            theta_val = rzphi_y_nodes[itheta]
            f, f1 = Spl.spline_eval(ff, theta_val, 1)

            rzphi_fs_nodes[ipsi, itheta, 1:3] = f[1:3]
            jac_term = (1.0 + f1[4]) * y_out[end, 2] * (2*pi) * psio
            rzphi_fs_nodes[ipsi, itheta, 4] = jac_term
        end

        # e. Store surface-averaged quantities for the `sq` spline
        sq_fs_nodes[ipsi, 1] = bf_start.f * (2pi) # 2pi*F
        sq_fs_nodes[ipsi, 2] = bf_start.p
        sq_fs_nodes[ipsi, 3] = y_out[end, 5] * (2pi) *psio # dV/d(psi)
        sq_fs_nodes[ipsi, 4] = y_out[end, 4] * bf_start.f / (2pi) # q-profile
    end
    println("...Loop over flux surfaces finished.")

    # 5. Finalize splines and perform q-profile revision if needed
    sq = Spl.spline_setup(sq_x_nodes, sq_fs_nodes, bctype=4)

    if equil_params.newq0 != 0.0
        println("Revising q-profile for newq0 = $(equil_params.newq0)...")
        f = Spl.spline_eval(sq, 0.0, 0)
        # q0_old = q(psi=0) = f[4] at x=0
        # f0_old = f[1] at x=0
        q0_old = f[4]
        f0_old = f[1]

        f0_fac_sq = f0_old^2 * ((equil_params.newq0 / q0_old)^2 - 1.0)

        for i in 1:(mpsi+1)
            f_current_sq = sq_fs_nodes[i, 1]^2
            ffac = sqrt(max(0.0, 1.0 + f0_fac_sq / f_current_sq)) * sign(equil_params.newq0/q0_old)
            sq_fs_nodes[i, 1] *= ffac      # F
            sq_fs_nodes[i, 4] *= ffac      # q
            rzphi_fs_nodes[i, :, 3] .*= ffac # Toroidal stream function
        end
        # Re-create the spline with the revised data
        sq = Spl.spline_setup(sq_x_nodes, sq_fs_nodes, bctype=4)
        println("...q-profile revision complete.")
    end

    # Create the final geometric spline `rzphi`. Periodic in theta (y-dimension)
    rzphi_y_nodes = range(0.0, 1.0, length=mtheta + 1)
    rzphi = Spl.bicube_setup(sq_x_nodes, collect(rzphi_y_nodes), rzphi_fs_nodes, bctypex=4, bctypey=2)
    println("Final geometric spline 'rzphi' is fitted.")

    # 6. Calculate final physics quantities (B-field, metric components, etc.)
    println("Calculating final physics quantities (B, g_ij)...")
    eqfun_fs_nodes = zeros(Float64, mpsi + 1, mtheta + 1, 3)
    v = zeros(Float64, 2, 3)

    for ipsi in 1:(mpsi+1), itheta in 1:(mtheta+1)
        psi_norm = sq_x_nodes[ipsi]
        theta_new = rzphi_y_nodes[itheta]

        f = Spl.spline_eval(sq, psi_norm, 0)
        q = f[4]
        f_val = f[1] 

        f, fx, fy = Spl.bicube_eval(rzphi, psi_norm, theta_new, 1)
        rfac_sq = max(0.0, f[1])
        rfac = sqrt(rfac_sq)
        eta = 2.0 * pi * (theta_new + f[2])
        r_coord = ro + rfac * cos(eta)
        jacfac = f[4]

 
        v[1,1] = (rfac > 0) ? fx[1] / (2.0 * rfac) : 0.0       # 1/(2rfac) * d(rfac)/d(psi_norm)
        v[1,2] = fx[2] * 2.0 * pi * rfac                       # 2pi*rfac * d(eta)/d(psi_norm)
        v[1,3] = fx[3] * r_coord                               # r_coord * d(phi_s)/d(psi_norm) 
        v[2,1] = (rfac > 0) ? fy[1] / (2.0 * rfac) : 0.0       # 1/(2rfac) d(rfac)/d(theta_new)
        v[2,2] = (1.0 + fy[2]) * 2.0 * pi * rfac               # 2pi*rfac * d(eta)/d(theta_new) 
        v[2,3] = fy[3] * r_coord                               # r_coord * d(phi_s)/d(theta_new)
        v33 = 2.0 * pi * r_coord                               # 2pi * r_coord 


        w11 = (1.0 + fy[2]) * (2.0*pi)^2 * rfac * r_coord / jacfac
        w12 = (jacfac*rfac != 0) ? -fy[1] * pi * r_coord / (rfac * jacfac) : 0.0
        delpsi_norm = sqrt(w11^2 + w12^2)

        b_sq = ((psio *2pi *delpsi_norm)^2 + (f_val/(2pi*r_coord))^2)
        eqfun_fs_nodes[ipsi, itheta, 1] = b_sq # B_total

        denom = jacfac * b_sq
        if abs(denom) > 1e-20
            # 2. Gyrokinetic coefficient C1 
            numerator_2 = dot(v[1,:], v[2,:]) + q * v33 * v[1,3] 
            eqfun_fs_nodes[ipsi, itheta, 2] = numerator_2 / denom
    
            # 3. Gyrokinetic coefficient C2 
            numerator_3 = v[2,3] * v33 + q * v33^2 
            eqfun_fs_nodes[ipsi, itheta, 3] = numerator_3 / denom
        else
            eqfun_fs_nodes[ipsi, itheta, 2] = 0.0
            eqfun_fs_nodes[ipsi, itheta, 3] = 0.0
        end
    end
    println("...done.")

    eqfun = Spl.bicube_setup(sq_x_nodes, collect(rzphi_y_nodes), eqfun_fs_nodes, bctypex=4, bctypey=2)

    println("--- Direct Equilibrium Processing Finished ---")

    return PlasmaEquilibrium(raw_profile.config,EquilibriumParameters(), sq, rzphi, eqfun, ro, zo, psio)
end

