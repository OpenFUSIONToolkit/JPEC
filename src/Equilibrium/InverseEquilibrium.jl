"""
Converts inverse equilibrium to straight-fieldline coordinates. Based on inverse.f
Of the entries we need to return: PlasmaEquilibrium(equil_params, sq_out, rzphi_out, eqfun_out, ro, zo, psio),
we only need to generate: sq_out, rzphi_out, eqfun_out. This is because we pass in InverseRunInput(equil_in, 
sq_in, rz_in, ro, zo, psio).

"""

# sq_in # spline_type object, sq_in.xs # 1D array: radial coordinate (normalized poloidal flux), sq_in.fs,
# 2D array: profiles vs. radius (F-function, pressure, q-profile), rz_in, # bicube_type object, rz_in.fs,
# 3D array: R and Z coordinates on the input grid, rz_in.mx, # Integer: Number of radial points in the input grid,
# rz_in.my, # Integer: Number of poloidal points in the input grid

function equilibrium_solver(raw_profile::InverseRunInput)
    println("--- Starting Inverse Equilibrium Processing ---")

    # 1.1 Unpack initial data
    equil_params = raw_profile.config.control # EquilInput object
    sq_in = raw_profile.sq_in # 1D spline input profile (e.g. F*Bt, Pressure, q)
    rz_in = raw_profile.rz_in # 2D bicubic spline input for (R,Z) geometry
    ro = raw_profile.ro # R magnetic axis location
    zo = raw_profile.zo # Z magnetic axis location
    psio = raw_profile.psio # Total flux difference |psi_axis - psi_boundary|
    newq0 = equil_params.newq0 # New target on-axis safety factor (q0)
    power_bp = equil_params.power_bp # Exponent for poloidal field weighting
    power_b = equil_params.power_b # Exponent for total magnetic field weighting
    power_r = equil_params.power_r # Exponent for major radius weighting

    # 1.2 Set flags
    power_flag = false #equil_params.power_flag
    direct_flag = false #equil_params.direct_flag
    out_eq_1d = false # Enable ASCII output for 1D input data
    bin_eq_1d = false # Enable binary output for 1D input data
    out_eq_2d = false # Enable ASCII output for 2D input data
    bin_eq_2d = false # Enable binary output for 2D input data
    input_only = false # If true, stop after writing input data
    interp = false # If true, interpolate data before writing output

    # 1.3 Other unimplemented parameters
    sp_pfac = 1 # (NEEDS IMPLEMENTATION) Packing factor for "mypack" grid_type
    jaceq = 1 # (NEEDS IMPLEMENTATION), DIMENSION(:,:), INTENT(IN): Jacobian from CHEASE, for inverse_chease4_run(jaceq)

    # 2. Setup the new output radial grid (`ψ_norm`)
    mtheta = equil_params.mtheta # Number of poloidal grid points
    grid_type = equil_params.grid_type # New radial grid type
    mpsi = equil_params.mpsi # Number of radial grid points
    psilow = equil_params.psilow # Lower bound of normalized poloidal flux for the new grid
    psihigh = equil_params.psihigh # Upper bound of normalized poloidal flux for the new grid

    sq_x_nodes = zeros(Float64, mpsi + 1)

    if grid_type == "ldp"
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

    # (Should make this function so both inverse and direct can use it instead of copy-pasting)
    # Helper function for robust Newton-Raphson search with restarts
    max_newton_iter = 200
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

    #rs1 = find_separatrix_crossing(ro, minimum(rz_in.fs[:, :, 1]), "inboard")
    #rs2 = find_separatrix_crossing(ro, maximum(rz_in.fs[:, :, 1]), "outboard")                            
                                
    normalized_profile = InverseRunInput(
        raw_profile.config,
        raw_profile.sq_in,
        raw_profile.rz_in,
        ro,
        zo,
        psio
    )

    # 3. Main integration loop over flux surfaces
    local rzphi::Spl.BicubicSplineType
    local eqfun::Spl.BicubicSplineType
    local rzphi_fs_nodes 

    println("Starting loop over flux surfaces...")
    # Loop from edge inwards (index mpsi+1 down to 1)
    # for ipsi in (mpsi+1):-1:1
    #     psi_norm_surf = sq_x_nodes[ipsi]
    #     @printf "--> Processing surface ipsi = %d / %d (ψ_norm = %.4f)\n" (ipsi-1) mpsi psi_norm_surf

    #     # a. Integrate along the field line for this surface
    #     y_out, bf_start = direct_fl_int(ipsi, psi_norm_surf, normalized_profile, ro, zo, rs2)
    #     #y_out
    #     #[1]: eta
    #     #[2]:∫(dl/Bp)
    #     #[3]: rfac (radial distance from magnetic axis)
    #     #[4]: ∫(dl/(R²Bp)) 
    #     #[5]: ∫(jac*dl/Bp) 

    #     # b. Process integration results into a temporary periodic spline `ff(θ_new)`
    #     theta_new_nodes = y_out[:, 2] ./ y_out[end, 2] #SFL angle θ_new
    #     ff_fs_nodes = hcat(
    #         y_out[:, 3].^2,                                            # 1: rfac²
    #         y_out[:, 1] / (2*pi) .- theta_new_nodes,                   # 2: η/(2π) - θ_new
    #         bf_start.f * (y_out[:, 4] .- theta_new_nodes .* y_out[end, 4]), # 3: Toroidal stream function term (term for calculate q)
    #         y_out[:, 5] ./ y_out[end, 5] .- theta_new_nodes            # 4: Jacobian-related term
    #     )
    #     ff = Spl.spline_setup(theta_new_nodes, ff_fs_nodes, bctype="periodic")

    #     # c. On first iteration, allocate the main output data array
    #     if ipsi == (mpsi+1)
    #         rzphi_fs_nodes = zeros(Float64, mpsi + 1, mtheta + 1, 4)
    #     end

    #     # d. Interpolate `ff` onto the uniform `theta` grid for `rzphi`
    #     rzphi_y_nodes = range(0.0, 1.0, length=mtheta + 1)
    #     for itheta in 1:(mtheta + 1)
    #         theta_val = rzphi_y_nodes[itheta]
    #         f, f1 = Spl.spline_eval(ff, theta_val, 1)

    #         rzphi_fs_nodes[ipsi, itheta, 1:3] = f[1:3]
    #         jac_term = (1.0 + f1[4]) * y_out[end, 5] * psio
    #         rzphi_fs_nodes[ipsi, itheta, 4] = jac_term
    #     end

    #     # e. Store surface-averaged quantities for the `sq` spline
    #     sq_fs_nodes[ipsi, 1] = bf_start.f * (2pi) # 2pi*F
    #     sq_fs_nodes[ipsi, 2] = bf_start.p
    #     sq_fs_nodes[ipsi, 3] = y_out[end, 5] * (2pi) *psio # dV/d(psi)
    #     sq_fs_nodes[ipsi, 4] = y_out[end, 4] * bf_start.f / (2pi) # q-profile
    # end
    # println("...Loop over flux surfaces finished.")


    # This is temporary code to replace the above loop (needs proper inverse implementation)
    rzphi_fs_nodes = zeros(Float64, mpsi + 1, mtheta + 1, 4) 
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

    println("--- Inverse Equilibrium Processing Complete ---")

    return PlasmaEquilibrium(raw_profile.config, sq, rzphi, eqfun, ro, zo, psio)

end