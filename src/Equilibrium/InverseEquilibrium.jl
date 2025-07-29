#=
This contains the logic for transforming an 'inverse' equilibrium
representation into a straight-fieldline coordinate system and associated quantities.
It uses bicubic splines for the geometry and supports several surface grid types.
=#


# --- Lagrange Polynomial Extrapolation ---
"""
    inverse_extrap(xx, ff, x)

Multi-point polynomial extrapolation of function values `ff` at `xx` to value at `x`.
Equivalent to Fortran's Lagrange extrapolation in inverse_extrap().
Assumes `xx` and `ff` are matrices with the first dimension being control points.

Returns a vector.
"""
function inverse_extrap(xx::AbstractMatrix, ff::AbstractMatrix, x::Number)
    m = size(xx,1)
    n = size(ff,2)
    f = zeros(eltype(ff), n)
    for i in 1:m
        term = view(ff, i, :)
        result = copy(term)
        for j in 1:m
            if j == i; continue; end
            denom = xx[i,:] .- xx[j,:]
            result .= result .* (x .- xx[j,:]) ./ denom
        end
        f .+= result
    end
    return f
end

# --- Diagnostics/Output Stubs (to be implemented) ---
function inverse_output(;diagnostics=false)
    # Add ASCII/Binary output, diagnostics etc. here if needed
    @info "inverse_output called (output stubs not implemented)"
end

"""
    equilibrium_solver(input::InverseRunInput)

Main routine for inverse equilibrium analysis.

Arguments:
    - input: InverseRunInput struct containing all splines and parameters.

Returns:
    - PlasmaEquilibrium object (as in direct_run).
"""
function equilibrium_solver(input::InverseRunInput)
    equil_params = deepcopy(input.equil_input)
    sq_in = deepcopy(input.sq_in)
    println(sq_in)
    rz_in = deepcopy(input.rz_in)
    ro = input.ro
    zo = input.zo
    psio = input.psio

    mpsi    = equil_params.mpsi
    psihigh = equil_params.psihigh
    psilow  = equil_params.psilow
    mtheta  = equil_params.mtheta
    grid_type = equil_params.grid_type
    power_bp = equil_params.power_bp
    power_b  = equil_params.power_b
    power_r  = equil_params.power_r


    # ---- Spline Preprocessing ----
    #!!!!!!NEED TO CHANGE THIS!!!!!! #sq_in.fs[:,4] = sqrt.(sq_in.xs)
    for name in fieldnames(typeof(sq_in))
        println("sq_in.$name = ", getfield(sq_in, name))
    end    
    #sq_in.name  = "  sq  "
    #sq_in.title = ["psifac", "  f   ", "mu0 p ", "  q   ", " rho  "]
    #Spl.spline_fit!(sq_in, "extrap")

    # Prepare rz_in grid
    rz_in.xs = sq_in.xs
    npsi = size(rz_in.fs, 1)
    ntheta = size(rz_in.fs, 2)
    rz_in.ys = range(0.0, stop=1.0, length=ntheta)

    rz_in.name = "  rz  "
    rz_in.xtitle = "psifac"
    rz_in.ytitle = "theta "
    rz_in.title = ["  r   ", "  z   "]

    inverse_output()  # Diagnostics (by default does nothing)

    # ---- Coordinates: Cartesian to Polar Transform ----
    x = rz_in.fs[:,:,1] .- ro   # R - Raxis
    y = rz_in.fs[:,:,2] .- zo   # Z - Zaxis
    r2 = x.^2 .+ y.^2

    # Calculate angle in units of 2π, care at r2==0 (origin)
    deta = zeros(size(r2))
    for ipsi in 1:npsi, itheta in 1:ntheta
        if r2[ipsi, itheta] == 0.0
            deta[ipsi, itheta] = 0.0
        else
            deta[ipsi, itheta] = atan(y[ipsi, itheta], x[ipsi, itheta]) / (2π)
        end
    end
    # Unwrapping: correct jumps over the periodicity boundary
    for ipsi in 1:npsi
        for itheta in 2:ntheta
            if deta[ipsi, itheta] - deta[ipsi, itheta-1] > 0.5
                deta[ipsi, itheta] -= 1.0
            elseif deta[ipsi, itheta] - deta[ipsi, itheta-1] < -0.5
                deta[ipsi, itheta] += 1.0
            end
        end
        for itheta in 1:ntheta
            if r2[ipsi, itheta] > 0
                deta[ipsi, itheta] -= rz_in.ys[itheta]
            end
        end
    end
    # Extrapolate deta at origin (first row) from next 3 rows
    me = 3
    if npsi >= me
        deta[1, :] .= inverse_extrap(r2[2:me+1, :], deta[2:me+1, :], 0.0)
    end

    # Store new coordinates back into rz_in for next phase
    rz_in.fs[:,:,1] .= r2
    rz_in.fs[:,:,2] .= deta

    Spl.bicube_fit!(rz_in, "extrap", "periodic")

    # ---- Set up Output Surface Grid ----
    sq_nodes = zeros(mpsi+1)
    if grid_type == "ldp"
        for i in 0:mpsi
            x_param = Float64(i) / mpsi
            sq_nodes[i+1] = psilow + (psihigh - psilow) * sin(x_param * π/2)^2
        end
    else
        error("Unsupported grid_type: $grid_type")
    end

    # Allocate output arrays
    if mtheta == 0; mtheta = ntheta; end
    rzphi_fs_nodes = zeros(mpsi+1, mtheta+1, 4)
    eqfun_fs_nodes = zeros(mpsi+1, mtheta+1, 3)
    sq_fs_nodes = zeros(mpsi+1, 4)

    rzphi_y_nodes = range(0.0, 1.0, length=mtheta+1)
    rzphi_x_nodes = sq_nodes

    # Local spline
    spl = Spl.alloc_spline(mtheta+1, 5)
    # ---- Loop Over Surfaces ----
    for ipsi in 1:(mpsi+1)
        psifac = rzphi_x_nodes[ipsi]
        Spl.spline_eval!(sq_in, psifac, 0)
        spl.xs .= rzphi_y_nodes

        # Interpolate rz_in and fill spl.fs
        for itheta in 1:(mtheta+1)
            theta = rzphi_y_nodes[itheta]
            Spl.bicube_eval!(rz_in, psifac, theta, 1)
            Spl.spline_eval!(sq_in, psifac, 0)
            if rz_in.f[1] < 0
                error("Invalid extrapolation near axis, rerun with larger psilow")
            end
            rfac = sqrt(rz_in.f[1])
            η = 2π * (theta + rz_in.f[2])
            r = ro + rfac * cos(η)
            jacfac = rz_in.fx[1] * (1+rz_in.fy[2]) - rz_in.fy[1]*rz_in.fx[2]
            w11 = (1+rz_in.fy[2]) * (2π)^2 * rfac / jacfac
            w12 = -rz_in.fy[1] * π / (rfac * jacfac)
            bp = psio * sqrt(w11^2 + w12^2) / r
            bt = sq_in.f[1] / r
            b = sqrt(bp^2 + bt^2)
            spl.fs[itheta, 1] = rz_in.f[1]
            spl.fs[itheta, 2] = rz_in.f[2]
            spl.fs[itheta, 3] = r * jacfac
            spl.fs[itheta, 4] = spl.fs[itheta,3] / (r*r)
            spl.fs[itheta, 5] = spl.fs[itheta,3] * bp^power_bp * b^power_b / r^power_r
        end

        Spl.spline_fit!(spl, "periodic")
        Spl.spline_int!(spl)

        # Normalize arc length and adjust
        spl.xs .= spl.fsi[:,5] ./ spl.fsi[end,5]
        spl.fs[:,2] .= spl.fs[:,2] .+ rzphi_y_nodes .- spl.xs
        spl.fs[:,4] .= (spl.fs[:,3] ./ spl.fsi[end,3]) ./ (spl.fs[:,5]./spl.fsi[end,5]) .* spl.fsi[end,3]*2π*π
        spl.fs[:,3] .= sq_in.f[1] * π / psio .* (spl.fsi[:,4] .- spl.fsi[end,4] .* spl.xs)
        Spl.spline_fit!(spl, "periodic")

        for itheta in 1:(mtheta+1)
            theta = rzphi_y_nodes[itheta]
            Spl.spline_eval!(spl, theta, 0)
            rzphi_fs_nodes[ipsi, itheta, :] .= spl.f[1:4]
        end

        sq_fs_nodes[ipsi,1] = sq_in.f[1] * 2π
        sq_fs_nodes[ipsi,2] = sq_in.f[2]
        sq_fs_nodes[ipsi,3] = spl.fsi[end,3] * 2π * π
        sq_fs_nodes[ipsi,4] = spl.fsi[end,4] * sq_fs_nodes[ipsi,1] / (2*2π*psio)
    end

    # Fit output spline(s)
    sq = Spl.spline_setup(sq_nodes, sq_fs_nodes, bctype=4)
    Spl.spline_fit!(sq, "extrap")

    # Q profile revision
    f = Spl.spline_eval(sq, 0.0, 0)
    q0 = f[4]
    if equil_params.newq0 == -1.0; equil_params.newq0 = -q0; end
    if equil_params.newq0 != 0.0
        f0 = f[1]
        f0_fac_sq = f0^2 * ((equil_params.newq0 / q0)^2 - 1.0)
        for i in 1:(mpsi+1)
            ffac = sqrt(1 + f0_fac_sq / sq_fs_nodes[i,1]^2) * sign(equil_params.newq0)
            sq_fs_nodes[i,1] *= ffac
            sq_fs_nodes[i,4] *= ffac
            rzphi_fs_nodes[i, :, 3] .*= ffac
        end
        sq = Spl.spline_setup(sq_nodes, sq_fs_nodes, bctype=4)
        Spl.spline_fit!(sq, "extrap")
    end

    rzphi = Spl.bicube_setup(rzphi_x_nodes, collect(rzphi_y_nodes), rzphi_fs_nodes, bctypex=4, bctypey=2)
    Spl.bicube_fit!(rzphi, "extrap", "periodic")

    # --- Compute eqfun (final physics quantities etc.) ---
    v = zeros(3,3)
    for ipsi in 1:(mpsi+1), itheta in 1:(mtheta+1)
        psi_norm = rzphi_x_nodes[ipsi]
        theta_new = rzphi_y_nodes[itheta]
        f = Spl.spline_eval(sq, psi_norm, 0)
        q = f[4]
        f, fx, fy = Spl.bicube_eval(rzphi, psi_norm, theta_new, 1)
        rfac = sqrt(f[1])
        eta = 2π * (theta_new + f[2])
        r = ro + rfac * cos(eta)
        jacfac = f[4]
        v[1,1] = fx[1]/(2*rfac)
        v[1,2] = fx[2]*2π*rfac
        v[1,3] = fx[3]*r
        v[2,1] = fy[1]/(2*rfac)
        v[2,2] = (1+fy[2])*2π*rfac
        v[2,3] = fy[3]*r
        v[3,3] = 2π*r
        w11 = (1+fy[2])* (2π)^2 * rfac * r / jacfac
        w12 = -fy[1]*π*r / (rfac*jacfac)
        delpsi = sqrt(w11^2 + w12^2)
        eqfun_fs_nodes[ipsi,itheta,1] = sqrt(((2π*psio*delpsi)^2 + f[1]^2)/(2π*r)^2)
        eqfun_fs_nodes[ipsi,itheta,2] = (dot(v[1,:], v[2,:]) + q*v[3,3]*v[1,3]) / (jacfac*eqfun_fs_nodes[ipsi,itheta,1]^2)
        eqfun_fs_nodes[ipsi,itheta,3] = (v[2,3]*v[3,3] + q*v[3,3]*v[3,3]) / (jacfac*eqfun_fs_nodes[ipsi,itheta,1]^2)
    end
    eqfun = Spl.bicube_setup(rzphi_x_nodes, collect(rzphi_y_nodes), eqfun_fs_nodes, bctypex=4, bctypey=2)
    Spl.bicube_fit!(eqfun, "extrap", "periodic")

    return PlasmaEquilibrium(equil_params, sq, rzphi, eqfun, ro, zo, psio)
end


