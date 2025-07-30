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


function update_rz_in_from_sq_in(rz_in::Spl.BicubicSplineType, sq_in::Spl.RealSplineType, ntheta::Int)
    # Create a new BicubicSplineType copying all fields, 
    # but replace _xs and _ys with updated arrays
    return Spl.BicubicSplineType(
        rz_in.handle,
        copy(sq_in._xs),                           # Replace _xs with copy from sq_in
        collect(range(0.0, stop=1.0, length=ntheta)),  # Replace _ys with range vector
        copy(rz_in._fs),                          # Copy fs arrays as needed
        rz_in.mx,
        ntheta,
        rz_in.nqty,
        rz_in.bctypex,
        rz_in.bctypey,
        copy(rz_in._fsx),
        copy(rz_in._fsy),
        copy(rz_in._fsxy)
    )
end


function bicube_fit!(
    bcs::Spl.BicubicSplineType,
    endmode1::String,
    endmode2::String;
    xperiodic::Bool = (endmode1 == "periodic"),
    yperiodic::Bool = (endmode2 == "periodic")
)
    fs = bcs._fs
    fsx = bcs._fsx
    fsy = bcs._fsy
    fsxy = bcs._fsxy

    mx, my, nqty = bcs.mx, bcs.my, bcs.nqty

    # --- Apply periodicity (boundary copy)
    if xperiodic
        fs[mx, :, :] .= fs[1, :, :]
    end
    if yperiodic
        fs[:, my, :] .= fs[:, 1, :]
    end

    # --- Y-derivatives
    spl = Spl.RealSplineType(
    Ptr{Cvoid}(0),
    Float64[],                # _xs
    zeros(0, 0),              # _fs
    0,                        # mx
    0,                        # nqty
    0,                        # bctype
    zeros(0, 0),              # _fsi
    zeros(0, 0)               # _fs1
)

    spline_alloc!(spl, my, mx + 1)
    spl._xs = bcs._ys
    for iqty in 1:nqty
        spl._fs .= transpose(fs[:, :, iqty])
        spline_fit!(spl, endmode2)
        fsy[:, :, iqty] .= transpose(spl._fs1)
    end
    spline_dealloc!(spl)

    # --- X-derivatives
    spline_alloc!(spl, mx, my + 1)
    spl._xs = bcs._xs
    for iqty in 1:nqty
        spl._fs .= fs[:, :, iqty]
        spline_fit!(spl, endmode1)
        fsx[:, :, iqty] .= spl._fs1
    end

    # --- Mixed derivatives
    for iqty in 1:nqty
        spl._fs .= fsy[:, :, iqty]
        spline_fit!(spl, endmode1)
        fsxy[:, :, iqty] .= spl._fs1
    end
    spline_dealloc!(spl)

    return nothing
end





# """
#     equilibrium_solver(input::InverseRunInput)

# Main routine for inverse equilibrium analysis.

# Arguments:
#     - input: InverseRunInput struct containing all splines and parameters.

# Returns:
#     - PlasmaEquilibrium object (as in direct_run).
# """

# function inverse_run(input::InverseRunInput)
#     equil_params = deepcopy(input.equil_input)
#     sq_in = deepcopy(input.sq_in)
#     println(sq_in)
#     rz_in = deepcopy(input.rz_in)
#     ro = input.ro
#     zo = input.zo
#     psio = input.psio

#     mpsi    = equil_params.mpsi
#     psihigh = equil_params.psihigh
#     psilow  = equil_params.psilow
#     mtheta  = equil_params.mtheta
#     grid_type = equil_params.grid_type
#     power_bp = equil_params.power_bp
#     power_b  = equil_params.power_b
#     power_r  = equil_params.power_r

#     # ---- Spline Preprocessing ----
#     # Here you may want to perform operations on sq_in.fs or related,
#     # but you must do so by creating new instances or mutating arrays inside fields
#     # (assuming arrays inside fields are mutable).
#     # For example, to replace sq_in._fs[:,4] with sqrt(sq_in._xs):
#     # sq_in._fs[:,4] .= sqrt.(sq_in._xs)

#     for name in fieldnames(typeof(sq_in))
#         println("sq_in.$name = ", getfield(sq_in, name))
#     end    

#     # Prepare rz_in grid with updated _xs and _ys arrays:
#     ntheta = size(rz_in._fs, 2)
#     rz_in = update_rz_in_from_sq_in(rz_in, sq_in, ntheta)

#     # rz_in.name = "  rz  "
#     # rz_in.xtitle = "psifac"
#     # rz_in.ytitle = "theta "
#     # rz_in.title = ["  r   ", "  z   "]

#     inverse_output()  # Diagnostics (by default does nothing)

#     println("rz_in._fs size = ", size(rz_in._fs))
#     (npsi, ntheta, _) = size(rz_in._fs)
#     if size(rz_in._fs, 3) < 2
#         rz_in._fs = zeros(npsi, ntheta, 2)
#     end

#     # Compute Cartesian coordinates relative to axis
#     x = rz_in._fs[:,:,1] .- ro   # This currently might be zeros, so maybe you want original R data here?
#     y = rz_in._fs[:,:,2] .- zo   # Same as above
#     r2 = x.^2 .+ y.^2

#     deta = zeros(size(r2))
#     for ipsi in 1:npsi, itheta in 1:ntheta
#         if r2[ipsi, itheta] == 0.0
#             deta[ipsi, itheta] = 0.0
#         else
#             deta[ipsi, itheta] = atan(y[ipsi, itheta], x[ipsi, itheta]) / (2π)
#         end
#     end

#     # Unwrapping and further deta corrections here...

#     # Now store back into rz_in._fs
#     rz_in._fs[:,:,1] .= r2
#     rz_in._fs[:,:,2] .= deta


#     # Unwrap deta to avoid discontinuities
#     for ipsi in 1:size(deta, 1)
#         for itheta in 2:size(deta, 2)
#             diff = deta[ipsi, itheta] - deta[ipsi, itheta-1]
#             if diff > 0.5
#                 deta[ipsi, itheta] -= 1.0
#             elseif diff < -0.5
#                 deta[ipsi, itheta] += 1.0
#             end
#         end
#         for itheta in 1:size(deta, 2)
#             if r2[ipsi, itheta] > 0
#                 deta[ipsi, itheta] -= rz_in._ys[itheta]
#             end
#         end
#     end

#     # Extrapolate deta at origin from next rows
#     me = 3
#     if size(deta, 1) >= me+1
#         deta[1, :] .= inverse_extrap(r2[2:me+1, :], deta[2:me+1, :], 0.0)
#     end

#     # Store new coordinates back into rz_in for next phase
#     rz_in._fs[:,:,1] .= r2
#     rz_in._fs[:,:,2] .= deta

#     bicube_fit!(rz_in, "extrap", "periodic")

#     # ---- Set up Output Surface Grid ----
#     sq_nodes = zeros(mpsi+1)
#     if grid_type == "ldp"
#         for i in 0:mpsi
#             x_param = Float64(i) / mpsi
#             sq_nodes[i+1] = psilow + (psihigh - psilow) * sin(x_param * π/2)^2
#         end
#     else
#         error("Unsupported grid_type: $grid_type")
#     end

#     # Allocate output arrays
#     if mtheta == 0; mtheta = ntheta; end
#     rzphi_fs_nodes = zeros(mpsi+1, mtheta+1, 4)
#     eqfun_fs_nodes = zeros(mpsi+1, mtheta+1, 3)
#     sq_fs_nodes = zeros(mpsi+1, 4)

#     rzphi_y_nodes = range(0.0, 1.0, length=mtheta+1)
#     rzphi_x_nodes = sq_nodes

#     # Local spline
#     spl = Spl.alloc_spline(mtheta+1, 5)

#     # ---- Loop Over Surfaces ----
#     for ipsi in 1:(mpsi+1)
#         psifac = rzphi_x_nodes[ipsi]
#         Spl.spline_eval!(sq_in, psifac, 0)
#         spl.xs .= rzphi_y_nodes

#         # Interpolate rz_in and fill spl.fs
#         for itheta in 1:(mtheta+1)
#             theta = rzphi_y_nodes[itheta]
#             Spl.bicube_eval!(rz_in, psifac, theta, 1)
#             Spl.spline_eval!(sq_in, psifac, 0)
#             if rz_in.f[1] < 0
#                 error("Invalid extrapolation near axis, rerun with larger psilow")
#             end
#             rfac = sqrt(rz_in.f[1])
#             η = 2π * (theta + rz_in.f[2])
#             r = ro + rfac * cos(η)
#             jacfac = rz_in.fx[1] * (1+rz_in.fy[2]) - rz_in.fy[1]*rz_in.fx[2]
#             w11 = (1+rz_in.fy[2]) * (2π)^2 * rfac / jacfac
#             w12 = -rz_in.fy[1] * π / (rfac * jacfac)
#             bp = psio * sqrt(w11^2 + w12^2) / r
#             bt = sq_in.f[1] / r
#             b = sqrt(bp^2 + bt^2)
#             spl.fs[itheta, 1] = rz_in.f[1]
#             spl.fs[itheta, 2] = rz_in.f[2]
#             spl.fs[itheta, 3] = r * jacfac
#             spl.fs[itheta, 4] = spl.fs[itheta,3] / (r*r)
#             spl.fs[itheta, 5] = spl.fs[itheta,3] * bp^power_bp * b^power_b / r^power_r
#         end

#         Spl.spline_fit!(spl, "periodic")
#         Spl.spline_int!(spl)

#         # Normalize arc length and adjust
#         spl.xs .= spl.fsi[:,5] ./ spl.fsi[end,5]
#         spl.fs[:,2] .= spl.fs[:,2] .+ rzphi_y_nodes .- spl.xs
#         spl.fs[:,4] .= (spl.fs[:,3] ./ spl.fsi[end,3]) ./ (spl.fs[:,5]./spl.fsi[end,5]) .* spl.fsi[end,3]*2π*π
#         spl.fs[:,3] .= sq_in.f[1] * π / psio .* (spl.fsi[:,4] .- spl.fsi[end,4] .* spl.xs)
#         Spl.spline_fit!(spl, "periodic")

#         for itheta in 1:(mtheta+1)
#             theta = rzphi_y_nodes[itheta]
#             Spl.spline_eval!(spl, theta, 0)
#             rzphi_fs_nodes[ipsi, itheta, :] .= spl.f[1:4]
#         end

#         sq_fs_nodes[ipsi,1] = sq_in.f[1] * 2π
#         sq_fs_nodes[ipsi,2] = sq_in.f[2]
#         sq_fs_nodes[ipsi,3] = spl.fsi[end,3] * 2π * π
#         sq_fs_nodes[ipsi,4] = spl.fsi[end,4] * sq_fs_nodes[ipsi,1] / (2*2π*psio)
#     end

#     # Fit output spline(s)
#     sq = Spl.spline_setup(sq_nodes, sq_fs_nodes, bctype=4)
#     Spl.spline_fit!(sq, "extrap")

#     # Q profile revision
#     f = Spl.spline_eval(sq, 0.0, 0)
#     q0 = f[4]
#     if equil_params.newq0 == -1.0; equil_params.newq0 = -q0; end
#     if equil_params.newq0 != 0.0
#         f0 = f[1]
#         f0_fac_sq = f0^2 * ((equil_params.newq0 / q0)^2 - 1.0)
#         for i in 1:(mpsi+1)
#             ffac = sqrt(1 + f0_fac_sq / sq_fs_nodes[i,1]^2) * sign(equil_params.newq0)
#             sq_fs_nodes[i,1] *= ffac
#             sq_fs_nodes[i,4] *= ffac
#             rzphi_fs_nodes[i, :, 3] .*= ffac
#         end
#         sq = Spl.spline_setup(sq_nodes, sq_fs_nodes, bctype=4)
#         Spl.spline_fit!(sq, "extrap")
#     end

#     rzphi = Spl.bicube_setup(rzphi_x_nodes, collect(rzphi_y_nodes), rzphi_fs_nodes, bctypex=4, bctypey=2)
#     Spl.bicube_fit!(rzphi, "extrap", "periodic")

#     # --- Compute eqfun (final physics quantities etc.) ---
#     v = zeros(3,3)
#     for ipsi in 1:(mpsi+1), itheta in 1:(mtheta+1)
#         psi_norm = rzphi_x_nodes[ipsi]
#         theta_new = rzphi_y_nodes[itheta]
#         f = Spl.spline_eval(sq, psi_norm, 0)
#         q = f[4]
#         f, fx, fy = Spl.bicube_eval(rzphi, psi_norm, theta_new, 1)
#         rfac = sqrt(f[1])
#         eta = 2π * (theta_new + f[2])
#         r = ro + rfac * cos(eta)
#         jacfac = f[4]
#         v[1,1] = fx[1]/(2*rfac)
#         v[1,2] = fx[2]*2π*rfac
#         v[1,3] = fx[3]*r
#         v[2,1] = fy[1]/(2*rfac)
#         v[2,2] = (1+fy[2])*2π*rfac
#         v[2,3] = fy[3]*r
#         v[3,3] = 2π*r
#         w11 = (1+fy[2])*(2π)^2 * rfac * r / jacfac
#         w12 = -fy[1]*π*r / (rfac*jacfac)
#         delpsi = sqrt(w11^2 + w12^2)
#         eqfun_fs_nodes[ipsi,itheta,1] = sqrt(((2π*psio*delpsi)^2 + f[1]^2)/(2π*r)^2)
#         eqfun_fs_nodes[ipsi,itheta,2] = (dot(v[1,:], v[2,:]) + q*v[3,3]*v[1,3]) / (jacfac*eqfun_fs_nodes[ipsi,itheta,1]^2)
#         eqfun_fs_nodes[ipsi,itheta,3] = (v[2,3]*v[3,3] + q*v[3,3]*v[3,3]) / (jacfac*eqfun_fs_nodes[ipsi,itheta,1]^2)
#     end

#     eqfun = Spl.bicube_setup(rzphi_x_nodes, collect(rzphi_y_nodes), eqfun_fs_nodes, bctypex=4, bctypey=2)
#     Spl.bicube_fit!(eqfun, "extrap", "periodic")

#     return PlasmaEquilibrium(equil_params, sq, rzphi, eqfun, ro, zo, psio)
# end


"""
    equilibrium_solver(raw_profile::InverseRunInput) -> PlasmaEquilibrium

Reconstructs a magnetic equilibrium from a CHEASE2 inverse input profile.
This routine fits input data to splines, performs transformation to
straight-fieldline coordinates, and evaluates magnetic field metrics
on the prescribed flux grid.

# Arguments
- `raw_profile::InverseRunInput`: Parsed input from CHEASE2 equilibrium file.

# Returns
- `PlasmaEquilibrium`: A fully processed equilibrium structure for simulation.
"""
function equilibrium_solver(raw_profile::InverseRunInput)::PlasmaEquilibrium
    println("--- Starting Inverse Equilibrium Processing ---")

    # --- Unpack input
    rz_in, sq_in, equil_params = raw_profile.rz_in, raw_profile.sq_in, raw_profile.params
    mpsi, mtheta = equil_params.mpsi, equil_params.mtheta
    psilow, psihigh = equil_params.psilow, equil_params.psihigh

    # --- Prepare normalized psi grid
    psi_norm_grid = range(0.0, 1.0; length=mpsi + 1)
    psi_grid = psilow .+ psi_norm_grid .* (psihigh - psilow)

    # --- Fit splines to input profiles
    println("--> Fitting 1D splines from CHEASE profiles...")
    bctype = 4  # Not-a-Knot
    sq_fs_nodes = zeros(mpsi + 1, 5)
    for i in 1:5
        sq_fs_nodes[:, i] .= sq_in.fs[:, i]
    end
    sq_in = Spl.spline_setup(collect(psi_grid), sq_fs_nodes, bctype=bctype)

    # --- Prepare rzphi grid and transform rz_in to straight-fieldline
    println("--> Transforming rz_in to straight-fieldline coordinates...")
    rzphi_fs_nodes = zeros(mpsi + 1, mtheta + 1, 3)
    deta = zeros(mpsi + 1, mtheta + 1)

    for ipsi in 0:mpsi
        for itheta in 0:mtheta
            rzphi_fs_nodes[ipsi + 1, itheta + 1, 1] = rz_in.fs[ipsi + 1, itheta + 1, 1]
            rzphi_fs_nodes[ipsi + 1, itheta + 1, 2] = rz_in.fs[ipsi + 1, itheta + 1, 2]
            rzphi_fs_nodes[ipsi + 1, itheta + 1, 3] = rz_in.fs[ipsi + 1, itheta + 1, 3]
            deta[ipsi + 1, itheta + 1] = rz_in.fs[ipsi + 1, itheta + 1, 3] / (2π)
        end
    end

    # Handle 2π wrapping in deta
    for ipsi in 1:(mpsi + 1)
        for itheta in 2:(mtheta + 1)
            if deta[ipsi, itheta] - deta[ipsi, itheta - 1] > 0.5
                deta[ipsi, itheta] -= 1.0
            end
        end
    end

    # Replace phi in rzphi with unwrapped theta
    for ipsi in 0:mpsi
        for itheta in 0:mtheta
            rzphi_fs_nodes[ipsi + 1, itheta + 1, 3] = 2π * deta[ipsi + 1, itheta + 1]
        end
    end

    # --- Fit bicubic spline to rzphi
    rzphi = Spl.bicube_setup(
        rzphi_fs_nodes,
        psi_norm_grid,
        range(0.0, 2π; length=mtheta + 1),
        bctypex=4,
        bctypey=1,
        extrap="extrap",
        periodic="periodic",
    )

    # --- Evaluate and normalize quantities
    println("--> Evaluating magnetic field and geometric metrics...")
    eqfun_fs_nodes = zeros(mpsi + 1, mtheta + 1, 4)
    ro = 0.0
    zo = 0.0
    psio = 1.0

    for ipsi in 0:mpsi
        psi_norm = psi_norm_grid[ipsi + 1]
        for itheta in 0:mtheta
            theta_val = 2π * itheta / mtheta

            r, dr_dpsi, dr_dtheta = Spl.bicube_eval(rzphi, psi_norm, theta_val, 1)
            z, dz_dpsi, dz_dtheta = Spl.bicube_eval(rzphi, psi_norm, theta_val, 2)
            phi, dphi_dpsi, dphi_dtheta = Spl.bicube_eval(rzphi, psi_norm, theta_val, 3)

            # Covariant basis vectors
            a = [dr_dpsi, dz_dpsi, dphi_dpsi]
            b = [dr_dtheta, dz_dtheta, dphi_dtheta]

            cross_ab = [
                a[2]*b[3] - a[3]*b[2],
                a[3]*b[1] - a[1]*b[3],
                a[1]*b[2] - a[2]*b[1],
            ]
            jac = dot([r, 0.0, 0.0], cross_ab)
            b_sq = sq_in(psi_grid[ipsi + 1], 2)^2 / r^2 + jac^2 / r^2

            # Metric terms
            eqfun_fs_nodes[ipsi + 1, itheta + 1, 1] = b_sq
            eqfun_fs_nodes[ipsi + 1, itheta + 1, 2] = cross_ab[1]
            eqfun_fs_nodes[ipsi + 1, itheta + 1, 3] = cross_ab[2]
            eqfun_fs_nodes[ipsi + 1, itheta + 1, 4] = cross_ab[3]
        end
    end

    # Fit final spline to metric functions
    eqfun = Spl.bicube_setup(
        eqfun_fs_nodes,
        psi_norm_grid,
        range(0.0, 2π; length=mtheta + 1),
        bctypex=4,
        bctypey=1,
        extrap="extrap",
        periodic="periodic",
    )

    # --- Return final equilibrium
    println("--- Inverse reconstruction complete ---")
    return PlasmaEquilibrium(equil_params, sq_in, rzphi, eqfun, ro, zo, psio)
end
