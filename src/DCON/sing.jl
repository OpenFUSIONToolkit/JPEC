"""
    sing_scan!(odet::OdeState, intr::DconInternal)

Scan all singular surfaces and calculate asymptotic vmat and mmat matrices and Mericer criterion.
"""
function sing_scan!(intr::DconInternal, ctrl::DconControl, equil::Equilibrium.PlasmaEquilibrium, ffit::FourFitVars)
    println("\n Singular Surfaces:")
    println("  i    psi      rho      q        q1      di0      di      err")
    # msol = intr.mpert # TODO: is msol used anywhere else? It's a global in Fortran and very confusing
    for ising in 1:intr.msing
        sing_vmat!(intr, ctrl, equil, ffit, ising)
    end
    println("  i    psi      rho      q        q1      di0      di      err")
end

"""
    sing_find!(ctrl::DconControl, equil::Equilibrium.PlasmaEquilibrium, intr::DconInternal; itmax=300)

Locate singular rational q-surfaces (q = m/nn) using a bisection method between extrema of the q-profile, and store their properties in `intr.sing`.

# Arguments
- 'itmax::Int`: Maximum number of iterations for the bisection method (default: 200)
"""
function sing_find!(ctrl::DconControl, equil::Equilibrium.PlasmaEquilibrium, intr::DconInternal; itmax=200)

    # Shorthand to evaluate q inside bisection search
    qval = psi -> Spl.spline_eval(equil.sq, psi, 0)[4]

    # Loop over extrema of q, find all rational values in between
    for iex in 2:equil.params.mextrema
        dq = equil.params.qextrema_q[iex] - equil.params.qextrema_q[iex-1]
        m = trunc(Int, ctrl.nn * equil.params.qextrema_q[iex-1])
        if dq > 0
            m += 1
        end
        dm = Int(sign(dq * ctrl.nn))

        # Loop over possible m's in interval
        while (m - ctrl.nn * equil.params.qextrema_q[iex-1]) * (m - ctrl.nn * equil.params.qextrema_q[iex]) <= 0
            it = 0
            psi0 = equil.params.qextrema_psi[iex-1]
            psi1 = equil.params.qextrema_psi[iex]
            psifac = (psi0 + psi1)/2 # initial guess for bisection

            # Bisection method to find singular surface
            while it < itmax
                it += 1
                psifac = (psi0 + psi1)/2
                singfac = (m - ctrl.nn * qval(psifac)) * dm
                if abs(singfac) <= 1e-8
                    break
                elseif singfac > 0
                    psi0 = psifac
                else
                    psi1 = psifac
                end
            end
            
            if it == itmax
                error("Bisection did not converge for m = $m")
            else
                push!(intr.sing, SingType(
                        m = m,
                        psifac = psifac,
                        rho = sqrt(psifac),
                        q = m / ctrl.nn,
                        q1 = Spl.spline_eval(equil.sq, psifac, 1)[2][4],
                ))
                intr.msing += 1
                if ctrl.verbose
                    println("Found singular surface: m=$(m), psifac=$(psifac), rho=$(sqrt(psifac)), q=$(m / ctrl.nn), q1=$(Spl.spline_eval(equil.sq, psifac, 1)[2][4])")
                end
            end
            m += dm
        end
    end
end

"""
    sing_lim!(intr::DconInternal, ctrl::DconControl, equil::Equilibrium.PlasmaEquilibrium, eqCtrl::Equilibrium.EquilibriumControl, eqPrm::Equilibrium.EquilibriumParameters)

Compute and set limiter values for the DCON analysis - handles cases where user
truncates before the last singular surface. 

# Arguments
- `itmax::Int`: The maximum number of iterations allowed for the psilim search (default: 50)
- `eps::Float64`: The convergence tolerance for the psilim search (default: 1e-10)
"""
# JMH - this has been checked against the DIID_Ideal example (with psiedge = 0.95 to enter the last if statement). 
# There are some discrepancies at the 4-5th digit between the two, not sure if this is numerical noise, 
# an issue in the spline evals, or a bug here. However, the differences in the main internal parameters are small
# so ignoring for now, can address in a unit test later.
function sing_lim!(intr::DconInternal, ctrl::DconControl, equil::Equilibrium.PlasmaEquilibrium; itmax=50, eps=1e-10)

    # Shorthand to evaluate q/q1 inside newton iteration
    qval = psi -> Spl.spline_eval(equil.sq, psi, 0)[4]
    q1val = psi -> Spl.spline_eval(equil.sq, psi, 1)[2][4]

    #compute and modify the DconInternal struct 
    intr.qlim   = min(equil.params.qmax, ctrl.qhigh)
    intr.q1lim  = equil.sq.fs1[equil.config.control.mpsi+1, 4]
    intr.psilim = equil.config.control.psihigh 
    #normalize dmlim to interval [0,1)
    if ctrl.sas_flag
        while ctrl.dmlim > 1.0
            ctrl.dmlim -= 1.0
        end
        while ctrl.dmlim < 0.0
            ctrl.dmlim += 1.0
        end
        #compute qlim
        intr.qlim = (trunc(Int, ctrl.nn * intr.qlim) + ctrl.dmlim) / ctrl.nn
        while intr.qlim > equil.params.qmax #could also be a while true with a break condition if (qlim <= qmax) like the Fortran code
            intr.qlim -= 1.0 / ctrl.nn
        end
    end

    #use newton iteration to find psilim
    if intr.qlim < equil.params.qmax
        # Find index jpsi that minimizes |equil.sq.fs[:,4] - qlim| and ensure within mpsi limits
        _, jpsi = findmin(abs.(equil.sq.fs[:,4] .- intr.qlim))
        jpsi = min(jpsi, equil.config.control.mpsi - 1)

        intr.psilim = equil.sq.xs[jpsi]
        it = 0
        while it < itmax
            it += 1
            dpsi = (intr.qlim - qval(intr.psilim)) / q1val(intr.psilim)
            intr.psilim += dpsi

            if abs(dpsi) < eps * abs(intr.psilim)
                intr.q1lim = q1val(intr.psilim)
                break
            end
        end

        # abort if not found
        if it > itmax
            error("Can't find psilim.")
        end

    else
        intr.qlim = equil.params.qmax
        intr.q1lim = equil.sq.fs1[equil.config.control.mpsi+1,4]
        intr.psilim = equil.config.control.psihigh 
    end

    #set up record for determining the peak in dW near the boundary.
    if ctrl.psiedge < intr.psilim
        qedgestart = trunc(Int, Spl.spline_eval(equil.sq, ctrl.psiedge, 0)[4])
        intr.size_edge = ceil(Int, (intr.qlim - qedgestart) * ctrl.nn * ctrl.nperq_edge)

        intr.dw_edge  = fill(-typemax(Float64) * (1 + im), intr.size_edge)
        intr.q_edge   = [qedgestart + i / (ctrl.nperq_edge * ctrl.nn) for i in 0:intr.size_edge-1]
        intr.psi_edge = zeros(intr.size_edge)

        # monitor some deeper points for an informative profile
        intr.pre_edge = 1
        for i in 1:intr.size_edge
            if intr.q_edge[i] < Spl.spline_eval(equil.sq, ctrl.psiedge, 0)[4] 
                intr.pre_edge += 1
            end
        end
    end
end

# Main differences: using 1 indexing, so zeroth order is the first index, whereas Fortran used 0 indexing for these arrays
function sing_vmat!(intr::DconInternal, ctrl::DconControl, equil::Equilibrium.PlasmaEquilibrium, ffit::FourFitVars, ising::Int)
    if ising < 1 || ising > intr.msing
        return
    end
    singp = intr.sing[ising]

    # TODO: move these to initialiation by passing in mpert and sing_order
    singp.vmat = zeros(ComplexF64, intr.mpert, 2*intr.mpert, 2, 2*ctrl.sing_order + 1)
    singp.mmat = zeros(ComplexF64, intr.mpert, 2*intr.mpert, 2, 2*ctrl.sing_order + 3)
    singp.power = zeros(ComplexF64, 2*intr.mpert)

    # Identify the resonant solution
    ipert0 = round(Int, ctrl.nn * singp.q, RoundFromZero) - intr.mlow + 1
    # TODO: is this needed if we initialize di to zero in the struct definition?
    if ipert0 <= 0 || intr.mlow > ctrl.nn * singp.q || intr.mhigh < ctrl.nn * singp.q
        singp.di = 0
        return
    end

    # Compute ranges
    singp.r1 = [ipert0]
    singp.r2 = [ipert0, ipert0 + intr.mpert]
    singp.n1 = [i for i in 1:intr.mpert if i != ipert0]
    singp.n2 = vcat(singp.n1, [i + intr.mpert for i in singp.n1])

    psifac = singp.psifac
    q = singp.q    
    # TODO: uncomment this once locstab (created in mercier.f) is working
    # spline_eval(locstab, singp.psifac, 0)
    # di0 = locstab.f[1] / singp.psifac
    di0 = 0.0 # placeholder
    q1 = singp.q1
    rho = singp.rho

    # Compute Mercier criterion and singular power
    sing_mmat!(intr, ctrl, equil, ffit, ising)
    singp.m0mat = transpose(singp.mmat[singp.r1[1], singp.r2, :, 1])
    # TODO: this is a little odd. di is complex since m0 is complex (but the imaginaries are negligible),
    # its then stored in singp where the type is real, and then converted to complex again for alpha and power.
    # There's gotta be a more clear way to do this? But this is the most faithful to the Fortran
    di = singp.m0mat[1,1]*singp.m0mat[2,2] - singp.m0mat[2,1]*singp.m0mat[1,2]
    singp.di = real(di)
    singp.alpha = sqrt(-complex(singp.di))
    singp.power[ipert0] = -singp.alpha
    singp.power[ipert0 + intr.mpert] = singp.alpha

    println(@sprintf("%3d %11.3e %11.3e %11.3e %11.3e %11.3e %11.3e %11.3e",
        ising, psifac, rho, q, q1, di0, singp.di, singp.di/di0-1))

    # Zeroth-order non-resonant solutions
    singp.vmat .= 0
    for ipert in 1:intr.mpert
        singp.vmat[ipert, ipert, 1, 1] = 1
        singp.vmat[ipert, ipert + intr.mpert, 2, 1] = 1
    end

    # Zeroth-order resonant solutions
    singp.vmat[ipert0, ipert0, 1, 1] = 1
    singp.vmat[ipert0, ipert0 + intr.mpert, 1, 1] = 1
    singp.vmat[ipert0, ipert0, 2, 1] = -(singp.m0mat[1,1] + singp.alpha) / singp.m0mat[1,2]
    singp.vmat[ipert0, ipert0 + intr.mpert, 2, 1] = -(singp.m0mat[1,1] - singp.alpha) / singp.m0mat[1,2]
    det = conj(singp.vmat[ipert0, ipert0, 1, 1]) * singp.vmat[ipert0, ipert0 + intr.mpert, 2, 1] -
          conj(singp.vmat[ipert0, ipert0 + intr.mpert, 1, 1]) * singp.vmat[ipert0, ipert0, 2, 1]
    singp.vmat[ipert0, :, :, 1] ./= sqrt(det)

    # Higher order solutions
    for k in 1:2*ctrl.sing_order
        sing_solve!(singp, intr, k)
    end
end

# Main differences: using 1 indexing, so zeroth order is the first index, whereas Fortran used 0 indexing for these arrays. Also adapted for
# dense f/g/k arrays here, so a lot of the LAPACK operations go under the hood.
function sing_mmat!(intr::DconInternal, ctrl::DconControl, equil::Equilibrium.PlasmaEquilibrium, ffit::FourFitVars, ising::Int)

    # ongoing TODO: list here
    # figure out if there's a direct way of doing the banded->dense matrix conversion to Julia (very messy right now)
    # avoided hardcoding of binomial coefficients for cleanup

    # Initial allocations
    msol = 2 * intr.mpert # TODO: make sure this doesn't get updated as a global elsewhere in the Fortran
    q = zeros(Float64, 4)
    singfac = zeros(Float64, intr.mpert, 4)
    f_interp = zeros(ComplexF64, intr.mpert, intr.mpert, 4)
    g_interp = zeros(ComplexF64, intr.mpert, intr.mpert, 4)
    k_interp = zeros(ComplexF64, intr.mpert, intr.mpert, 4)
    f = zeros(ComplexF64, intr.mpert, intr.mpert, ctrl.sing_order + 1)
    f0 = zeros(ComplexF64, intr.mpert, intr.mpert)
    ff = zeros(ComplexF64, intr.mpert, intr.mpert, ctrl.sing_order + 1)
    g = zeros(ComplexF64, intr.mpert, intr.mpert, ctrl.sing_order + 1)
    k = zeros(ComplexF64, intr.mpert, intr.mpert, ctrl.sing_order + 1)
    v = zeros(ComplexF64, intr.mpert, msol, 2)
    x = zeros(ComplexF64, intr.mpert, msol, 2, ctrl.sing_order + 1)

    singp = intr.sing[ising]

    # Evaluate cubic splines
    # TODO: third derivative has some error, but only included via sing_fac[ipert0] for sing_order < 3. Tests with solovev ideal indicate little sensitivity
    # TODO: this is an annoying way to have to take apart this tuple of vectors, I think
    # this is a planned fix already (i.e. separating cubic splines)
    q .= getindex.(Spl.spline_eval(equil.sq, singp.psifac, 3), 4)
    f_interp[:, :, 1], f_interp[:, :, 2], f_interp[:, :, 3], f_interp[:, :, 4] = Spl.spline_eval(ffit.fmats, singp.psifac, 3)
    g_interp[:, :, 1], g_interp[:, :, 2], g_interp[:, :, 3], g_interp[:, :, 4] = Spl.spline_eval(ffit.gmats, singp.psifac, 3)
    k_interp[:, :, 1], k_interp[:, :, 2], k_interp[:, :, 3], k_interp[:, :, 4] = Spl.spline_eval(ffit.kmats, singp.psifac, 3)

    # Evaluate singfac and its derivatives
    ipert0 = singp.m - intr.mlow + 1
    mvec = intr.mlow:intr.mhigh
    singfac[:, 1] .= mvec .- ctrl.nn * q[1]
    singfac[:, 2] .= -ctrl.nn * q[2]
    singfac[:, 3] .= -ctrl.nn * q[3]
    singfac[:, 4] .= -ctrl.nn * q[4]
    singfac[ipert0, 1] = -ctrl.nn * q[2]
    singfac[ipert0, 2] = -ctrl.nn * q[3] / 2
    singfac[ipert0, 3] = -ctrl.nn * q[4] / 3
    singfac[ipert0, 4] = 0

    # Compute NON-factored Hermitian F and its derivatives
    # TODO: only supports dense matrices for now, implement Banded?
    # TODO: this can be simplified to extracting binomial coefficients (see ChatGPT)
    for jpert in 1:intr.mpert
        for ipert in jpert:min(intr.mpert, jpert + intr.mband)
            f[ipert, jpert, 1] = singfac[ipert, 1] * f_interp[ipert, jpert, 1]
            if ctrl.sing_order < 1 continue end
            f[ipert, jpert, 2] = singfac[ipert, 1] * f_interp[ipert, jpert, 2] +
                                 singfac[ipert, 2] * f_interp[ipert, jpert, 1]
            if ctrl.sing_order < 2 continue end
            f[ipert, jpert, 3] = singfac[ipert, 1] * f_interp[ipert, jpert, 3] + 
                                 2 * singfac[ipert, 2] * f_interp[ipert, jpert, 2] + 
                                 singfac[ipert, 3] * f_interp[ipert, jpert, 1]
            if ctrl.sing_order < 3 continue end
            f[ipert, jpert, 4] = singfac[ipert, 1] * f_interp[ipert, jpert, 4] + 
                                 3 * singfac[ipert, 2] * f_interp[ipert, jpert, 3] + 
                                 3 * singfac[ipert, 3] * f_interp[ipert, jpert, 2] + 
                                 singfac[ipert, 4] * f_interp[ipert, jpert, 1]
            if ctrl.sing_order < 4 continue end
            f[ipert, jpert, 5] = 4 * singfac[ipert, 2] * f_interp[ipert, jpert, 4] +
                                 6 * singfac[ipert, 3] * f_interp[ipert, jpert, 3] +
                                 4 * singfac[ipert, 4] * f_interp[ipert, jpert, 2]
            if ctrl.sing_order < 5 continue end
            f[ipert, jpert, 6] = 10 * singfac[ipert, 3] * f_interp[ipert, jpert, 4] +
                                 10 * singfac[ipert, 4] * f_interp[ipert, jpert, 3]
            if ctrl.sing_order < 6 continue end
            f[ipert, jpert, 7] = 20 * singfac[ipert, 4] * f_interp[ipert, jpert, 4]
        end
    end
    f0 .= f[:, :, 1]

    # Compute product of Hermitian matrix F
    fac0 = 1
    for n in 0:ctrl.sing_order
        fac1 = 1
        for j in 0:n
            for jpert in 1:intr.mpert
                for ipert in jpert:min(intr.mpert, jpert + intr.mband)
                    for kpert in max(1, ipert - intr.mband):jpert
                        ff[ipert, jpert, n + 1] += fac1 * f[ipert, kpert, j + 1] * conj(f[jpert, kpert, n - j + 1])
                    end
                end
            end
            fac1 *= (n - j) / (j + 1)
        end
        ff[:, :, n + 1] ./= fac0
        fac0 *= (n + 1)
    end

    # NOTE: In Fortran, f, g, k are stored as banded Hermitian matrices (lower triangle only), 
    # with f already factored, this is passed into LAPACK which handles all of this under the hood.
    # However, for now we are using dense matrices in Julia, so ff is stored as a full matrix (both triangles), 
    # and we need to fill in the upper triangle here.
    # Use Hermitian property to fill the upper triangle
    for n in 1:ctrl.sing_order
        for i in 1:intr.mpert
            for j in 1:i-1
                ff[j, i, n + 1] = conj(ff[i, j, n + 1])
            end
        end
    end

    # Compute non-Hermitian matrix K
    for jpert in 1:intr.mpert
        for ipert in max(1, jpert - intr.mband):min(intr.mpert, jpert + intr.mband)
            k[ipert, jpert, 1] = singfac[ipert, 1] * k_interp[ipert, jpert, 1]
            if ctrl.sing_order < 1 continue end
            k[ipert, jpert, 2] = singfac[ipert, 1] * k_interp[ipert, jpert, 2] +
                                 singfac[ipert, 2] * k_interp[ipert, jpert, 1]
            if ctrl.sing_order < 2 continue end
            k[ipert, jpert, 3] = singfac[ipert, 1] * k_interp[ipert, jpert, 3] / 2 +
                                 singfac[ipert, 2] * k_interp[ipert, jpert, 2] +
                                 singfac[ipert, 3] * k_interp[ipert, jpert, 1] / 2
            if ctrl.sing_order < 3 continue end
            k[ipert, jpert, 4] = singfac[ipert, 1] * k_interp[ipert, jpert, 4] / 6 +
                                 singfac[ipert, 2] * k_interp[ipert, jpert, 3] / 2 +
                                 singfac[ipert, 3] * k_interp[ipert, jpert, 2] / 2 +
                                 singfac[ipert, 4] * k_interp[ipert, jpert, 1] / 6
            if ctrl.sing_order < 4 continue end
            k[ipert, jpert, 5] = singfac[ipert, 2] * k_interp[ipert, jpert, 4] / 6 +
                                 singfac[ipert, 3] * k_interp[ipert, jpert, 3] / 4 +
                                 singfac[ipert, 4] * k_interp[ipert, jpert, 2] / 6
            if ctrl.sing_order < 5 continue end
            k[ipert, jpert, 6] = singfac[ipert, 3] * k_interp[ipert, jpert, 4] / 12 +
                                 singfac[ipert, 4] * k_interp[ipert, jpert, 3] / 12
            if ctrl.sing_order < 6 continue end
            k[ipert, jpert, 7] = singfac[ipert, 4] * k_interp[ipert, jpert, 4] / 36
        end
    end
    
    # Compute Hermitian matrix G
    for jpert in 1:intr.mpert
        for ipert in jpert:min(intr.mpert, jpert + intr.mband)
            g[ipert, jpert, 1] = g_interp[ipert, jpert, 1]
            if ctrl.sing_order < 1 continue end
            g[ipert, jpert, 2] = g_interp[ipert, jpert, 2]
            if ctrl.sing_order < 2 continue end
            g[ipert, jpert, 3] = g_interp[ipert, jpert, 3] / 2
            if ctrl.sing_order < 3 continue end
            g[ipert, jpert, 4] = g_interp[ipert, jpert, 3] / 6
        end
    end

    # Compute identity
    for ipert in 1:intr.mpert
        v[ipert, ipert, 1] = 1
        v[ipert, ipert + intr.mpert, 2] = 1
    end

    # Compute zeroth-order x1
    for isol=1:msol
        x[:, isol, 1, 1] .= v[:, isol, 2] .- k[:, :, 1] * v[:, isol, 1]
    end
    # f is prefactorized so can just use this calculation to get F⁻¹x
    x[:, :, 1, 1] = UpperTriangular(f0') \ (LowerTriangular(f0) \ x[:, :, 1, 1])

    # Compute higher-order x1
    for i=1:ctrl.sing_order
        for isol=1:msol
            for j=1:i
                x[:, isol, 1, i + 1] .-= ff[:, :, j + 1] * x[:, isol, 1, i - j + 1]
            end
            x[:, isol, 1, i + 1] .-= k[:, :, i + 1] * v[:, isol, 1]
        end
        x[:, :, 1, i + 1] = UpperTriangular(f0') \ (LowerTriangular(f0) \ x[:, :, 1, i + 1])
    end

    # Use Hermitian property to fill the upper triangle of g
    # TODO: very hacky and unclear. Eventually, completely remove all mention of banding/Hermitian storage or make it cleaner
    # Lots of stuff hidden in LAPACK calls in Fortran that we are doing manually here
    for n in 0:ctrl.sing_order
        for i in 1:intr.mpert
            for j in 1:i-1
                g[j, i, n + 1] = conj(g[i, j, n + 1])
            end
        end
    end

    # Compute x2
    for i=0:ctrl.sing_order
        for isol=1:msol
            for j = 0:i
                x[:, isol, 2, i + 1] .+= k[:, :, j + 1]' * x[:, isol, 1, i - j + 1]
            end
            x[:, isol, 2, i + 1] .+= g[:, :, i + 1] * v[:, isol, 1]
        end
    end

    # Principal terms of mmat
    singp.mmat .= 0
    r1 = singp.r1
    r2 = singp.r2
    n1 = singp.n1
    n2 = singp.n2
    j = 0
    for i = 0:ctrl.sing_order
        singp.mmat[r1, r2, :, j+1] .= x[r1, r2, :, i+1]
        singp.mmat[r1, n2, :, j+2] .= x[r1, n2, :, i+1]
        singp.mmat[n1, r2, :, j+2] .= x[n1, r2, :, i+1]
        singp.mmat[n1, n2, :, j+3] .= x[n1, n2, :, i+1]
        j += 2
    end

    # Shearing terms
    singp.mmat[r1, r2[1], 1, 1] .+= 0.5
    singp.mmat[r1, r2[2], 2, 1] .-= 0.5
end

"""
    sing_solve!(singp::SingType, k::Int)

Solves iteratively for the next order in the power series `singp.vmat`. See equation 47 in the Glass 2016 DCON paper.

Arguments
---------
- `singp::SingType`: The singularity data structure containing all relevant matrices and parameters.
- `k::Int`: The current order in the power series expansion.
"""
function sing_solve!(singp::SingType, intr::DconInternal, k::Int)
    for l in 1:k
        singp.vmat[:,:,:,k + 1] .+= sing_matmul(singp.mmat[:,:,:,l + 1], singp.vmat[:,:,:,k - l + 1])
    end
    for isol in 1:2*intr.mpert
        a = copy(singp.m0mat)
        a[1,1] -= k/2.0 + singp.power[isol]
        a[2,2] -= k/2.0 + singp.power[isol]
        det = a[1,1] * a[2,2] - a[1,2] * a[2,1]
        x = -singp.vmat[singp.r1[1], isol, :, k+1]
        singp.vmat[singp.r1[1], isol, 1, k+1] = (a[2,2]*x[1] - a[1,2]*x[2]) / det
        singp.vmat[singp.r1[1], isol, 2, k+1] = (a[1,1]*x[2] - a[2,1]*x[1]) / det
        singp.vmat[singp.n1, isol, :, k+1] ./= (singp.power[isol] + k/2.0)
    end
end

"""
    sing_matmul(a, b) -> c

Matrix multiplication specific to singular matrices.

Arguments
---------
- `a::Array{ComplexF64,3}`: shape (mpert, 2 * mpert, 2)
- `b::Array{ComplexF64,3}`: shape (mmpert, 2 * mpert, 2)

Returns
-------
- `c::Array{ComplexF64,3}`: shape (mpert, 2 * mpert, 2)
"""
function sing_matmul(a::Array{ComplexF64,3}, b::Array{ComplexF64,3})
    m = size(b,1)
    n = size(b,2)

    # consistency check
    if size(a,2) != 2*m
        error("Sing_matmul: size(a,2) = $(size(a,2)) != 2*size(b,1) = $(2*m)")
    end

    c = zeros(ComplexF64, size(a,1), n, 2)

    # main computation
    for i in 1:n
        for j in 1:2
            c[:, i, j] .+= a[:, 1:m, j] * b[:, i, 1] +
                           a[:, m+1:2*m, j] * b[:, i, 2]
        end
    end

    return c
end

"""
    sing_get_ua!(odet::OdeState, intr::DconInternal, ctrl::DconControl)

Compute the asymptotic series solution for a given singularity.

- `odet`: Current ODE state (contains psifac, ua, ising, etc.)
- `intr`: DCON internal state (contains singularity data, mpert, etc.)
- `ctrl`: DCON control parameters (contains sing_order, etc.)

Fills `ua` with the asymptotic solution vmat for the specified singular surface and psi value.
"""
function sing_get_ua!(odet::OdeState, intr::DconInternal, ctrl::DconControl)

    singp = intr.sing[odet.ising]
    r1 = singp.r1
    r2 = singp.r2

    # Compute distance from singular surface
    dpsi = odet.psifac - singp.psifac
    sqrtfac = sqrt(dpsi)
    pfac = abs(dpsi)^singp.alpha

    # Compute power series via Horner's method
    odet.ua .= singp.vmat[:,:,:,2 * ctrl.sing_order + 1]
    for iorder in (2 * ctrl.sing_order - 1):-1:0
        odet.ua .= odet.ua .* sqrtfac .+ singp.vmat[:,:,:,iorder+1]
    end

    # Restore powers
    odet.ua[r1, :, 1] ./= sqrtfac
    odet.ua[r1, :, 2] .*= sqrtfac
    odet.ua[:, r2[1], :] ./= pfac
    odet.ua[:, r2[2], :] .*= pfac

    # Renormalize
    if odet.psifac < singp.psifac
        odet.ua[:, r2[1], :] .*= abs(odet.ua[r1[1], r2[1], 1]) / odet.ua[r1[1], r2[1], 1]
        odet.ua[:, r2[2], :] .*= abs(odet.ua[r1[1], r2[2], 1]) / odet.ua[r1[1], r2[2], 1]
    end
end

"""
    sing_get_ca!(ca::AbstractArray{ComplexF64,3}, intr::DconInternal, ising::Int, psifac::Real, u::AbstractArray{ComplexF64,3})

Compute the asymptotic expansion coefficients for a given singularity and solution array.

- `ca`: Output 3D complex array (mutated in-place), shape (mpert, msol, 2)
- `intr`: DCON internal state (contains singularity data, mpert, etc.)
- `ising`: Index of the singularity (integer)
- `psifac`: Input psi factor (real scalar)
- `u`: Input 3D complex array (solution array), shape (mpert, msol, 2)

Fills `ca` with the expansion coefficients for the specified singularity and psi value.
"""
#TODO: Are these comments correct?
#Compute asymptotic coefficients 
# Inputs:
#   ising is just passed through if needed externally
#   psifac: input psi factor
#   u: input array
# Output:
#   ca::Array{ComplexF64,3}
function sing_get_ca!(ca::AbstractArray{ComplexF64,3},
    intr::DconInternal,
    ising::Int,
    psifac::Real,
    u::AbstractArray{ComplexF64,3}) #TODO: ca is being modified so from Julia conventions, it should be listed first. In Fortran, it is listed last

    # number of solutions
    msol = size(u,2)
    mpert = intr.mpert

    # call asymptotic solution generator
    ua = similar(u)
    ua = sing_get_ua!(ua,intr, ising, psifac)

    # build system matrix temp1 (2*mpert × 2*mpert)
    temp1 = Array{ComplexF64}(undef, 2*mpert, 2*mpert)
    temp1[1:mpert, :] .= ua[:, :, 1]
    temp1[mpert+1:2*mpert, :] .= ua[:, :, 2]

    # build right-hand side temp2 (2*mpert × msol)
    temp2 = Array{ComplexF64}(undef, 2*mpert, msol)
    temp2[1:mpert, :] .= u[:, :, 1]
    temp2[mpert+1:2*mpert, :] .= u[:, :, 2]

    # LU factorization and solve
    F = lu(temp1)
    temp2 .= F \ temp2

    # output coefficients ca (mpert × msol × 2)
    #ca = Array{ComplexF64}(undef, mpert, msol, 2) #don't think we need this
    ca[:, 1:msol, 1] .= temp2[1:mpert, :]
    ca[:, 1:msol, 2] .= temp2[mpert+1:2*mpert, :]

    return ca
end

"""
    sing_der!(
        du::Array{ComplexF64,3},
        u::Array{ComplexF64,3},
        params::Tuple{DconControl, Equilibrium.PlasmaEquilibrium, DconInternal, FourFitVars},
        psieval::Float64
    )

Evaluate the derivative of the Euler-Lagrange equations, i.e. u' in equation 24 of Glasser 2016.
This follows the Julia DifferentialEquations package format for in place updating.

    ode_function!(du, u, p, t)

# "Defining your ODE function to be in-place updating can have performance benefits. 
# What this means is that, instead of writing a function which outputs its solution, 
# you write a function which updates a vector that is designated to hold the solution. 
# By doing this, DifferentialEquations.jl's solver packages are able to reduce the 
# amount of array allocations and achieve better performance."

Wherever possible, in-place operations on pre-allocated arrays are used to minimize memory allocations.
All LAPACK operations are handled under the hood by Julia's LinearAlgebra package, so we can obtain a much
more simplistic code with similar performance.
"""
function sing_der!(du::Array{ComplexF64, 3}, u::Array{ComplexF64, 3},
                   params::Tuple{DconControl, Equilibrium.PlasmaEquilibrium,
                                 DconInternal, OdeState, FourFitVars},
                   psieval::Float64)
    # Unpack structs
    ctrl, equil, intr, odet, ffit = params

    # Spline evaluation
    odet.q = Spl.spline_eval(equil.sq, psieval, 0)[4]
    odet.singfac_vec .= 1.0 ./ (intr.mlow .- ctrl.nn * odet.q .+ (0:intr.mpert-1))
    chi1 = 2π * equil.psio

    # kinetic stuff - skip for now
    if false #(TODO: kin_flag)
        error("kin_flag not implemented yet")
    else
        # Evaluate splines at psieval and reshape avoiding new allocations
        # TODO: is this actually more efficient?
        Spl.spline_eval!(odet.amat, ffit.amats, psieval; derivs=0)
        amat = reshape(odet.amat, intr.mpert, intr.mpert)
        Spl.spline_eval!(odet.bmat, ffit.bmats, psieval; derivs=0)
        bmat = reshape(odet.bmat, intr.mpert, intr.mpert)
        Spl.spline_eval!(odet.cmat, ffit.cmats, psieval; derivs=0)
        cmat = reshape(odet.cmat, intr.mpert, intr.mpert)
        Spl.spline_eval!(odet.fmat, ffit.fmats, psieval; derivs=0)
        fmat = reshape(odet.fmat, intr.mpert, intr.mpert)
        Spl.spline_eval!(odet.kmat, ffit.kmats, psieval; derivs=0)
        kmat = reshape(odet.kmat, intr.mpert, intr.mpert)
        Spl.spline_eval!(odet.gmat, ffit.gmats, psieval; derivs=0)
        gmat = reshape(odet.gmat, intr.mpert, intr.mpert)

        # ldiv!(A,B): Compute A \ B in-place and overwriting B to store the result,
        # where A is a factorization object.
        odet.Afact = cholesky(Hermitian(amat))
        # Solve Afact \ bmat
        ldiv!(odet.Afact, bmat)
        # Solve Afact \ cmat
        ldiv!(odet.Afact, cmat)

        # TODO: banded matrix calculations would go here
    end
    
    # Compute du
    if false #(TODO: kin_flag)
        error("kin_flag not implemented yet")
    else
        # See equation 23 in Glasser 2016 DCON paper for derivation
        # TODO: these comments are a work in progress and should be checked, but I think they're close
        # main thing is I'm not sure where singfac comes from

        # du[1] = - K * u[1] + u[2]
        du[:, :, 1] .= u[:, :, 2] .* odet.singfac_vec .- kmat * u[:, :, 1]

        # du[1] = - F⁻¹ * K * u[1] + F⁻¹ * u[2] (remember F is already stored in factored form)
        du[:, :, 1] .= UpperTriangular(fmat') \ (LowerTriangular(fmat) \ du[:, :, 1])

        # du[2] = G * u[1] + K' * du[1] = G * u[1] - K^† * F⁻¹ * K * u[1] + K^† * F⁻¹ * u[2]
        du[:, :, 2] .= gmat * u[:, :, 1] .+ adjoint(kmat) * du[:, :, 1]
        du[:, :, 1] .*= odet.singfac_vec
    end
    # TODO: this is used in GPEC, so will need to dump to file later
    # ud[:,:,1] .= du[:,:,1]
    # ud[:,:,2] .= -bmat * du[:,:,1] - cmat * u[:,:,1]
end
