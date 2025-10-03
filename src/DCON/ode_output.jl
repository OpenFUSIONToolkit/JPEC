"""
    ode_output_init(ctrl, equil, outp, intr, odet)

    Write header info to output files at the start of the integration.
"""
function ode_output_init(ctrl::DconControl, equil::Equilibrium.PlasmaEquilibrium, outp::DconOutput, intr::DconInternal, odet::OdeState)
    
    # TODO: mess with this to condense the number of write_output calls? Maybe allow it to pass in dicts
    # Write euler.bin header info
    if outp.write_euler_h5
        # Write dcon run parameters
        write_output(outp, :euler_h5, intr.mpert; dsetname="info/mpert")
        write_output(outp, :euler_h5, intr.mband; dsetname="info/mband")
        write_output(outp, :euler_h5, intr.mlow; dsetname="info/mlow")
        write_output(outp, :euler_h5, intr.mhigh; dsetname="info/mhigh")
        write_output(outp, :euler_h5, ctrl.nn; dsetname="info/nn")
        write_output(outp, :euler_h5, ctrl.singfac_min; dsetname="info/singfac_min")
        write_output(outp, :euler_h5, ctrl.kin_flag; dsetname="info/kin_flag")
        write_output(outp, :euler_h5, ctrl.con_flag; dsetname="info/con_flag")
        write_output(outp, :euler_h5, ctrl.mthvac; dsetname="info/mthvac")
        write_output(outp, :euler_h5, intr.psilim; dsetname="info/psilim")
        write_output(outp, :euler_h5, intr.qlim; dsetname="info/qlim")
        write_output(outp, :euler_h5, outp.mthsurf0; dsetname="info/mthsurf0") #TODO: mthsurf0 is deprecated

        # Write equilibrium parameters
        write_output(outp, :euler_h5, length(equil.rzphi.xs); dsetname="equil/nr") # TODO: equil save mpsi as really mpsi - 1, fix this
        write_output(outp, :euler_h5, length(equil.rzphi.ys); dsetname="equil/nz")
        write_output(outp, :euler_h5, equil.ro; dsetname="equil/ro")
        write_output(outp, :euler_h5, equil.zo; dsetname="equil/zo")
        write_output(outp, :euler_h5, equil.params.amean; dsetname="equil/amean")
        write_output(outp, :euler_h5, equil.params.rmean; dsetname="equil/rmean")
        write_output(outp, :euler_h5, equil.params.aratio; dsetname="equil/aratio")
        write_output(outp, :euler_h5, equil.params.kappa; dsetname="equil/kappa")
        write_output(outp, :euler_h5, equil.params.delta1; dsetname="equil/delta1")
        write_output(outp, :euler_h5, equil.params.delta2; dsetname="equil/delta2")
        write_output(outp, :euler_h5, equil.params.li1; dsetname="equil/li1")
        write_output(outp, :euler_h5, equil.params.li2; dsetname="equil/li2")
        write_output(outp, :euler_h5, equil.params.li3; dsetname="equil/li3")
        write_output(outp, :euler_h5, equil.params.betap1; dsetname="equil/betap1")
        write_output(outp, :euler_h5, equil.params.betap2; dsetname="equil/betap2")
        write_output(outp, :euler_h5, equil.params.betap3; dsetname="equil/betap3")
        write_output(outp, :euler_h5, equil.params.betat; dsetname="equil/betat")
        write_output(outp, :euler_h5, equil.params.betan; dsetname="equil/betan")
        write_output(outp, :euler_h5, equil.params.bt0; dsetname="equil/bt0")
        write_output(outp, :euler_h5, equil.params.q0; dsetname="equil/q0")
        write_output(outp, :euler_h5, equil.params.q95; dsetname="equil/q95")
        write_output(outp, :euler_h5, equil.params.qmin; dsetname="equil/qmin")
        write_output(outp, :euler_h5, equil.params.qmax; dsetname="equil/qmax")
        write_output(outp, :euler_h5, equil.params.qa; dsetname="equil/qa")
        write_output(outp, :euler_h5, equil.params.crnt; dsetname="equil/crnt")
        write_output(outp, :euler_h5, equil.params.psio; dsetname="equil/psio")
        write_output(outp, :euler_h5, equil.params.psilow; dsetname="equil/psilow")
        write_output(outp, :euler_h5, equil.config.control.power_b; dsetname="equil/power_b")
        write_output(outp, :euler_h5, equil.config.control.power_r; dsetname="equil/power_r")
        write_output(outp, :euler_h5, equil.config.control.power_bp; dsetname="equil/power_bp")
        write_output(outp, :euler_h5, 0; dsetname="equil/shotnum") # TODO: equil.params.shotnum)
        write_output(outp, :euler_h5, 0; dsetname="equil/shottime") # TODO: equil.params.shottime)

        # Write spline arrays
        write_output(outp, :euler_h5, equil.sq.xs; dsetname="splines/sq/xs")
        write_output(outp, :euler_h5, equil.sq.fs; dsetname="splines/sq/fs")
        write_output(outp, :euler_h5, equil.sq.fs1; dsetname="splines/sq/fs1")
        write_output(outp, :euler_h5, 0; dsetname="splines/sq/xpower") # TODO: equil.sq.xpower
        write_output(outp, :euler_h5, equil.rzphi.xs; dsetname="splines/rzphi/xs")
        write_output(outp, :euler_h5, equil.rzphi.ys; dsetname="splines/rzphi/ys")
        write_output(outp, :euler_h5, equil.rzphi.fs; dsetname="splines/rzphi/fs")
        write_output(outp, :euler_h5, equil.rzphi.fs1; dsetname="splines/rzphi/fsx")
        write_output(outp, :euler_h5, equil.rzphi.fs2; dsetname="splines/rzphi/fsy")
        write_output(outp, :euler_h5, equil.rzphi.fs12; dsetname="splines/rzphi/fsxy")
        write_output(outp, :euler_h5, 0; dsetname="splines/rzphi/x0") # TODO: equil.rzphi.x0
        write_output(outp, :euler_h5, 0; dsetname="splines/rzphi/y0") # TODO: equil.rzphi.y0
        write_output(outp, :euler_h5, 0; dsetname="splines/rzphi/xpower") # TODO: equil.rzphi.xpower
        write_output(outp, :euler_h5, 0; dsetname="splines/rzphi/fpower") # TODO: equil.rzphi.fpower
    end

    # Write crit.out header
    if odet.ising > 0 && odet.ising <= intr.msing
        singp = intr.sing[odet.ising]
        write_output(outp, :crit_out, @sprintf("   ising   psi         q          di      re alpha   im alpha\n"))
        write_output(outp, :crit_out, @sprintf("%6d%11.3e%11.3e%11.3e%11.3e%11.3e\n", odet.ising, singp.psifac, singp.q, singp.di, real(singp.alpha), imag(singp.alpha)))
    end
    write_output(outp, :crit_out, "    psifac      dpsi        q       singfac     eval1\n")
end

"""
    ode_output_step(odet, intr, ctrl, equil; force=false)

Performs output and monitoring tasks at each integration step.

This function calls `ode_output_monitor` to track critical eigenvalue behavior
and handle any diagnostics or logging associated with the current step.
Additional output (e.g., eigenvalue dumps, binary solution logging) may be added later.
"""
#TODO: depending on how we restructure our outputs, this function might be uncessary.
# It currently just calls ode_output_monitor! and does not write any outputs.
function ode_output_step!(odet::OdeState, intr::DconInternal, ctrl::DconControl, equil::Equilibrium.PlasmaEquilibrium, ffit::FourFitVars, outp::DconOutput; force::Bool=false)

    # Compute and print critical data for each time step
    ode_output_monitor!(odet, intr, ctrl, equil, outp)
    # if dout.out_evals || dout.bin_evals
    #     ode_output_get_evals(intr, ctrl, dout, fNames, equil, odet) # this is just outputs
    # end

    # Write solutions to euler.bin
    # TODO: for now, always writing to euler.bin every step. Can add option to write at all (through outp data struct)
    # or multiples of steps through eulers_stride. However, we don't have an istep in this version, so might need to do
    # multiples of psi instead of integer steps?
    # if bin_euler && (mod(odet.istep, euler_stride) == 0 || force)
    if outp.write_euler_h5 # TODO: update this logic for set number of psi vals
        du = zeros(ComplexF64, intr.mpert, intr.mpert, 2)
        params = (ctrl, equil, intr, odet, ffit, outp)
        sing_der!(du, odet.u, params, odet.psifac)
        write_to!(outp, :euler_h5, odet.psifac, odet.q, odet.msol)
        write_to!(outp, :euler_h5, odet.u)
        write_to!(outp, :euler_h5, odet.ud)
    end

    # Output solutions components for each time step
    #if bin_sol
    #    ode_output_sol(psifac, u, unorm) # again, more outputs
    #end

    return
end

# TODO: I don't think this function is essentially for the minimal working version - can convert later if needed
# function ode_output_get_evals(intr::DconInternal, ctrl::DconControl, dout::DconOutput, fNames::DconFileNames, equil::Equilibrium.PlasmaEquilibrium, odet::OdeState)
#     # Access all variables from structs, not globals.
#     # Pretty much just giving them all aliases so we don't have to type `intr.` and `ctrl.` every time.
#     u = odet.u #TODO: is ud the same as u
#     mpert = intr.mpert
#     out_evals = dout.out_evals
#     bin_evals = dout.bin_evals
#     nn = intr.nn
#     evals_out_unit = fNames.evals_out_unit
#     evals_bin_unit = fNames.evals_bin_unit
#     sq = equil.sq
#     psifac = intr.sing.psifac
#     q = equil.sq.f[4]

#     # Compute plasma response matrix
#     temp = conj.(transpose(u[:, 1:mpert, 1]))
#     wp = u[:, 1:mpert, 2]
#     wp = conj.(transpose(wp))

#     ipiv = zeros(Int, mpert)
#     info = Ref{Int}(0)

#     temp_lapack = copy(temp) # copy just renames a function pretty much
#     LAPACK.zgetrf!(temp_lapack, ipiv)
#     wp_lapack = copy(wp)
#     LAPACK.zgetrs!('N', temp_lapack, ipiv, wp_lapack)
#     wp = (wp_lapack + conj.(transpose(wp_lapack))) / 2

#     # Compute and sort eigenvalues
#     evals = zeros(Float64, mpert)
#     work = zeros(ComplexF64, 2*mpert-1)
#     rwork = zeros(Float64, 3*mpert-2)
#     evals, _ = LAPACK.zheev!('N', 'U', wp)
#     index = sortperm(-abs.(1 ./ evals), rev=true) # TODO: Should this still be rev

#     # Compute and sort inverse eigenvalues
#     evalsi, _ = LAPACK.zheev!('N', 'U', wp)
#     indexi = sortperm(-abs.(evalsi), rev = true) # TODO: Should this still be rev

#     # Write ascii eigenvalues
#     if out_evals && odet.istep > 0
#         @printf(evals_out_unit, "\nistep = %d, psifac = %.5e, q = %.5e\n", odet.istep, psifac, q)
#         @printf(evals_out_unit, "\n     i     eval       evali      err\n")
#         for ipert in 1:mpert
#             @printf(evals_out_unit, "%6d %11.3e %11.3e %11.3e\n",
#                 ipert,
#                 evalsi[indexi[ipert]],
#                 evals[index[ipert]],
#                 evalsi[indexi[ipert]] * evals[index[ipert]] - 1)
#         end
#         @printf(evals_out_unit, "\n     i     eval       evali      err\n")
#     end

#     # Write binary eigenvalues
#     if bin_evals && psifac > 0.1
#         spline_eval(sq, psifac, 0) #TODO: Where does sq come from?
#         q = sq.f[4]
#         singfac = abs(odet.m1 - nn * q)
#         logpsi1 = log10(psifac)
#         logpsi2 = log10(singfac)
#         for ipert in 2:mpert
#             write(evals_bin_unit, Float32(ipert), Float32(psifac),
#                 Float32(logpsi1), Float32(logpsi2),
#                 Float32(q), Float32(1 / evals[index[ipert]]))
#         end
#         write(io, evals_bin_unit) # End record #TODO: Looks like the function call should be write(io, stuff). Will it be written to rhe right place?
#     end

#     return
# end

"""
    ode_output_monitor!(odet, intr, ctrl, equil, outp)

Monitor the evolution of a critical eigenvalue (`crit`) during ODE integration and detect zero crossings, which indicate resonant or singular behavior.
The function evaluates `crit` using `ode_output_get_crit`, and if a sign change is found, it estimates the crossing point via linear interpolation.
If the crossing satisfies sharpness and consistency conditions, it's logged and `nzero` is incremented.
A termination flag (`termbycross_flag`) set within DconControl can be used to stop integration upon detection.

### Notes
- Zero-crossing detection is based on changes in the sign of `crit`.
- `u_med` is constructed as a linear interpolation between current and saved `u`.
- This function updates psi_prev, crit_prev, u_prev, and nzero in the `odet` state.

### TODO
- Restore or redirect output to appropriate logging units.
- Replace `error(...)` with graceful shutdown if zero crossing is an intended exit condition.
"""
function ode_output_monitor!(odet::OdeState, intr::DconInternal, ctrl::DconControl, equil::Equilibrium.PlasmaEquilibrium, outp::DconOutput)

    # Compute new crit
    dpsi = odet.psifac - odet.psi_prev
    q, singfac, logpsi1, logpsi2, crit = ode_output_get_crit(odet.psifac, odet.u, intr.mpert, odet.m1, ctrl.nn, equil.sq)

    # Check for zero crossing
    if crit * odet.crit_prev < 0 # if crit changes sign, zero crossing has occurred
        fac = crit / (crit - odet.crit_prev)
        psi_med = odet.psifac - fac * (odet.psifac - odet.psi_prev)
        dpsi = psi_med - odet.psi_prev
        u_med = odet.u .- fac .* (odet.u .- odet.u_prev)
        q_med, singfac_med, logpsi1_med, logpsi2_med, crit_med = ode_output_get_crit(psi_med, u_med, intr.mpert, odet.m1, ctrl.nn, equil.sq)

        if (crit_med - crit) * (crit_med - odet.crit_prev) < 0 &&
           abs(crit_med) < 0.5 * min(abs(crit), abs(odet.crit_prev))
            println("Zero crossing detected at psi = $psi_med, q = $q_med")
            write_output(outp, :dcon_out, @sprintf("Zero crossing at psi = %10.3e, q = %10.3e", psi_med, q_med))
            write_output(outp, :crit_out, @sprintf("Zero crossing at psi = %10.3e, q = %10.3e", psi_med, q_med))
            write_output(outp, :crit_out, @sprintf("%11.3e%11.3e%11.3e%11.3e%11.3e", psi_med,dpsi,q_med,singfac_med,crit_med))
            odet.nzero += 1
        end
        if ctrl.termbycross_flag
            #TODO: not sure if this is really an error, or just a user specified termination condition. So the program shoudln't error out, just stop integrating
            # In fortran, these were both handled by program_stop()
            error("Terminated by zero crossing.")
        end
    end

    # Write new crit
    write_output(outp, :crit_out, @sprintf("%11.3e%11.3e%11.3e%11.3e%11.3e", odet.psifac,dpsi,odet.q,singfac,crit))

    # Update saved values
    odet.psi_prev = odet.psifac
    odet.crit_prev = crit
    odet.u_prev .= odet.u
end

"""
    ode_output_get_crit(psi, u, mpert, m1, nn, sq) -> (q, singfac, logpsi1, logpsi2, crit)

Compute critical quantities at a given flux surface by constructing and inverting the plasma response matrix from the complex 
array `u`, symmetrizing it, computing its Hermitian eigenvalues, and using the smallest (in magnitude) 
inverse eigenvalue in combination with the equilibrium profiles to form `crit`.

This uses Julia’s built-in linear algebra:
`lu(temp) = zgetrf`, `temp_fact \\ wp = zgetrs`
Under the hood, Julia calls optimized LAPACK routines via the LinearAlgebra standard library.

### Main Arguments (others defined in main structs)
- `psi::Float64`: The flux surface at which to evaluate.
- `u::Array{ComplexF64, 3}`: Complex response data array of shape `(mpert, mpert, 2)`.

### Returns
A tuple `(q, singfac, logpsi1, logpsi2, crit)` of critical data for the time step.
"""
# TODO: on self-contained functions like this, it feels silly to pass in these large structs when we only need one variable
# Should we just pass in the variables we need? This seems more readable to me personally
#function ode_output_get_crit(psi::Float64, u::Array{ComplexF64, 3}, intr::DconInternal, odet::OdeState, ctrl::DconControl, equil::Equilibrium.PlasmaEquilibrium) # (alternate version)
function ode_output_get_crit(psi::Float64, u::Array{ComplexF64, 3}, mpert::Int, m1::Int, nn::Int, sq::Spl.CubicSpline)

    # Compute inverse plasma response matrix
    uu = u[:, 1:mpert, :]
    wp = adjoint(uu[:, :, 1])         # adjoint = conjugate transpose
    temp = adjoint(uu[:, :, 2])

    # Compute wp using LU decomposition
    temp_fact = lu(temp)
    wp = temp_fact \ wp

    # Symmetrize
    wp = (wp + adjoint(wp)) / 2

    # Compute and sort inverse eigenvalues
    evalsi = eigen(Hermitian(wp)).values
    indexi = sortperm(abs.(evalsi))  # bubble in Fortran sorts in descending order of -|evalsi|, we just do ascending order of |evalsi|
    
    # Compute critical data for each time step
    profiles = Spl.spline_eval(sq, psi, 0)
    q = profiles[4]
    singfac = abs(m1 - nn * profiles[4])
    logpsi1 = log10(psi)
    logpsi2 = log10(singfac) # TODO: why is this called logpsi2 and not logsingfac?
    crit = evalsi[indexi[1]] # * profiles[3]^2 # TODO: appears to be a bug in profiles[3]

    return q, singfac, logpsi1, logpsi2, crit
end