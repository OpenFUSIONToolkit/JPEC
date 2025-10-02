"""
    ode_output_init(ctrl, equil, fnames, intr, odet)

    Write header info to output files at the start of the integration.
"""
function ode_output_init(ctrl::DconControl, equil::Equilibrium.PlasmaEquilibrium, fnames::DconFileNames, intr::DconInternal, odet::OdeState)
    
    # Write euler.bin header

    # Write equilibrium grid and parameters
    write_to!(fnames, :euler_bin_unit, intr.mlow, intr.mhigh, ctrl.nn, length(equil.rzphi.xs), length(equil.rzphi.ys), equil.ro, equil.zo)
    write_to!(fnames, :euler_bin_unit, intr.mband, NaN, ctrl.mthvac, equil.psio, equil.config.control.psilow, intr.psilim, intr.qlim, ctrl.singfac_min) #TODO: mthsurf0
    write_to!(fnames, :euler_bin_unit, equil.config.control.power_b, equil.config.control.power_r, equil.config.control.power_bp)
    write_to!(fnames, :euler_bin_unit, ctrl.kin_flag, ctrl.con_flag)

    # Write equilibrium summary parameters
    write_to!(fnames, :euler_bin_unit, equil.params.amean, equil.params.rmean, equil.params.aratio, equil.params.kappa, equil.params.delta1, equil.params.delta2,
        equil.params.li1, equil.params.li2, equil.params.li3, equil.params.betap1, equil.params.betap2, equil.params.betap3, equil.params.betat, equil.params.betan, equil.params.bt0,
        equil.params.q0, equil.params.qmin, equil.params.qmax, equil.params.qa, equil.params.crnt, equil.params.q95, NaN, NaN) # TODO: , equil.params.shotnum, equil.params.shottime)

    # Write spline and grid arrays
    write_to!(fnames, :euler_bin_unit, equil.sq.xs, equil.sq.fs, equil.sq.fs1, 0) # TODO:, equil.sq.xpower)
    write_to!(fnames, :euler_bin_unit, equil.rzphi.xs, equil.rzphi.ys,
        equil.rzphi.fs, equil.rzphi.fsx, equil.rzphi.fsy, equil.rzphi.fsxy,
        NaN, NaN, NaN, NaN) # TODO:equil.rzphi.x0, equil.rzphi.y0,, equil.rzphi.xpower, equil.rzphi.ypower)

    # Write crit.out header
    if odet.ising > 0 && odet.ising <= intr.msing 
        write_to!(fnames, :crit_out_unit, @sprintf(
            "   ising   psi         q          di      re alpha   im alpha\n\n%6d%11.3e%11.3e%11.3e%11.3e%11.3e\n",
            odet.ising,
            intr.sing[odet.ising].psifac,
            intr.sing[odet.ising].q,
            intr.sing[odet.ising].di,
            real(intr.sing[odet.ising].alpha),
            imag(intr.sing[odet.ising].alpha))
        )
    end
    write_to!(fnames, :crit_out_unit, "    psifac      dpsi        q       singfac     eval1\n")
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
function ode_output_step!(odet::OdeState, intr::DconInternal, ctrl::DconControl, equil::Equilibrium.PlasmaEquilibrium, ffit::FourFitVars, fnames::DconFileNames; force::Bool=false)

    # Compute and print critical data for each time step
    ode_output_monitor!(odet, intr, ctrl, equil, fnames)
    # if dout.out_evals || dout.bin_evals
    #     ode_output_get_evals(intr, ctrl, dout, fNames, equil, odet) # this is just outputs
    # end

    # Write solutions to euler.bin
    # TODO: for now, always writing to euler.bin every step. Can add option to write at all (through outp data struct)
    # or multiples of steps through eulers_stride. However, we don't have an istep in this version, so might need to do
    # multiples of psi instead of integer steps?
    # if bin_euler && (mod(odet.istep, euler_stride) == 0 || force)
    du = zeros(ComplexF64, intr.mpert, intr.mpert, 2)
    params = (ctrl, equil, intr, odet, ffit, fnames)
    sing_der!(du, odet.u, params, odet.psifac)
    write_to!(fnames, :euler_bin_unit, 1)
    write_to!(fnames, :euler_bin_unit, odet.psifac, odet.q, odet.msol)
    write_to!(fnames, :euler_bin_unit, odet.u)
    write_to!(fnames, :euler_bin_unit, odet.ud)
    # end

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
    ode_output_monitor!(odet, intr, ctrl, equil, fnames)

Monitor the evolution of a critical eigenvalue (`crit`) during ODE integration and detect zero crossings, which indicate resonant or singular behavior.
The function evaluates `crit` using `ode_output_get_crit`, and if a sign change is found, it estimates the crossing point via linear interpolation.
If the crossing satisfies sharpness and consistency conditions, it's logged and `nzero` is incremented.
A termination flag (`termbycross_flag`) set within DconControl can be used to stop integration upon detection.

### Notes
- Zero-crossing detection is based on changes in the sign of `crit`.
- `u_med` is constructed as a linear interpolation between current and saved `u`.
- This function updates psi_save, crit_save, u_save, and nzero in the `odet` state.

### TODO
- Restore or redirect output to appropriate logging units.
- Replace `error(...)` with graceful shutdown if zero crossing is an intended exit condition.
"""
#function ode_output_monitor!(odet::OdeState, intr::DconInternal, ctrl::DconControl, fNames::DconFileNames, equil::Equilibrium.PlasmaEquilibrium)#, sVars::SingVars)
function ode_output_monitor!(odet::OdeState, intr::DconInternal, ctrl::DconControl, equil::Equilibrium.PlasmaEquilibrium, fnames::DconFileNames)

    # Compute new crit
    dpsi = odet.psifac - odet.psi_save
    q, singfac, logpsi1, logpsi2, crit = ode_output_get_crit(odet.psifac, odet.u, intr.mpert, odet.m1, ctrl.nn, equil.sq)

    # Check for zero crossing
    if crit * odet.crit_save < 0 # if crit changes sign, zero crossing has occurred
        fac = crit / (crit - odet.crit_save)
        psi_med = odet.psifac - fac * (odet.psifac - odet.psi_save)
        dpsi = psi_med - odet.psi_save
        u_med = odet.u .- fac .* (odet.u .- odet.u_save)
        q_med, singfac_med, logpsi1_med, logpsi2_med, crit_med = ode_output_get_crit(psi_med, u_med, intr.mpert, odet.m1, ctrl.nn, equil.sq)

        if (crit_med - crit) * (crit_med - odet.crit_save) < 0 &&
           abs(crit_med) < 0.5 * min(abs(crit), abs(odet.crit_save))
            # Very basic terminal printing for now
            println("Zero crossing detected at psi = $psi_med, q = $q_med")
            # println(term_unit, "Zero crossing at psi = $psi_med, q = $q_med")
            write_to!(fnames, :dcon_unit, @sprintf("Zero crossing at psi = %10.3e, q = %10.3e", psi_med, q_med))
            write_to!(fnames, :crit_out_unit, @sprintf("Zero crossing at psi = %10.3e, q = %10.3e", psi_med, q_med))
            write_to!(fnames, :crit_out_unit, @sprintf("%11.3e%11.3e%11.3e%11.3e%11.3e", psi_med,dpsi,q_med,singfac_med,crit_med))
            # write(crit_bin_unit, Float32(psi_med), Float32(logpsi1_med),
            #       Float32(logpsi2_med), Float32(q_med), Float32(crit_med))
            odet.nzero += 1
        end
        if ctrl.termbycross_flag
            #TODO: not sure if this is really an error, or just a user specified termination condition. So the program shoudln't error out, just stop integrating
            # In fortran, these were both handled by program_stop()
            error("Terminated by zero crossing.")
        end
    end

    # Write new crit
    write_to!(fnames, :crit_out_unit, @sprintf("%11.3e%11.3e%11.3e%11.3e%11.3e", odet.psifac,dpsi,odet.q,singfac,crit))
    # write(crit_bin_unit, Float32(psifac), Float32(logpsi1),
    #       Float32(logpsi2), Float32(q), Float32(crit))

    # Update saved values
    odet.psi_save = odet.psifac
    odet.crit_save = crit
    odet.u_save .= odet.u
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