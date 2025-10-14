"""
    ode_output_init(ctrl, equil, outp, intr, odet)

Write header info to output files at the start of the integration.
Performs similar output writing as the Fortran `ode_output_open`,
except we no longer need to open binary files.

### TODOs

Remove deprecated outputs
Combine spline unpacking if possible, too many extra lines
Replace `println` statements with logging to appropriate output files
Remove `mthsurf0` if deprecated
Place outputs in a Julia do loop for automatic file closing
"""
function ode_output_init(ctrl::DconControl, equil::Equilibrium.PlasmaEquilibrium, intr::DconInternal, odet::OdeState, outp::DconOutput)

    # TODO: mess with this to condense the number of write_output calls? Maybe allow it to pass in dicts
    # Write euler.h5 header info
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
        write_output(outp, :euler_h5, equil.psio; dsetname="equil/psio")
        write_output(outp, :euler_h5, equil.config.control.psilow; dsetname="equil/psilow")
        write_output(outp, :euler_h5, equil.config.control.power_b; dsetname="equil/power_b")
        write_output(outp, :euler_h5, equil.config.control.power_r; dsetname="equil/power_r")
        write_output(outp, :euler_h5, equil.config.control.power_bp; dsetname="equil/power_bp")
        write_output(outp, :euler_h5, 0; dsetname="equil/shotnum") # TODO: equil.params.shotnum)
        write_output(outp, :euler_h5, 0; dsetname="equil/shottime") # TODO: equil.params.shottime)

        # Write spline arrays
        write_output(outp, :euler_h5, Vector(equil.sq.xs); dsetname="splines/sq/xs")
        # TODO: getting errors when trying to dump just fs, so splitting for now, which adds so many lines
        #  This should be fixed if we separate these like Nik mentioned
        write_output(outp, :euler_h5, equil.sq.fs[:, 1]; dsetname="splines/sq/fs/2piF")
        write_output(outp, :euler_h5, equil.sq.fs[:, 2]; dsetname="splines/sq/fs/mu0p")
        write_output(outp, :euler_h5, equil.sq.fs[:, 3]; dsetname="splines/sq/fs/dVdpsi")
        write_output(outp, :euler_h5, equil.sq.fs[:, 4]; dsetname="splines/sq/fs/q")
        write_output(outp, :euler_h5, equil.sq.fs1[:, 1]; dsetname="splines/sq/fs1/2piF")
        write_output(outp, :euler_h5, equil.sq.fs1[:, 2]; dsetname="splines/sq/fs1/mu0p")
        write_output(outp, :euler_h5, equil.sq.fs1[:, 3]; dsetname="splines/sq/fs1/dVdpsi")
        write_output(outp, :euler_h5, equil.sq.fs1[:, 4]; dsetname="splines/sq/fs1/q")
        write_output(outp, :euler_h5, 0; dsetname="splines/sq/xpower") # TODO: equil.sq.xpower
        write_output(outp, :euler_h5, Vector(equil.rzphi.xs); dsetname="splines/rzphi/xs")
        write_output(outp, :euler_h5, Vector(equil.rzphi.ys); dsetname="splines/rzphi/ys")
        write_output(outp, :euler_h5, equil.rzphi.fs[:, 1]; dsetname="splines/rzphi/fs/rcoords")
        write_output(outp, :euler_h5, equil.rzphi.fs[:, 2]; dsetname="splines/rzphi/fs/offset")
        write_output(outp, :euler_h5, equil.rzphi.fs[:, 3]; dsetname="splines/rzphi/fs/nu")
        write_output(outp, :euler_h5, equil.rzphi.fs[:, 4]; dsetname="splines/rzphi/fs/jac")
        write_output(outp, :euler_h5, equil.rzphi.fsx[:, 1]; dsetname="splines/rzphi/fsx/rcoords")
        write_output(outp, :euler_h5, equil.rzphi.fsx[:, 2]; dsetname="splines/rzphi/fsx/offset")
        write_output(outp, :euler_h5, equil.rzphi.fsx[:, 3]; dsetname="splines/rzphi/fsx/nu")
        write_output(outp, :euler_h5, equil.rzphi.fsx[:, 4]; dsetname="splines/rzphi/fsx/jac")
        write_output(outp, :euler_h5, equil.rzphi.fsy[:, 1]; dsetname="splines/rzphi/fsy/rcoords")
        write_output(outp, :euler_h5, equil.rzphi.fsy[:, 2]; dsetname="splines/rzphi/fsy/offset")
        write_output(outp, :euler_h5, equil.rzphi.fsy[:, 3]; dsetname="splines/rzphi/fsy/nu")
        write_output(outp, :euler_h5, equil.rzphi.fsy[:, 4]; dsetname="splines/rzphi/fsy/jac")
        write_output(outp, :euler_h5, equil.rzphi.fsxy[:, 1]; dsetname="splines/rzphi/fsxy/rcoords")
        write_output(outp, :euler_h5, equil.rzphi.fsxy[:, 2]; dsetname="splines/rzphi/fsxy/offset")
        write_output(outp, :euler_h5, equil.rzphi.fsxy[:, 3]; dsetname="splines/rzphi/fsxy/nu")
        write_output(outp, :euler_h5, equil.rzphi.fsxy[:, 4]; dsetname="splines/rzphi/fsxy/jac")
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
    ode_output_step(ctrl::DconControl, equil::Equilibrium.PlasmaEquilibrium, intr::DconInternal, odet::OdeState, outp::DconOutput)

Performs output and monitoring tasks at each integration step by calling
`ode_output_monitor`. to track critical eigenvalue behavior. Unlike Fortran,
we don't implement any output file dumping here, as Julia's in-memory storage
structure handles this differently.

### TODOs

Depending on output needs, may need to implement `ode_output_get_evals` and `ode_output_sol`
Determine if this function is even necessary of if we should just call `ode_output_monitor` directly
"""
#TODO: depending on how we restructure our outputs, this function might be uncessary
# (i.e. if we don't need an ode_output_get_evals or ode_output_sol, can just replace calls with ode_output_monitor)
function ode_output_step(ctrl::DconControl, equil::Equilibrium.PlasmaEquilibrium, intr::DconInternal, odet::OdeState, outp::DconOutput)

    # Compute and print critical data for each time step
    ode_output_monitor(ctrl, equil, intr, odet, outp)
    # if dout.out_evals || dout.bin_evals
    #     ode_output_get_evals(intr, ctrl, dout, fNames, equil, odet) # this is just outputs
    # end

    # All bin_euler logic performed after integration

    # if outp.bin_sol
    #     ode_output_sol()
    # end
end

# TODO: I don't think this function is essentially for the minimal working version - can convert later if needed
# function ode_output_get_evals(intr::DconInternal, ctrl::DconControl, dout::DconOutput, fNames::DconFileNames, equil::Equilibrium.PlasmaEquilibrium, odet::OdeState)
#     return
# end

"""
    ode_output_monitor(ctrl::DconControl, equil::Equilibrium.PlasmaEquilibrium, intr::DconInternal, odet::OdeState, outp::DconOutput)

Monitor the evolution of a critical eigenvalue (`crit`) during integration
using `ode_output_get_crit` and detect zero crossings, which indicate instability.
If a sign change is found, it estimates the crossing point via linear interpolation.
If the crossing satisfies sharpness and consistency conditions, it's logged and
`nzero` is incremented. Performs the same function as `ode_output_monitor` in the
Fortran code, with small differences in output handling.

### TODO

Restore or redirect output to appropriate logging units.
Replace `error(...)` with graceful shutdown if zero crossing is an intended exit condition.
All the _prev variables can probably be removed and the logic can be simplified to just take the odet.step - 1 values when needed
"""
function ode_output_monitor(ctrl::DconControl, equil::Equilibrium.PlasmaEquilibrium, intr::DconInternal, odet::OdeState, outp::DconOutput)

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
            write_output(outp, :crit_out, @sprintf("%11.3e%11.3e%11.3e%11.3e%11.3e", psi_med, dpsi, q_med, singfac_med, crit_med))
            odet.nzero += 1
        end
        if ctrl.termbycross_flag
            #TODO: not sure if this is really an error, or just a user specified termination condition. So the program shoudln't error out, just stop integrating
            # In fortran, these were both handled by program_stop()
            error("Terminated by zero crossing.")
        end
    end

    # Write new crit
    write_output(outp, :crit_out, @sprintf("%11.3e%11.3e%11.3e%11.3e%11.3e", odet.psifac, dpsi, odet.q, singfac, crit))

    # Update saved values
    odet.psi_prev = odet.psifac
    odet.crit_prev = crit
    odet.u_prev .= odet.u
end

"""
    ode_output_get_crit(psi, u, mpert, m1, nn, sq) -> (q, singfac, logpsi1, logpsi2, crit)

Compute critical quantities at a given flux surface and using the smallest (in magnitude)
inverse eigenvalue in combination with the equilibrium profiles to form `crit`.
Performs the same function as `ode_output_get_crit` in the Fortran code.

### Arguments

  - `psi::Float64`: The flux surface at which to evaluate.
  - `u::Array{ComplexF64, 3}`: Solution matrix at `psi`.
  - `mpert::Int`: Number of poloidal modes considered.
  - `m1::Int`: Poloidal mode number of the perturbation.
  - `nn::Int`: Toroidal mode number of the perturbation.
  - `sq::Spl.CubicSpline`: Spline object containing equilibrium profiles.

### Returns

A tuple `(q, singfac, logpsi1, logpsi2, crit)` of critical data for the time step.

### TODOs

Add dVdpsi multiplication back once equilibrium bug is fixed
Decide if we want to just pass in the relevant quantities instead of structs for functions like this
"""
function ode_output_get_crit(psi::Float64, u::Array{ComplexF64,3}, mpert::Int, m1::Int, nn::Int, sq::Spl.CubicSpline)

    # Compute inverse plasma response matrix
    uu = u[:, 1:mpert, :]
    wp = adjoint(uu[:, :, 1])
    temp = adjoint(uu[:, :, 2])

    # Compute wp using LU decomposition
    temp_fact = lu(temp)
    wp = temp_fact \ wp
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