function free_run(odet::OdeState, ctrl::DconControl, intr::DconInternal, equil::Equilibrium.PlasmaEquilibrium, ffit::FourFitVars, outp::DconOutput; op_netcdf_out::Bool=false)

    # TODO: it looks like vac_memory is always true - remove all ahg things and just assume true?
    vac_memory = true
    # TODO: this is always true in fortran - just get rid of it?
    normalize = true
    # Flags used within VACUUM
    complex_flag = true
    wall_flag = false
    ahg_file = "ahg2msc_dcon.out" # Deprecated

    # Allocations
    star = fill(' ', intr.mpert, odet.msol)
    ep = zeros(ComplexF64, intr.mpert)
    ev = zeros(ComplexF64, intr.mpert)
    et = zeros(ComplexF64, intr.mpert)
    tt = zeros(ComplexF64, intr.mpert)
    wv = zeros(ComplexF64, intr.mpert, intr.mpert)
    wt = zeros(ComplexF64, intr.mpert, intr.mpert)
    wt0 = zeros(ComplexF64, intr.mpert, intr.mpert)
    wp = zeros(ComplexF64, intr.mpert, intr.mpert)
    temp = zeros(ComplexF64, intr.mpert, intr.mpert)
    wpt = zeros(ComplexF64, intr.mpert, intr.mpert)
    wvt = zeros(ComplexF64, intr.mpert, intr.mpert)

    mtheta = length(equil.rzphi.ys)

    # Evaluate dV/dpsi at the plasma edge
    v1 = Spl.spline_eval(equil.sq, intr.psilim, 0)[3]

    # Compute plasma response matrix.
    if ctrl.ode_flag
        temp = adjoint(odet.u[:, 1:intr.mpert, 1])
        wp .= adjoint(odet.u[:, 1:intr.mpert, 2])
        # Compute wp using LU decomposition
        temp_fact = lu(temp)
        wp .= temp_fact \ wp
        wp .= adjoint(wp) / equil.psio^2
    end

    # Write file for mscvac
    # TODO: can likely remove last two arguments, ahgstr_op is deprecated
    # TODO: actually, can probably remove this function entirely and just call set_dcon_params directly
    free_write_msc(intr.psilim, ctrl, intr, equil; inmemory_op = vac_memory, ahgstr_op = ahg_file)

    # TODO: there is some ahb_flag logic here that has a comment "must be false if using GPEC"
    # I am assuming this means we don't have to implement any of it
    # if ahb_flag
    #     free_ahb_prep(wp, nmat, smat, asmat, bsmat, csmat, ipiva)
    # end

    # Compute vacuum response matrix.
    grri = Array{Float64}(undef, 2*(ctrl.mthvac+5), intr.mpert*2)
    xzpts = Array{Float64}(undef, ctrl.mthvac+5, 4)

    farwal_flag = true
    kernelsignin = -1.0
    # TODO: make this a ! function
    VacuumMod.mscvac(wv, intr.mpert, mtheta, ctrl.mthvac, complex_flag, kernelsignin,
           wall_flag, farwal_flag, grri, xzpts, ahg_file, intr.dir_path)

    # TODO: assuming bin_vac is deprecated for good, will remove all calls later after checking
    # if bin_vac
    #     @warn "!! WARNING: Use of vacuum.bin is deprecated in GPEC. Set bin_vac = false in dcon.in to reduce file IO."
    #     bin_open(vac_unit, "vacuum.bin", "UNKNOWN", "REWIND", "none")
    #     write(vac_unit, grri)
    # end

    kernelsignin = 1.0
    VacuumMod.mscvac(wv, intr.mpert, mtheta, ctrl.mthvac, complex_flag, kernelsignin,
           wall_flag, farwal_flag, grri, xzpts, ahg_file, intr.dir_path)
    # if bin_vac
    #     write(vac_unit, grri)
    # end
    if ctrl.wv_farwall_flag
        temp .= wv
    end

    farwal_flag = false
    kernelsignin = -1.0
    VacuumMod.mscvac(wv, intr.mpert, mtheta, ctrl.mthvac, complex_flag, kernelsignin,
           wall_flag, farwal_flag, grri, xzpts, ahg_file, intr.dir_path)
    # if bin_vac
    #     write(vac_unit, grri)
    # end

    kernelsignin = 1.0
    VacuumMod.mscvac(wv, intr.mpert, mtheta, ctrl.mthvac, complex_flag, kernelsignin,
           wall_flag, farwal_flag, grri, xzpts, ahg_file, intr.dir_path)
    # if bin_vac
    #     write(vac_unit, grri)
    #     write(vac_unit, xzpts)
    #     bin_close(vac_unit)
    # end

    if ctrl.wv_farwall_flag
        wv .= temp
    end

    # Scale vacuum matrix by singfac = (m - nn*qlim)
    singfac = (intr.mlow .- ctrl.nn .* intr.qlim) .+ collect(0:intr.mpert-1)
    for ipert in 1:intr.mpert
        wv[ipert, :] .*= singfac[ipert]
        wv[:, ipert] .*= singfac[ipert]
    end

    # Compute complex energy eigenvalues.
    wt .= wp .+ wv
    wt0 .= wt

    # use Eigen decomposition for general complex matrix (get left & right eigenvectors)
    # For general complex: use eigen(wt) which returns eigenvalues and eigenvectors (right)
    # For left eigenvectors we can compute eigen(wt') and conj-transpose as needed,
    # but Fortran uses zgeev('V','V',...) returning both vl and vr.
    Ev = eigen(wt)  # Ev.values, Ev.vectors (columns are eigenvectors)
    # Ev.vectors are right eigenvectors; we need to reorder by magnitude using bubble on real parts of eigenvalues
    et .= Ev.values
    eindex = sortperm(real.(et), rev=true)

    tt .= et
    # rearrange wt columns to correspond to eigenvector reordering similar to Fortran
    for ipert in 1:intr.mpert
        wt[:, ipert] .= Ev.vectors[:, eindex[intr.mpert+1-ipert]]
        et[ipert] = tt[eindex[intr.mpert+1-ipert]]
    end

    # Normalize eigenfunction and energy.
    if normalize
        for isol in 1:intr.mpert
            norm = 0.0 + 0.0im
            for ipert in 1:intr.mpert, jpert in 1:intr.mpert
                norm += ffit.jmat[jpert - ipert + intr.mband + 1] * wt[ipert, isol] * conj(wt[jpert, isol])
            end
            norm /= v1
            wt[:, isol] ./= sqrt(norm)
            et[isol] /= norm
        end
    end

    # Normalize phase and label largest component.
    imax = 0
    for isol in 1:intr.mpert
        # get index of largest absolute component (first occurrence)
        imax = argmax(abs.(wt[:, isol]))
        phase = abs(wt[imax, isol]) / wt[imax, isol]
        wt[:, isol] .*= phase
        # Mark largest component with *, all others initalized to ' '
        star[imax, isol] = '*'
    end

    # Compute plasma and vacuum contributions.
    # wpt = wt' * wp * wt  ; wvt = wt' * wv * wt
    wpt .= adjoint(wt) * (wp * wt)
    wvt .= adjoint(wt) * (wv * wt)

    for ipert in 1:intr.mpert
        ep[ipert] = wpt[ipert, ipert]
        ev[ipert] = wvt[ipert, ipert]
    end

    plasma1 = ComplexF64(real(ep[1]), 0.0)
    vacuum1 = ComplexF64(real(ev[1]), 0.0)
    total1  = ComplexF64(real(et[1]), 0.0)

    # Write data for ahb and deallocate.
    # if ahb_flag
    #     free_ahb_write(nmat, smat, wt, et)
    #     # Fortran did DEALLOCATE(r,z,theta,dphi,thetas,project) - we assume they are module arrays
    #     # and will be GC'd or freed by dcon_dealloc below
    # end

    # Write to euler.h5
    if outp.write_euler_h5
        write_output(outp, :euler_h5, 3)
        write_output(outp, :euler_h5, ep)
        write_output(outp, :euler_h5, ev)
        write_output(outp, :euler_h5, et)
        write_output(outp, :euler_h5, wt)
        write_output(outp, :euler_h5, wt0)
        write_output(outp, :euler_h5, ctrl.wv_farwall_flag)
    end

    # Write to screen and copy to output.
    if ctrl.verbose
        println("Energies: plasma = ", real(ep[1]), ", vacuum = ", real(ev[1]),
                ", real = ", real(et[1]), ", imaginary = ", imag(et[1]))
    end

    # Write eigenvalues to file
    write_output(outp, :dcon_out, "\nTotal Energy Eigenvalues:")
    write_output(outp, :dcon_out, "\n   isol   plasma      vacuum   re total   im total\n")
    for isol in 1:intr.mpert
        write_output(outp, :dcon_out, @sprintf("%6d %11.3e %11.3e %11.3e %11.3e",
                    isol, real(ep[isol]), real(ev[isol]), real(et[isol]), imag(et[isol])))
    end
    write_output(outp, :dcon_out, "\n   isol   plasma      vacuum   re total   im total\n")

    # Write eigenvectors to file
    write_output(outp, :dcon_out, "Total Energy Eigenvectors:")
    m = intr.mlow .+ collect(0:intr.mpert-1)
    for isol in 1:intr.mpert
        write_output(outp, :dcon_out, "\n   isol   imax   plasma      vacuum   re total   im total\n")
        write_output(outp, :dcon_out, @sprintf("%6d %6d %11.3e %11.3e %11.3e %11.3e",
            isol, imax, real(ep[isol]), real(ev[isol]), real(et[isol]), imag(et[isol])))
        write_output(outp, :dcon_out, "\n  ipert     m      re wt      im wt      abs wt\n")
        for ipert in 1:intr.mpert
            write_output(outp, :dcon_out,
            @sprintf("%6d %6d %11.3e %11.3e %11.3e %s",
                ipert, m[ipert], real(wt[ipert, isol]), imag(wt[ipert, isol]), abs(wt[ipert, isol]), star[ipert, isol]))
        end
       write_output(outp, :dcon_out, "\n  ipert     m      re wt      im wt      abs wt\n")
    end

    # Write the plasma matrix.
    write_output(outp, :dcon_out, "Plasma Energy Matrix:\n")
    for isol in 1:intr.mpert
        write_output(outp, :dcon_out, "isol = $(isol), m = $(m[isol])")
        write_output(outp, :dcon_out, "\n  i     re wp        im wp        abs wp\n")
        for ipert in 1:intr.mpert
            write_output(outp, :dcon_out, @sprintf("%3d%13.5e%13.5e%13.5e",
                ipert, real(wp[ipert, isol]), imag(wp[ipert, isol]), abs(wp[ipert, isol])))
        end
        write_output(outp, :dcon_out, "\n  i     re wp        im wp        abs wp\n")
    end

    # Compute separate plasma and vacuum eigenvalues.
    # Right eigenvalues/vectors of wp
    Ev_wp = eigen(wp)
    ep .= Ev_wp.values
    eindex = sortperm(real.(ep), rev=true)
    tt .= ep
    for ipert in 1:intr.mpert
        wp[:, ipert] .= Ev_wp.vectors[:, eindex[intr.mpert+1-ipert]]
        ep[ipert] = tt[eindex[intr.mpert+1-ipert]]
    end
    # For vacuum (Hermitian) use eigen of Hermitian wv
    Ev_wv = eigen(Hermitian(wv, :U))
    ev .= Ev_wv.values

    # Optionally write netcdf file
    # TODO: move this to HDF5 or something else later, make it an input flag?
    # if op_netcdf_out
    #     dcon_netcdf_out(wp, wv, wt, wt0, ep, ev, et)
    # end

    return plasma1, vacuum1, total1
end

function free_write_msc(psifac::Float64, ctrl::DconControl, intr::DconInternal, equil::Equilibrium.PlasmaEquilibrium; inmemory_op::Union{Bool,Nothing}=nothing,
                        ahgstr_op::Union{String,Nothing}=nothing)

    # Defaults for optional arguments
    inmemory = isnothing(inmemory_op) ? false : inmemory_op
    ahgstr   = isnothing(ahgstr_op)   ? "ahg2msc_dcon.out" : ahgstr_op
    inmemory = true # TODO: remove the above, and modify the code logic so VACUUM is always in memory

    # Allocations
    theta_norm = Vector(equil.rzphi.ys)
    mtheta = length(theta_norm)
    angle = zeros(Float64, mtheta)
    r     = zeros(Float64, mtheta)
    z     = zeros(Float64, mtheta)
    delta = zeros(Float64, mtheta)
    rfac  = zeros(Float64, mtheta)

    # Compute output
    qa = Spl.spline_eval(equil.sq, psifac, 0)[4] # TODO: this had a deriv = 1 in Fortran, but not used?
    for itheta in 1:mtheta
        f = Spl.bicube_eval(equil.rzphi, psifac, theta_norm[itheta], 0)
        rfac[itheta]  = sqrt(f[1])
        angle[itheta] = 2π * (theta_norm[itheta] + f[2])
        delta[itheta] = -f[3] / qa
    end
    r .= equil.ro .+ rfac .* cos.(angle)
    z .= equil.zo .+ rfac .* sin.(angle)

    # Invert values for nn < 0
    n = ctrl.nn
    if ctrl.nn < 0
        qa    = -qa
        delta .= -delta
        n     = -n
    end

    # Pass all required values to VACUUM
    if inmemory
        VacuumMod.set_dcon_params(mtheta, intr.mlow, intr.mhigh, n, qa,
                        reverse(r), reverse(z), reverse(delta))
    else
        # TODO: this section contains ahg2msc file writing, which is deprecated, just remove
    end
end