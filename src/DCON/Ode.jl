

using LinearAlgebra

using JPEC

# Placeholder for dependencies
# using Fourfit
# using OdeOutput
# using Debug
# using Free

# For now, we set up the equilibrium spline locally, but this should be replaced with a proper import
plasma_eq = JPEC.Equilibrium.setup_equilibrium(equil_input)
sq = plasma_eq.sq

const diagnose_ca = false
const eps = 1e-10

# Global variables (should be refactored into structs for safety)
new = true
next = ""
flag_count = 0

# Example types for variables (adjust as needed)
# These should be replaced with proper structs and types
@kwdef struct OdeState(mpert::Int, msol::Int)
    istate::Int=1
    istep::Int=0
    liw::Int=20
    lrw::Int=22+64*mpert*msol
    iopt::Int=1
    itask::Int=5
    itol::Int=2
    mf::Int=10
    ix::Int=0
    iwork::Vector{Int}=zeros(Int, 20)
    rwork::Vector{Float64}=zeros(Float64, 22+64*mpert*msol)
    atol::Array{ComplexF64,3}=zeros(ComplexF64, mpert, msol, 2)
    index::Vector{Int}=collect(1:mpert)
    unorm::Vector{Float64}=zeros(Float64, 2*mpert)
    unorm0::Vector{Float64}=zeros(Float64, 2*mpert)
    fixfac::Array{ComplexF64,2}=zeros(ComplexF64, msol, msol)
end

function init_ode_state(mpert, msol)
    return OdeState(mpert, msol)
end


function ode_run(ctrl::DconControl, equil::DconEquilibrium, intr::DconInternal)
    # Initialization
    odet = init_ode_state(intr.mpert, intr.msol)

    if ctrl.sing_start <= 0
        ode_axis_init(ctrl, equil, intr)
    elseif ctrl.sing_start <= intr.msing
        @warn "sing_start = 0 not implemented yet!"
        ode_sing_init()
    else
        message = "sing_start = $(ctrl.sing_start) > msing = $(intr.msing)"
        program_stop(message)
    end
    flag_count = 0
    #ode_output_open() # TODO: have to handle io
    if diagnose_ca
        ascii_open(ca_unit, "ca.out", "UNKNOWN")
    end

    # Integration loop
    while true
        first = true
        while true
            if odet.istep > 0
                ode_unorm(false)
            end
            # Always record the first and last point in an inter-rational region
            test = ode_test()
            force_output = first || test
            ode_output_step(unorm; op_force=force_output)
            ode_record_edge()
            if test
                break
            end
            ode_step()
            first = false
        end

        # Re-initialize
        if ising == ksing
            break
        end
        if next == "cross"
            if res_flag
                ode_resist_cross()
            elseif kin_flag
                ode_kin_cross()
            else
                ode_ideal_cross()
            end
        else
            break
        end
        flag_count = 0
    end

    # Finalize
    ode_output_close()
    # Deallocate arrays (in Julia, set to nothing or let GC handle)
    rwork = nothing
    atol = nothing
    unorm0 = nothing
    unorm = nothing
    index = nothing
    fixfac = nothing
    if diagnose_ca
        ascii_close(ca_unit)
    end
end


# Example stub for axis initialization
"""
    ode_axis_init(sing_surf_data, sq; nn = 1, ψlim = 0.9936, ψlow = 0.01, mlow = -12, mhigh = 20, singfac_min = 1e-5, qlow = 0.0, sort_type = "absm")

Initializes the ODE system near the magnetic axis for MHD or kinetic-MHD stability
analysis. This routine prepares solution vectors and sorting indices, locates
the relevant singular surfaces in the plasma equilibrium based on input data,
and sets integration boundaries and normalization factors for subsequent ODE integration
through singularities.

# Arguments
- `sing_surf_data`: Array of structures, each containing information about singular surfaces,
    including fields like `.ψfac`, `.q`, and `.q1`. Filled from sing_find function
- `sq`: Spline object controlling the equilibrium safety factor profile (typically contains equilibrium profile splines).

# Keyword Arguments (these are likely to be replaced by the user input file structure)
- `nn`: Integer. Toroidal mode number.
- `ψlim`: Float. Upper bound for normalized flux coordinate ψ.
- `ψlow`: Float. Lower bound for normalized flux coordinate ψ.
- `mlow`: Integer. Lowest poloidal mode number to consider.
- `mhigh`: Integer. Highest poloidal mode number to consider.
- `singfac_min`: Float. Minimum separation from a singular surface.
- `qlow`: Float. Lower bound for safety factor q. (Currently not implemented if > 0.)
- `sort_type`: String. Sorting criteria for mode numbers (`"absm"` for |m|, or
    `"sing"` for proximity to the first singular surface's m).

# Details
This function does the following:
- Determines relevant singular surfaces from the provided equilibrium/singularity data.
- Sets up arrays for the ODE solution and related workspaces for all poloidal harmonics in the range [`mlow`, `mhigh`].
- Sorts the initial harmonics according to the specified `sort_type`.
- Initializes normalization factors and solution arrays for use in the ODE integrator.
- Computes the initial offset (`singfac`) from the next singular surface.

Several features for kinetic MHD (indicated by `kin_flag`) or for `qlow > 0` are noted but not yet implemented.

"""

#function ode_axis_init(sing_surf_data, plasma_eq; nn = 1, ψlim=0.9936, ψlow = 0.01, mlow = -12, mhigh = 20, singfac_min = 1e-5, qlow = 0.0, sort_type = "absm")
function ode_axis_init(ctrl::DconControl, equil::DconEquilibrium, intr::DconInternal)

    # JMH - I think all default values here will be supplied by inputs
    # Right now, there are a lot of globals - not sure how we want to handle these
    # But at least the major outputs here match fortran dcon

    # This might be wrong - double check once equil is more defined. 
    qval(psi) = JPEC.SplinesMod.spline_eval(equil.sq, psi, 0)[4]
    q1val(psi) = JPEC.SplinesMod.spline_eval(equil.sq, psi, 1)[2][4]

    # Preliminary computations
    psifac = equil.psilow  # TODO: this was Fortran sq%xs(0)? - confirm this is psilow? Don't know how to access xs 

    # Use Newton iteration to find starting psi if qlow is above q0
    # TODO: add this. For now, assume qlow = 0.0 and start at the first singular surface
    if ctrl.qlow > 0.0
        @warn "qlow > 0.0 is not implemented in ode_axis_init yet"
    end

    # Find inner singular surface
    ising = 0
    if false #(TODO: kin_flag)
        # for ising = 1:kmsing
        #     if kinsing[ising].psifac > psifac
        #         break
        #     end
        # end
    else
        for ising in 1:intr.msing
            if intr.sing[ising].psifac > psifac
                break
            end
        end
    end
    ising = max(0, ising - 1)

    # Find next singular surface
    if false # TODO: (kin_flag)
        # while true
        #     ising += 1
        #     if ising > kmsing
        #         break
        #     end
        #     if psilim < kinsing[ising].psifac
        #         break
        #     end
        #     q = kinsing[ising].q
        #     if mlow <= nn * q && mhigh >= nn * q
        #         break
        #     end
        # end
        # if ising > kmsing || singfac_min == 0
        #     psimax = psilim * (1 - eps)
        #     next = "finish"
        # elseif psilim < kinsing[ising].psifac
        #     psimax = psilim * (1 - eps)
        #     next = "finish"
        # else
        #     psimax = kinsing[ising].psifac - singfac_min / abs(nn * kinsing[ising].q1)
        #     next = "cross"
        # end
    else
        while true
            ising += 1
            if ising > intr.msing || psilim < intr.sing[min(ising, msing)].psifac
                break
            end
            q = intr.sing[ising].q
            if intr.mlow <= ctrl.nn * q && intr.mhigh >= ctrl.nn * q
                break
            end
        end
        if ising > intr.msing || psilim < intr.sing[min(ising, intr.msing)].psifac || ctrl.singfac_min == 0
            psimax = psilim * (1 - eps)
            next = "finish"
        else
            psimax = intr.sing[ising].psifac - ctrl.singfac_min / abs(ctrl.nn * intr.sing[ising].q1)
            next = "cross"
        end
    end

    # Allocate and sort solutions by increasing value of |m-ms1|
    # TODO: THESE PROBABLY SHOULDN'T BE GLOBAL - WHERE DO THEY GO? 
    global u = zeros(ComplexF64, intr.mpert, intr.mpert, 2)
    global du = zeros(ComplexF64, intr.mpert, intr.mpert, 2)
    global u_save = zeros(ComplexF64, intr.mpert, intr.mpert, 2)
    global unorm0 = zeros(Float64, 2 * intr.mpert)
    global unorm = zeros(Float64, 2 * intr.mpert)
    global index = collect(1:intr.mpert)
    m = intr.mlow - 1 .+ index
    if ctrl.sort_type == "absm"
        key = abs.(m)
    elseif ctrl.sort_type == "sing"
        key = m
        if msing > 0
            key = key .- sing[1].m
        end
        key = -abs.(key)
    else
        error("Cannot recognize sort_type = $ctrl.sort_type")
    end
    bubble!(key, index, 1, intr.mpert) # in original Fortran: bubble(key, index, 1, mpert)

    # Initialize solutions
    for ipert = 1:mpert
        u[index[ipert], ipert, 2] = 1
    end

    # TODO: I THINK THESE ARE IN THE SINGVARS STRUCT - NEED TO ADD AT SOME POINT
    global msol = intr.mpert
    global neq = 4 * intr.mpert * msol
    u_save .= u
    global psi_save = psifac

    # Initialize integrator parameters (JMH - not sure if these are needed in Julia)
    # global istep = 0
    # global iopt = 1
    # global itask = 5
    # global itol = 2
    # global mf = 10
    # global istate = 1

    # Compute conditions at next singular surface
    q = qval(equil.psilow) # Fortran: q=sq%fs(0,4) # IS THIS RIGHT? 
    if false #TODO: (kin_flag)
        # if kmsing > 0
        #     m1 = round(Int, nn * kinsing[ising].q)
        # else
        #     m1 = round(Int, nn * qlim) + sign(one, nn * sq.fs1[mpsi, 5])
        # end
    else
        if intr.msing > 0
            m1 = round(Int, ctrl.nn * intr.sing[ising].q)
        else
            m1 = round(Int, ctrl.nn * intr.qlim) + sign(ctrl.nn * q1val[end])
        end
    end
    # TODO: not sure where this is used
    global singfac = abs(m1 - ctrl.nn * q)

    # Set up work arrays (JMH - not sure if these are needed in Julia)
    # global liw = length(iwork)
    # global lrw = 22 + 64 * mpert * msol
    # global rwork = zeros(Float64, lrw)
    # global atol = zeros(ComplexF64, mpert, msol, 2)
    # global fixfac = zeros(ComplexF64, msol, msol)
    # iwork = 0
    # rwork = 0
    # rwork[1] = psimax
    # rwork[5] = psifac * 1e-3
    # rwork[11] = rwork[5]
end

function ode_sing_init()
    # Declare and initialize local variables
    star = fill(' ', mpert)
    # local variables
    ua = Array{Complex{r8}}(undef, mpert, 2*mpert, 2)
    dpsi = 0.0
    new = true

    ising = sing_start
    dpsi = singfac_min/abs(nn*sing[ising].q1)*10
    psifac = sing[ising].psifac + dpsi
    q = sing[ising].q + dpsi*sing[ising].q1

    # Allocate and initialize solution arrays
    u      = zeros(Complex{r8}, mpert, mpert, 2)
    du     = zeros(Complex{r8}, mpert, mpert, 2)
    u_save = zeros(Complex{r8}, mpert, mpert, 2)
    unorm0 = zeros(r8, 2*mpert)
    unorm  = zeros(r8, 2*mpert)
    index  = zeros(Int, 2*mpert)

    if old_init
        u .= 0.0
        for ipert ∈ 1:mpert
            u[ipert, ipert, 2] = 1.0
        end
    else
        sing_get_ua(ising, psifac, ua)
        # Big slice: u = ua[:, mpert+1:2*mpert,:]
        for i = 1:mpert, j = 1:mpert, k = 1:2
            u[i, j, k] = ua[i, mpert+j, k]
        end
    end
    u_save .= u
    psi_save = psifac
    msol = mpert
    neq = 4 * mpert * msol

    # Diagnose output: for demonstration, I'll use console prints only
    if diagnose
        sing_der(neq, psifac, u, du)
        file = ascii_open("init.out", "w")
        println(file, "Output from ode_sing_init")
        println(file, "mlow mhigh mpert q psifac dpsi order")
        println(file, mlow, " ", mhigh, " ", mpert, " ", q, " ", psifac, " ", dpsi, " ", sing_order)
        ipert0 = sing[ising].m - mlow + 1
        star = fill(' ', mpert)
        star[ipert0] = '*'
        m = mlow
        for isol in 1:msol
            println(file, "isol = $isol, m = $m", star[isol])
            # The original Fortran code prints a header and then u and du for each ipert.
            for ipert in 1:mpert
                println(file, "$ipert ", u[ipert, isol, 1], " ", u[ipert, isol, 2],
                        " ", du[ipert, isol, 1], " ", du[ipert, isol, 2], " ",
                        star[ipert])
            end
            m += 1
        end
        ascii_close(file)
        program_stop("Termination by ode_sing_init.")
    end

    # Compute conditions at next singular surface
    while true
        ising += 1
        if ising > msing || psilim < sing[ising ≤ msing ? ising : msing].psifac
            break
        end
        q = sing[ising].q
        if mlow ≤ nn*q && mhigh ≥ nn*q
            break
        end
    end

    if ising > msing || psilim < sing[ising ≤ msing ? ising : msing].psifac
        m1 = round(Int, nn*qlim) + round(Int, sign(one, nn*sq.fs1[mpsi, 4]))
        psimax = psilim * (1-eps)
        next_ = "finish"
    else
        m1 = round(Int, nn*sing[ising].q)
        psimax = sing[ising].psifac - singfac_min/abs(nn*sing[ising].q1)
        next_ = "cross"
    end

    # Set up integrator parameters
    istep = 0
    istate = 1
    iopt = 1
    itask = 5
    itol = 2
    mf = 10
    istate = 1

    # Set up work arrays
    liw = length(iwork)    # iwork must exist!
    lrw = 22 + 64*mpert*msol
    rwork = zeros(r8, lrw)
    atol  = zeros(r8, mpert, msol, 2)
    fixfac = zeros(r8, msol, msol)
    iwork = zeros(Int, liw)  # or readjust as needed
    rwork[1] = psimax
    rwork[5] = dpsi
    rwork[11] = rwork[5]

    # Terminate, or in Julia just return (no need for RETURN)
    return nothing  # or could return a struct with all these values, for a more Julian approach
end

function ode_ideal_cross()
    # ...existing code...

    # Fixup solution at singular surface
    if verbose
        println("psi = $(sing[ising].psifac), q = $(sing[ising].q)")
    end
    ode_unorm(true)

    # Diagnose solution before reinitialization
    if diagnose
        sing_get_ca(ising, psifac, u, ca)
        init_out_unit = ascii_open("reinit.out", "UNKNOWN")
        println(init_out_unit, "Output from ode_ideal_cross, singfac = $singfac")
        println(init_out_unit, "Asymptotic coefficients matrix before reinit:")
        for jsol in 1:msol
            println(init_out_unit, "jsol = $jsol, m = $(mlow + jsol - 1)")
            # Print ca values for each ipert
            for ipert in 1:mpert
                println(init_out_unit, ca[ipert, jsol, 1], abs(ca[ipert, jsol, 1]),
                        ca[ipert, jsol, 2], abs(ca[ipert, jsol, 2]))
            end
        end
    end

    # Write asymptotic coefficients before reinit
    if bin_euler
        sing_get_ca(ising, psifac, u, ca)
        write(euler_bin_unit, 4)
        write(euler_bin_unit, sing[ising].psifac, sing[ising].q, sing[ising].q1)
        write(euler_bin_unit, msol)
        write(euler_bin_unit, ca)
    end

    # Re-initialize
    psi_old = psifac
    ipert0 = round(Int, nn * sing[ising].q) - mlow + 1
    dpsi = sing[ising].psifac - psifac
    psifac = sing[ising].psifac + dpsi
    sing_get_ua(ising, psifac, ua)
    if !con_flag
        u[:, index[1], :] .= 0  # originally u(ipert0,:,:) = 0
    end
    sing_der(neq, psi_old, u, du1)
    sing_der(neq, psifac, u, du2)
    u .= u .+ (du1 .+ du2) .* dpsi
    if !con_flag
        u[ipert0, :, :] .= 0
        u[:, index[1], :] .= ua[:, ipert0 + mpert, :]
    end

    # Diagnose solution after reinitialization
    if diagnose
        sing_get_ca(ising, psifac, u, ca)
        println(init_out_unit, "Asymptotic coefficients matrix after reinit:")
        for jsol in 1:msol
            println(init_out_unit, "jsol = $jsol, m = $(mlow + jsol - 1)")
            for ipert in 1:mpert
                println(init_out_unit, ca[ipert, jsol, 1], abs(ca[ipert, jsol, 1]),
                        ca[ipert, jsol, 2], abs(ca[ipert, jsol, 2]))
            end
        end
        ascii_close(init_out_unit)
        program_stop("Termination by ode_ideal_cross.")
    end

    # Write asymptotic coefficients after reinit
    if bin_euler
        sing_get_ca(ising, psifac, u, ca)
        write(euler_bin_unit, msol)
        write(euler_bin_unit, ca)
        write(euler_bin_unit, sing[ising].restype.e, sing[ising].restype.f,
              sing[ising].restype.h, sing[ising].restype.m,
              sing[ising].restype.g, sing[ising].restype.k,
              sing[ising].restype.eta, sing[ising].restype.rho,
              sing[ising].restype.taua, sing[ising].restype.taur)
    end

    # Find next ising
    while true
        ising += 1
        if ising > msing || psilim < sing[min(ising, msing)].psifac
            break
        end
        q = sing[ising].q
        if mlow <= nn * q && mhigh >= nn * q
            break
        end
    end

    # Compute conditions at next singular surface
    if ising > msing || psilim < sing[min(ising, msing)].psifac
        psimax = psilim * (1 - eps)
        m1 = round(Int, nn * qlim) + sign(one, nn * sq.fs1[mpsi, 4])
        next = "finish"
    else
        psimax = sing[ising].psifac - singfac_min / abs(nn * sing[ising].q1)
        m1 = round(Int, nn * sing[ising].q)
        println(crit_out_unit, "ising = $ising, psifac = $(sing[ising].psifac), q = $(sing[ising].q), di = $(sing[ising].di), alpha = $(sing[ising].alpha)")
    end

    # Restart ode solver
    istate = 1
    istep += 1
    rwork[1] = psimax
    new = true
    u_save .= u
    psi_save = psifac

    # Write to files
    println(crit_out_unit, "step info")
    if crit_break
        write(crit_bin_unit)
    end

    # ...existing code...
end

# Example stub for kinetic crossing
function ode_kin_cross()
    # Implement kinetic crossing logic here
    return
end

# Example stub for resistive crossing
function ode_resist_cross()
    # Implement resistive crossing logic here
    return
end

function ode_step()
    # ODE integrator step
    # Compute relative tolerances
    singfac_local = typemax(Float64)
    if kin_flag
        if ising == 1 && kmsing >= 1
            singfac_local = abs(psifac - kinsing[ising].psifac) /
                            (kinsing[ising].psifac - psilow)
        elseif ising <= kmsing
            singfac_local = min(
                abs(psifac - kinsing[ising].psifac),
                abs(psifac - kinsing[ising-1].psifac)
            ) / abs(kinsing[ising].psifac - kinsing[ising-1].psifac)
        end
    else
        if ising <= msing
            singfac_local = abs(sing[ising].m - nn * q)
        end
        if ising > 1
            singfac_local = min(
                singfac_local,
                abs(sing[ising-1].m - nn * q)
            )
        end
    end
    tol = singfac_local < crossover ? tol_r : tol_nr
    rtol = tol

    # Compute absolute tolerances
    for ieq in 1:2, isol in 1:msol
        atol0 = maximum(abs.(u[:, isol, ieq])) * tol
        if atol0 == 0
            atol0 = typemax(Float64)
        end
        atol[:, isol, ieq] .= ComplexF64(atol0, atol0)
    end

    # Choose psiout
    if node_flag
        while psifac < sq.xs[ix]
            ix += 1
        end
        psiout = sq.xs[ix]
        psiout = min(psiout, psimax)
        rwork[1] = psiout
    else
        psiout = psimax
    end

    # Advance differential equations
    global istep
    istep += 1

    # Use DifferentialEquations.jl for general ODE solving in Julia
    using DifferentialEquations

    # Define the ODE function in the DifferentialEquations.jl format
    function ode_func!(du, u, p, t)
        sing_der(neq, t, u, du)
    end

    # Set up the problem
    u0 = copy(u)
    tspan = (psifac, psiout)
    prob = ODEProblem(ode_func!, u0, tspan)

    # Set tolerances
    abstol = maximum(abs.(atol))
    reltol = rtol

    # Solve the ODE
    sol = solve(prob, abstol=abstol, reltol=reltol)

    # Update u and psifac with the solution at the end of the interval
    u .= sol.u[end]
    psifac = sol.t[end]

    # Diagnose error
    if rwork[11] < 1e-14 && diagnose
        dt = rwork[11]
        sing_der(neq, psifac, u, du)
        dewset(neq, itol, rtol, atol, u, ewt)
        err = du .* dt ./ ewt
        errabs = abs.(err)
        errmax = maximum(errabs)
        errloc = ind2sub(size(errabs), argmax(errabs))
        ipert, isol, ieq = errloc
        ewtmax = ewt[ipert, isol, ieq]
        msg = """
        Termination by ode_step
        ipert = $ipert, ieq = $ieq, isol = $isol, msol = $msol
        errmax = $(errmax), ewt = $(ewtmax), atol = $(abs(atol[ipert, isol, ieq]))
        """
        program_stop(msg)
    end

end

function ode_unorm(sing_flag::Bool)
    # Compute norms of first solution vectors, abort if any are zero
    unorm[1:mpert] .= sqrt.(sum(abs.(u[:, 1:mpert, 1]).^2, dims=1)[:])
    unorm[mpert+1:msol] .= sqrt.(sum(abs.(u[:, mpert+1:msol, 2]).^2, dims=1)[:])
    if minimum(unorm[1:msol]) == 0
        jmax = argmin(unorm[1:msol])
        message = "_unorm: unorm(1,$jmax) = 0"
        program_stop(message)
    end

    # Normalize unorm and perform Gaussian reduction if required
    if new
        new = false
        unorm0[1:msol] .= unorm[1:msol]
    else
        unorm[1:msol] .= unorm[1:msol] ./ unorm0[1:msol]
        uratio = maximum(unorm[1:msol]) / minimum(unorm[1:msol])
        if uratio > ucrit || sing_flag
            ode_fixup(sing_flag, false)
            if diagnose
                ode_test_fixup()
            end
            new = true
        end
    end

    return
end

function ode_fixup(sing_flag::Bool, test::Bool)
    # ...existing code...

    diagnose = false
    secondary = false
    # new is a persistent variable in Fortran; use a global or Ref in Julia if needed
    global new
    mask = trues(2, msol)
    jmax = zeros(Int, 1)

    # Open output file and write initial values
    if diagnose
        if new
            ascii_open(init_out_unit, "fixup.out", "UNKNOWN")
            new = false
        end
        println(init_out_unit, "mlow = $mlow, mhigh = $mhigh, mpert = $mpert, msol = $msol, psifac = $psifac, q = $q")
        println(init_out_unit, "input values:")
        for isol in 1:msol
            println(init_out_unit, "isol = $isol, m = $(mlow + isol - 1)")
            # Print first two components
            for ipert in 1:2
                println(init_out_unit, "re u($ipert), im u($ipert)")
            end
            # Normalize and print
            norm_val = maximum(abs.(u[:, isol, 1]))
            for ipert in 1:mpert
                println(init_out_unit, ipert, u[ipert, isol, :] / norm_val)
            end
            for ipert in 1:2
                println(init_out_unit, "re u($ipert), im u($ipert)")
            end
        end
    end

    # Initial output
    println(crit_out_unit, "Gaussian Reduction at istep = $istep, psi = $psifac, q = $q")
    if !sing_flag
        println(crit_out_unit, "Gaussian Reduction at istep = $istep, psi = $psifac, q = $q")
    end
    istate = 1
    flag_count = 0

    # Initialize fixfac
    if !test
        fixfac .= 0
        for isol in 1:msol
            fixfac[isol, isol] = 1
        end
        # Sort unorm
        index[1:msol] .= collect(1:msol)
        index[1:mpert] .= sortperm_subrange(unorm, 1:mpert) # in original Fortran: bubble(unorm, index, 1, mpert)
        index[mpert+1:msol] .= sortperm_subrange(unorm, mpert+1:msol) # in original Fortran: bubble(unorm, index, mpert + 1, msol)
    end

    # Triangularize primary solutions
    mask .= true
    for isol in 1:mpert
        ksol = index[isol]
        mask[2, ksol] = false
        if !test
            # Find max location
            absvals = abs.(u[:, ksol, 1])
            masked = absvals .* mask[1, 1:mpert]
            kpert = argmax(masked)
            mask[1, kpert] = false
        end
        for jsol in 1:msol
            if mask[2, jsol]
                if !test
                    fixfac[ksol, jsol] = -u[kpert, jsol, 1] / u[kpert, ksol, 1]
                end
                u[:, jsol, :] .= u[:, jsol, :] .+ u[:, ksol, :] .* fixfac[ksol, jsol]
                if !test
                    u[kpert, jsol, 1] = 0
                end
            end
        end
    end

    # Triangularize secondary solutions
    if msol > mpert && secondary
        mask .= true
        for isol in mpert + 1:msol
            ksol = index[isol]
            mask[2, ksol] = false
            if !test
                absvals = abs.(u[:, ksol, 2])
                masked = absvals .* mask[1, 1:mpert]
                kpert = argmax(masked)
                mask[1, kpert] = false
            end
            for jsol in mpert + 1:msol
                if mask[2, jsol]
                    if !test
                        fixfac[ksol, jsol] = -u[kpert, jsol, 2] / u[kpert, ksol, 2]
                    end
                    u[:, jsol, :] .= u[:, jsol, :] .+ u[:, ksol, :] .* fixfac[ksol, jsol]
                    if !test
                        u[kpert, jsol, 2] = 0
                    end
                end
            end
        end
    end

    # Save fixfac to file
    if bin_euler && !test
        write(euler_bin_unit, 2)
        write(euler_bin_unit, sing_flag, msol)
        write(euler_bin_unit, fixfac, index[1:msol])
    end

    # Write output values and close output file
    if diagnose
        println(init_out_unit, "output values:")
        for isol in 1:msol
            println(init_out_unit, "isol = $isol, m = $(mlow + isol - 1)")
            for ipert in 1:2
                println(init_out_unit, "re u($ipert), im u($ipert)")
            end
            for ipert in 1:mpert
                println(init_out_unit, ipert, u[ipert, isol, :])
            end
            for ipert in 1:2
                println(init_out_unit, "re u($ipert), im u($ipert)")
            end
        end
    end

    # ...existing code...
    return
end

function ode_test()
    # Returns true/false according to stopping/flag criteria

    # All globals are assumed available (see your input Fortran)
    # psifac, psimax, istep, nstep, istate, res_flag
    # ising, msing, singfac, singfac_max
    # ca_old, ca, msol, mpert, u, flag_count, ksing, err_unit
    # Also: sing_get_ca, etc

    # Static storage: emulated with `global`
    global singfac_old, powmax
    if !isdefined(@__MODULE__, :singfac_old)
        singfac_old = singfac
    end
    if !isdefined(@__MODULE__, :powmax)
        powmax = 0.0
    end

    flag = psifac == psimax || istep == nstep || istate < 0
    if !res_flag || flag || ising > msing || singfac > singfac_max
        return flag
    end

    ca_old .= ca           # Save previous surface ca
    sing_get_ca(ising, psifac, u, ca)    # Updates global ca
    dca = ca .- ca_old
    dsingfac = abs((singfac - singfac_old)/singfac)
    singfac_old = singfac

    power = zeros(Float64, msol)
    powmax_old = powmax
    for isol in 1:msol
        # ca[:,isol,:]::matrix of size (mpert,2)
        # Flatten the relevant slice for norm calculation
        ca_isol = ca[:, isol, :]           # mpert × 2
        dca_isol = dca[:, isol, :]
        norm = abs(sum(conj.(ca_isol) .* ca_isol))
        dnorm = abs(sum(conj.(ca_isol) .* dca_isol))
        power[isol] = dnorm / (norm * dsingfac)
    end
    powmax = maximum(power)

    global flag_count
    flag_count += 1
    if flag_count < 3
        return flag
    end
    flag = flag || ((singfac < singfac_max) && (powmax > powmax_old))

    if ising == ksing
        println(err_unit, float(log10(singfac)), "  ", float(log10(powmax)))
    end

    return flag
end

# Example stub for test fixup
function ode_test_fixup()
    # Implement test fixup logic here
    return
end

function ode_record_edge()
    debug = false
    # calc_number is persistent in Fortran; use a global or Ref in Julia if needed
    global calc_number
    total1 = 0.0 + 0.0im
    vacuum1 = 0.0 + 0.0im
    plasma1 = 0.0 + 0.0im

    #spline_eval(sq, psifac, 0)
    q = JPEC.SplinesMod.spline_eval(sq, psifac, 0)
    if size_edge > 0
        if sq.f[4] >= q_edge[i_edge] && psifac >= psiedge
            free_test(plasma1, vacuum1, total1, psifac)
            if debug
                println("$(calc_number) $(i_edge) $(psifac) $(sq.f[4]) $(q_edge[i_edge]) $(real(total1)) $(real(vacuum1)) $(real(plasma1))")
            end
            calc_number += 1
            dw_edge[i_edge] = total1
            q_edge[i_edge] = sq.f[4]
            psi_edge[i_edge] = psifac
            i_edge = min(i_edge + 1, size_edge)  # just to be extra safe
        end
    end

    return
end
