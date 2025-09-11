# Placeholder for dependencies
# using Fourfit
# using OdeOutput
# using Debug
# using Free
# using DifferentialEquations

const eps = 1e-10 # TODO: this shouldn't be here, right?

"""
OdeState

A mutable struct to hold the state of the ODE solver for DCON.
This struct contains all necessary fields to manage the ODE integration process,
including solution vectors, tolerances, and flags for the integration process.
"""
#TODO: variable description review by DCON expert (once finished), evaluate if all variables need to be part of the struct
# TODO: will msol ever be different than mpert? If so, need to create a separate tmp for mpert x msol
@kwdef mutable struct OdeState
    mpert::Int                  # poloidal mode number count
    msol::Int                   # number of solutions
    psi_save::Float64 = 0.0     # last saved psi value
    istep::Int= 0               # integration step count
    ix::Int= 0                  # index for psiout in spline #TODO: does this really need to be part of the struct?
    singfac::Float64 = 0.0      # separation from singular surface
    psifac::Float64 = 0.0       # normalized flux coordinate
    next::String = ""           # next action ("cross" or "finish")
    ising::Int = 0               # index of next singular surface
    flag_count::Int = 0         # count of flags raised
    new::Bool = true            # flag for new solution
    first::Bool = true          # flag for first step
    u::Array{ComplexF64, 3} = zeros(ComplexF64, mpert, mpert, 2)            # solution vectors
    u_save::Array{ComplexF64, 3} = zeros(ComplexF64, mpert, mpert, 2)       # saved solution vectors
    index::Vector{Int} = collect(1:mpert)                                   # indices for sorting solutions
    unorm::Vector{Float64} = zeros(Float64, 2*mpert)                        # norms of solution vectors
    unorm0::Vector{Float64} = zeros(Float64, 2*mpert)                       # initial norms of solution vectors
    fixfac::Array{ComplexF64,2} = zeros(ComplexF64, msol, msol)             # fixup factors for Gaussian reduction
    crit_save::Float64 = 0.0    # saved critical value for zero crossing detection
    nzero::Int = 0              # count of zero crossings detected
    m1::Int = 0                 # poloidal mode number for the first singular surface (?)
    q::Float64 = 0.0            # q value at the next singular surface
    psimax::Float64 = 0.0         # maximum psi value for the integrator (TODO: is this correct?)
    ua::Array{ComplexF64, 3} = zeros(ComplexF64, mpert, 2*mpert, 2) # auxiliary array for singular surface calculations
    # Temporary matrices for derivative calculations
    du_temp::Array{ComplexF64, 3} = zeros(ComplexF64, mpert, msol, 2)    
    amat::Vector{ComplexF64} = Vector{ComplexF64}(undef, mpert^2)
    bmat::Vector{ComplexF64} = Vector{ComplexF64}(undef, mpert^2)
    cmat::Vector{ComplexF64} = Vector{ComplexF64}(undef, mpert^2)
    fmat::Vector{ComplexF64} = Vector{ComplexF64}(undef, mpert^2)
    kmat::Vector{ComplexF64} = Vector{ComplexF64}(undef, mpert^2)
    gmat::Vector{ComplexF64} = Vector{ComplexF64}(undef, mpert^2)
    tmp::Matrix{ComplexF64} = Matrix{ComplexF64}(undef, mpert, mpert)
    Afact::Union{Cholesky{ComplexF64,Matrix{ComplexF64}}, Nothing} = nothing
    Ffact::Union{Cholesky{ComplexF64,Matrix{ComplexF64}}, Nothing} = nothing
    singfac_vec::Vector{Float64} = zeros(Float64, mpert)
end

OdeState(mpert::Int, msol::Int) = OdeState(; mpert, msol)

"""
    `ode_run(ctrl::DconControl, equil::Equilibrium.PlasmaEquilibrium, intr::DconInternal)`

Main driver for integrating plasma equilibrium and detecting singular surfaces.

Initializes state and iterates through flux surfaces, calling appropriate update
routines and recording output. Terminates when a target singularity is reached
or integration ends. Support for `res_flag` and `kin_flag` is not yet implemented.

### Returns
nzero: number of zero crossings of the critical determinant detected during the integration.
"""
# TODO: We pass the same structs into most functions; however, do to our use of multiple ! functions that modify the structs
# and the convention that the modified struct comes first, the struct order is never consistent. Is there a better way? 
function ode_run(ctrl::DconControl, equil::Equilibrium.PlasmaEquilibrium, intr::DconInternal, ffit::FourFitVars)
    # Initialization
    odet = OdeState(intr.mpert, intr.mpert)

    if ctrl.sing_start <= 0
        ode_axis_init!(ctrl, equil, intr, odet)
    elseif ctrl.sing_start <= intr.msing
        error("sing_start > 0 not implemented yet!")
        # ode_sing_init!(ctrl, equil, intr, odet)
    else
        error("Invalid value for sing_start: $(ctrl.sing_start) > msing = $(intr.msing)")
    end

    # Integration loop
    if ctrl.verbose # mimicing an output from ode_output_open
        println("   ψ=$(odet.psifac), q=$(SplinesMod.spline_eval(equil.sq, odet.psifac, 0)[4])")
    end
    # TODO: Eventually figure out a way to get rid of these - use callbacks in the integrator for the second one? maybe just do while ising != ksing and next != cross for the first?
    while true # inside plasma (break out when at the edge)
        odet.first = true
        while true # inside sing surface (break out when close to singular surface)
            # if odet.istep > 0
            #     ode_unorm!(odet, intr, ctrl, false)
            # end
            # Always record the first and last point in an inter-rational region
            # these are important for resonant quantities
            # recording the last point is critical for matching the nominal edge
            test = ode_test(odet, intr, ctrl)
            force_output = odet.first || test
            #ode_output_step(unorm, intr, ctrl, DconFileNames(), equil; force=force_output)
            ode_output_step(odet, intr, ctrl, equil; force=force_output)
            ode_record_edge!(intr, odet, ctrl, equil)
            test && break # break out of loop if ode_test returns true
            ode_step!(odet, equil, intr, ctrl, ffit)
            odet.first = false
        end

        # Re-initialize
        odet.ising == ctrl.ksing && break
        if odet.next == "cross"
            if ctrl.res_flag
                error("res_flag = true not implemented yet!")
                #ode_resist_cross()
            elseif ctrl.kin_flag
                error("kin_flag = true not implemented yet!")
                #ode_kin_cross()
            else
                ode_ideal_cross!(ctrl, equil, intr, odet)
            end
        else
            break
        end
        odet.flag_count = 0
    end
    return odet.nzero
end


"""
    ode_axis_init!(ctrl::DconControl, equil::Equilibrium.PlasmaEquilibrium,
                   intr::DconInternal, odet::OdeState; itmax = 50) -> Int

Initialize the ODE state for the case of sing_start = 0.

Finds the starting flux surface where `q = qlow` using Newton iteration, locates
nearby singular surfaces, and sets integration bounds. Sorts mode indices and
initializes solution arrays for perturbation equations.

Returns the index of the next relevant singular surface.
"""

function ode_axis_init!(ctrl::DconControl, equil::Equilibrium.PlasmaEquilibrium, intr::DconInternal, odet::OdeState; itmax = 50)

    # Shorthand to evaluate q/q1 inside newton iteration
    qval = psi -> SplinesMod.spline_eval(equil.sq, psi, 0)[4]
    q1val = psi -> SplinesMod.spline_eval(equil.sq, psi, 1)[2][4]

    # Preliminary computations
    odet.psifac = equil.sq.xs[1]

    # Use Newton iteration to find starting psi if qlow is above q0
    if ctrl.qlow > equil.sq.fs[1, 4]
        # Start check from the edge for robustness in reverse shear cores
        for jpsi = equil.sq.mx:-1:2  # avoid starting iteration on endpoints
            if equil.sq.fs[jpsi - 1, 4] < ctrl.qlow
                odet.psifac = equil.sq.xs[jpsi]
                break
            end
        end
        it = 0
        while it ≤ itmax
            it += 1
            dpsi = (ctrl.qlow - qval(odet.psifac)) / q1val(odet.psifac)
            odet.psifac += dpsi
            if abs(dpsi) < eps * abs(odet.psifac)
                break
            end
        end
    end

    # Find inner singular surface
    if false #(TODO: kin_flag)
        # for ising = 1:kmsing
        #     if kinsing[ising].psifac > psifac
        #         break
        #     end
        # end
    else
        odet.ising = 0
        for i in 1:intr.msing
            if intr.sing[i].psifac > odet.psifac
                odet.ising = max(0, i - 1)
                break
            end
        end
    end

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
        #     odet.q = kinsing[ising].q
        #     if mlow <= nn * odet.q && mhigh >= nn * odet.q
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
            odet.ising += 1
            if odet.ising > intr.msing || intr.psilim < intr.sing[min(odet.ising, intr.msing)].psifac
                break
            end
            odet.q = intr.sing[odet.ising].q
            if intr.mlow <= ctrl.nn * odet.q && intr.mhigh >= ctrl.nn * odet.q
                break
            end
        end
        if odet.ising > intr.msing || intr.psilim < intr.sing[min(odet.ising, intr.msing)].psifac || ctrl.singfac_min == 0
            odet.psimax = intr.psilim * (1 - eps)
            odet.next = "finish"
        else
            odet.psimax = intr.sing[odet.ising].psifac - ctrl.singfac_min / abs(ctrl.nn * intr.sing[odet.ising].q1)
            odet.next = "cross"
        end
    end

    # Allocate and sort solutions by increasing value of |m-ms1|
    m = intr.mlow - 1 .+ collect(1:intr.mpert)
    if ctrl.sort_type == "absm"
        key = abs.(m)
    elseif ctrl.sort_type == "sing"
        key = m
        if intr.msing > 0
            key .-= intr.sing[1].m
        end
        @. key = -abs(key)
    else
        error("Cannot recognize sort_type = $(ctrl.sort_type)")
    end
    odet.index = sortperm(key, rev = true) # in original Fortran: bubble(key, index, 1, mpert)

    # Initialize solutions
    for ipert = 1:intr.mpert
        odet.u[odet.index[ipert], ipert, 2] = 1
    end
    odet.msol = intr.mpert
    odet.u_save .= odet.u
    odet.psi_save = odet.psifac

    # Compute conditions at next singular surface
    if false #TODO: (kin_flag)
        # if kmsing > 0
        #     m1 = round(Int, nn * kinsing[ising].q)
        # else
        #     m1 = round(Int, nn * qlim) + sign(one, nn * sq.fs1[mpsi, 5])
        # end
    else
        # note: Julia's default round does Banker's rounding, to match NINT in fortran we need to specify RoundFromZero
        if intr.msing > 0
            odet.m1 = round(Int, ctrl.nn * intr.sing[odet.ising].q, RoundFromZero)
        else
            odet.m1 = round(Int, ctrl.nn * intr.qlim, RoundFromZero) + sign(ctrl.nn * equil.sq.fs1[end, 4])
        end
    end
    odet.singfac = abs(odet.m1 - ctrl.nn * equil.sq.fs[1, 4]) # Fortran: q=sq%fs(0,4)
end

# TODO: NOT IMPLEMENTED YET! (low priority, just make sure sing_start = 0 in dcon.toml)
function ode_sing_init()
    return
end
#     # Declare and initialize local variables
#     star = fill(' ', mpert)
#     # local variables
#     ua = Array{Complex{r8}}(undef, mpert, 2*mpert, 2)
#     dpsi = 0.0
#     new = true

#     ising = sing_start
#     dpsi = singfac_min/abs(nn*sing[ising].q1)*10
#     odet.psifac = sing[ising].psifac + dpsi
#     odet.q = sing[ising].q + dpsi*sing[ising].q1

#     # Allocate and initialize solution arrays
#     u      = zeros(Complex{r8}, mpert, mpert, 2)
#     du     = zeros(Complex{r8}, mpert, mpert, 2)
#     u_save = zeros(Complex{r8}, mpert, mpert, 2)
#     unorm0 = zeros(r8, 2*mpert)
#     unorm  = zeros(r8, 2*mpert)
#     index  = zeros(Int, 2*mpert)

#     if old_init
#         u .= 0.0
#         for ipert ∈ 1:mpert
#             u[ipert, ipert, 2] = 1.0
#         end
#     else
#         sing_get_ua(ising, odet.psifac, ua)
#         # Big slice: u = ua[:, mpert+1:2*mpert,:]
#         for i = 1:mpert, j = 1:mpert, k = 1:2
#             u[i, j, k] = ua[i, mpert+j, k]
#         end
#     end
#     u_save .= u
#     psi_save = psifac
#     msol = mpert
#     neq = 4 * mpert * msol

#     # Compute conditions at next singular surface
#     while true
#         ising += 1
#         if ising > msing || psilim < sing[ising ≤ msing ? ising : msing].psifac
#             break
#         end
#         odet.q = sing[ising].q
#         if mlow ≤ nn*q && mhigh ≥ nn*q
#             break
#         end
#     end
# This needs to be fixed up
#     if ising > msing || psilim < sing[ising ≤ msing ? ising : msing].psifac
#         m1 = round(Int, ctrl.nn*intr.qlim) + round(Int, sign(one, ctrl.nn*equil.sq.fs1[mpsi, 4]))
#         psimax = psilim * (1-eps)
#         next_ = "finish"
#     else
#         m1 = round(Int, nn*sing[ising].q)
#         psimax = sing[ising].psifac - singfac_min/abs(nn*sing[ising].q1)
#         next_ = "cross"
#     end

#     # Set up integrator parameters
#     istep = 0
#     istate = 1
#     iopt = 1
#     itask = 5
#     itol = 2
#     mf = 10
#     istate = 1

#     # Set up work arrays
#     liw = length(iwork)    # iwork must exist!
#     lrw = 22 + 64*mpert*msol
#     rwork = zeros(r8, lrw)
#     atol  = zeros(r8, mpert, msol, 2)
#     fixfac = zeros(r8, msol, msol)
#     iwork = zeros(Int, liw)  # or readjust as needed
#     rwork[1] = psimax
#     rwork[5] = dpsi
#     rwork[11] = rwork[5]

#     # Terminate, or in Julia just return (no need for RETURN)
#     return nothing  # or could return a struct with all these values, for a more Julian approach
# end

"""
    ode_ideal_cross!(odet::OdeState, equil::Equilibrium.PlasmaEquilibrium, 
                     intr::DconInternal, ctrl::DconControl)

Handle the crossing of an ideal MHD singular surface during ODE integration 
in a plasma equilibrium calculation.

# Arguments
- `odet::OdeState`: The current state of the ODE system (mutable, updated in-place).
- `equil::Equilibrium.PlasmaEquilibrium`: Object representing plasma equilibrium parameters.
- `intr::DconInternal`: Structure containing internal data, including singular surface information.
- `ctrl::DconControl`: Holds control flags and integration parameters.

# Description
This function updates the ODE solution as it crosses a (resonant) singular surface in the 
ideal MHD calculation. It normalizes and reinitializes the solution vector at the singularity, 
handles singular asymptotics, manages control flags, and updates relevant state variables 
in `odet` for continued integration. It also determines the location and parameters of 
the next singular surface, updates tolerance/stepping logic for the ODE solver, and may 
write auxiliary diagnostic output or coefficients as configured.

# Notes
- Handles resetting or reinitializing the solver after singular surface crossing.
- Assumes certain global state (such as `u`, `ua`, `index`, etc.) is properly managed.
- Intended for ideal MHD regions; kinetic surface handling is not included in this function.
"""
function ode_ideal_cross!(odet::OdeState, equil::Equilibrium.PlasmaEquilibrium, intr::DconInternal, ctrl::DconControl)

    # Fixup solution at singular surface
    if ctrl.verbose
        println("psi = $(intr.sing[odet.ising].psifac), q = $(intr.sing[odet.ising].q)")
    end
    ode_unorm(true)

    # Write asymptotic coefficients before reinit
    if ctrl.bin_euler
   #     sing_get_ca(ising, odet.psifac, u, ca)
   #     write(euler_bin_unit, 4)
   #     write(euler_bin_unit, intr.sing[ising].psifac, intr.sing[ising].q, intr.sing[ising].q1)
   #     write(euler_bin_unit, msol)
   #     write(euler_bin_unit, ca)
    end

    # Re-initialize
    psi_old = odet.psifac
    ipert0 = round(Int, ctrl.nn * intr.sing[odet.ising].q, RoundFromZero) - intr.mlow + 1
    dpsi = intr.sing[odet.ising].psifac - odet.psifac
    odet.psifac = intr.sing[odet.ising].psifac + dpsi
    # TODO: uncomment this once sing_get_ua! is implemented
    # sing_get_ua!(ising, odet.psifac, ua)
    if !ctrl.con_flag
        u[:, index[1], :] .= 0 #TODO: need to add index to odet
    end

    # Update solution vectors
    du1 = zeros(ComplexF64, intr.mpert, intr.mpert, 2)
    du2 = zeros(ComplexF64, intr.mpert, intr.mpert, 2)
    params = (ctrl, equil, intr, ffit)
    sing_der!(du1, odet.u, params, psi_old) #TODO: where does du1 come from?
    sing_der!(du2, odet.u, params, odet.psifac)
    odet.u .= odet.u .+ (du1 .+ du2) .* dpsi
    if !ctrl.con_flag
        error("con_flag = false not implemented yet!")
        odet.u[ipert0, :, :] .= 0
        odet.u[:, index[1], :] .= odet.ua[:, ipert0 + mpert, :]
    end

    # Write asymptotic coefficients after reinit
    if ctrl.bin_euler
  #      sing_get_ca(ising, odetpsifac, u, ca)
  #      write(euler_bin_unit, msol)
  #      write(euler_bin_unit, ca)
  #      write(euler_bin_unit, sing[ising].restype.e, sing[ising].restype.f,
  #            sing[ising].restype.h, sing[ising].restype.m,
  #            sing[ising].restype.g, sing[ising].restype.k,
  #            sing[ising].restype.eta, sing[ising].restype.rho,
  #            sing[ising].restype.taua, sing[ising].restype.taur)
    end

    # Find next ising
    while true
        odet.ising += 1
        if odet.ising > intr.msing || intr.psilim < intr.sing[min(odet.ising, intr.msing)].psifac
            break
        end
        odet.q = intr.sing[odet.ising].q
        if intr.mlow <= ctrl.nn * odet.q && intr.mhigh >= ctrl.nn * odet.q
            break
        end
    end
    
    # Compute conditions at next singular surface
    if odet.ising > intr.msing || intr.psilim < intr.sing[min(odet.ising, intr.msing)].psifac
        odet.psimax = intr.psilim * (1 - eps)
        odet.m1 = round(Int, ctrl.nn * intr.qlim, RoundFromZero) + sign(ctrl.nn * equil.sq.fs1[mpsi, 4])
        odet.next = "finish"
    else
        odet.psimax = intr.sing[odet.ising].psifac - ctrl.singfac_min / abs(ctrl.nn * intr.sing[odet.ising].q1)
        odet.m1 = round(Int,ctrl.nn * intr.sing[ising].q, RoundFromZero)
        # println(crit_out_unit, "ising = $ising, psifac = $(intr.sing[ising].psifac), q = $(intr.sing[ising].q), di = $(intr.sing[ising].di), alpha = $(intr.sing[ising].alpha)")

    end

    # Restart ode solver
    odet.istep += 1
    odet.new = true
    odet.u_save .= u
    odet.psi_save = odet.psifac

    # Write to files
    # println(crit_out_unit, "step info")
    # if crit_break
    #     write(crit_bin_unit)
    # end
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

"""
    ode_step!(odet::OdeState, equil::Equilibrium.PlasmaEquilibrium, 
             intr::DconInternal, ctrl::DconControl, ffit::FourFitVars)

Advance one adaptive integration step of the ODE system for a plasma equilibrium calculation.

# Arguments
- `odet::OdeState`: Structure storing the current ODE state and solution arrays.
- `equil::Equilibrium.PlasmaEquilibrium`: Plasma equilibrium parameters and profiles.
- `intr::DconInternal`: Internal data, including singular surfaces.
- `ctrl::DconControl`: Control and solver settings (e.g., tolerances, flags).
- `ffit::FourFitVars`: Additional fitting variables or Fourier fit coefficients.

# Description
This function computes and sets appropriate relative and absolute tolerances
for the ODE integration depending on proximity to singular surfaces, 
chooses the next integration endpoint, and advances the solution using 
an adaptive ODE solver. The state in `odet` is updated in-place with
the solution at the new point. Designed for use in stability or equilibrium 
scanning in plasma physics and MHD codes.

# Notes
- Tolerance logic distinguishes between "near-singular" and "regular" regions.
- ODE integration is performed using the DifferentialEquations.jl interface.
- Handles both node-based and maximal-psi stepping, depending on control flags.
"""
function ode_step!(odet::OdeState, equil::Equilibrium.PlasmaEquilibrium, intr::DconInternal, ctrl::DconControl, ffit::FourFitVars)

    # TODO: is there better relative/absolute tolerance handling for Julia? This might be specific to lsode
    
#     # Compute relative tolerances
#     singfac_local = typemax(Float64)
#     if false #(TODO: kin_flag)
#  #       if ising == 1 && intr.kmsing >= 1
#  #           singfac_local = abs(odet.psifac - kinsing[ising].psifac) /
#  #                           (kinsing[ising].psifac - psilow)
#  #       elseif ising <= intr.kmsing
#  #           singfac_local = min(
#  #               abs(odet.psifac - kinsing[ising].psifac),
#  #               abs(odet.psifac - kinsing[ising-1].psifac)
#  #           ) / abs(kinsing[ising].psifac - kinsing[ising-1].psifac)
#  #       end
#     else
#         if odet.ising <= intr.msing
#             singfac_local = abs(intr.sing[odet.ising].m - ctrl.nn * odet.q)
#         end
#         if odet.ising > 1
#             singfac_local = min(singfac_local, abs(intr.sing[odet.ising-1].m - ctrl.nn * odet.q))
#         end
#     end
#     rtol = tol = singfac_local < ctrl.crossover ? ctrl.tol_r : ctrl.tol_nr

#     # Compute absolute tolerances
#     atol = zeros(ComplexF64, intr.mpert, odet.msol, 2)
#     atol0 = 0.0
#     for ieq in 1:2, isol in 1:odet.msol
#         atol0 = maximum(abs.(odet.u[:, isol, ieq])) * tol
#         atol0 = atol0 == 0 ? typemax(Float64) : atol0
#         atol[:, isol, ieq] .= ComplexF64(atol0, atol0)
#     end

    # Choose psiout
    if ctrl.node_flag
        while odet.psifac < equil.sq.xs[ix]
            ix += 1
        end
        psiout = equil.sq.xs[ix]
        psiout = min(psiout, odet.psimax)
    else
        psiout = odet.psimax
    end

    function compute_tols(intr, odet, ctrl)
        singfac_local = typemax(Float64)
        # Relative tolerance
        if false  # kin_flag (not implemented)
            # Insert kin_flag branch if needed
        else
            if odet.ising <= intr.msing
                singfac_local = abs(intr.sing[odet.ising].m - ctrl.nn * odet.q)
            end
            if odet.ising > 1
                singfac_local = min(singfac_local,
                                    abs(intr.sing[odet.ising-1].m - ctrl.nn * odet.q))
            end
        end
        rtol = tol = singfac_local < ctrl.crossover ? ctrl.tol_r : ctrl.tol_nr
        # Absolute tolerances
        atol = similar(odet.u, Float64) # same shape as u
        for ieq in 1:size(odet.u,3), isol in 1:size(odet.u,2)
            atol0 = maximum(abs.(odet.u[:, isol, ieq])) * tol
            atol0 = atol0 == 0 ? typemax(Float64) : atol0
            atol[:, isol, ieq] .= atol0
        end
        return rtol, atol
    end    

    function update_tols!(integrator)
        integrator.p[4].u .= integrator.u
        rtol, atol = compute_tols(integrator.p[3], integrator.p[4], integrator.p[1])
        integrator.opts.reltol = rtol
        integrator.opts.abstol = atol
    end 
    tol_cb = DiscreteCallback((u,t,integrator)->true, (integrator)->update_tols!(integrator)) # update every step

    # Print out every tenth of the interval
    psi_start = odet.psifac
    Δψ = (odet.psimax - psi_start) / 10
    last_psi = Ref(0.0)

    print!(integrator) = begin
        u = integrator.u
        # params = integrator.p
        # params is a tuple: (ctrl, equil, intr, odet, ffit)
        # du = similar(u)
        # sing_der!(du, u, params, integrator.t)
        fmt(x) = @sprintf("%.5e", x)
        println("       ψ=", fmt(integrator.t), 
            ", max|u|=", fmt(maximum(abs.(u))))
            # ", max|du|=", fmt(maximum(abs.(du))))
        last_psi[] = integrator.t
    end
    print_cb = DiscreteCallback((u, t, integrator) -> (t >= last_psi[] + Δψ), print!)
    
    # Callback to normalize solution after each step if t > 0
    function unorm!(integrator)
        # Extract odet from parameters tuple
        params = integrator.p
        odet = params[4]
        intr = params[3]
        ctrl = params[1]
        odet.u .= integrator.u
        odet.psifac = integrator.t
        # TODO: I don't love how we double copy from integrator.u to odet.u then back, is there a better way?
        ode_unorm!(odet, intr, ctrl, false)
        # Need to update integrator.u with normalized odet.u!
        integrator.u .= odet.u
    end
    unorm_cb = DiscreteCallback((u, t, integrator) -> (t > odet.psifac), unorm!) # call on all but the first step

    # Existing callback functionality to match lsode/Fortran logic:
    # print_cb: print status every tenth of the interval
    # unorm_cb: unorming/fixup
    # tol_cb: update tolerances
    cb = CallbackSet(print_cb, unorm_cb, tol_cb)

    # Advance differential equation
    odet.istep += 1
    rtol, atol = compute_tols(intr, odet, ctrl)
    prob = ODEProblem(sing_der!, odet.u, (odet.psifac, psiout), (ctrl, equil, intr, odet, ffit))
    # TODO: check absolute tolerances, check how sensitive outputs are to tolerances
    # TODO: check if relative is updating correctly
    sol = solve(prob, Tsit5(), reltol=rtol, abstol=atol, callback=cb)

    # Update u and psifac with the solution at the end of the interval
    odet.u .= sol.u[end]
    odet.psifac = sol.t[end]

    println("Final psifac = $(odet.psifac), max u = $(maximum(abs.(odet.u)))")

    # dump_matrix("/Users/jakehalpern/Github/JPEC/notebooks/u_integration1.dat", odet.u[:, :, 1])
    # dump_matrix("/Users/jakehalpern/Github/JPEC/notebooks/u_integration2.dat", odet.u[:, :, 2])

    # error("stopping after ode step with psi = $(odet.psifac)")
end

"""
    ode_unorm!(odet::OdeState, intr::DconInternal,  sing_flag::Bool)

Computes norms of the solution vectors, normalizes them
relative to initial values, and applies Gaussian reduction via `ode_fixup!`
if the variation exceeds a threshold or if `sing_flag` is true.
Throws an error if any vector norm is zero.
"""
#function ode_unorm(intr::DconInternal, odet::OdeState, dout::DconOutput, fNames::DconFileNames, sing_flag::Bool)
# odet should be first, its modified. Can add output structs back later, removing for simplicity now
function ode_unorm!(odet::OdeState, intr::DconInternal, ctrl::DconControl, sing_flag::Bool)
    # Compute norms of first solution vectors, abort if any are zero
    odet.unorm[1:intr.mpert] .= [norm(odet.u[:, j, 1]) for j in 1:intr.mpert]
    odet.unorm[intr.mpert+1:odet.msol] .= [norm(odet.u[:, j, 2]) for j in intr.mpert+1:odet.msol]
    if minimum(odet.unorm[1:odet.msol]) == 0
        jmax = argmin(odet.unorm[1:odet.msol])
        error("One of the first solution vector norms unorm(1,$jmax) = 0")
    end

    # Normalize unorm and perform Gaussian reduction if required
    if odet.new
        odet.new = false
        odet.unorm0[1:odet.msol] .= odet.unorm[1:odet.msol]
    else
        odet.unorm[1:odet.msol] .= odet.unorm[1:odet.msol] ./ odet.unorm0[1:odet.msol]
        uratio = maximum(odet.unorm[1:odet.msol]) / minimum(odet.unorm[1:odet.msol])
        if uratio > ctrl.ucrit || sing_flag
            # ode_fixup(intr, odet, dout, fNames, sing_flag, false)
            ode_fixup!(odet, intr, sing_flag, false)
            odet.new = true
        end
    end
end

"""
    ode_fixup!(odet::OdeState, intr::DconInternal, sing_flag::Bool, test::Bool)

Applies Gaussian reduction to orthogonalize solution vectors in `odet.u`.

Used when the spread in norms exceeds a threshold or a singularity is detected.
Sorts solutions by norm, eliminates dependencies using forward elimination,
and optionally computes the transformation matrix `fixfac`.
Behavior can be suppressed for testing via the `test` flag.
"""

#function ode_fixup(intr::DconInternal, odet::OdeState, dout::DconOutput, fNames::DconFileNames, sing_flag::Bool, test::Bool)
#TODO: I think sing_flag is only needed for writing to binary in fortran DCON, might be able to remove
function ode_fixup!(odet::OdeState, intr::DconInternal, sing_flag::Bool, test::Bool)

    # TODO: seems like secondary is always false in fortran DCON (unless manually changed). is this needed?
    secondary = false

    # TODO: can test be removed if we aren't implementing ode_test_fixup?

    # TODO: add ctrl to use verbose flag for determining whether to print

    # Initial output
    println("Gaussian Reduction at psi = $(odet.psifac)")
    # println(fNames.crit_out_unit, "Gaussian Reduction at istep = $(odet.istep), psi = $(odet.psifac), q = $q")
    # if !sing_flag
    #     println(fNames.crit_out_unit, "Gaussian Reduction at istep = $(odet.istep), psi = $(odet.psifac), q = $q")
    # end
    odet.flag_count = 0

    # Initialize fixfac
    if !test
        odet.fixfac .= 0
        for isol in 1:odet.msol
            odet.fixfac[isol, isol] = 1
        end
        # Sort unorm
        odet.index[1:odet.msol] .= collect(1:odet.msol)
        odet.index[1:intr.mpert] .= sortperm_subrange(odet.unorm, 1:intr.mpert) # in original Fortran: bubble(unorm, index, 1, mpert)
        odet.index[intr.mpert+1:odet.msol] .= sortperm_subrange(odet.unorm, intr.mpert+1:odet.msol) # in original Fortran: bubble(unorm, index, mpert + 1, msol)
    end

    # Triangularize primary solutions
    mask = trues(2, odet.msol)
    for isol in 1:intr.mpert
        ksol = odet.index[isol]
        mask[2, ksol] = false
        if !test
            # Find max location
            absvals = abs.(odet.u[:, ksol, 1])
            masked = absvals .* mask[1, 1:intr.mpert]
            kpert = argmax(masked)
            #kpert = findmax(abs.(odet.u[:, ksol, 1]) .* (mask[1, 1:intr.mpert] .|> Int))[2] #TODO: ChatGPT suggests this might be better than the above three lines
            mask[1, kpert] = false
        end
        for jsol in 1:odet.msol
            if mask[2, jsol]
                if !test
                    odet.fixfac[ksol, jsol] = -odet.u[kpert, jsol, 1] / odet.u[kpert, ksol, 1]
                end
                odet.u[:, jsol, :] .= odet.u[:, jsol, :] .+ odet.u[:, ksol, :] .* odet.fixfac[ksol, jsol]
                if !test
                    odet.u[kpert, jsol, 1] = 0
                end
            end
        end
    end

    # Triangularize secondary solutions
    if odet.msol > intr.mpert && secondary
        mask = trues(2, odet.msol)
        for isol in intr.mpert + 1:odet.msol
            ksol = odet.index[isol]
            mask[2, ksol] = false
            if !test
                absvals = abs.(odet.u[:, ksol, 2])
                masked = absvals .* mask[1, 1:intr.mpert]
                kpert = argmax(masked)
                mask[1, kpert] = false
            end
            for jsol in intr.mpert + 1:odet.msol
                if mask[2, jsol]
                    if !test
                        odet.fixfac[ksol, jsol] = -odet.u[kpert, jsol, 2] / odet.u[kpert, ksol, 2]
                    end
                    odet.u[:, jsol, :] .= odet.u[:, jsol, :] .+ odet.u[:, ksol, :] .* odet.fixfac[ksol, jsol]
                    if !test
                        odet.u[kpert, jsol, 2] = 0
                    end
                end
            end
        end
    end

    # # Save fixfac to file
    # if dout.bin_euler && !test
    #     write(dout.euler_bin_unit, 2)
    #     write(fNames.euler_bin_unit, sing_flag, msol)
    #     write(fNames.euler_bin_unit, fixfac, index[1:msol])
    # end
end

"""
    ode_test(odet::OdeState, intr::DconInternal, ctrl::DconControl)::Bool

Checks whether integration should terminate based on step count,
singularity index, or diagnostic thresholds.

Returns `true` if integration is complete or if any stopping criteria
are met.
If `res_flag` is enabled, also computes changes in mode structure
and evaluates a power growth metric to trigger early stopping (TODO: not sure if this is true).
"""
function ode_test(odet::OdeState, intr::DconInternal, ctrl::DconControl)::Bool

    # check if we are at end of integration
    flag = odet.psifac == odet.psimax || odet.istep == ctrl.nstep

    # if not running with res_flag = true this function will exit and return flag here)
    if !ctrl.res_flag || flag || odet.ising > intr.msing || odet.singfac > ctrl.singfac_max
        return flag
    end

    @error "ode_test with res_flag = true not implemented yet!"
    # TODO: none of this has been checked, wait until we implement res_flag = true
    # ca_old .= ca           # Save previous surface ca
    # sing_get_ca(ising, psifac, u, ca)    # Updates global ca
    # dca = ca .- ca_old
    # dsingfac = abs((singfac - singfac_old)/singfac)
    # singfac_old = singfac

    # power = zeros(Float64, msol)
    # powmax_old = powmax
    # for isol in 1:msol
    #     # ca[:,isol,:]::matrix of size (mpert,2)
    #     # Flatten the relevant slice for norm calculation
    #     ca_isol = ca[:, isol, :]           # mpert × 2
    #     dca_isol = dca[:, isol, :]
    #     norm = abs(sum(conj.(ca_isol) .* ca_isol))
    #     dnorm = abs(sum(conj.(ca_isol) .* dca_isol))
    #     power[isol] = dnorm / (norm * dsingfac)
    # end
    # powmax = maximum(power)

    # odet.flag_count += 1
    # if odet.flag_count < 3
    #     return flag
    # end
    # flag = flag || ((singfac < singfac_max) && (powmax > powmax_old))

    # if ising == ksing
    #     println(err_unit, float(log10(singfac)), "  ", float(log10(powmax)))
    # end

    # return flag
end

"""
    ode_record_edge!(intr::DconInternal, odet::OdeState, ctrl::DconControl, equil::Equilibrium.PlasmaEquilibrium)

Records edge quantities at the current integration point if edge conditions are met.

Checks whether `psifac` and `q` exceed edge thresholds, and if so,
stores relevant values in `intr`. Currently, vacuum calculations are
not implemented; this function raises an error if invoked in a vacuum region.
"""
# TODO: this function requires integration the vacuum code for free_test
# for now, going to make it look ok anmd error out if we get to free_test
function ode_record_edge!(intr::DconInternal, odet::OdeState, ctrl::DconControl, equil::Equilibrium.PlasmaEquilibrium)

    total1 = ComplexF64(0.0, 0.0)
    vacuum1 = ComplexF64(0.0, 0.0)
    plasma1 = ComplexF64(0.0, 0.0)

    q_psifac = SplinesMod.spline_eval(equil.sq, odet.psifac, 0)
    if intr.size_edge > 0
        #TODO: WTH? fortran has both a psiedge and psi_edge and they appear to be different.
        # someone smarter than me please double check this
        if q_psifac >= intr.q_edge[intr.i_edge] && odet.psifac >= ctrl.psiedge
            @error "Vacuum code not yet integrated. This run has size_edge = $(intr.size_edge) > 0 and integrator passed psiedge/q_edge"
            #free_test(plasma1, vacuum1, total1, odet.psifac)
            intr.dw_edge[intr.i_edge] = total1
            intr.q_edge[intr.i_edge] = q_psifac
            intr.psi_edge[intr.i_edge] = odet.psifac
            intr.i_edge = min(intr.i_edge + 1, intr.size_edge)  # just to be extra safe
        end
    end
end