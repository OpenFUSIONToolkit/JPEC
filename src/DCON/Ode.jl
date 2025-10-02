# This is used for various small tolerances throughout this file
const eps = 1e-10

"""
OdeState

A mutable struct to hold the state of the ODE solver for DCON.
This struct contains all necessary fields to manage the ODE integration process,
including solution vectors, tolerances, and flags for the integration process.
"""
#TODO: variable description review by DCON expert (once finished), evaluate if all variables need to be part of the struct
# TODO: remove msol variable, I think it just makes things more confusing. Either use mpert or 2*mpert in loops
@kwdef mutable struct OdeState
    mpert::Int                  # poloidal mode number count
    msol::Int                   # number of solutions
    psi_save::Float64 = 0.0     # last saved psi value
    singfac::Float64 = 0.0      # separation from singular surface in terms of m - nq
    psifac::Float64 = 0.0       # normalized flux coordinate
    q::Float64 = 0.0            # q value at psifac
    next::String = ""           # next action ("cross" or "finish")
    ising::Int = 0               # index of next singular surface
    flag_count::Int = 0         # count of flags raised
    new::Bool = true            # flag for new solution
    u::Array{ComplexF64, 3} = zeros(ComplexF64, mpert, mpert, 2)            # solution vectors
    u_save::Array{ComplexF64, 3} = zeros(ComplexF64, mpert, mpert, 2)       # saved solution vectors
    ud::Array{ComplexF64, 3} = zeros(ComplexF64, mpert, mpert, 2)           # derivative of solution vectors used in GPEC
    index::Vector{Int} = collect(1:mpert)                                   # indices for sorting solutions
    unorm::Vector{Float64} = zeros(Float64, 2*mpert)                        # norms of solution vectors
    unorm0::Vector{Float64} = zeros(Float64, 2*mpert)                       # initial norms of solution vectors
    fixfac::Array{ComplexF64,2} = zeros(ComplexF64, mpert, mpert)             # fixup factors for Gaussian reduction
    crit_save::Float64 = 0.0    # saved critical value for zero crossing detection
    nzero::Int = 0              # count of zero crossings detected
    m1::Int = 0                 # poloidal mode number for the first singular surface (?)
    psimax::Float64 = 0.0         # maximum psi value for the integrator
    ua::Array{ComplexF64, 3} = zeros(ComplexF64, mpert, 2*mpert, 2) # asympototic solutions at singular surface
    ca::Array{ComplexF64, 3} = zeros(ComplexF64, mpert, 2*mpert, 2) # asymptotic coefficients at singular surface
    # Temporary matrices for sing_der calculations
    amat::Vector{ComplexF64} = Vector{ComplexF64}(undef, mpert^2)
    bmat::Vector{ComplexF64} = Vector{ComplexF64}(undef, mpert^2)
    cmat::Vector{ComplexF64} = Vector{ComplexF64}(undef, mpert^2)
    fmat::Vector{ComplexF64} = Vector{ComplexF64}(undef, mpert^2)
    kmat::Vector{ComplexF64} = Vector{ComplexF64}(undef, mpert^2)
    gmat::Vector{ComplexF64} = Vector{ComplexF64}(undef, mpert^2)
    Afact::Union{Cholesky{ComplexF64,Matrix{ComplexF64}}, Nothing} = nothing
    singfac_vec::Vector{Float64} = zeros(Float64, mpert)
end

OdeState(mpert::Int, msol::Int) = OdeState(; mpert, msol)

"""
    `ode_run(ctrl::DconControl, equil::Equilibrium.PlasmaEquilibrium, intr::DconInternal)`

Main driver for integrating plasma equilibrium and detecting singular surfaces.

Initializes state and iterates through flux surfaces, calling appropriate update
routines and recording output. Terminates when a target singularity is reached
or integration ends. Support for `res_flag` and `kin_flag` is not yet implemented.

# TODO: add a detailed description of how the integration loop varies from the Fortran version

### Returns
nzero: number of zero crossings of the critical determinant detected during the integration.
"""
# TODO: We pass the same structs into most functions; however, do to our use of multiple ! functions that modify the structs
# and the convention that the modified struct comes first, the struct order is never consistent. Is there a better way? 
function ode_run(ctrl::DconControl, equil::Equilibrium.PlasmaEquilibrium, intr::DconInternal, ffit::FourFitVars, fnames::DconFileNames)
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

    if ctrl.verbose # mimicing an output from ode_output_open
        println("   ψ=$(odet.psifac), q=$(Spl.spline_eval(equil.sq, odet.psifac, 0)[4])")
    end

    # Write header data to files
    ode_output_init(ctrl, equil, fnames, intr, odet)

    # TODO: does this exactly reproduce all functionality of the Fortran in all cases?
    # Might be good for ideal case, but maybe things I haven't considered break this
    # Always integrate once, even if no rational surfaces are crossed
    ode_step!(odet, equil, intr, ctrl, ffit, fnames)

    # If at a rational surface, do the appropriate crossing routine, then integrate again
    while odet.ising != ctrl.ksing && odet.next == "cross"
        if ctrl.res_flag
            error("res_flag = true not implemented yet!")
        elseif ctrl.kin_flag
            error("kin_flag = true not implemented yet!")
        else
            ode_ideal_cross!(odet, equil, intr, ctrl, ffit, fnames)
        end

        odet.flag_count = 0
        ode_step!(odet, equil, intr, ctrl, ffit, fnames)
    end

    # This is a rough sketch of the original Fortran logic for comparison, TODO: get rid of this eventually
    # while true # inside plasma (break out when at the edge)
        # Rough sketch of Fortran logic for comparison
        # while true # inside sing surface (break out when close to singular surface)
            # if odet.istep > 0
            #     ode_unorm!(odet, intr, ctrl, false)
            # end
            # Always record the first and last point in an inter-rational region
            # these are important for resonant quantities
            # recording the last point is critical for matching the nominal edge
            # test = ode_test(odet, intr, ctrl)
            # force_output = odet.first || test
            # ode_output_step(odet, intr, ctrl, equil; force=force_output)
            # ode_record_edge!(intr, odet, ctrl, equil)
            # test && break # break out of loop if ode_test returns true
            # ode_step!(odet, equil, intr, ctrl, ffit)
            # odet.first = false
        # end

    #     # Re-initialize
    #     odet.ising == ctrl.ksing && break
    #     if odet.next == "cross"
    #         if ctrl.res_flag
    #             error("res_flag = true not implemented yet!")
    #             #ode_resist_cross()
    #         elseif ctrl.kin_flag
    #             error("kin_flag = true not implemented yet!")
    #             #ode_kin_cross()
    #         else
    #             ode_ideal_cross!(odet, equil, intr, ctrl, ffit)
    #         end
    #     else
    #         break
    #     end
    #     odet.flag_count = 0
    # end
    return odet
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
    qval = psi -> Spl.spline_eval(equil.sq, psi, 0)[4]
    q1val = psi -> Spl.spline_eval(equil.sq, psi, 1)[2][4]

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
            if intr.mlow <= ctrl.nn * intr.sing[odet.ising].q && intr.mhigh >= ctrl.nn * intr.sing[odet.ising].q
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
function ode_ideal_cross!(odet::OdeState, equil::Equilibrium.PlasmaEquilibrium, intr::DconInternal, ctrl::DconControl, ffit::FourFitVars, fnames::DconFileNames)

    # Fixup solution at singular surface
    if ctrl.verbose
        println("psi = $(intr.sing[odet.ising].psifac), q = $(intr.sing[odet.ising].q)")
    end
    ode_unorm!(odet, intr, ctrl, fnames, true)

    # Write asymptotic coefficients before reinit
    # if ctrl.bin_euler
    sing_get_ca!(odet, intr, ctrl)
    write_to!(fnames, :euler_bin_unit, 4)
    write_to!(fnames, :euler_bin_unit, intr.sing[odet.ising].psifac, intr.sing[odet.ising].q, intr.sing[odet.ising].q1)
    write_to!(fnames, :euler_bin_unit, odet.msol)
    write_to!(fnames, :euler_bin_unit, odet.ca)
    # end

    # Re-initialize on opposite side of rational surface
    psi_old = odet.psifac
    ipert0 = round(Int, ctrl.nn * intr.sing[odet.ising].q, RoundFromZero) - intr.mlow + 1
    dpsi = intr.sing[odet.ising].psifac - odet.psifac
    odet.psifac = intr.sing[odet.ising].psifac + dpsi
    sing_get_ua!(odet, intr, ctrl)
    if !ctrl.con_flag
        odet.u[:, odet.index[1], :] .= 0
    end

    # Update solution vectors
    du1 = zeros(ComplexF64, intr.mpert, intr.mpert, 2)
    du2 = zeros(ComplexF64, intr.mpert, intr.mpert, 2)
    params = (ctrl, equil, intr, odet, ffit, fnames)
    sing_der!(du1, odet.u, params, psi_old)
    sing_der!(du2, odet.u, params, odet.psifac)
    odet.u .+= (du1 .+ du2) .* dpsi
    if !ctrl.con_flag
        odet.u[ipert0, :, :] .= 0
        odet.u[:, odet.index[1], :] .= odet.ua[:, ipert0 + intr.mpert, :]
    end

    # Write asymptotic coefficients after reinit
    # if ctrl.bin_euler
    sing_get_ca!(odet, intr, ctrl)
    write_to!(fnames, :euler_bin_unit, odet.msol)
    write_to!(fnames, :euler_bin_unit, odet.ca)
    write_to!(fnames, :euler_bin_unit, intr.sing[odet.ising].restype.e, intr.sing[odet.ising].restype.f,
            intr.sing[odet.ising].restype.h, intr.sing[odet.ising].restype.m,
            intr.sing[odet.ising].restype.g, intr.sing[odet.ising].restype.k,
            intr.sing[odet.ising].restype.eta, intr.sing[odet.ising].restype.rho,
            intr.sing[odet.ising].restype.taua, intr.sing[odet.ising].restype.taur)
    # end

    # Find next ising
    while true
        odet.ising += 1
        if odet.ising > intr.msing || intr.psilim < intr.sing[min(odet.ising, intr.msing)].psifac
            break
        end
        if intr.mlow <= ctrl.nn * intr.sing[odet.ising].q && intr.mhigh >= ctrl.nn * intr.sing[odet.ising].q
            break
        end
    end
    
    # Compute conditions at next singular surface
    if odet.ising > intr.msing || intr.psilim < intr.sing[min(odet.ising, intr.msing)].psifac
        odet.psimax = intr.psilim * (1 - eps)
        odet.m1 = round(Int, ctrl.nn * intr.qlim, RoundFromZero) + sign(ctrl.nn * equil.sq.fs1[end, 4])
        odet.next = "finish"
    else
        odet.psimax = intr.sing[odet.ising].psifac - ctrl.singfac_min / abs(ctrl.nn * intr.sing[odet.ising].q1)
        odet.m1 = round(Int,ctrl.nn * intr.sing[odet.ising].q, RoundFromZero)
        write_to!(fnames, :crit_out_unit, @sprintf("   ising   psi         q          di      re alpha   im alpha\n\n%6d%11.3e%11.3e%11.3e%11.3e%11.3e\n",
            odet.ising, intr.sing[odet.ising].psifac, intr.sing[odet.ising].q, intr.sing[odet.ising].di, real(intr.sing[odet.ising].alpha), imag(intr.sing[odet.ising].alpha)))
    end

    # Restart ode solver
    odet.new = true
    odet.u_save .= odet.u
    odet.psi_save = odet.psifac

    # Write next header before continuing integration
    write_to!(fnames, :crit_out_unit, "    psifac      dpsi        q       singfac     eval1")
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

Integrate the Euler-Lagrange equations to the next rational surface or edge.

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
the solution at the new point.

# Notes
- Tolerance logic distinguishes between "near-singular" and "regular" regions.
- ODE integration is performed using the DifferentialEquations.jl interface.
- Handles both node-based and maximal-psi stepping, depending on control flags.
"""
function ode_step!(odet::OdeState, equil::Equilibrium.PlasmaEquilibrium, intr::DconInternal, ctrl::DconControl, ffit::FourFitVars, fnames::DconFileNames)

    # print!(integrator) = begin
    #     u = integrator.u
    #     # params = integrator.p
    #     # params is a tuple: (ctrl, equil, intr, odet, ffit)
    #     # du = similar(u)
    #     # sing_der!(du, u, params, integrator.t)
    #     fmt(x) = @sprintf("%.5e", x)
    #     println("       ψ=", fmt(integrator.t), 
    #         ", max|u|=", fmt(maximum(abs.(u))))
    #         # ", max|du|=", fmt(maximum(abs.(du))))
    #     last_psi[] = integrator.t
    # end

    # This callback handles the logic that was previously in a DO loop within ode_run
    # and called every step by running LSODE in one step mode
    function integrator_callback!(integrator)
        ctrl, equil, intr, odet, ffit, fnames = integrator.p
        odet.u .= integrator.u
        odet.psifac = integrator.t
        
        # Adaptively update integration tolerances
        rtol, atol = compute_tols(intr, odet, ctrl)
        integrator.opts.reltol = rtol
        # integrator.opts.abstol = atol

        # TODO: Print status every Δψ, maybe get rid of this at some point? Good for debugging rn
        # if integrator.t >= last_psi[] + Δψ
        #     print!(integrator)
        # end

        # Note: no need for istep > 0 condition, since this is called after the first integration step always        
        ode_unorm!(odet, intr, ctrl, fnames, false)
        # Need to update integrator.u with normalized odet.u!
        integrator.u .= odet.u

        # TODO: ode_test effectively controls whether we print at the end of integration or if a number of integration steps is reached
        # this then controls the output of ode_output_step and then has a break condition below
        # Given that the Julia DiffEq solver will handle step sizes and errors, and we could just add a call to ode_output_step after the Integration
        # to force output, is this logic needed here? Might depend on how we do outputs in Julia
        force_output = ode_test(odet, intr, ctrl)
        ode_output_step!(odet, intr, ctrl, equil, ffit, fnames; force=force_output)
        ode_record_edge!(intr, odet, ctrl, equil)
    end
    cb = DiscreteCallback((u,t,integrator)->true, integrator_callback!) # perform callback at every step

    # Choose integration bounds
    if ctrl.node_flag # TODO: I can't find any mention of what this flag is for, can we get rid of it?
        while odet.psifac >= equil.sq.xs[ix]
            ix += 1
        end
        psiout = equil.sq.xs[ix]
        psiout = min(psiout, odet.psimax)
    else
        psiout = odet.psimax
    end

    # Print out every tenth of the inter-rational interval
    # Δψ = (odet.psimax - odet.psifac) / 10
    # last_psi = Ref(0.0)

    # Always record the first and last point in an inter-rational region these are important
    # for resonant quantities. Recording the last point is critical for matching the nominal edge
    # We don't use the first variable, and instead call the functions directly here before integration
    rtol, atol = compute_tols(intr, odet, ctrl) # initial tolerances
    ode_output_step!(odet, intr, ctrl, equil, ffit, fnames; force=true) # always output initial condition

    # Advance differential equation to next singular surface or edge
    prob = ODEProblem(sing_der!, odet.u, (odet.psifac, psiout), (ctrl, equil, intr, odet, ffit, fnames))
    sol = solve(prob, Tsit5(), reltol=rtol, abstol=atol, callback=cb)
    # TODO: check absolute tolerances, check how sensitive outputs are to tolerances

    # Update u and psifac with the solution at the end of the interval
    odet.u .= sol.u[end]
    odet.psifac = sol.t[end]
    ode_output_step!(odet, intr, ctrl, equil, ffit, fnames; force=true) # always output final condition

    println("   ψ = $(odet.psifac), max u = $(maximum(abs.(odet.u)))")
end

# Compute tolerances based on distance to singular surface during integration
function compute_tols(intr, odet, ctrl)
    singfac_local = typemax(Float64)
    # Relative tolerance
    if false  # kin_flag (not implemented)
        # Insert kin_flag branch if needed
    else
        # Note: odet.q is updated within the derivative calculation
        if odet.ising <= intr.msing
            singfac_local = abs(intr.sing[odet.ising].m - ctrl.nn * odet.q)
        end
        # If in between singular surfaces, check distance to both
        if odet.ising > 1
            singfac_local = min(singfac_local, abs(intr.sing[odet.ising-1].m - ctrl.nn * odet.q))
        end
    end
    rtol = tol = singfac_local < ctrl.crossover ? ctrl.tol_r : ctrl.tol_nr
    # Absolute tolerances
    atol = similar(odet.u, Float64)
    for ieq in 1:size(odet.u,3), isol in 1:size(odet.u,2)
        atol0 = maximum(abs.(odet.u[:, isol, ieq])) * tol
        if (atol0 == 0) atol0 = typemax(Float64) end
        atol[:, isol, ieq] .= atol0
    end
    return rtol, atol
end

"""
    ode_unorm!(odet::OdeState, intr::DconInternal,  sing_flag::Bool)

Computes norms of the solution vectors, normalizes them
relative to initial values, and applies Gaussian reduction via `ode_fixup!`
if the variation exceeds a threshold or if `sing_flag` is true.
Throws an error if any vector norm is zero.
"""
function ode_unorm!(odet::OdeState, intr::DconInternal, ctrl::DconControl, fnames::DconFileNames, sing_flag::Bool)
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
            ode_fixup!(odet, intr, fnames, sing_flag, false)
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
#TODO: I think sing_flag is only needed for writing to binary in fortran DCON, might be able to remove
function ode_fixup!(odet::OdeState, intr::DconInternal, fnames::DconFileNames, sing_flag::Bool, test::Bool)

    # TODO: seems like secondary is always false in fortran DCON (unless manually changed). is this needed?
    secondary = false

    # TODO: can test be removed if we aren't implementing ode_test_fixup?

    # Write to output
    write_to!(fnames, :crit_out_unit, "\n   psifac      dpsi        q       singfac     eval1")
    write_to!(fnames, :crit_out_unit, @sprintf("\nGaussian Reduction at psi = %10.3e, q = %6.3f\n", odet.psifac, odet.q))
    if !sing_flag
        write_to!(fnames, :crit_out_unit, "   psifac      dpsi        q       singfac     eval1\n")
    end
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
            kpert = argmax(abs.(odet.u[:, ksol, 1]) .* mask[1, 1:intr.mpert])
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
    write_to!(fnames, :euler_bin_unit, 2)
    write_to!(fnames, :euler_bin_unit, sing_flag, odet.msol)
    write_to!(fnames, :euler_bin_unit, odet.fixfac, odet.index[1:odet.msol])
    # end
end

"""
    ode_test(odet::OdeState, intr::DconInternal, ctrl::DconControl)::Bool

Returns `true` if integration is complete or if any stopping criteria are met.

TODO: implement res_flag = true logic
"""
function ode_test(odet::OdeState, intr::DconInternal, ctrl::DconControl)::Bool

    # check if we are at end of integration
    flag = odet.psifac == odet.psimax

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
# for now, going to make it look ok and error out if we get to free_test
function ode_record_edge!(intr::DconInternal, odet::OdeState, ctrl::DconControl, equil::Equilibrium.PlasmaEquilibrium)

    total1 = ComplexF64(0.0, 0.0)
    vacuum1 = ComplexF64(0.0, 0.0)
    plasma1 = ComplexF64(0.0, 0.0)

    q_psifac = Spl.spline_eval(equil.sq, odet.psifac, 0)
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