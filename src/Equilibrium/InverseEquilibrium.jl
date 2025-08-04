"""
Converts inverse equilibrium to straight-fieldline coordinates. Based on inverse.f
Of the entries we need to return: PlasmaEquilibrium(equil_params, sq_out, rzphi_out, eqfun_out, ro, zo, psio),
we only need to generate: sq_out, rzphi_out, eqfun_out. This is because we pass in InverseRunInput(equil_in, 
sq_in, rz_in, ro, zo, psio).

"""

"""
    inverse_extrap(xx::Matrix{Float64}, ff::Matrix{Float64}, x::Float64) -> Vector{Float64}

Performs component-wise Lagrange extrapolation for a vector-valued function.

## Arguments:
- `xx`: A (m × n) matrix where each row contains the x-values for each component.
- `ff`: A (m × n) matrix where each row contains function values at the corresponding `xx`.
- `x`: A scalar Float64 value at which to extrapolate.

## Returns:
- A vector of length n representing the extrapolated function values at `x`.
"""
function inverse_extrap(xx::Matrix{Float64}, ff::Matrix{Float64}, x::Float64)::Vector{Float64}
    m, n = size(ff)             # m = number of data points, n = number of components
    f = zeros(Float64, n)       # Output vector
    
    for i in 1:m
        term = copy(ff[i, :])   # Start with f_i (vector)
        for j in 1:m
            if j == i
                continue
            end
            term .= term .* ((x .- xx[j, :]) ./ (xx[i, :] .- xx[j, :]))
        end
        f .+= term              # Accumulate to output
    end
    
    return f
end



function equilibrium_solver(input::InverseRunInput)
    println("--- Starting Inverse Equilibrium Processing ---")
    # Extract input parameters

    config = input.config
    rz_in = input.rz_in
    sq_in = input.sq_in
    ro = input.ro
    zo = input.zo
    psio = input.psio

    grid_type = config.control.grid_type
    mpsi = config.control.mpsi
    mtheta = config.control.mtheta
    psilow = config.control.psilow
    psihigh = config.control.psihigh
    newq0 = config.control.newq0

    me = 3
    interp = false
    diagnose_rz_in = false
    diagnose_rzphi = false

    # c-----------------------------------------------------------------------
    # c     allocate and define local arrays.
    # c-----------------------------------------------------------------------
    # sq_in._fs[:, 3] .= sqrt.(sq_in._xs)
    rz_in._xs = sq_in._xs
    rz_in._ys = collect(0:rz_in.my) ./ rz_in.my

    mx = rz_in.mx
    my = rz_in.my

    x = rz_in.fs[:, :, 1] .- ro
    y = rz_in.fs[:, :, 2] .- zo
    r2 = x.^2 .+ y.^2

    twopi = 2 * pi

    deta = zeros(Float64, mx+1, my+1)
    for ipsi in 0:mx, itheta in 0:my
        if r2[ipsi+1, itheta+1] == 0.0
            deta[ipsi+1, itheta+1] = 0.0
        else
            deta[ipsi+1, itheta+1] = atan(y[ipsi+1, itheta+1], x[ipsi+1, itheta+1]) / twopi
        end
    end




end