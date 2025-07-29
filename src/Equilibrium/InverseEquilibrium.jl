# src/Equilibrium/InverseEquilibrium.jl
module InverseEquilibrium

using ..Types: InverseRunInput, PlasmaEquilibrium
using ..Spl: bicube_setup, bicube_eval, spline_setup, spline_eval
using Printf

export inverse_run

"""
    inverse_run(input::InverseRunInput)

Processes the inverse equilibrium to transform it to straight-fieldline coordinates.

## Arguments:
- `input`: An `InverseRunInput` object containing all setup parameters and splines.
## Returns:
- A `PlasmaEquilibrium` object containing the transformed equilibrium.
"""
function inverse_run(input::InverseRunInput)
    println("Starting inverse equilibrium process...")

    # Unpack input fields
    equil_input, sq_in, rz_in, ro, zo, psio = input.equil_input, input.sq_in, input.rz_in, input.ro, input.zo, input.psio

    # Allocate arrays
    mx, my = rz_in.mx, rz_in.my
    x, y, r2, deta = rz_in.fs[:,:,1] .- ro, rz_in.fs[:,:,2] .- zo, Float64[], Float64[]

    # Transform input coordinates from cartesian to polar
    function transform_coordinates(mx, my)
        r2 = x .^ 2 .+ y .^ 2
        deta = atan.(y, x) ./ (2π)
        for ipsi in 1:mx, itheta in 2:my
            if deta[ipsi, itheta] - deta[ipsi, itheta - 1] > 0.5
                deta[ipsi, itheta] -= 1
            elseif deta[ipsi, itheta] - deta[ipsi, itheta - 1] < -0.5
                deta[ipsi, itheta] += 1
            end
        end
        return r2, deta
    end

    r2, deta = transform_coordinates(mx, my)

    # Fit rz_in to bicubic splines
    rz_in.fs[:,:,1] = r2
    rz_in.fs[:,:,2] = deta
    bicube_fit(rz_in, "extrap", "periodic")

    # Prepare new spline type for surface quantities
    mpsi = sq_in.mx
    sq = spline_setup(rz_in.xs, rz_in.fs, bctype=4)

    # Set up radial grid
    sq_xs = collect(0.0:1.0/mpsi:1.0)
    sq_xs .= eq_grid_choice(sq_xs, equil_input.grid_type)

    # Prepare new bicube type for coordinates
    rzphi_fs_nodes = Array{Float64}(undef, mpsi + 1, my + 1, 4)
    rzphi = bicube_setup(sq_xs, collect(0.0:1.0/my:1.0), rzphi_fs_nodes, bctypex=4, bctypey=2)

    # Process new surfaces
    for ipsi in 0:mpsi
        psifac = rzphi.xs[ipsi + 1]
        spl = spline_setup(rzphi.ys, rzphi.fs[:,:,1], bctype="periodic")

        # Interpolate to new surface
        for itheta in 0:my
            theta = rzphi.ys[itheta + 1]
            f, fx, fy = bicube_eval(rz_in, psifac, theta, 1)
            rfac = sqrt(max(f[1], 0))
            r = ro + rfac * cos(2π * (theta + f[2]))
            # Prepare spl values for storage
            # ...
        end

        spline_fit(spl, "periodic")

        # Interpolate to final storage
        for itheta in 0:my
            theta = rzphi.ys[itheta + 1]
            f = spline_eval(spl, theta, 0)
            rzphi.fs[ipsi + 1, itheta + 1, :] = f[1:4]
        end
    end

    println("Finished inverse equilibrium process.")
    return PlasmaEquilibrium(equil_input, sq, rzphi, nothing, ro, zo, psio)
end

end # module Inverse

# Function to choose grid based on type
function eq_grid_choice(xs, grid_type)
    # Implement function to choose correct grid transformation
    return xs # Placeholder, actual transformation logic needed
end