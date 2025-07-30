"""
    mercier_scan()

Evaluates Mercier criterion and related quantities for MHD stability analysis.

# Global Variables Used:
- `r8`: Real precision type (likely Float64) - from DconMod
- `mtheta`, `mpsi`: Grid dimensions - from DconMod  
- `twopi`: 2π constant - from DconMod
- `psio`: Normalization constant - from DconMod
- `ro`: Reference radius - from DconMod
- `rzphi`: 2D spline structure for coordinates - from DconMod
- `sq`: 1D spline structure for surface quantities - from DconMod
- `locstab`: Output spline structure for stability results - from DconMod

# External Methods Used:
- `spline_alloc()`: Allocates spline structure - from DconMod
- `spline_dealloc()`: Deallocates spline structure - from DconMod
- `bicube_eval()`: Evaluates bicubic spline - from DconMod
- `spline_fit()`: Fits spline with periodic boundary conditions - from DconMod
- `spline_int()`: Integrates splined functions - from DconMod
"""
function mercier_scan(plasma_eq)
    
    # Local variables
    #local bsq::Float64, chi1::Float64, di::Float64, dpsisq::Float64
    #local eta::Float64, h::Float64, jac::Float64, p1::Float64
    #local psifac::Float64, q::Float64, q1::Float64, r::Float64
    #local rfac::Float64, term::Float64, theta::Float64, twopif::Float64
    #local v1::Float64, v2::Float64, v21::Float64, v22::Float64
    #local v23::Float64, v33::Float64
    
    sq = plasma_eq.sq
    rzphi = plasma_eq.rzphi

    # Prepare spline types
    psis = sq.xs
    ff_fs = zeros(mpsi, 5)
    locstab = zeros(mpsi, 3)
    #ff.xs = rzphi.ys  # Copy coordinate array
    
    # Compute surface quantities
    for ipsi in 1:mpsi
        psifac = sq.xs[ipsi]
        twopif = sq.fs[ipsi, 1]   # Surface quantity 1
        p1 = sq.fs1[ipsi, 2]     # Derivative of surface quantity 2
        v1 = sq.fs[ipsi, 3]      # Surface quantity 3
        v2 = sq.fs1[ipsi, 3]     # Derivative of surface quantity 3
        q = sq.fs[ipsi, 4]       # Safety factor
        q1 = sq.fs1[ipsi, 4]     # Derivative of safety factor
        chi1 = 2π * plasma_eq.psio         # Normalized flux
        
        # Evaluate coordinates and jacobian
        for itheta in 1:mtheta
            # Evaluate bicubic spline at grid point
            Spl.bicube_eval(rzphi, rzphi.xs[ipsi], rzphi.ys[itheta], 1)
            
            theta = rzphi.ys[itheta]
            rfac = sqrt(rzphi.f[1])
            eta = twopi * (theta + rzphi.f[2])
            r = ro + rfac * cos(eta)
            jac = rzphi.f[4]  # Jacobian
            
            # Evaluate other local quantities (metric components)
            v21 = rzphi.fy[1] / (2 * rfac * jac)
            v22 = (1 + rzphi.fy[2]) * twopi * rfac / jac
            v23 = rzphi.fy[3] * r / jac
            v33 = twopi * r / jac
            
            # Magnetic field squared and flux derivative squared
            bsq = chi1^2 * (v21^2 + v22^2 + (v23 + q * v33)^2)
            dpsisq = (twopi * r)^2 * (v21^2 + v22^2)
            
            # Evaluate integrands for Mercier analysis
            ff_fs[itheta, 1] = bsq / dpsisq
            ff_fs[itheta, 2] = 1 / dpsisq
            ff_fs[itheta, 3] = 1 / bsq
            ff_fs[itheta, 4] = 1 / (bsq * dpsisq)
            ff_fs[itheta, 5] = bsq
            
            # Weight by jacobian and volume element
            ff_fs[itheta, :] .*= jac / v1
        end
            
        ff = Spl.spline_setup(psis, ff_fs; bctype=2) # bctype=2 is periodic

        # Integrate quantities with respect to theta using periodic splines
      #  spline_fit(ff, "periodic")  # Fit with periodic boundary conditions
      #  spline_int(ff)              # Integrate the fitted splines

        avg = ff.fsi[mtheta, :]   # Get integrated values
        
        # Evaluate Mercier criterion and related quantities
        term = twopif * p1 * v1 / (q1 * chi1^3) * avg[2]
        
        # Mercier stability parameter
        di = -0.25 + term * (1 - term) + 
             p1 * (v1 / (q1 * chi1^2))^2 * avg[1] * 
             (p1 * (avg[3] + (twopif / chi1)^2 * avg[4]) - v2 / v1)
        
        # Additional stability parameter
        h = twopif * p1 * v1 / (q1 * chi1^3) * (avg[2] - avg[1] / avg[5])
        
        # Store results in output spline structure
        locstab_fs[ipsi, 1] = di * psis # is this the right psis?
        locstab_fs[ipsi, 2] = (di + (h - 0.5)^2) * psis # is this the right psis?
        locstab_fs[ipsi, 3] = h
    end
    
    return locstab
end

end # module MercierMod