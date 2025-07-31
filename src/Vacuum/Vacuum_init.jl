"""
    build_vacuum_globals(
        mthvac::Int,                # Number of poloidal grid points (like nths0 in Fortran)
        mpert::Int,                 # Number of Fourier harmonics (like nfm)
        settings::VacuumSettingsType,
        input::VacuumInputType
    ) -> VacuumGlobalsType

Constructs a VacuumGlobalsType, mimicking the Fortran workflow.

WIP WIP WIP WIP WIP WIP WIP WIP WIP WIP WIP WIP WIP WIP WIP WIP 

"""

using Interpolations

function build_vacuum_globals(
    mthvac::Int,
    mpert::Int,
    wall::Bool,
    farwal::Bool,
    kernelsign::Float64,
    settings::VacuumSettingsType,
    input::VacuumInputType,
)

    #-------
    # Fortran's defglo and cardmo
    #-------
    # These are the canonical Fortran variables:
    mth = mthvac
    mth1 = mth + 1
    mth2 = mth1 + 1
    nfm = mpert
    mtot = mpert
    dth = 2π / mth


    # 
    # Mode index arrays (lmin/lmax): get from input or compute
    lmin = zeros(1)
    lmax = zeros(1)

    lmin[1] = input.mlow
    lmax[1] = input.mhigh

    n = input.n
    thgr = range(0, stop=2π, length=mth1)

    # Arrays from input
    xinf = copy(input.r)
    zinf = copy(input.z)
    delta = copy(input.delta)

    # Physical derived parameters
    qa1 = input.qa
    ga1 = 1.0  # Placeholder, fill with correct logic as needed
    fa1 = 1.0  # Placeholder, fill with correct logic as needed

    # Wall arrays
    xwal = zeros(mth1)
    zwal = zeros(mth1)
    xwalp = zeros(mth1)
    zwalp = zeros(mth1)

    if !farwal
        if settings.Shape.a >= 10.
            farwal = true
        end
    end

    delx, delz, xwal, zwal, xwalp, zwalp,
        cnqd, snqd, sinlt, coslt, snlth, 
        cslth, xplap, zplap = setuparrays!(
            mth, dth, thgr, lmin, lmax,
            qa1, n, settings.Shape.delfac,
            delta, input.xpla, input.zpla
        )


    
    

    return VacuumGlobalsType(
        mth=mth,
        mth1=mth1,
        mth2=mth2,
        nfm=nfm,
        mtot=mtot,
        lmin=lmin,
        lmax=lmax,
        xinf=xinf,
        zinf=zinf,
        delta=delta,
        xplap=xplap,
        zplap=zplap,
        qa1=qa1,
        ga1=ga1,
        fa1=fa1,
        delx=delx,
        delz=delz,
        xwal=xwal,
        zwal=zwal,
        xwalp=xwalp,
        zwalp=zwalp,
        dth=dth,
        wall=wall,
        farwal=farwal,
        kernelsign=kernelsign
    )
end

function setuparrays!(
    mth::Int,
    dth::Float64,
    thgr::Vector{Float64},
    lmin::Vector{Int},
    lmax::Vector{Int},
    q::Float64,
    n::Int,
    delfac::Float64,
    delta::Vector{Float64},
    xpla::Vector{Float64},
    zpla::Vector{Float64},
)
    # Sizes
    mth1 = mth + 1
    jmax1 = lmax[1] - lmin[1] + 1
    nq = n * q

    # Compute geometry bounds
    xmin = minimum(xpla)
    xmax = maximum(xpla)
    # zmin = minimum(zpla)
    # zmax = maximum(zpla)

    plrad = 0.5 * (xmax - xmin)
    delx = plrad * delfac
    delz = plrad * delfac

    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # Wall shape (placeholder: straight line wall matching plasma)
    # This is a placeholder for the wall shape, which should be replaced
    # with the actual wall shape logic in wwall.
    xwal = copy(xpla[1:mth1])
    zwal = copy(zpla[1:mth1])
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    # Derivatives (simple finite difference instead of spline)
    # Will convert to spline once equil branch is merged
    # xwalp = [ (xwal[mod1(i+1,mth)] - xwal[i]) / dth for i in 1:mth ]
    # zwalp = [ (zwal[mod1(i+1,mth)] - zwal[i]) / dth for i in 1:mth ]

    itp_x = Interpolations.interpolate((theta,), xwal, Interpolations.BSpline(Interpolations.Cubic(Interpolations.Periodic(OnGrid()))))
    itp_z = Interpolations.interpolate((theta,), zwal, Interpolations.BSpline(Interpolations.Cubic(Interpolations.Periodic(OnGrid()))))

    xwalp = [Interpolations.derivative(itp_x, 1, theta) for theta in thgr[:mth1]]
    zwalp = [Interpolations.derivative(itp_z, 1, theta) for theta in thgr[:mth1]]

    itp_x = Interpolations.interpolate((theta,), xpla, Interpolations.BSpline(Interpolations.Cubic(Interpolations.Periodic(OnGrid()))))
    itp_z = Interpolations.interpolate((theta,), zpla, Interpolations.BSpline(Interpolations.Cubic(Interpolations.Periodic(OnGrid()))))

    xplap = [Interpolations.derivative(itp_x, 1, theta) for theta in thgr[:mth1]]
    zplap = [Interpolations.derivative(itp_z, 1, theta) for theta in thgr[:mth1]]

    push!(xplap, xwalp[1])
    push!(zplap, zwalp[1])

    # Allocate arrays
    cnqd = zeros(Float64, mth1)
    snqd = zeros(Float64, mth1)
    sinlt = zeros(Float64, mth1, jmax1)
    coslt = zeros(Float64, mth1, jmax1)
    snlth = zeros(Float64, mth1, jmax1)
    cslth = zeros(Float64, mth1, jmax1)

    # Trigonometric basis arrays
    for is in 1:mth1
        theta = (is-1) * dth
        znqd = nq * delta[is]
        cnqd[is] = cos(znqd)
        snqd[is] = sin(znqd)
        for l1 in 1:jmax1
            ll = lmin[1] - 1 + l1
            elth = ll * theta
            elthnq = ll * theta + znqd
            sinlt[is,l1] = sin(elth)
            coslt[is,l1] = cos(elth)
            snlth[is,l1] = sin(elthnq)
            cslth[is,l1] = cos(elthnq)
        end
    end

    return (
        delx, delz, xwal, zwal, xwalp, zwalp,
        cnqd, snqd, sinlt, coslt, snlth, cslth, xplap, zplap
    )
end
