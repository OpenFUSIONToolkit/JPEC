module VacuumMod

export mscvac, set_dcon_params


const libdir = joinpath(@__DIR__, "..", "..", "deps")
const libvac = joinpath(libdir, "libvac")

function set_dcon_params(mthin::Int32, lmin::Int32, lmax::Int32, nnin::Int32,
                        qa1in::Float64,
                        xin::Vector{Float64}, zin::Vector{Float64}, deltain::Vector{Float64})

    ccall((:set_dcon_params_, libvac),
          Nothing,
          (Ref{Cint}, Ref{Cint}, Ref{Cint}, Ref{Cint},
           Ref{Cdouble},
           Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
          Ref(mthin), Ref(lmin), Ref(lmax), Ref(nnin),
          Ref(qa1in),
          pointer(xin), pointer(zin), pointer(deltain))
end



function mscvac(
    wv::Array{ComplexF64,2},
    mpert::Int32,
    mtheta::Int32,
    mthvac::Int32,
    complex_flag::Bool,
    kernelsignin::Float64,
    wall_flag::Bool,
    farwal_flag::Bool,
    grrio::Array{Float64,2},
    xzptso::Array{Float64,2},
    op_ahgfile::Union{Nothing,String}=nothing
)

    ahgfile_ptr = if op_ahgfile === nothing
        Ptr{UInt8}(C_NULL)
    else
        Base.cconvert(Ptr{UInt8}, op_ahgfile)
    end

    ccall((:__vacuum_mod_MOD_mscvac, libvac),
        Nothing,
        (Ptr{ComplexF64},       # wv(mpert,mpert)
         Ref{Cint},            # mpert
         Ref{Cint},            # mtheta
         Ref{Cint},            # mthvac
         Ref{Cint},            # complex_flag (logical)
         Ref{Cdouble},         # kernelsignin
         Ref{Cint},            # wall_flag (logical)
         Ref{Cint},            # farwal_flag (logical)
         Ptr{Cdouble},         # grrio(:,:)
         Ptr{Cdouble},         # xzptso(:,:)
         Ptr{UInt8}),          # op_ahgfile (optional)
        pointer(wv),
        Ref(mpert),
        Ref(mtheta),
        Ref(mthvac),
        Ref(Int32(complex_flag)),
        Ref(kernelsignin),
        Ref(Int32(wall_flag)),
        Ref(Int32(farwal_flag)),
        pointer(grrio),
        pointer(xzptso),
        ahgfile_ptr
    )

    return wv
end

end