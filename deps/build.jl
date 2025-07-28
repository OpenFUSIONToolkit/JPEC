



# How to modify build.jl:
# 1) Declare build functions for each Fortran codes as 'build_*_fortran()'
# 2) Add these functions to the results array in build_fortran()
# 3) Export the functions you want to be accessible from outside


const parent_dir = joinpath(@__DIR__, "..", "src")

export build_fortran
export build_spline_fortran, build_vacuum_fortran




function build_fortran()

    results = [
        # build_jpec_fortran() add here
        build_spline_fortran(),
        build_vacuum_fortran(),
    ]

    if all(results)
        @info "All Fortran builds succeeded!"
    else
        @error "Some Fortran builds failed."
    end

end

function build_spline_fortran()
    dir = joinpath(parent_dir, "Splines", "fortran")
    try
        run(pipeline(`make -C $dir`))
        @info "Splines-fortran compiled well"
        return true
    catch e
        @error "Failed to build Splines-fortran: $e"
        return false
    end
end

function build_vacuum_fortran()
    dir = joinpath(parent_dir, "Vacuum", "fortran")
    try
        run(pipeline(`make -C $dir`))
        @info "Vacuum-fortran compiled well"
        return true
    catch e
        @error "Failed to build Vacuum-fortran: $e"
        return false
    end

end

# function build_jpec_fortran()
#     dir = joinpath(parent_dir, "Gpec", "fortran")
#     try
#         run(pipeline(`make -C $dir`, stdout=DevNull, stderr=DevNull))
#         println("Gpec-fortran compiled well")
#         return true
#     catch e
#         println("Failed to build Gpec-fortran: $e")
#         return false
#     end
# end


