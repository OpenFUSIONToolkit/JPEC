


using Test
using Pkg

# Activate the project environment one level up
Pkg.activate(joinpath(@__DIR__, ".."))
using JPEC

@testset "Make fortran files" begin
    try
        # Run the build function
        result = JPEC.build_fortran()
        @test true
        @info "Fortran builds completed successfully"
    catch e
        @test false
        @warn "build_spline_fortran() failed: $e"
    end
end


@testset "splines" begin
    # This lines are  test for splines
end

@testset "fourfit" begin
    # This lines are test for fourfit
end

@testset "vacuum" begin
    # This lines are test for vacuum code
end
