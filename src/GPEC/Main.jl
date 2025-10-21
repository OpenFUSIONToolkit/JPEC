function Main(path::String, euler_solution::Union(Nothing,DCON.DconSolution) = nothing)

    if euler_solution !== nothing
        euler = euler_solution
    else
        euler = DCON.Main(path)
    end

    println("IPEC START")
    println("----------------------------------")
    start_time = time()

    # Read input files
    inputtoml = TOML.parsefile(joinpath(path, "gpec.toml"))
    input = GpecInputParameters(; (Symbol(k) => v for (k, v) in inputtoml["GPEC_INPUT"])...)
    ctrl = GpecControlParameters(; (Symbol(k) => v for (k, v) in inputtoml["GPEC_CONTROL"])...)
    outp = GpecOutputParameters(; (Symbol(k) => v for (k, v) in inputtoml["GPEC_OUTPUT"])...)

    # TODO: Call field_bs_psi to compute the field at the control surface



end