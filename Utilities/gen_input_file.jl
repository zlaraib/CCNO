
# This file extracts the initial condition from the test script and generates an input file

function generate_inputs_file(directory, filename, data)
    filepath = joinpath(directory, filename)
    
    # Open the file in write mode
    open(filepath, "w") do file
        # Write data to the file
        for line in data
            println(file, line)
        end
    end
end

function extract_initial_conditions()
    # Putting all my initial conditions into a dictionary
    initial_conditions = Dict(
        "N_sites" => N_sites,
        "N_sites_eachflavor" => N_sites_eachflavor,
        "τ" => τ,
        "ttotal" => ttotal,
        "tolerance" => tolerance,
        "Δm²" => Δm²,
        "maxdim" => maxdim,
        "cutoff" => cutoff,
        "x" => x,
        "p" => p,
        "ψ₀" => ψ₀,
        "L" => L,
        "Δx" => Δx,
        "n_νₑ" => n_νₑ,
        "n_νₑ̄" => n_νₑ̄,
        "Eνₑ" => Eνₑ,
        "Eνₑ̄" => Eνₑ̄,
        "B" => B,
        "N" => N,
        "shape_name" => shape_name,
        "Δp" => Δp,
        "periodic" => periodic
    )
    

    # Convert the dictionary to a string and write it into the input file
    input_data = []
    for (key, value) in initial_conditions
        push!(input_data, "$key = $value")
    end

    return input_data
end
