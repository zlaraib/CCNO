using ITensors
using ITensorMPS

# pass in the simulationstate to have access to both lists of indices
# brute force search through data to output the values for each site based on their original ordering
function write_ordered_array(state::CCNO.SimulationState, filename::String, t::Float64, data)
    Nsites = length(state.s)
    @assert Nsites == length(data)

    # create the ouput array
    data_reordered = Vector{Float64}(undef, Nsites)
    
    # loop through the output array.
    # for each element of the output array, find the data corresponding to the index value stored in s0
    for i in 1:Nsites
        for j in 1:Nsites
            if state.s0[i] == state.s[j]
                data_reordered[i] = data[j]
            end
        end
    end

    # attach time time value and write transposed array as a single row
    open(filename, "a") do f
        writedlm(f, [t transpose(data_reordered)])
    end
end

function store_data(datadir::String, t::Float64, state::CCNO.SimulationState)
    #================================#
    # calculate necessary quantities #
    #================================#

    # compute the avg expectation value of Sz at all sites
    sz_arr = expect(state.ψ, "Sz")  # Compute Sz for each site and store the values in sz_arr
    
    # compute expectation value of sy and sx (inbuilt operator in ITensors library) at all sites on the chain
    sy_arr = expect(complex(state.ψ), "Sy")
    sx_arr = expect(state.ψ, "Sx")
    
    #survival probability for all sites (neutrino) to be found in its initial flavor state
    #prob_surv_arr = 0.5 * (1 .- 2 .* sz_arr)

    # recall that in our code sigma_z = 2*Sz so make sure these expressions are consistent with "Sz in ITensors" 
    ρₑₑ_arr = ((2 .* sz_arr) .+ 1) ./ 2
    
    ρ_μμ_arr = ((-2 .* sz_arr) .+ 1) ./ 2
    
    ρₑμ_arr = sqrt.(sx_arr.^2 .+ sy_arr.^2)
    
    #===============#
    # save the data #
    #===============#
    
    # Writing data to files with corresponding headers
    # 1. Sz arrays
    write_ordered_array(state,
                        joinpath(datadir, "t_<Sz>.dat"),
                        t,
                        sz_arr)

    # 2. Sy arrays
    write_ordered_array(state,
                        joinpath(datadir, "t_<Sy>.dat"),
                        t,
                        sy_arr)

    # 3. Sx arrays
    write_ordered_array(state,
                        joinpath(datadir, "t_<Sx>.dat"),
                        t,
                        sx_arr)

    # 5. xsite values
    write_ordered_array(state,
                        joinpath(datadir, "t_xsiteval.dat"),
                        t,
                        state.xyz[:,1])
    write_ordered_array(state,
                        joinpath(datadir, "t_ysiteval.dat"),
                        t,
                        state.xyz[:,2])
    write_ordered_array(state,
                        joinpath(datadir, "t_zsiteval.dat"),
                        t,
                        state.xyz[:,3])

    # 6. pxsite values
    write_ordered_array(state,
                        joinpath(datadir, "t_pxsiteval.dat"),
                        t,
                        state.p[:,1])
    write_ordered_array(state,
                        joinpath(datadir, "t_pysiteval.dat"),
                        t,
                        state.p[:,2])
    write_ordered_array(state,
                        joinpath(datadir, "t_pzsiteval.dat"),
                        t,
                        state.p[:,3])
    
    # 7. ρₑₑ arrays
    write_ordered_array(state,
                        joinpath(datadir, "t_ρₑₑ.dat"),
                        t,
                        ρₑₑ_arr)

    # 8. ρ_μμ arrays
    write_ordered_array(state,
                        joinpath(datadir, "t_ρ_μμ.dat"),
                        t,
                        ρ_μμ_arr)

    # 9. ρₑμ arrays
    write_ordered_array(state,
                        joinpath(datadir, "t_ρₑμ.dat"),
                        t,
                        ρₑμ_arr)
    
end 
