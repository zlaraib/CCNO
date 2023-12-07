include("main_self_interaction.jl")

using BenchmarkTools
using Measures
using DelimitedFiles

# Define the parameters for N range
start_N = 4
end_N = 16
increment = 2

# Specify the relative directory path
datadir = joinpath(@__DIR__, "misc","datafiles","w 2tmine", "par_"*string(start_N), "_"*string(end_N))
# check if a directory exists, and if it doesn't, create it using mkpath
isdir(datadir) || mkpath(datadir)
plotdir = joinpath(@__DIR__, "misc","plots","w 2tmine", "par_"*string(start_N), "_"*string(end_N))
# check if a directory exists, and if it doesn't, create it using mkpath
isdir(plotdir) || mkpath(plotdir)
function run_main_loop(start_N, end_N, increment)
    times_no_threaded = Float64[]
    times_with_threaded = Float64[]
    N_values = Int[]
    
    for N in start_N:increment:end_N
        println("Running main for N = ", N)
        
        time_no_threaded = @time main(N; use_threaded_blocksparse=false)
        push!(times_no_threaded, time_no_threaded)
        
        time_with_threaded = @time main(N; use_threaded_blocksparse=true)
        push!(times_with_threaded, time_with_threaded)
        
        push!(N_values, N)
    end
    # Writing data to files with corresponding headers
    fname1 = joinpath(datadir, "N_tnothreaded_twiththreaded.dat")
    writedlm(fname1, [N_values times_no_threaded times_with_threaded])

    plot(N_values, times_no_threaded, label="Without threaded block sparse", xlabel="N", ylabel="Computation Time", legend=:topleft, left_margin=15mm, bottom_margin=15mm)
    plot!(N_values, times_with_threaded, label="With threaded block sparse")
    savefig(joinpath(plotdir,"Benchmark performance analysis w 4tmine.pdf"))
end

# Call the function to run main in a loop for different N values
run_main_loop(start_N, end_N, increment)

