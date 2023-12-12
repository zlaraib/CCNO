include("main_self_interaction (test_copy_for_running_in_loop_of_N).jl")

using Measures
using DelimitedFiles
using ITensors
using LinearAlgebra
# Define the parameters for N range
start_N = 4
end_N = 12
increment = 2

# Specify the relative directory path
datadir = joinpath(@__DIR__, "misc","datafiles", "par_"*string(start_N), "_"*string(end_N))
# check if a directory exists, and if it doesn't, create it using mkpath
isdir(datadir) || mkpath(datadir)
plotdir = joinpath(@__DIR__, "misc","plots", "par_"*string(start_N), "_"*string(end_N))
# check if a directory exists, and if it doesn't, create it using mkpath
isdir(plotdir) || mkpath(plotdir)
function run_main_loop(start_N, end_N, increment)
    times_serial = Float64[]
    times_no_ITensors_threaded = Float64[]
    times_with_threaded = Float64[]
    N_values = Int[]
    
    for N in start_N:increment:end_N
        println("Running main function for N = ", N)
        
        
        time_serial=@time main(N; use_splitblocks = true,nsweeps=10, blas_num_threads=1,
        strided_num_threads=1, use_threaded_blocksparse=false, outputlevel=0)
        push!(times_serial, time_serial)

       
        time_no_ITensors_threaded = @time main(N; use_splitblocks = true,nsweeps=10, blas_num_threads=128,
                            strided_num_threads=128, use_threaded_blocksparse=false, outputlevel=1)
        push!(times_no_ITensors_threaded, time_no_ITensors_threaded)
        
        
        time_with_threaded = @time main(N; use_splitblocks = true,nsweeps=10, blas_num_threads=1,
                            strided_num_threads=1, use_threaded_blocksparse=true, outputlevel=1)
        push!(times_with_threaded, time_with_threaded)
        
        push!(N_values, N)
    end
    # Writing data to files with corresponding headers
    fname1 = joinpath(datadir, "N_t_Serial_tnoIthreaded_twithIthreaded.dat")
    writedlm(fname1, [N_values times_serial times_no_ITensors_threaded times_with_threaded])

    
    plot(N_values, times_serial, label="Serial", xlabel="N", ylabel="Computation Time(s)", legend=:topleft, left_margin=15mm, bottom_margin=15mm)
    plot!(N_values, times_no_ITensors_threaded, label="Without threaded block sparse(my_MT)", xlabel="N", ylabel="Computation Time(s)", legend=:topleft, left_margin=15mm, bottom_margin=15mm)
    plot!(N_values, times_with_threaded, label="With threaded block sparse(I_MT)")
    savefig(joinpath(plotdir,"Benchmark performance analysis serial vs myMT vs_ITensorsMT.pdf"))
end

# Call the function to run main in a loop for different N values
run_main_loop(start_N, end_N, increment)

