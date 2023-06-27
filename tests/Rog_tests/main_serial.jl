using ITensors
using Plots
using Measures
#include("src/gates_function.jl")
include("src/expect.jl")

function main()
    N = 5 # number of sites NEED TO GO TILL 96
    cutoff = 1E-14 # specifies a truncation threshold for the SVD in MPS representation THE SMALLER THE BETTER BECUASE SMALL CUTOFF MEANS MORE ENTANGLEMENT
    tau = 0.05 # time step NEED TO BE 0.05
    ttotal = 10 # total time of evolution NEED TO GO TILL 50
    tolerance  = 5E-1 
    s = siteinds("S=1/2", N; conserve_qns=true)  

    # Constants for the curve
    a_t = 0
    b_t = 2.10
    c_t = 0
    
    # # Specify the directory path
    # #directory_path = "/home/zohalaraib/Test_rep/tests/Rog_tests"
    # directory_path = joinpath(@__DIR__)

    # # Create the file path within the specified directory
    # datafile_path = joinpath(directory_path, "datafiles", string(N) * "(par)_" * string(ttotal) * "(ttotal)final.txt")
    
    # # Open the file for writing
    # datafile = open(datafile_path, "w")
    
    Sz_array, prob_surv_array = calc_expect(s, tau, N, cutoff, ttotal)

    # close(datafile)  # Close the file

    i_min = argmin(prob_surv_array)
    t_min = tau * i_min - tau

    t_p_Rog = a_t*log(N) + b_t * sqrt(N) + c_t
    println("t_p_Rog= ",t_p_Rog)
    println("i_min =", i_min)
    println("t_min= ", t_min)
    @assert abs(t_min - t_p_Rog) <  tau + tolerance

    # # Plotting P_surv vs t
    # plot(0.0:tau:tau*(length(prob_surv_array)-1), prob_surv_array, xlabel = "t", ylabel = "prob_surv", legend = false, size=(800, 600), aspect_ratio=:auto,margin= 10mm) 


    # # Save the plot in the same directory
    # plot_path = joinpath(directory_path, "plots", string(N) * "(par)_" * string(ttotal) * "(ttotal)final.pdf")

    # savefig(plot_path)
end 

@time main()
