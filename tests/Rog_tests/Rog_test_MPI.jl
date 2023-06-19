

# ENV["MKL_NUM_THREADS"] = "64"
# ENV["OPENBLAS_NUM_THREADS"] = "64"
# ENV["OMP_NUM_THREADS"] = "64"

using ITensors
using Plots
using Measures
using LinearAlgebra

function main(; blas_num_threads = 1, strided_num_threads = 1, use_threaded_blocksparse = true)
    N = 10
    cutoff = 1E-8
    tau = 0.05
    ttotal = 10
    tolerance = 5E-1

    # Constants for the curve
    a_t = 0
    b_t = 2.105
    c_t = 0
    
    BLAS.set_num_threads(128)  # Set the number of threads to 128

    ITensors.Strided.set_num_threads(strided_num_threads)
    BLAS.set_num_threads(blas_num_threads)
    ITensors.enable_threaded_blocksparse(use_threaded_blocksparse)

    s = siteinds("S=1/2", N; conserve_qns=true)

    gates = ITensor[]

    for i in 1:(N-1)
        for j in i+1:N
            s1 = s[i]
            s2 = s[j]
            
            hj = 2.0/N * (op("Sz", s1) * op("Sz", s2) +
                         1/2 * op("S+", s1) * op("S-", s2) +
                         1/2 * op("S-", s1) * op("S+", s2))
            Gj = exp(-im * tau/2 * hj)
            push!(gates, Gj)
        end
    end

    append!(gates, reverse(gates))

    global psi = productMPS(s, n -> isodd(n) ? "Dn" : "Up")

    c = div(N, 2)

    Sz_array = Float64[]
    prob_surv_array = Float64[]

    # directory_path = "/home/zohalaraib/Test_rep/tests/Rog_tests"
    # datafile_path = joinpath(directory_path, "datafiles", string(N) * "(par)_" * string(ttotal) * "(ttotal)MPI.txt")

    # datafile = open(datafile_path, "w")

    for t in 0.0:tau:ttotal
        sz = expect(psi, "Sz"; sites=1)
        push!(Sz_array, sz)

        t ≈ ttotal && break

        global psi = apply(gates, psi; cutoff)
        normalize!(psi)

        prob_surv = 0.5 * (1 - 2 * sz)
        push!(prob_surv_array, prob_surv)
        println("$t $prob_surv")

        # println(datafile, "$t $prob_surv")
        # flush(datafile)
    end

    #close(datafile)
    i_min = argmin(prob_surv_array)
    t_min = tau * i_min - tau
    
    t_p_Rog = a_t*log(N) + b_t * sqrt(N) + c_t
    println("t_p_Rog= ",t_p_Rog)
    println("i_min =", i_min)
    println("t_min= ", t_min)
    @assert abs(t_min - t_p_Rog) <  tau + tolerance 
    

    # plot(0.0:tau:tau * (length(prob_surv_array) - 1), prob_surv_array,
    #     xlabel = "t", ylabel = "prob_surv", legend = false, size=(800, 600), aspect_ratio=:auto, margin=10mm)

    # plot_path = joinpath(directory_path, "plots", string(N) * "(par)_" * string(ttotal) * "(ttotal)MPI.png")
    # savefig(plot_path)
end


# Run the main function
@time main(use_threaded_blocksparse = true)
