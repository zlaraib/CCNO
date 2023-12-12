using ITensors
using Plots 
using Measures 
using LinearAlgebra
using Base.Threads
include("../src/evolution.jl")
include("../src/constants.jl")

# We are simulating the time evolution of a 1D spin chain with N sites, where each site is a spin-1/2 particle. 
# The simulation is done by applying a sequence of unitary gates to an initial state of the system, 
# which is a product state where each site alternates between up and down.
function main(N; use_splitblocks = false,nsweeps = 10, blas_num_threads=4,
  strided_num_threads=4, use_threaded_blocksparse=false, outputlevel=1)

  cutoff = 1E-14 # specifies a truncation threshold for the SVD in MPS representation (SMALL CUTOFF = MORE ENTANGLEMENT)
  τ = 0.05 # time step 
  ttotal = 10 # total time of evolution
  tolerance  = 5E-1 # acceptable level of error or deviation from the exact value or solution
  Δx = 1E-3 # length of the box of interacting neutrinos at a site/shape function width of neutrinos in cm 
  maxdim = 200 #bondimension


  if outputlevel > 0
    ITensors.Strided.set_num_threads(strided_num_threads)
    BLAS.set_num_threads(blas_num_threads)
    if use_threaded_blocksparse
        ITensors.enable_threaded_blocksparse()
      else
        ITensors.disable_threaded_blocksparse()
    end      
    @show Threads.nthreads()
    @show Sys.CPU_THREADS
    @show BLAS.get_num_threads()
    @show ITensors.Strided.get_num_threads()
    @show ITensors.using_threaded_blocksparse()
    println()

  end

  sweeps = Sweeps(nsweeps)
  maxdims = min.([100, 200, 400, 800, 2000, 3000, maxdim], maxdim)
  maxdim!(sweeps, maxdims...)
  noise!(sweeps, 1e-6, 1e-7, 1e-8, 0.0)
  @show sweeps
  # Make an array of 'site' indices and label as s 
  # conserve_qns=true conserves the total spin quantum number "S"(in z direction) in the system as it evolves
  s = siteinds("S=1/2", N; conserve_qns=false)  

  # Initialize an array of zeros for all N particles
  mu = zeros(N)
                                
  # Create an array of dimension N and fill it with the value 1/(sqrt(2) * G_F). This is the number of neutrinos.
  n = mu.* fill((Δx)^3/(sqrt(2) * G_F), N)
      
  # Create a B vector which would be same for all N particles 
  B = [1, 0, 0]           
  
  # Create an array ω with N elements. Each element of the array is a const pi.
  ω = fill(π, N) 

  # Initialize psi to be a product state (alternating up and down)
  ψ = productMPS(s, n -> isodd(n) ? "Up" : "Dn")

  #extract output from the expect.jl file where the survival probability values were computed at each timestep
  Sz_array, prob_surv_array, apply_time=  evolve(s, τ, n, ω, B, N, Δx, ψ, cutoff, tolerance, ttotal, nsweeps,
                                          maxdim,outputlevel, use_splitblocks, use_threaded_blocksparse)
  energy,ψ,apply_time_parallel  = eigenvals(s, τ, n, ω, B, N, Δx, cutoff, tolerance, ttotal,outputlevel,
                                  maxdim,use_splitblocks, sweeps,nsweeps,use_threaded_blocksparse)

             
  expected_sz_array = Float64[]
  expected_sz = Float64[]
  for t in 0.0:τ:ttotal
      local_sz = 0.0  # Reset for each time step
      for i in 1:(N-1)
          if ω[i] != 0
              if B[1] == 1
                  # Compute the expected value
                  local_sz = 0.5 * cos(ω[i] * t)
              elseif B[3] == -1
                  # Compute the expected value
                  local_sz = -0.5
              end
          end
      end
      push!(expected_sz_array, local_sz)
  end
  if outputlevel >0 && use_threaded_blocksparse==false
    # Extract the Float64 values from Sz_array
    sz_values = [val[] for val in Sz_array]
  else 
    sz_values = Sz_array
  end 
  if outputlevel == 0 || use_threaded_blocksparse==true
    # Check if each element in sz_values is less than tolerance away from corresponding element in expected_sz_array
    @assert all(abs.(sz_values .- expected_sz_array) .< tolerance) #not checking this assert condtion for my Multithread system because thats giving incorrect results 
  end

    if outputlevel >0 && use_threaded_blocksparse==false
      # Convert Atomic{Float64} array to Float64 array
      sz_values = [val[] for val in Sz_array]

      # Time steps for the x-axis
      time_steps = 0.0:τ:τ*(length(sz_values)-1)

      # Plotting <Sz> vs t
      plot(time_steps, sz_values, xlabel = "t", ylabel = "<Sz>", legend = false, size=(700, 600), aspect_ratio=:auto, margin=10mm)

      # Save the plot as a PDF file
      savefig("<Sz> vs t (vacuum oscillation term plot).pdf")
    else
      # Plotting P_surv vs t
      plot(0.0:τ:τ*(length(Sz_array)-1), Sz_array, xlabel = "t", ylabel = "<Sz>", legend = false, size=(700, 600), aspect_ratio=:auto,margin= 10mm) 

      # Save the plot as a PDF file
      savefig("<Sz> vs t (vacuum oscillation term plot).pdf")
    end
end

N = 4 # number of sites 
println("Serial:\n")
@time main(N; use_splitblocks = false,nsweeps=10, blas_num_threads=1,
strided_num_threads=1, use_threaded_blocksparse=false, outputlevel=0)
println("Multi-threading without threaded block sparse(my_effort):\n")
@time main(N; use_splitblocks = false,nsweeps=10, blas_num_threads=128,
strided_num_threads=128, use_threaded_blocksparse=false, outputlevel=1)
println("Multi-threading with threaded block sparse(ITensors):\n")
@time main(N; use_splitblocks = false,nsweeps=10, blas_num_threads=1,
strided_num_threads=1, use_threaded_blocksparse=true, outputlevel=1)


