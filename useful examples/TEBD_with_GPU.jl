using MKL
using Pkg
using ITensors
using CUDA
using CSV
using DataFrames
using DelimitedFiles
using LinearAlgebra
using SpecialFunctions
using IterTools
ITensors.Strided.disable_threads()
using ITensors.HDF5
using Printf
using NDTensors


# define custom Hilbert space
import ITensors.op
import ITensors.state

"""
    space(::SiteType"mySpin";
          dim = 2,
          conserve_qns = false,
          conserve_number = false,
          qnname_number = "Number")

Create the Hilbert space for a site of type "mySpin".

Optionally specify the conserved symmetries and their quantum number labels.
"""

function ITensors.space(
  ::SiteType"mySpin";
  dim=2,
)
  return dim
end

function ITensors.op(::OpName"Id", ::SiteType"mySpin", ds::Int...)
  d = prod(ds)
  return Matrix(1.0I, d, d)
end
op(on::OpName"I", st::SiteType"mySpin", ds::Int...) = op(OpName"Id"(), st, ds...)

function op(::OpName"J-", ::SiteType"mySpin", d::Int)
  mat = zeros(d, d)
  for k in 1:(d - 1)
    mat[k+1,k]=sqrt(k*(d-k))
  end
  return mat
end


function op(::OpName"J+", ::SiteType"mySpin", d::Int)
  mat = zeros(d, d)
  for k in 1:(d - 1)
    mat[k,k+1]=sqrt(k*(d-k))
  end
  return mat
end


function op(::OpName"Jz2", ::SiteType"mySpin", d::Int)
    J=(d-1)/2;
    mat=diagm(J:-1:-J);
  return mat*mat
end


# interface
function op(on::OpName, st::SiteType"mySpin", s1::Index, s_tail::Index...; kwargs...)
  rs = reverse((s1, s_tail...))
  ds = dim.(rs)
  opmat = op(on, st, ds...; kwargs...)
  return itensor(opmat, prime.(rs)..., dag.(rs)...)
end

function op(on::OpName, st::SiteType"mySpin"; kwargs...)
  return error("`op` can't be called without indices or dimensions.")
end


function HamLattice_OBC(N, sites, lambda, v)

    # construct Hamiltonian

    os = OpSum()
    for j=1:N
        os += 1.0 , "Jz2", j
    end

    for j=1:N
        os += lambda, "J+", j
        os += lambda', "J-", j
    end

    for j=1:N-1
        os += - v, "J+", j, "J-", j+1
    end

    for j=1:N-1
        os += - v, "J-", j, "J+", j+1
    end
    
    H = MPO(os, sites)

    return H
    
end

function Gates_diag(N, sites, lambda, v, tau)

    # define single-site gates
    
    gates_diag = ITensor[]
    
    for j in 1:N
        s= sites[j]
        h = op("Jz2", s)
        push!(gates_diag, exp(-im * tau * h))
    end
    
    return gates_diag
                
end
            
function GatesOdd_offdiag(N, sites, lambda, v, tau)

    # define two-site gates on odd bonds

    gates_odd_offdiag = ITensor[]
    
    for j in 1:2:(N-1)
        sl= sites[j]
        sr = sites[j+1]
        
        if j==1
            h = lambda * op("J+", sl)  * op("I", sr)
            h += lambda' * op("J-", sl)  * op("I", sr)
        else
            h = lambda * op("J+", sl) * op("I", sr)
            h += lambda' * op("J-", sl)  * op("I", sr)
        end
        
        if j==N-1 && iseven(N)
            h += lambda * op("I", sl) *  op("J+", sr) 
            h += lambda' * op("I", sl) *  op("J-", sr)
        else
            h += lambda * op("I", sl) *  op("J+", sr) 
            h += lambda' * op("I", sl) *  op("J-", sr)
        end
        
        h += - v * op("J+", sl) * op("J-", sr)
        h += - v * op("J-", sl) * op("J+", sr)
        
        push!(gates_odd_offdiag, exp(-im * tau * h))
    end
    
    return gates_odd_offdiag
                
end
        
function GatesEven_offdiag(N, sites, lambda, v, tau)

    # define two-site gates on even bonds

    gates_even_offdiag = ITensor[]
    
    
    for j in 2:2:(N-1)
        
        sl= sites[j]
        sr = sites[j+1]
        
        h = lambda * op("J+", sl) * op("I", sr)
        h += lambda' * op("J-", sl) * op("I", sr)
        
        if j==N-1 && isodd(N)
            h += lambda  * op("I", sl) * op("J+", sr) 
            h += lambda'  * op("I", sl) * op("J-", sr)
        else
            h += lambda  * op("I", sl) *  op("J+", sr)
            h += lambda' * op("I", sl) *  op("J-", sr) 
        end
        
        h += -v * op("J+", sl) * op("J-", sr)
        h += -v * op("J-", sl) * op("J+", sr)
        
        push!(gates_even_offdiag, exp(-im * tau * h))
    
    end
    
    return gates_even_offdiag
                
end


function TEBD_obs_GPU(psi_start, N, sites, H, tau, T, lambda, v; restart_time=0, cutoff=1e-10, maxdim=128, psi_0::MPS=psi_start)

    # here I am evolving with TEBD and using GPU

    @show cutoff
    @show maxdim
    
    steps=Int64(ceil(abs(T)/abs(tau)))
        
    start_step=1
    # psi_init=cu(psi_start)    
    psi=cu(psi_start)

    gates_odd_offdiag=GatesOdd_offdiag(N, sites, lambda, v, tau)
    gates_even_offdiag=GatesEven_offdiag(N, sites, lambda, v, tau)
    gates_diag=Gates_diag(N, sites, lambda, v, tau/2)

    gates_odd_offdiag=cu.(gates_odd_offdiag)
    gates_even_offdiag=cu.(gates_even_offdiag)
    gates_diag=cu.(gates_diag)
        
    for i in start_step:steps
        
        idx_time=round(i*tau, digits=2)

        time_elapsed = @elapsed begin
            psi=apply(gates_diag, psi; cutoff=cutoff, maxdim=maxdim)
            normalize!(psi)

            psi=apply(gates_odd_offdiag, psi; cutoff=cutoff, maxdim=maxdim)
            normalize!(psi) 

            psi=apply(gates_even_offdiag, psi; cutoff=cutoff, maxdim=maxdim)
            normalize!(psi) 

            psi=apply(gates_diag, psi; cutoff=cutoff, maxdim=maxdim)
            normalize!(psi)
        end
        
        energy=inner(psi', H, psi)
        
        println("\nStep n. $(i), \t time=$(idx_time), \t energy=$(energy), \t maxlinkdim=$(maxlinkdim(psi)), \t elapsed=$(time_elapsed)")

        flush(stdout)

        GC.gc()

    end

    return psi
        
end

function TEBD_obs_CPU(psi_start, N, sites, H, tau, T, lambda, v; restart_time=0, cutoff=1e-10, maxdim=128, psi_0::MPS=psi_start)
    
    #here I am evolving with TEBD and using standard CPU

    @show cutoff
    @show maxdim
    
    steps=Int64(ceil(abs(T)/abs(tau)))
        
    start_step=1
    psi=psi_start

    gates_odd_offdiag=GatesOdd_offdiag(N, sites, lambda, v, tau)
    gates_even_offdiag=GatesEven_offdiag(N, sites, lambda, v, tau)
    gates_diag=Gates_diag(N, sites, lambda, v, tau/2)
        
    for i in start_step:steps
        
        idx_time=round(i*tau, digits=2)

        time_elapsed = @elapsed begin
            psi=apply(gates_diag, psi; cutoff=cutoff, maxdim=maxdim)
            normalize!(psi)

            psi=apply(gates_odd_offdiag, psi; cutoff=cutoff, maxdim=maxdim)
            normalize!(psi) 

            psi=apply(gates_even_offdiag, psi; cutoff=cutoff, maxdim=maxdim)
            normalize!(psi) 

            psi=apply(gates_diag, psi; cutoff=cutoff, maxdim=maxdim)
            normalize!(psi)
        end
        
        energy=inner(psi', H, psi)
        
        println("\nStep n. $(i), \t time=$(idx_time), \t energy=$(energy), \t maxlinkdim=$(maxlinkdim(psi)), \t elapsed=$(time_elapsed)")

        flush(stdout)

        GC.gc()

    end

    return psi
        
end


N=11
J=6
tau=0.01
T=1

dimension=Int(2*J+1)
sites=siteinds("mySpin", dim=dimension,N);

psi_init = randomMPS(sites,10);

v1=4*13^2/pi^2;
v=v1/(J*(J+1));
l1=0.3^2*v1/2;
l=-l1/(2*(J*(J+1))^(1/2));

#construct Hamiltonian with open boundary conditions
H=HamLattice_OBC(N, sites, l, v)


# ground state minimization with DMRG
obs = DMRGObserver(; energy_tol = 1E-7, minsweeps = 2)
energy, psi_init = dmrg(H,psi_init, nsweeps=25, maxdim=100, cutoff=1e-8; observer=obs)
@show energy


# TEBD with CPU
psi=TEBD_obs_CPU(psi_init, N, sites, H, tau, T, l, v; cutoff=1e-9, maxdim=128)


# # TEBD with GPU 
# H= NDTensors.cu(H)
# psi=TEBD_obs_GPU(psi_init, N, sites, H, tau, T, l, v; cutoff=1e-9, maxdim=128)