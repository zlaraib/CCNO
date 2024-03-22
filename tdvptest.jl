using ITensors
using ITensorTDVP

function run_tdvp()
  N = 100
  cutoff = 1E-8
  tau = 0.1
  ttotal = 5.0

  # Make an array of 'site' indices
  s = siteinds("S=1/2", N; conserve_qns=true)

  # Define Hamiltonian as an MPO
  ampo = AutoMPO()
  for j in 1:(N - 1)
    add!(ampo, 0.5, "S+", j, "S-", j+1)
    add!(ampo, 0.5, "S-", j, "S+", j+1)
    add!(ampo, "Sz", j, "Sz", j+1)
  end
  H = MPO(ampo, s)

  # Initialize psi to be a product state (alternating up and down)
  psi = MPS(s, n -> isodd(n) ? "Up" : "Dn")

  c = div(N, 2) # center site

  # Time evolution using TDVP
  for t in 0.0:tau:ttotal
    Sz = expect(psi, "Sz"; sites=c)
    println("$t $Sz")

    t≈ttotal && break

    # Evolve the state using TDVP
    psi = tdvp(H, psi, tau, cutoff=cutoff)
  end

  return
end

run_tdvp()
