module Mc
  function measure_energy(A, J, h)
      energy=0
      for ix in 1:size(A)[1], iy in 1:size(A)[2]
          energy+=-h*A[ix,iy]
          energy += -J/2.0*A[ix,iy]*calc_S(A, ix, iy)
      end
      return energy
  end

  function calc_S(A, ix, iy)
    Lx = size(A)[1]
    Ly = size(A)[2]
    S = 0.0
    S+=A[mod1(ix-1,Lx), iy]
    S+=A[mod1(ix+1,Lx), iy]
    S+=A[ix, mod1(iy-1,Ly)]
    S+=A[ix, mod1(iy+1,Ly)]
    return S
  end

  function calc_ΔE(A, J, h, ix, iy)
    ΔE = 0.0
    ΔE += 2.0*J*A[ix,iy]*calc_S(A, ix, iy)
    ΔE += 2.0*h*A[ix, iy]
    return ΔE
  end

  function metropolis(ΔE, T)
    is_accepted=ifelse(rand()<= exp(-ΔE/T), true, false)
    return is_accepted
  end

  function local_metropolis_update(A, J, h, T, ix, iy)
    ΔE = Mc.calc_ΔE(A, J, h, ix, iy)
    is_accepted=ifelse(metropolis(ΔE, T), -A[ix, iy], A[ix,iy])
    A[ix,iy]=ifelse(metropolis(ΔE, T), -A[ix, iy], A[ix,iy])
    return A, is_accepted
  end
end