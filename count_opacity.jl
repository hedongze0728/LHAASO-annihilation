using QuadGK

const melectron = 9.109e-31 # kg
const relectron = 2.81794e-15 # m
const light_speed = 2.99792458e8 # m / s
const GeVtokg = 1.60217e-10 # GeV -> kg m^2 / s^-2
const alpha = 1/137
const lambda = 3.86e-11 # cm
const Tcmb = 2.72548 # K
const k = 1.38064852e-23 # m^2 kg s^-2 K^-1

function sigma(s::Real)
  beta = (1 - 1/s)^0.5
  return 0.5 * pi * relectron^2 * (1-beta^2) * ((3-beta^4) * log((1+beta)/(1-beta)) - 2 * beta * (2-beta^2))
end

function phi_bar(s0::Real)
  if s0 <= 0
    return 0
  elseif s0 > 5
    return 2 * s0 * (log(4 * s0) - 2) + log(4 * s0) * (log(4 * s0) - 2) - (pi^2 - 9) / 3 + s0^(-1) * (log(4 * s0) + 9 / 8)
  end

  phi, err = quadgk(s->s * 2 * sigma(s) / pi / relectron^2, 1, s0, rtol=1e-8)
  return phi
end

function f(nu::Real)
  #if (nu > 2)
  #  return (pi * nu / 4)^0.5 * exp(-nu) * (1 + 75/8 * nu)
  #end

  result, err = quadgk(eps->(exp(eps)-1)^(-1) * phi_bar(eps/nu), nu, max(100, 10 * nu); rtol=1e-8)
  return nu^2 * result
end

function dtaudx(E::Real)
  return alpha^2 / pi / lambda * (k * Tcmb / melectron / light_speed / light_speed)^3 * f(melectron^2 * light_speed^4 / (E * GeVtokg) / k / Tcmb)
end

function dtaudx_isrf(E::Real)
  electron_mass = melectron * light_speed^2 / GeVtokg
  Ethred = electron_mass^2 / E

  n(eps) = 1

  inte, err = quadgk(eps->n(eps) * phi_bar(eps * E / electron_mass^2) / eps^2)

  return pi * relectron^2 * Ethred^2 * inte
end
