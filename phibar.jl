using QuadGK

function phibar_compare(s0::Real)
  beta = (1-1/s0)^0.5

  w0 = (1 + beta) / (1 - beta)
  L, err = quadgk(w->log(w+1)/w, 1, w0, rtol=1e-8)

  result = 0
  result += (1 + beta^2) / (1 - beta^2) * log(w0)
  result += - beta^2 * log(w0)
  result += - (log(w0))^2
  result += 4 * log(w0) * log(w0 + 1)
  result += 2 * beta - 4 * beta / (1-beta^2)
  result += - 4 * L

  return result
end
