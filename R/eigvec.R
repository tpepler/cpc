eigvec <- function(U)                        
{
  if(U[1, 2] != 0){ratio <- (U[2, 2] - U[1, 1]) / U[1, 2]}
  else{ratio = 0} #nie seker of hierdie reg is nie?
  discr <- sqrt(ratio^2 + 4)
  T1 <- (ratio + discr) / 2
  T2 <- (ratio - discr) / 2
  if(abs(T1) > abs(T2)){T.waarde <- T2}
  else{T.waarde <- T1}
  C.waarde <- 1 / sqrt(1 + T.waarde^2)
  S <- T.waarde * C.waarde
  J <- matrix(NA, nrow = 2, ncol = 2)
  J[1, 1] <- C.waarde
  J[2, 1] <- S
  J[1, 2] <- -S
  J[2, 2] <- C.waarde
  return(J)
}
