flury.phi <- function(datamat)
{
  # Calculates phi measure of goodness of diagonalisation as proposed by Flury (1988)
  # Interpretation: phi smaller value --> better diagonalisation
  
  # datamat: square symmetric matrix
  
  return(det(diag(diag(datamat))) / det(datamat))
}
