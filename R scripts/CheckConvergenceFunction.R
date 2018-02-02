chkconv <- function(model) {
  relgrad <- with(model@optinfo$derivs,solve(Hessian,gradient))
  print("If convergence <0.001, model is acceptable.")
  max(abs(relgrad))
}

