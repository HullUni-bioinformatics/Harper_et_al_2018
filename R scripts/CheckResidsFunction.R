chkres <- function(model) {
  require(RVAideMemoire)
  sresid <- resid(model, type = "deviance")
  hist(sresid)
  fitted.glmm <- fitted(model, level=1)        # Extract the fitted (predicted) values
  plot(sresid ~ fitted.glmm)                   # Check for homoscedasticity
  plot(model)                                  # gives a heteroscedasticity plot #fans out, not good
  require(arm)
  binnedplot(fitted(model),resid(model))
}




chkres.PQL <- function(model) {
  require(RVAideMemoire)
  plot(model)
  presid <- resid(model, type = "pearson")
  hist(presid)
  fitted.glmm <- fitted(model, level=1)        # Extract the fitted (predicted) values
  plot(presid ~ fitted.glmm)                   # Check for homoscedasticity
  
}


chkres.zi <- function(model) {
  require(RVAideMemoire)
  presid <- residuals(model)
  hist(presid)
  fitted.glmm <- fitted(model, level=1)        # Extract the fitted (predicted) values
  plot(presid ~ fitted.glmm)                   # Check for homoscedasticity
  require(arm)
  binnedplot(fitted.glmm,presid)
}
