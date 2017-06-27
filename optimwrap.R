optimwrap <- function(fn,par,lower,upper,control=list(),
                      ...) {
  if (is.null(control$method)) stop("must specify method in optCtrl")
  method <- control$method
  control$method <- NULL
  ## "L-BFGS-B" requires finite upper and lower values, set the lower values to competition
  if (method=="L-BFGS-B")
  par <- rep(1, times=3)#fixef(mlist[[i]]) # use starting values from previous models to avoid overfitting
  upper <- c(Inf, 0.0001, Inf)
  control=list(maxit=1000000)
  res <- optim(par=par, fn=fn, lower=lower,upper=upper,
               control=control,method=method,...)
  with(res, list(par = par,
                 fval = value,
                 feval= counts[1],
                 conv = convergence,
                 message = message))
}