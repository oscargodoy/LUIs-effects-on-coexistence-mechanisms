glmmConstr <- function(devfun,mod,mc,beta.lwr,beta.upr,debug=FALSE) {
  opt0 <- minqa::bobyqa(fn=devfun,par=1,lower=0,upper=Inf)
  if (debug) cat("initial theta: ",opt0$par,"\n")
  opt1 <- optimizeGlmer(devfun)
  if (debug) cat("next theta: ",opt1$par,"\n")
  rho <- environment(devfun)
  nbeta <- ncol(rho$pp$X)
  theta <- rho$pp$theta
  theta.lwr <- rho$lower  ## this changes after update!
  devfun <- updateGlmerDevfun(devfun, mod$reTrms)
  opt2 <- nloptwrap(par=c(theta,rep(0,nbeta)),
                    fn=devfun,
                    lower=c(theta.lwr,beta.lwr),
                    upper=c(rep(Inf,length(theta.lwr)),beta.upr))
  if (debug) cat("final theta: ",opt2$par,"\n")
  mkMerMod(rho, opt2, mod$reTrms, mod$fr, mc=mc)
}