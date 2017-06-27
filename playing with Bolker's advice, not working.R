fm0  <- lmer(as.formula(paste(yy[i], "~ LUI + LUI*(", paste(top20.short, collapse="+"),"+Rest)","+(0+Yeart|Plot)+(1|Plot)+(1|Year_change)")), data= pchange.all2, REML=FALSE)
            

fmod <- glFormula(as.formula(paste(yy[i], "~ LUI + LUI*(", paste(top20.short, collapse="+"),"+Rest)","+(0+Yeart|Plot)+(1|Plot)+(1|Year_change)")), data= pchange.all2, family=gaussian,REML=FALSE)
## note we need family=gaussian() here -- an actual family object,
## not a string ("gaussian") or a family function (gaussian)
ldevfun <- do.call(mkGlmerDevfun, c(fmod, list(family=gaussian())))
source("glmmConstr.R")
fm1 <- glmmConstr(ldevfun,fmod,getCall(fm0),rep(1,2),rep(1,2), debug=TRUE)
