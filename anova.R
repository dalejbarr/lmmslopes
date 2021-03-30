fitanova <- function(xd) {
  xd$C <- factor(xd$A)
  xd$S <- factor(xd$subj_id)
  xd$I <- factor(xd$item_id)
  F1.all <- summary(aov(Y ~ A + Error(subj_id/A), data=xd))[[2]][[1]]
  F1 <- F1.all$`F value`[1]
  F1.df.denom <- F1.all$`Df`[2]
  p1 <- F1.all$`Pr(>F)`[1]
  F2.all <- summary(aov(Y ~ A + Error(item_id/A), data=xd))[[2]][[1]]
  F2 <- F2.all$`F value`[1]
  F2.df.denom <- F2.all$`Df`[2]
  p2 <- F2.all$`Pr(>F)`[1]
  minF <- (F1*F2) / (F1 + F2)
  minF.df.denom <- ( (F1+F2)^2 ) / ( (F1^2)/F2.df.denom + (F2^2)/F1.df.denom ) 
  v1 <- c(F1=F1, F2=F2, minF=minF,
          p1=p1, p2=p2,  pmf=pf(minF, 1, minF.df.denom, lower.tail=FALSE),
          pmax=max(c(p1,p2)))
  return(v1)
}
