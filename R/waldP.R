# wald.test return p value
waldP <- function(X,Y){
  df <- lm(Y~X)
  P <- wald.test(Sigma = vcov(df),b = coef(df),Terms = 2)$result$chi2['P']
  return(c(P))
}
