slurm_singleTrait_likelihood = function(par){
  theta=unlist(par)
  test = optim(theta, singleTrait_likelihood,
               betX=betX, pi1=pi1, sig1=sig1, weights=weights,
               m0=m0, nX=nX, bn=2^7, bins=10,
               method = "Nelder-Mead",
               control = list(maxit = 5e3))

  test.res=c(test$value,test$par,test$convergence)
  names(test.res)=c("mLL","piX", "h2X","iX","conv")
  return(test.res)
}
