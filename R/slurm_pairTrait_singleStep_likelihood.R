#' Title
#'
#' @param par
#'
#' @return
#' @export
#'
#' @examples
slurm_pairTrait_singleStep_likelihood = function(par){
  theta=unlist(par)

  test = optim(theta, pairTrait_singleStep_likelihood,
               betXY=betXY, pi1=pi1, sig1=sig1, weights=weights, pi_U=pi_U,
               i_X=i_X, i_Y=i_Y,
               m0=m0, nX=nX, nY=nY, bn=bn, bins=bins, model=param,
               method = "Nelder-Mead",
               control = list(maxit = 5e3,
                              parscale = parscale))

  test.res=c(test$value,test$par,test$convergence)
  cnames_res = c("mLL","piX","piY","h2X","h2Y","tX","tY","alp","bet","iXY","conv")

  if(param=="comp"){
    names(test.res)=cnames_res
  }else if(param=="U"){
    param_temp = c("tX","tY")
    names(test.res)=cnames_res[!(cnames_res %in% param_temp)]
  }else{
    names(test.res)=cnames_res[!(cnames_res %in% param)]
  }

  return(list(res = test.res))
}
