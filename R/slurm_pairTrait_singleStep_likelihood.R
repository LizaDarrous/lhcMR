#' Title
#'
#' @param par
#'
#' @return
#' NOT EXPORTED @export
#' @keywords internal
#'
#' @examples
slurm_pairTrait_singleStep_likelihood = function(par){
  theta=unlist(par)

  test = optim(theta, pairTrait_singleStep_likelihood,
               betXY=betXY, pi1=pi1, sig1=sig1, w8s=w8s, pi_U=piU,
               i_X=iX, i_Y=iY,
               m0=m0, M=M, nX=nX, nY=nY, bn=bn, bins=bins, model=param,
               method = "Nelder-Mead",
               control = list(maxit = 5e3,
                              parscale = parscale1))

  test.res=c(test$value,test$par,test$convergence)
  cnames_res = c("mLL","piX","piY","h2X","h2Y","tX","tY","axy","ayx","iXY","conv")

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
