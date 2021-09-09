#' Title
#'
#' @param par pp
#' @param start_ind
#' @param end_ind
#'
#' @return
#' @export
#'
#' @examples
slurm_blockJK_twoStep_likelihood = function(par,start_ind, end_ind){
  theta=unlist(par)

  test = optim(theta, pairTrait_twoStep_likelihood,
               betXY=betXY[-(start_ind:end_ind),], pi1=pi1[-(start_ind:end_ind)], sig1=sig1[-(start_ind:end_ind)],
               w8s=w8s[-(start_ind:end_ind)], pi_U=piU,
               pi_X=piX, pi_Y=piY, i_X=iX, i_Y=iY,
               m0=m0, M=M, nX=nX, nY=nY, bn=bn, bins=bins, model=param,
               method = "Nelder-Mead",
               control = list(maxit = 5e3,
                              parscale = parscale2))

  test.res=c(test$value,test$par,test$convergence,start_ind, end_ind)
  cnames_res = c("mLL","h2X","h2Y","tX","tY","axy","ayx","iXY","conv","start_ind", "end_ind")

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
