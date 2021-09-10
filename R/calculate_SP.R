#' Title
#'
#' @param input.df
#' @param trait.names
#' @param log.file
#' @param run_ldsc
#' @param run_MR
#' @param hm3
#' @param ld
#' @param jobs
#'
#' @importFrom dplyr slice select rename
#' @importFrom stats optim
#' @importFrom parallel mclapply
#' @return
#' @export
#'
#' @examples
calculate_SP <- function(input.df,trait.names,log.file=NA,run_ldsc=TRUE,run_MR=TRUE,hm3,ld,StepNum=2,
                         SP_single=30,SP_pair=100,SNP_filter=10){

  if(StepNum>2 || StepNum<1){
    cat(print("Please choose 1 or 2 for the number of analysis steps"))
    stop()
  }

  EXP = trait.names[1]
  OUT = trait.names[2]

  ## get estimates of starting points from LDSC and standard MR method (IVW)
  SP = gettingSP_ldscMR(input.df,trait.names,log.file,run_ldsc,run_MR,hm3,ld)
  i_XY = as.numeric(SP[[1]])
  axy_MR = as.numeric(SP[[2]])
  ayx_MR = as.numeric(SP[[3]])

  ## for LHC single step analysis we reduce the number of SNPs for faster computation
  input.df %>%
    slice(seq(1, nrow(input.df), by=SNP_filter)) -> input.df_filtered
  nrow(input.df_filtered)
  # 468,993

  nX = mean(input.df_filtered$N.x)  #Get sample size for trait X
  nY = mean(input.df_filtered$N.y)  #Get sample size for trait Y
  m0 = mean(input.df_filtered$M_LOCAL) ##Get number of SNPs used to calculate the local LD

  bX = input.df_filtered$TSTAT.x/sqrt(nX)   #Get standardised beta for trait X
  bY = input.df_filtered$TSTAT.y/sqrt(nY)   #Get standardised beta for trait Y

  ld = input.df_filtered$LDSC
  w8s = input.df_filtered$WEIGHT
  pi1 = input.df_filtered$PIK
  sig1 = input.df_filtered$SIGK

  sp_piX = runif(SP_single,0,0.01)
  sp_h2X = runif(SP_single,0,0.5)
  sp_iX = runif(SP_single,0.5,1.5)

  ## Get piX, iX (and total_h2) in a separate step - single trait analysis for X and Y
    para=cbind(sp_piX,sp_h2X,sp_iX)
    sp_mat=matrix(unlist(para), ncol=3, byrow = FALSE)
    colnames(sp_mat)=colnames(para)
    par.df = data.frame(par=I(apply(sp_mat,1,as.list))) #generate a dataframe of lists for each row of parameters - input for rslurm/lapply

    test.exp = matrix(NA,nrow=SP_single,ncol=5)
    test.out = matrix(NA,nrow=SP_single,ncol=5)

    test.exp <- parallel::mclapply(par.df[[1]], function(x) {
      theta = unlist(x)
      test1 = optim(theta, singleTrait_likelihood,
                    betX=bX, pi1=pi1, sig1=sig1, w8s=w8s,
                    m0=m0, nX=nX, bn=2^7, bins=10,
                    method = "Nelder-Mead",
                    control = list(maxit = 5e3))

     list("mLL"=test1$value,"par"=test1$par,"conv"=test1$convergence)
      })

    test.exp = as.data.frame(t(matrix(unlist(test.exp), nrow=length(unlist(test.exp[1])))))

    test.out <- parallel::mclapply(par.df[[1]], function(x) {
      theta = unlist(x)
      test1 = optim(theta, singleTrait_likelihood,
                    betX=bY, pi1=pi1, sig1=sig1, w8s=w8s,
                    m0=m0, nX=nX, bn=2^7, bins=10,
                    method = "Nelder-Mead",
                    control = list(maxit = 5e3))

      list("mLL"=test1$value,"par"=test1$par,"conv"=test1$convergence)
    })

    test.out = as.data.frame(t(matrix(unlist(test.out), nrow=length(unlist(test.out[1])))))

    colnames(test.exp)=c("mLL","piX", "h2X","iX","conv")
    colnames(test.out)=c("mLL","piX", "h2X","iX","conv")
    res_exp_min = test.exp[which(test.exp$mLL == min(test.exp$mLL)), ]
    res_exp = abs(res_exp_min[2:4])
    res_out_min = test.out[which(test.out$mLL == min(test.out$mLL)), ]
    res_out = abs(res_out_min[2:4])

  ## get fixed points due to single trait analysis estimate
  pi_X = as.numeric(res_exp[1])
  pi_Y = as.numeric(res_out[1])
  h2_x = as.numeric(res_exp[2]) #total heritability, not used
  h2_y = as.numeric(res_out[2]) #total heritability, not used
  i_X = as.numeric(res_exp[3])
  i_Y = as.numeric(res_out[3])


  ## generate the rest of the starting points depending on TwoStep vs SingleStep
  if(StepNum==2){
    sp_tX = runif(SP_pair,0,0.5)
    sp_tY = runif(SP_pair,-0.5,0.5)
    sp_h2X = max(0,h2_x-(sp_tX^2))
    sp_h2Y = max(0,h2_y-(sp_tY^2))
    sp_axy = replicate(SP_pair, (axy_MR+runif(1,-0.1,0.1)))
    sp_ayx = replicate(SP_pair, (ayx_MR+runif(1,-0.1,0.1)))
    sp_iXY = rep(i_XY,SP_pair)

    para=cbind(sp_h2X,sp_h2Y,sp_tX,sp_tY,sp_axy,sp_ayx,sp_iXY)
    sp_mat1=matrix(unlist(para), ncol=7, byrow = FALSE)
    colnames(sp_mat1)=c("sp_h2X","sp_h2Y","sp_tX","sp_tY","sp_axy","sp_ayx","sp_iXY")
    #sp_mat1 = cbind("SP"=c(1:nrow(sp_mat1)),sp_mat1)
    write.csv(sp_mat1,"StartingPoints.csv", row.names=F)
    return(list("iX"=i_X,"iY"=i_Y,"piX"=pi_X,"piY"=pi_Y,"input.df_filtered"=input.df_filtered,"sp_mat"=sp_mat1))
  }

  if(StepNum==1){
    sp_piX = runif(SP_pair,0,1e-4) #rep(pi_X,SP_pair)
    sp_piY = runif(SP_pair,0,1e-4) #rep(pi_Y,SP_pair)
    sp_tX = runif(SP_pair,0,0.5)
    sp_tY = runif(SP_pair,-0.5,0.5)
    sp_h2X = h2_x-(sp_tX^2)
    sp_h2Y = h2_y-(sp_tY^2)
    sp_axy = replicate(SP_pair, (axy_MR+runif(1,-0.1,0.1)))
    sp_ayx = replicate(SP_pair, (ayx_MR+runif(1,-0.1,0.1)))
    sp_iXY = rep(i_XY,SP_pair)

    para=cbind(sp_piX,sp_piY,sp_h2X,sp_h2Y,sp_tX,sp_tY,sp_axy,sp_ayx,sp_iXY)
    sp_mat1=matrix(unlist(para), ncol=9, byrow = FALSE)
    colnames(sp_mat1)=c("sp_piX","sp_piY","sp_h2X","sp_h2Y","sp_tX","sp_tY","sp_axy","sp_ayx","sp_iXY")
    write.csv(sp_mat1,"StartingPoints.csv", row.names=F) #SP_pair
    return(list("iX"=i_X,"iY"=i_Y,"input.df_filtered"=input.df_filtered,"sp_mat"=sp_mat1))

  }

}
