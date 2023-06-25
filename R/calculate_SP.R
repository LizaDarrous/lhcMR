#' Calculate starting points to be used in the likelihood function optimisation
#'
#' @param input.df The resulting data frame from merge_sumstats(), where the effect size, SE, RSID and other columns are present, in addition to columns representing LD scores, weights and local LD structure
#' @param trait.names Vector containing the trait names in the order they were used in merge_sumstats(): Exposure, Outcome
#' @param log.file
#' @param run_ldsc Boolean. Whether GenomicSEM::ldsc should be run to obtain the cross trait-intercept (i_XY). If FALSE, a random value will be generated. Default value = TRUE
#' @param run_MR Boolean. Whether TwoSampleMR::mr should be run to obtain the bidirectional causal effects (axy_MR, ayx_MR). If FALSE, random values will be generated. Default value = TRUE
#' @param saveRFiles Boolean, whether to write the results of GenomicSEM::ldsc,TwoSampleMR::mr, and the single trait analysis of LHC-MR (returns trait intercept and polygenicity) Default value = TRUE
#' @param hm3 Path to the input file (HAPMAP3 SNPs) required by GenomicSEM::ldsc
#' @param ld Path to the input file (LD scores) required by GenomicSEM::ldsc
#' @param nStep Can take two numerical values: 1 or 2. Represents the number of steps the lhcMR analysis will undertake. One single step estimates all 9 parameters simultaneously while fixing only
#' the traits' intercepts iX and iY, while two steps estimates 7 parameters after having estimated traits' intercepts and polygenicity (iX, piX, iY, piY) from the single trait analysis and fixed
#' their values in the likelihood optimisation and parameter estimation
#' @param SP_single Numerical value indicating how many starting points should the single trait analysis use in the likelihood optimisation. Best to range between 3-5, default value = 3
#' @param SP_pair Numerical value indicating how many starting points should the pair trait analysis use in the likelihood optimisation. Best to range between 50-100, default value = 50
#' @param SNP_filter Numerical value indicating the filtering of every nth SNP to reduce large datasets and speed up analysis. Default value = 10
#' @param SNP_filter_ldsc Numerical value indicating the filtering of every nth SNP to reduce large datasets and speed up the LDSC analysis. Set to 1 if no filtering is needed, otherwise default = 10
#' @param nCores Numerical value indicating number of cores to be used in 'mclapply' to parallelise the analysis. If set to NA, then it will be calculated as 2/3 of the available cores, default value = 1 to avoid parallelisation
#' @param M Numerical value indicating the number of SNPs used to calculate the LD reported in the LD file (for genotyped SNPs). Default value = 1e7

#'
#' @importFrom dplyr slice select rename
#' @importFrom stats optim
#' @importFrom parallel mclapply detectCores
#' @return Returns a list containing the filtered dataset (by every `SNP_filter`th SNP), the starting points to be used in the pair trait optimisation, the traits' intercepts,
#' the traits' polygenicity if nStep = 2, as well as some extra parameters like the cross-trait intercept and bidirectional causal effect estimated by IVW
#' @export
#'
#' @examples
calculate_SP <- function(input.df,trait.names,run_ldsc=TRUE,run_MR=TRUE,saveRFiles=TRUE,hm3=NA,ld=NA,nStep=2,
                         SP_single=3,SP_pair=50,SNP_filter=10,SNP_filter_ldsc=NA,nCores=1,M=1e7, b_file=NULL, plink_path=NULL){

  if(nStep>2 || nStep<1){
    cat(print("Please choose 1 or 2 for the number of analysis steps\n"))
    stop()
  }

  if(is.na(nCores)){
    nCores = max(1,floor((parallel::detectCores())/3))
  } else {
    if(nCores > parallel::detectCores()){
      cat(print("Core number chosen is greater than cores available\n"))
      stop()
    }
  }

  if(is.na(SNP_filter_ldsc)){SNP_filter_ldsc = SNP_filter}
  if(SNP_filter_ldsc<1 | SNP_filter <1){
    cat(print("Please choose a value equal to or greater than 1 for the thinning of every nth SNP\n"))
    stop()
  }

  EXP = trait.names[1]
  OUT = trait.names[2]

  # Get estimates of starting points from LDSC and standard MR method (IVW)
  SP = gettingSP_ldscMR(input.df,trait.names,run_ldsc,run_MR,saveRFiles,SNP_filter_ldsc,hm3,ld, b_file = b_file, plink_path = plink_path)
  i_XY = as.numeric(SP[[1]])
  axy_MR = as.numeric(SP[[2]])
  ayx_MR = as.numeric(SP[[3]])

  # For LHC single step analysis we reduce the number of SNPs for faster computation
  input.df %>%
    slice(seq(1, nrow(input.df), by=SNP_filter)) -> input.df_filtered
  #nrow(input.df_filtered)

  nX = mean(input.df_filtered$N.x)  #get sample size for trait X
  nY = mean(input.df_filtered$N.y)  #get sample size for trait Y
  m0 = mean(input.df_filtered$M_LOCAL) #get number of SNPs used to calculate the local LD

  bX = input.df_filtered$TSTAT.x/sqrt(nX)   #get standardised beta for trait X
  bY = input.df_filtered$TSTAT.y/sqrt(nY)   #get standardised beta for trait Y

  ld = input.df_filtered$LDSC
  w8s = input.df_filtered$WEIGHT
  pi1 = input.df_filtered$PIK
  sig1 = input.df_filtered$SIGK

  sp_piX = runif(SP_single,0,0.01)
  sp_h2X = runif(SP_single,0,0.5)
  sp_iX = runif(SP_single,0.5,1.5)

  # Get piX, iX (and total_h2) in a separate step - single trait analysis for X and Y
    para=cbind(sp_piX,sp_h2X,sp_iX)
    sp_mat=matrix(unlist(para), ncol=3, byrow = FALSE)
    colnames(sp_mat)=colnames(para)
    par.df = data.frame(par=I(apply(sp_mat,1,as.list))) #generate a dataframe of lists for each row of parameters - input for rslurm/lapply

    #test.exp = matrix(NA,nrow=SP_single,ncol=5)
    #test.out = matrix(NA,nrow=SP_single,ncol=5)

    test.exp <- parallel::mclapply(par.df[[1]], function(x) {
      theta = unlist(x)
      test1 = optim(theta, singleTrait_likelihood,
                    betX=bX, pi1=pi1, sig1=sig1, w8s=w8s, M=M,
                    m0=m0, nX=nX, bn=2^7, bins=10,
                    method = "Nelder-Mead",
                    control = list(maxit = 5e3))

     list("mLL"=test1$value,"par"=test1$par,"conv"=test1$convergence)
      }, mc.cores = nCores)

    test.exp = as.data.frame(t(matrix(unlist(test.exp), nrow=length(unlist(test.exp[1])))))

    test.out <- parallel::mclapply(par.df[[1]], function(x) {
      theta = unlist(x)
      test2 = optim(theta, singleTrait_likelihood,
                    betX=bY, pi1=pi1, sig1=sig1, w8s=w8s, M=M,
                    m0=m0, nX=nY, bn=2^7, bins=10,
                    method = "Nelder-Mead",
                    control = list(maxit = 5e3))

      list("mLL"=test2$value,"par"=test2$par,"conv"=test2$convergence)
    }, mc.cores = nCores)

    test.out = as.data.frame(t(matrix(unlist(test.out), nrow=length(unlist(test.out[1])))))

    colnames(test.exp)=c("mLL","piX", "h2X","iX","conv")
    colnames(test.out)=c("mLL","piX", "h2X","iX","conv")
    res_exp_min = test.exp[which(test.exp$mLL == min(test.exp$mLL)), ]
    res_exp = abs(res_exp_min[2:4])
    res_out_min = test.out[which(test.out$mLL == min(test.out$mLL)), ]
    res_out = abs(res_out_min[2:4])

    if(saveRFiles){
      res_singleTrait = rbind(res_exp_min, res_out_min)
      res_singleTrait = cbind("Trait"=c(EXP,OUT), res_singleTrait)
      write.csv(res_singleTrait, paste0("SingleTraitAnalysis_",EXP,"-",OUT,".csv"), row.names = F)
    }

  # Get fixed points due to single trait analysis estimate
  pi_X = as.numeric(res_exp[1])
  pi_Y = as.numeric(res_out[1])
  h2_x = as.numeric(res_exp[2]) #total heritability, not used
  h2_y = as.numeric(res_out[2]) #total heritability, not used
  i_X = as.numeric(res_exp[3])
  i_Y = as.numeric(res_out[3])


  # Generate the rest of the starting points depending on TwoStep vs SingleStep
  if(nStep==2){
    sp_tX = runif(SP_pair,0,0.5)
    sp_tY = runif(SP_pair,-0.5,0.5)
    sp_h2X = max(0,h2_x-(sp_tX^2))
    sp_h2Y = max(0,h2_y-(sp_tY^2))
    sp_axy = replicate(SP_pair, (axy_MR+runif(1,-0.1,0.1)))
    sp_ayx = replicate(SP_pair, (ayx_MR+runif(1,-0.1,0.1)))
    sp_iXY = replicate(SP_pair, (i_XY+runif(1,-0.05,0.05))) #rep(i_XY,SP_pair)

    para=cbind(sp_h2X,sp_h2Y,sp_tX,sp_tY,sp_axy,sp_ayx,sp_iXY)
    sp_mat1=matrix(unlist(para), ncol=7, byrow = FALSE)
    colnames(sp_mat1)=c("sp_h2X","sp_h2Y","sp_tX","sp_tY","sp_axy","sp_ayx","sp_iXY")
    #sp_mat1 = cbind("SP"=c(1:nrow(sp_mat1)),sp_mat1)
    if(saveRFiles){write.csv(sp_mat1,"StartingPoints.csv", row.names=F)}
    return(list("iX"=i_X,"iY"=i_Y,"piX"=pi_X,"piY"=pi_Y,"i_XY"=i_XY,"axy_MR"=axy_MR,"ayx_MR"=ayx_MR,
                "input.df_filtered"=input.df_filtered,"sp_mat"=sp_mat1))
  }

  if(nStep==1){
    sp_piX = runif(SP_pair,0,1e-4) #rep(pi_X,SP_pair)
    sp_piY = runif(SP_pair,0,1e-4) #rep(pi_Y,SP_pair)
    sp_tX = runif(SP_pair,0,0.5)
    sp_tY = runif(SP_pair,-0.5,0.5)
    sp_h2X = max(0,h2_x-(sp_tX^2))
    sp_h2Y = max(0,h2_y-(sp_tY^2))
    sp_axy = replicate(SP_pair, (axy_MR+runif(1,-0.1,0.1)))
    sp_ayx = replicate(SP_pair, (ayx_MR+runif(1,-0.1,0.1)))
    sp_iXY = replicate(SP_pair, (i_XY+runif(1,-0.05,0.05)))

    para=cbind(sp_piX,sp_piY,sp_h2X,sp_h2Y,sp_tX,sp_tY,sp_axy,sp_ayx,sp_iXY)
    sp_mat1=matrix(unlist(para), ncol=9, byrow = FALSE)
    colnames(sp_mat1)=c("sp_piX","sp_piY","sp_h2X","sp_h2Y","sp_tX","sp_tY","sp_axy","sp_ayx","sp_iXY")
    if(saveRFiles){write.csv(sp_mat1,"StartingPoints.csv", row.names=F)}
    return(list("iX"=i_X,"iY"=i_Y,"i_XY"=i_XY,"axy_MR"=axy_MR,"ayx_MR"=ayx_MR,"input.df_filtered"=input.df_filtered,"sp_mat"=sp_mat1))

  }

}
