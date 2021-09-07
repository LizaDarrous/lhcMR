#' Title
#'
#' @param input.df
#' @param trait.names
#' @param log.file
#' @param run_ldsc
#' @param run_MR
#' @param hm3
#' @param ld
#'
#' @importFrom stats fft filter qnorm runif sd
#' @importFrom utils write.csv write.table
#' @importFrom dplyr slice
#' @return
# #' @export

gettingSP_ldscMR = function(input.df,trait.names,log.file,run_ldsc=TRUE,run_MR=TRUE,hm3,ld){

  EXP = trait.names[1]
  OUT = trait.names[2]
  nX = mean(input.df$N.x)  #Get sample size for trait X
  nY = mean(input.df$N.y)  #Get sample size for trait Y

  bX = input.df$TSTAT.x/sqrt(nX)   #Get standardised beta for trait X
  bY = input.df$TSTAT.y/sqrt(nY)   #Get standardised beta for trait Y

  X = select(input.df, RSID, CHR, A1, A2, BETA.x, SE.x, PVAL.x, TSTAT.x, N.x) %>% rename(unstdb = BETA.x, sderr = SE.x, pval = PVAL.x, TSTAT=TSTAT.x, N = N.x)
  Y = select(input.df, RSID, CHR, A1, A2, BETA.y, SE.y, PVAL.y, TSTAT.y, N.y) %>% rename(unstdb = BETA.y, sderr = SE.y, pval = PVAL.y, TSTAT=TSTAT.y, N = N.y)
  X$BETA = bX
  Y$BETA = bY
  X$SE = 1/sqrt(nX)
  Y$SE = 1/sqrt(nY)

  if(run_ldsc){
    # slice, every 10th SNP for faster computation
    X %>%
      slice(seq(1, nrow(X), by=10)) -> X_filtered
    #nrow(X_filtered)
    Y %>%
      slice(seq(1, nrow(Y), by=10)) -> Y_filtered
    #nrow(Y_filtered)

    write.table(X_filtered, file = paste0(EXP, "_GWAS.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names=TRUE)
    write.table(Y_filtered, file = paste0(OUT, "_GWAS.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names=TRUE)

    cat("Munging exposure and outcome data FOR LDSC: \n")
    invisible(utils::capture.output(GenomicSEM::munge( paste0(EXP, "_GWAS.txt"),
                                                       hm3,
                                                       trait.names=EXP)))
    invisible(utils::capture.output(GenomicSEM::munge( paste0(OUT, "_GWAS.txt"),
                                                       hm3,
                                                       trait.names=OUT)))

    cat("Please check the log files", paste0(EXP, "_munge.log and ", OUT, "_munge.log" ),  "to ensure that all columns were interpreted correctly and no warnings were issued for any of the summary statistics files\n")

    traits = c(paste0(EXP, ".sumstats.gz"), paste0(OUT, ".sumstats.gz"))
    sample.prev <- c(NA,NA) # continuous traits
    population.prev <- c(NA,NA) # continuous traits

    trait.names<-c(EXP, OUT)
    invisible(utils::capture.output(LDSCoutput <- GenomicSEM::ldsc(traits,
                                                                   sample.prev,
                                                                   population.prev,
                                                                   ld,
                                                                   trait.names)))
    #save(LDSCoutput, file="Pfactor.RData")

    cat("Please check the log file", paste0(c(traits, "ldsc.log"), collapse="_"),  "for in-depth results of the cross-trait LDSC analysis\n")

    # Clean up generated files, but save log file
    file.remove(paste0(EXP, c("_GWAS.txt", ".sumstats.gz")))  #need to silent the result
    file.remove(paste0(OUT, c("_GWAS.txt", ".sumstats.gz"))) #need to silent the result
    #i_X = as.numeric(LDSCoutput$I[1,1])
    #i_Y = as.numeric(LDSCoutput$I[2,2])
    i_XY = as.numeric(LDSCoutput$I[1,2])
    #h2_x_ldsc = as.numeric(LDSCoutput$S[1,1])
    #h2_y_ldsc = as.numeric(LDSCoutput$S[2,2])
  } else {i_XY = runif(1,-0.2,0.2)}


  ### Running standard MR
  if(run_MR){
    MR_output = paste0(EXP,"-",OUT,"_MRresults.csv")
    ## Get significant SNPs above certain Z-statistic corresponding to set p-value
    prune_X = function(zX,p_limit=1e-5){
      zX=zX
      z_limit=abs(qnorm(0.5*p_limit))
      ind_keep=which(abs(zX)>z_limit)
      ind_keep=unique(ind_keep)
      ind_keep=list(ind_keep)
      return(ind_keep)
    }

    ## Taken from Jonathan Sulc to create bins that fit a maximum of 50k SNPs (max threshold for clumping)
    snp_bin  =  function( snp_ranks,
                          chunk_size = 50000 ){
      if (nrow( snp_ranks ) == 0) {
        return()
      }

      max_chr  =  snp_ranks$chr %>%
        table %>%
        cumsum %>%
        (function(x) x < chunk_size) %>%
        (function(x) names(x)[ max(which(x)) ] ) %>%
        as.numeric
      if (is.na( max_chr )) {
        max_chr = min( snp_ranks$chr )
      }

      bin = snp_ranks %>%
        dplyr::filter( chr <= max_chr ) %>%
        list
      return( c( bin,
                 snp_bin( snp_ranks[ snp_ranks$chr > max_chr, ],
                          chunk_size ) ) )
    }

    # Set values
    pval=2*pnorm(-abs(5.45))
    pval1=2*pnorm(-abs(4))
    reverse_t_threshold  =  qnorm( 5e-2 )

    ### Forward
    mr_dataX = cbind.data.frame(SNP = X$RSID, beta = X$BETA, se = X$SE, effect_allele = X$A1, other_allele = X$A2, chr=X$CHR, Phenotype=EXP, tstat=X$TSTAT )
    mr_dataY = cbind.data.frame(SNP = Y$RSID, beta = Y$BETA, se = Y$SE, effect_allele = Y$A1, other_allele = Y$A2, chr=Y$CHR, Phenotype=OUT, tstat=Y$TSTAT )

    mr_ind=unlist(prune_X(mr_dataX$tstat,pval))
    print(length(mr_ind))
    if(length(mr_ind)==0){
      mr_ind=unlist(prune_X(mr_dataX$tstat,pval1))
      print(length(mr_ind))
    }
    mr_dataX = mr_dataX[mr_ind,]
    mr_dataY = mr_dataY[mr_ind,]

    # remove SNPs that are more strongly associated with the outcome than the exposure
    #filter( ( beta.exposure - beta.outcome ) / sqrt( se.exposure^2 + se.outcome^2 ) > reverse_t_threshold ) %>%
    ind_keep=which((abs(mr_dataX$beta)-abs(mr_dataY$beta))/sqrt(mr_dataX$se^2+mr_dataY$se^2) > reverse_t_threshold)
    print(length(ind_keep))
    mr_dataX = mr_dataX[ind_keep,]
    mr_dataY = mr_dataY[ind_keep,]

    exp_dat <- TwoSampleMR::format_data(mr_dataX, type="exposure")  ##same rows as mr_dataX
    clump_bin = snp_bin(mr_dataX,50000)

    #exp_dat2 <- clump_data(exp_dat)
    exp_data = c()
    for (x in 1:length(clump_bin)) {
      temp = exp_dat[exp_dat$SNP %in% clump_bin[[x]]$SNP,]
      temp1 = TwoSampleMR::clump_data(temp)
      exp_data=rbind(exp_data,temp1)
    }

    dups=which(duplicated(exp_data$SNP)==TRUE)
    if(length(dups)>0){
      exp_dat2 = exp_data[-dups,]
    }else{
      exp_dat2 = exp_data
    }

    out_dat <- TwoSampleMR::format_data(mr_dataY, type="outcome")
    out_dat2=out_dat[out_dat$SNP %in% exp_dat2$SNP,]

    exp_dat2=exp_dat2[order(exp_dat2$SNP),]
    out_dat2=out_dat2[order(out_dat2$SNP),]

    if(all(exp_dat2$`effect_allele.exposure` == out_dat2$`effect_allele.outcome`)){
      print("action=1")
      action = 1
    } else {
      print("action=2/3")
      aligned = which(exp_dat2$`effect_allele.exposure` == out_dat2$`effect_allele.outcome` &
                        exp_dat2$`other_allele.exposure` == out_dat2$`other_allele.outcome`)
      swapped = which(exp_dat2$`effect_allele.exposure` == out_dat2$`other_allele.outcome` &
                        exp_dat2$`other_allele.exposure` == out_dat2$`effect_allele.outcome`)
      exp_dat2[swapped,'beta.exposure']=exp_dat2[swapped,'beta.exposure']*-1
      exp_dat2 = exp_dat2[c(aligned,swapped),]
      out_dat2 = out_dat2[c(aligned,swapped),]
      action = 1  ## made sure all strands are okay
    }
    # Harmonise the exposure and outcome data
    dat <- TwoSampleMR::harmonise_data(
      exposure_dat = exp_dat2,
      outcome_dat = out_dat2, action = action
    )

    #Sensitivity - Q-test
    het <- TwoSampleMR::mr_heterogeneity(dat)
    het$I2 = ((het$Q-het$Q_df)/het$Q)*100
    plei <- TwoSampleMR::mr_pleiotropy_test(dat)

    smaller=FALSE
    tryCatch( {res1 <- TwoSampleMR::mr(dat); print("Bigger MR list") }, error = function(e) {smaller<<-TRUE})
    if(smaller){
      print("Smaller MR list")
      res1 <- TwoSampleMR::mr(dat, method_list=c("mr_egger_regression","mr_weighted_median","mr_ivw"))
    }else{
      res1 <- TwoSampleMR::mr(dat)
    }

    if(nrow(res1)<2){
      axy_MR = res1[,'b']
    }else{
      axy_MR = res1[which(res1$method=="Inverse variance weighted"),'b']
    }
    write.table(as.data.frame(res1), file=MR_output, sep = ",", append = TRUE, row.names = FALSE)
    write.table(as.data.frame(het), file=MR_output, sep = ",", append = TRUE, row.names = FALSE)
    write.table(as.data.frame(plei), file=MR_output, sep = ",", append = TRUE, row.names = FALSE)
    write.table("*", MR_output, sep = ",", quote = FALSE, col.names = FALSE, row.names = FALSE, append = TRUE)

    ## reverse MR
    rm(mr_dataX,mr_dataY,exp_dat,exp_data,exp_dat2,out_dat,out_dat2,dups,ind_keep,mr_ind,clump_bin, action, temp, temp1,dat,het,plei,smaller)

    #reverse the exposure and outcome to Y - X, nothing else besides this needs to change
    mr_dataX = cbind.data.frame(SNP = Y$RSID, beta = Y$BETA, se = Y$SE, effect_allele = Y$A1, other_allele = Y$A2, chr=Y$CHR, Phenotype=OUT, tstat = Y$TSTAT )
    mr_dataY = cbind.data.frame(SNP = X$RSID, beta = X$BETA, se = X$SE, effect_allele = X$A1, other_allele = X$A2, chr=X$CHR, Phenotype=EXP, tstat = X$TSTAT )

    mr_ind=unlist(prune_X(mr_dataX$tstat,pval))
    print(length(mr_ind))
    if(length(mr_ind)==0){
      mr_ind=unlist(prune_X(mr_dataY$tstat,pval1))
      print(length(mr_ind))
    }
    mr_dataX = mr_dataX[mr_ind,]
    mr_dataY = mr_dataY[mr_ind,]

    #filter( ( beta.exposure - beta.outcome ) / sqrt( se.exposure^2 + se.outcome^2 ) > reverse_t_threshold ) %>%  ind_keep=which((mr_dataX$beta-mr_dataY$beta)/sqrt(mr_dataX$se^2+mr_dataY$se^2) > reverse_t_threshold)
    #ind_keep=which((mr_dataX$beta-mr_dataY$beta)/sqrt(mr_dataX$se^2+mr_dataY$se^2) > reverse_t_threshold)
    ind_keep=which((abs(mr_dataX$beta)-abs(mr_dataY$beta))/sqrt(mr_dataX$se^2+mr_dataY$se^2) > reverse_t_threshold)
    print(length(ind_keep))
    mr_dataX = mr_dataX[ind_keep,]
    mr_dataY = mr_dataY[ind_keep,]

    exp_dat <- TwoSampleMR::format_data(mr_dataX, type="exposure")  ##same rows as mr_dataX
    clump_bin = snp_bin(mr_dataX,50000)

    #exp_dat2 <- clump_data(exp_dat)
    exp_data = c()
    for (x in 1:length(clump_bin)) {
      temp = exp_dat[exp_dat$SNP %in% clump_bin[[x]]$SNP,]
      temp1 = TwoSampleMR::clump_data(temp)
      exp_data=rbind(exp_data,temp1)
    }

    dups=which(duplicated(exp_data$SNP)==TRUE)
    if(length(dups)>0){
      exp_dat2 = exp_data[-dups,]
    }else{
      exp_dat2 = exp_data
    }

    out_dat <- TwoSampleMR::format_data(mr_dataY, type="outcome")
    out_dat2=out_dat[out_dat$SNP %in% exp_dat2$SNP,]

    exp_dat2=exp_dat2[order(exp_dat2$SNP),]
    out_dat2=out_dat2[order(out_dat2$SNP),]

    if(all(exp_dat2$`effect_allele.exposure` == out_dat2$`effect_allele.outcome`)){
      print("action=1")
      action = 1
    } else {
      print("action=2/3")
      aligned = which(exp_dat2$`effect_allele.exposure` == out_dat2$`effect_allele.outcome` &
                        exp_dat2$`other_allele.exposure` == out_dat2$`other_allele.outcome`)
      swapped = which(exp_dat2$`effect_allele.exposure` == out_dat2$`other_allele.outcome` &
                        exp_dat2$`other_allele.exposure` == out_dat2$`effect_allele.outcome`)
      exp_dat2[swapped,'beta.exposure']=exp_dat2[swapped,'beta.exposure']*-1
      exp_dat2 = exp_dat2[c(aligned,swapped),]
      out_dat2 = out_dat2[c(aligned,swapped),]
      action = 1  ## made sure all strands are okay
    }

    # Harmonise the exposure and outcome data
    dat <- TwoSampleMR::harmonise_data(
      exposure_dat = exp_dat2,
      outcome_dat = out_dat2, action = action
    )

    #Sensitivity - Q-test
    het <- TwoSampleMR::mr_heterogeneity(dat)
    het$I2 = ((het$Q-het$Q_df)/het$Q)*100
    plei <- TwoSampleMR::mr_pleiotropy_test(dat)

    smaller=FALSE
    tryCatch( {res2 <- TwoSampleMR::mr(dat); print("Bigger MR list") }, error = function(e) {smaller<<-TRUE})
    if(smaller){
      print("Smaller MR list")
      res2 <- TwoSampleMR::mr(dat, method_list=c("mr_egger_regression","mr_weighted_median","mr_ivw"))
    }else{
      res2 <- TwoSampleMR::mr(dat)
    }

    if(nrow(res2)<2){
      ayx_MR = res2[,'b']
    }else{
      ayx_MR = res2[which(res2$method=="Inverse variance weighted"),'b']
    }

    write.table(as.data.frame(res2), file=MR_output, sep = ",", append = TRUE, row.names = FALSE)
    write.table(as.data.frame(het), file=MR_output, sep = ",", append = TRUE, row.names = FALSE)
    write.table(as.data.frame(plei), file=MR_output, sep = ",", append = TRUE, row.names = FALSE)

  } else{axy_MR = runif(1,-0.5,0.5);ayx_MR = runif(1,-0.5,0.5);}

  return(list(i_XY, axy_MR, ayx_MR))

}
