#' Title - Read the columns from the summary statistics files and edit them as required
#'
#' @param input.files  - list of data frames, where each data frame contains the summary statistics of a trait to use
#' @param trait.names - vector containing the trait names in the order they're found in 'input files'
#' @return - returns the same list of data frames, but with updated column names, begin.time and log.file
#' @importFrom stats pnorm
#' @keywords internal

# NOT EXPORTED @export
# @usage munge_sumstats(list(X,Y),c("X","Y"))

# HEAVILY INSPIRED by genomicSEM::munge @ https://github.com/GenomicSEM/GenomicSEM/blob/master/R/munge.R
munge_sumstats <- function(input.files,trait.names){
  log.name = paste(trait.names,collapse="_")
  log.file = file(paste0(log.name, "-munging.log"),open="wt")
  begin.time = Sys.time()

  cat(print(paste0("Munging ", length(trait.names), " summary statistics files: ",log.name, ", started at ",begin.time), sep = ""),file=log.file,sep="\n",append=TRUE)
  cat(paste0("     "),file=log.file,sep="\n",append=TRUE)

  # for loop over the trait names (i)
  for(i in 1:length(trait.names)){
    file = input.files[[i]]
    Xcols = toupper(names(file))
    namesX = Xcols
    #RSID
    if("RSID" %in% Xcols){cat(print(paste0("Interpreting the RSID column as the SNP column in ",trait.names[i],".")),file=log.file,sep="\n",append=TRUE)
    }else{
      Xcols[Xcols %in% c("SNP","SNPID","RS_NUMBER","RS_NUMBERS", "MARKERNAME", "ID","PREDICTOR","SNP_ID","RS")] <- "RSID"
      if(length(base::setdiff(namesX,Xcols)) > 0) cat(print(paste0("Interpreting the ", setdiff(namesX, Xcols), " column in ", trait.names[i] ," as the SNP column.")),file=log.file,sep="\n",append=TRUE)
    }
    namesX = Xcols
    #beta
    if("BETA" %in% Xcols){cat(print(paste0("Interpreting the BETA column as the effect column in ",trait.names[i],".")),file=log.file,sep="\n",append=TRUE)
    }else{
      Xcols[Xcols %in% c("B","BETA","LOG_ODDS","EFFECTS","EFFECT","SIGNED_SUMSTAT","EST", "BETA1", "LOGOR","BETA_WF")] <- "BETA"
      if(length(base::setdiff(namesX,Xcols)) > 0) cat(print(paste0("Interpreting the ", setdiff(namesX, Xcols), " column in ", trait.names[i] ," as the effect column.")),file=log.file,sep="\n",append=TRUE)
    }
    namesX = Xcols
    #se
    if("SE" %in% Xcols){cat(print(paste0("Interpreting the SE column as the SE column in ",trait.names[i],".")),file=log.file,sep="\n",append=TRUE)
    }else{
      Xcols[Xcols %in% c("STANDARD_ERROR","STD_ERR","STDERR","SE_BETA_WF","SE_BETA")] <- "SE"
      if(length(base::setdiff(namesX,Xcols)) > 0) cat(print(paste0("Interpreting the ", setdiff(namesX, Xcols), " column in ", trait.names[i] ," as the SE column.")),file=log.file,sep="\n",append=TRUE)
    }
    namesX = Xcols
    #TSTAT - MAKE COLUMN IF MISSING
    if("TSTAT" %in% Xcols){
      cat(print(paste0("Interpreting the TSTAT column as the standardised effect column in ",trait.names[i],".")),file=log.file,sep="\n",append=TRUE)
    }else if(sum(Xcols %in% c("TSTAT","Z","ZSCORE","ZSTAT","ZSTATISTIC")) > 0){
      Xcols[Xcols %in% c("TSTAT","Z","ZSCORE","ZSTAT","ZSTATISTIC")] <- "TSTAT"
      if(length(base::setdiff(namesX,Xcols)) > 0) cat(print(paste0("Interpreting the ", setdiff(namesX, Xcols), " column in ", trait.names[i] ," as the standardised effect column.")),file=log.file,sep="\n",append=TRUE)
    }else{
      input.files[[i]]$TSTAT = input.files[[i]]$BETA/input.files[[i]]$SE
      Xcols = c(Xcols,"TSTAT")
      cat(print(paste0("Calculated the standardised effect column in ", trait.names[i] ,".")),file=log.file,sep="\n",append=TRUE)
    }
    namesX = Xcols
    #pval - MAKE COLUMN IF MISSING
    if("PVAL" %in% Xcols){cat(print(paste0("Interpreting the PVAL column as the P-value column in ",trait.names[i],".")),file=log.file,sep="\n",append=TRUE)
    }else if(sum(Xcols %in% c("P","PVALUE","PVAL","P_VALUE","P-VALUE","P.VALUE","P_VAL","GC_PVALUE","WALD_P","P_BETA_WF")) > 0){
      Xcols[Xcols %in% c("P","PVALUE","PVAL","P_VALUE","P-VALUE","P.VALUE","P_VAL","GC_PVALUE","WALD_P","P_BETA_WF")] <- "PVAL"
      if(length(base::setdiff(namesX,Xcols)) > 0) cat(print(paste0("Interpreting the ", setdiff(namesX, Xcols), " column in ", trait.names[i] ," as the P-value column.")),file=log.file,sep="\n",append=TRUE)
    }else{
      input.files[[i]]$PVAL = 2*pnorm(-abs(input.files[[i]]$TSTAT))
      Xcols = c(Xcols,"PVAL")
      cat(print(paste0("Calculated the p-value column in ", trait.names[i] ,".")),file=log.file,sep="\n",append=TRUE)
    }
    namesX = Xcols
    #alt
    if("A1" %in% Xcols){cat(print(paste0("Interpreting the A1 column as the effect allele column in ",trait.names[i],".")),file=log.file,sep="\n",append=TRUE)
    }else{
      Xcols[Xcols %in% c("A1", "ALLELE1","EFFECT_ALLELE","ALT","EA")] <- "A1"
      if(length(base::setdiff(namesX,Xcols)) > 0) cat(print(paste0("Interpreting the ", setdiff(namesX, Xcols), " column in ", trait.names[i] ," as the effect allele column.")),file=log.file,sep="\n",append=TRUE)
    }
    namesX = Xcols
    #ref
    if("A2" %in% Xcols){cat(print(paste0("Interpreting the A2 column as the alternate allele column in ",trait.names[i],".")),file=log.file,sep="\n",append=TRUE)
    }else{
      Xcols[Xcols %in% c("A2","ALLELE2","ALLELE0","OTHER_ALLELE","REF","NON_EFFECT_ALLELE","DEC_ALLELE","OA","NEA", "ALT")] <- "A2"
      if(length(base::setdiff(namesX,Xcols)) > 0) cat(print(paste0("Interpreting the ", setdiff(namesX, Xcols), " column in ", trait.names[i] ," as the alternate allele column.")),file=log.file,sep="\n",append=TRUE)
    }
    namesX = Xcols
    #N
    if("N" %in% Xcols){cat(print(paste0("Interpreting the N column as the sample size column in ",trait.names[i],".")),file=log.file,sep="\n",append=TRUE)
    }else{
      Xcols[Xcols %in% c("N","NCOMPLETESAMPLES", "TOTALSAMPLESIZE", "TOTALN", "TOTAL_N","N_COMPLETE_SAMPLES", "SAMPLESIZE","N_REG")] <- "N"
      if(length(base::setdiff(namesX,Xcols)) > 0) cat(print(paste0("Interpreting the ", setdiff(namesX, Xcols), " column in ", trait.names[i] ," as the sample size column.")),file=log.file,sep="\n",append=TRUE)
    }
    namesX = Xcols
    #CHR
    if("CHR" %in% Xcols){cat(print(paste0("Interpreting the CHR column as the chromosome column in ",trait.names[i],".")),file=log.file,sep="\n",append=TRUE)
    }else{
      Xcols[Xcols %in% c("HG18CHR","CHR", "CHROM", "CHROMOSOME")] <- "CHR"
      if(length(base::setdiff(namesX,Xcols)) > 0) cat(print(paste0("Interpreting the ", setdiff(namesX, Xcols), " column in ", trait.names[i] ," as the chromosome column.")),file=log.file,sep="\n",append=TRUE)
    }
    namesX = Xcols
    #BP
    if("POS" %in% Xcols){cat(print(paste0("Interpreting the POS column as the base pair/position column in ",trait.names[i],".")),file=log.file,sep="\n",append=TRUE)
    }else{
      Xcols[Xcols %in% c("BP","POS","POSITION")] <- "POS"
      if(length(base::setdiff(namesX,Xcols)) > 0) cat(print(paste0("Interpreting the ", setdiff(namesX, Xcols), " column in ", trait.names[i] ," as the base pair/position column.")),file=log.file,sep="\n",append=TRUE)
    }
    namesX = Xcols

    # Print a message for missing RSID, BETA, Pvalue, effect or other allele columns and sample size
    if(sum(Xcols %in% "RSID") == 0) cat(print(paste0('Cannot find an \'rsid\' column, try renaming it to RSID in the summary statistics file for:',trait.names[i])),file=log.file,sep="\n",append=TRUE)
    if(sum(Xcols %in% "BETA") == 0) cat(print(paste0('Cannot find an \'effect\' column, try renaming it to BETA in the summary statistics file for:',trait.names[i])),file=log.file,sep="\n",append=TRUE)
    if(sum(Xcols %in% "PVAL") == 0) cat(print(paste0('Cannot find a \'Pvalue\' column, try renaming it to PVAL in the summary statistics file for:',trait.names[i])),file=log.file,sep="\n",append=TRUE)
    if(sum(Xcols %in% "A1") == 0) cat(print(paste0('Cannot find an \'effect allele\' column, try renaming it to A1 in the summary statistics file for:',trait.names[i])),file=log.file,sep="\n",append=TRUE)
    if(sum(Xcols %in% "A2") == 0) cat(print(paste0('Cannot find an \'other allele\' column, try renaming it to A2 in the summary statistics file for:',trait.names[i])),file=log.file,sep="\n",append=TRUE)
    if(sum(Xcols %in% "N") == 0) cat(print(paste0('Cannot find a \'sample size\' column, try renaming it to N in the summary statistics file for:',trait.names[i])),file=log.file,sep="\n",append=TRUE)
    if(sum(Xcols %in% "CHR") == 0) cat(print(paste0('Cannot find a \'chromosome\' column, try renaming it to CHR in the summary statistics file for:',trait.names[i])),file=log.file,sep="\n",append=TRUE)
    if(sum(Xcols %in% "POS") == 0) cat(print(paste0('Cannot find a \'position\' column, try renaming it to POS in the summary statistics file for:',trait.names[i])),file=log.file,sep="\n",append=TRUE)

    # Throw warnings for missing RSID, BETA, Pvalue, effect or other allele columns and sample size
    if(sum(Xcols %in% "RSID") == 0) warning(paste0('Cannot find an \'rsid\' column, try renaming it to RSID in the summary statistics file for:',trait.names[i]))
    if(sum(Xcols %in% "BETA") == 0) warning(paste0('Cannot find an \'effect\' column, try renaming it to BETA in the summary statistics file for:',trait.names[i]))
    if(sum(Xcols %in% "PVAL") == 0) warning(paste0('Cannot find a \'Pvalue\' column, try renaming it to PVAL in the summary statistics file for:',trait.names[i]))
    if(sum(Xcols %in% "A1") == 0) warning(paste0('Cannot find an \'effect allele\' column, try renaming it to A1 in the summary statistics file for:',trait.names[i]))
    if(sum(Xcols %in% "A2") == 0) warning(paste0('Cannot find an \'other allele\' column, try renaming it to A2 in the summary statistics file for:',trait.names[i]))
    if(sum(Xcols %in% "N") == 0) warning(paste0('Cannot find a \'sample size\' column, try renaming it to N in the summary statistics file for:',trait.names[i]))
    if(sum(Xcols %in% "CHR") == 0) warning(paste0('Cannot find a \'chromosome\' column, try renaming it to CHR in the summary statistics file for:',trait.names[i]))
    if(sum(Xcols %in% "POS") == 0) warning(paste0('Cannot find a \'position\' column, try renaming it to POS in the summary statistics file for:',trait.names[i]))

    # Print a warning message when multiple columns are interpreted as RSID, BETA, Pvalue, effect or other allele columns and sample size
    if(sum(Xcols %in% "RSID") > 1) cat(print(paste0('Multiple columns are being interpreted as the rsid column. Try renaming the column you dont want interpreted as the rsid to RSID2 for:',trait.names[i])),file=log.file,sep="\n",append=TRUE)
    if(sum(Xcols %in% "BETA") > 1) cat(print(paste0('Multiple columns are being interpreted as the beta or effect column. Try renaming the column you dont want interpreted as the beta or effect column to BETA2 for:',trait.names[i])),file=log.file,sep="\n",append=TRUE)
    if(sum(Xcols %in% "PVAL") > 1) cat(print(paste0('Multiple columns are being interpreted as the Pvalue column. Try renaming the column you dont want interpreted as PVAL to PVAL2 for:',trait.names[i])),file=log.file,sep="\n",append=TRUE)
    if(sum(Xcols %in% "A1") > 1) cat(print(paste0('Multiple columns are being interpreted as the effect allele column. Try renaming the column you dont want interpreted as effect allele column to A1_2 for:',trait.names[i])),file=log.file,sep="\n",append=TRUE)
    if(sum(Xcols %in% "A2") > 1) cat(print(paste0('Multiple columns are being interpreted as the other allele column. Try renaming the column you dont want interpreted as the other allele column to A2_2 for:',trait.names[i])),file=log.file,sep="\n",append=TRUE)
    if(sum(Xcols %in% "N") > 1) cat(print(paste0('Multiple columns are being interpreted as the sample size column. Try renaming the column you dont want interpreted as the sample size column to N2 for:',trait.names[i])),file=log.file,sep="\n",append=TRUE)
    if(sum(Xcols %in% "CHR") > 1) cat(print(paste0('Multiple columns are being interpreted as the chromosome column. Try renaming the column you dont want interpreted as the chromosome column to CHR2 for:',trait.names[i])),file=log.file,sep="\n",append=TRUE)
    if(sum(Xcols %in% "POS") > 1) cat(print(paste0('Multiple columns are being interpreted as the position column. Try renaming the column you dont want interpreted as the position column to POS2 for:',trait.names[i])),file=log.file,sep="\n",append=TRUE)

    # Replace the original column names
    names(input.files[[i]]) <- Xcols

  }
  close(log.file)
  return(list(input.files, begin.time,log.name))
}
