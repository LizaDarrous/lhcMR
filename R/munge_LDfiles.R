#' Title
#'
#' @param input.files  - list of data frames, where each data frame contains the summary statistics of a trait to use
#' @param trait.names - vector containing the trait names in the order they're found in 'input files'
#' @param log.name - the name of the log.file connection to append log descriptions
#' @return - returns the same list of data frames, but with updated column names, and log.name
#' @importFrom stats pnorm

# NOT EXPORTED @export
# @usage munge_sumstats(list(X,Y),c("X","Y"))

# HEAVILY INSPIRED by genomicSEM::munge @ https://github.com/GenomicSEM/GenomicSEM/blob/master/R/munge.R
munge_LDfiles <- function(input.files,trait.names,log.name){
  log.file = file(paste0(log.name, "-munging.log"), open="a")
  begin.time = Sys.time()

  cat(print(paste0("Munging LD score files, started at ",begin.time), sep = ""),file=log.file,sep="\n",append=TRUE)
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
    if(sum(Xcols %in% "CHR") == 0) cat(print(paste0('Cannot find a \'chromosome\' column, try renaming it to CHR in the summary statistics file for:',trait.names[i])),file=log.file,sep="\n",append=TRUE)
    if(sum(Xcols %in% "POS") == 0) cat(print(paste0('Cannot find a \'position\' column, try renaming it to POS in the summary statistics file for:',trait.names[i])),file=log.file,sep="\n",append=TRUE)

    # Throw warnings for missing RSID, BETA, Pvalue, effect or other allele columns and sample size
    if(sum(Xcols %in% "RSID") == 0) warning(paste0('Cannot find an \'rsid\' column, try renaming it to RSID in the summary statistics file for:',trait.names[i]))
    if(sum(Xcols %in% "CHR") == 0) warning(paste0('Cannot find a \'chromosome\' column, try renaming it to XHR in the summary statistics file for:',trait.names[i]))
    if(sum(Xcols %in% "POS") == 0) warning(paste0('Cannot find a \'position\' column, try renaming it to POS in the summary statistics file for:',trait.names[i]))

    # Print a warning message when multiple columns are interpreted as for missing RSID, BETA, Pvalue, effect or other allele columns and sample size
    if(sum(Xcols %in% "RSID") > 1) cat(print(paste0('Multiple columns are being interpreted as the rsid column. Try renaming the column you dont want interpreted as the rsid to RSID2 for:',trait.names[i])),file=log.file,sep="\n",append=TRUE)
    if(sum(Xcols %in% "CHR") > 1) cat(print(paste0('Multiple columns are being interpreted as the chromosome column. Try renaming the column you dont want interpreted as the chromosome column to CHR2 for:',trait.names[i])),file=log.file,sep="\n",append=TRUE)
    if(sum(Xcols %in% "POS") > 1) cat(print(paste0('Multiple columns are being interpreted as the position column. Try renaming the column you dont want interpreted as the position column to POS2 for:',trait.names[i])),file=log.file,sep="\n",append=TRUE)

    if(grepl(trait.names[i],"LD")){
      #LDSC
      if("LDSC" %in% Xcols){cat(print(paste0("Interpreting the LDSC column as the LD score column in ",trait.names[i],".")),file=log.file,sep="\n",append=TRUE)
      }else{
        Xcols[Xcols %in% c("LDSC","LD", "L2")] <- "LDSC"
        if(length(base::setdiff(namesX,Xcols)) > 0) cat(print(paste0("Interpreting the ", setdiff(namesX, Xcols), " column in ", trait.names[i] ," as the LD score column.")),file=log.file,sep="\n",append=TRUE)
      }
      namesX = Xcols
      #weight
      if("WEIGHT" %in% Xcols){cat(print(paste0("Interpreting the WEIGHT column as the weight column in ",trait.names[i],".")),file=log.file,sep="\n",append=TRUE)
      }
      namesX = Xcols

      # Print a message for missing LDSC or WEIGHT
      if(sum(Xcols %in% "LDSC") == 0) cat(print(paste0('Cannot find an \'LDSC\' column, try renaming it to LDSC in the file for:',trait.names[i])),file=log.file,sep="\n",append=TRUE)
      if(sum(Xcols %in% "WEIGHT") == 0) cat(print(paste0('Cannot find a \'weight\' column, try renaming it to WEIGHT in the file for:',trait.names[i])),file=log.file,sep="\n",append=TRUE)

      # Throw warnings for missing LDSC or WEIGHT
      if(sum(Xcols %in% "LDSC") == 0) warning(paste0('Cannot find an \'LDSC\' column, try renaming it to LDSC in the file for:',trait.names[i]))
      if(sum(Xcols %in% "WEIGHT") == 0) warning(paste0('Cannot find a \'weight\' column, try renaming it to WEIGHT in the file for:',trait.names[i]))

      # Print a warning message when multiple columns are interpreted as LDSC or WEIGHT
      if(sum(Xcols %in% "LDSC") > 1) cat(print(paste0('Multiple columns are being interpreted as the LDSC column. Try renaming the column you dont want interpreted as the LDSC to LDSC2 for:',trait.names[i])),file=log.file,sep="\n",append=TRUE)
      if(sum(Xcols %in% "WEIGHT") > 1) cat(print(paste0('Multiple columns are being interpreted as the weight column. Try renaming the column you dont want interpreted as the weight column to weight2 for:',trait.names[i])),file=log.file,sep="\n",append=TRUE)
    }

    if(grepl(trait.names[i],"rho")){
      #LDSC
      if("PIK" %in% Xcols){cat(print(paste0("Interpreting the piK column as the proportion of effective SNPs in the local-LD distribution column in ",trait.names[i],".")),file=log.file,sep="\n",append=TRUE)
      }
      namesX = Xcols
      #weight
      if("SIGK" %in% Xcols){cat(print(paste0("Interpreting the sigK column as the variance of the local-LD distribution column in ",trait.names[i],".")),file=log.file,sep="\n",append=TRUE)
      }
      namesX = Xcols
      #local_m
      if("LOCAL_M" %in% Xcols){cat(print(paste0("Interpreting the local_m column as the number of SNPs used to calculate the local-LD distribution column in ",trait.names[i],".")),file=log.file,sep="\n",append=TRUE)
      }
      namesX = Xcols

      # Print a message for missing LDSC or WEIGHT or LOCAL_M
      if(sum(Xcols %in% "PIK") == 0) cat(print(paste0('Cannot find a \'piK\' column, try renaming it to piK in the file for:',trait.names[i])),file=log.file,sep="\n",append=TRUE)
      if(sum(Xcols %in% "SIGK") == 0) cat(print(paste0('Cannot find a \'sigK\' column, try renaming it to sigK in the file for:',trait.names[i])),file=log.file,sep="\n",append=TRUE)
      if(sum(Xcols %in% "LOCAL_M") == 0) cat(print(paste0('Cannot find a \'local_m\' column, try renaming it to local_m in the file for:',trait.names[i])),file=log.file,sep="\n",append=TRUE)

      # Throw warnings for missing LDSC or WEIGHT or LOCAL_M
      if(sum(Xcols %in% "PIK") == 0) warning(paste0('Cannot find an \'piK\' column, try renaming it to piK in the file for:',trait.names[i]))
      if(sum(Xcols %in% "SIGK") == 0) warning(paste0('Cannot find a \'sigK\' column, try renaming it to sigK in the file for:',trait.names[i]))
      if(sum(Xcols %in% "LOCAL_M") == 0) warning(paste0('Cannot find a \'local_m\' column, try renaming it to local_m in the file for:',trait.names[i]))

      # Print a warning message when multiple columns are interpreted as LDSC or WEIGHT or LOCAL_M
      if(sum(Xcols %in% "PIK") > 1) cat(print(paste0('Multiple columns are being interpreted as the piK column. Try renaming the column you dont want interpreted as the piK to piK2 for:',trait.names[i])),file=log.file,sep="\n",append=TRUE)
      if(sum(Xcols %in% "SIGK") > 1) cat(print(paste0('Multiple columns are being interpreted as the sigK column. Try renaming the column you dont want interpreted as the sigK column to sigK2 for:',trait.names[i])),file=log.file,sep="\n",append=TRUE)
      if(sum(Xcols %in% "LOCAL_M") > 1) cat(print(paste0('Multiple columns are being interpreted as the local_m column. Try renaming the column you dont want interpreted as the local_m column to local_m2 for:',trait.names[i])),file=log.file,sep="\n",append=TRUE)
    }

    # Replace the original column names
    names(input.files[[i]]) <- Xcols
  }
  close(log.file)
  return(list(input.files, log.name))
}
