#' Merge summary statistics into a single input data frame
#'
#' @param input.files - list of data frames, where each data frame contains the summary statistics of a trait to use in the order of Exposure - Outcome
#' @param trait.names - Vector containing the trait names in the order they're found in 'input files'
#' @param LD.file - LD scores file, either obtained from Alkes group (1000G) or the one provided in the github (UK10K)
#' @param rho.file - Genotyped SNP-specific (local) LD scores
#' @param mafT - Minor allele frequency threshold of selection, to be used if a MAF column is found in the summary statistics file. Default value = 0.005
#' @param infoT - SNP imputation quality threshold, to be used if an INFO column is found in the summary statistics file. Default value = 0.99
#'
#' @return Returns a data frame where the summary statistics file, the LD file, and the SNP-specific LD file are merged
#' @importFrom data.table fread
#' @importFrom dplyr arrange
#' @importFrom dplyr inner_join
#' @export
#'
#' @examples
merge_sumstats <- function(input.files,trait.names,LD.filepath,rho.filepath,mafT=0.005,infoT=0.99){

  # Obtain the correct column names for the input files
  input1 = munge_sumstats(input.files,trait.names)
  begin.time = input1[[2]]
  log.file = input1[[3]]
  Xfile = input1[[1]][[1]]
  Yfile = input1[[1]][[2]]

  # Obtain the correct column names for the LD files
  LD.file = fread(LD.filepath)
  rho.file = fread(rho.filepath)

  input2 = munge_LDfiles(list(LD.file,rho.file),c("LDfile","rhofile"),log.file)
  LDfile = input2[[1]][[1]]
  RHOfile = input2[[1]][[2]]

  # Remove the HLA region due to highly associated SNPs
  LDfile_ind = which(!(LDfile$CHR==6 & LDfile$POS>=28.5e6 & LDfile$POS<=33.5e6))
  LDfile = LDfile[LDfile_ind]

  # Filter for MAF/info if columns are present + order by chr/pos before slicing
  if("MAF" %in% colnames(Xfile)){
    Xfile_ind = which(Xfile$MAF>mafT)
    Xfile = Xfile[Xfile_ind]
  }
  if("MAF" %in% colnames(Yfile)){
    Yfile_ind = which(Yfile$MAF>mafT)
    Yfile = Yfile[Yfile_ind]
  }
  if("INFO" %in% colnames(Xfile)){
    Xfile_ind = which(Xfile$INFO>infoT)
    Xfile = Xfile[Xfile_ind]
  }
  if("INFO" %in% colnames(Yfile)){
    Yfile_ind = which(Yfile$INFO>infoT)
    Yfile = Yfile[Yfile_ind]
  }

  # Join the exposure and outcome files
  Data = inner_join(Xfile, Yfile,
                    by = c("CHR", "POS", "RSID"))

  # Join with rho file
  Data = inner_join(Data, RHOfile) #based on rsid / chr / pos
  #nrow(Data)

  # Join with LD file
  LDfile$INFO = NULL #lots of SNPs have different info (UKBB vs UK10K?), needed otherwise only 207,914 SNPs left
  Data = inner_join(Data, LDfile) #based on rsid / chr / pos
  #nrow(Data)

  # Reorder based on chromosome and position
  Data %>% arrange(CHR, POS) -> Data

  aligned = which(Data$A1.x==Data$A1.y &
                    Data$A2.x==Data$A2.y)
  swapped = which(Data$A1.x==Data$A2.y &
                    Data$A2.x==Data$A1.y)
  # Correct the effect of swapped alleles as well as the t-stat for one of the two traits
  Data[swapped,'TSTAT.x']=Data[swapped,'TSTAT.x']*-1
  Data[swapped,'BETA.x']=Data[swapped,'BETA.x']*-1
  temp_alt=Data$A1.x
  Data[swapped,'A1.x']=Data[swapped,'A2.x']
  Data[swapped,'A2.x']=temp_alt[swapped]

  Data1=Data[c(aligned,swapped),]

  # Make sure all alleles are upper case for ease of comparison
  Data1$A1.x <- factor(toupper(Data1$A1.x), c("A", "C", "G", "T"))
  Data1$A2.x <- factor(toupper(Data1$A2.x), c("A", "C", "G", "T"))
  Data1$A1.y <- factor(toupper(Data1$A1.y), c("A", "C", "G", "T"))
  Data1$A2.y <- factor(toupper(Data1$A2.y), c("A", "C", "G", "T"))

  # Test swapping
  all(Data1$A1.x==Data1$A1.y)
  all(Data1$A2.x==Data1$A2.y)
  Data = Data1
  Data$A1 = Data$A1.x
  Data$A2 = Data$A2.x
  #nrow(Data)


  # Filter out any NA in the effects
  Data.ind = which(is.na(Data$BETA.x)==F & is.na(Data$BETA.y)==F)
  Data = Data[Data.ind,]

  return(Data)

}
