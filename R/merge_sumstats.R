#' Title
#'
#' @param input.files - list of data frames, where each data frame contains the summary statistics of a trait to use
#' @param trait.names - vector containing the trait names in the order they're found in 'input files'
#' @param LD.file - LD scores file, either obtained from Alkes group (1000G) or the one provided (UK10K)
#' @param rho.file - Genotyped SNP-specific LD scores
#' @param mafT - Minor allele frequency threshold of selection
#' @param infoT - SNP imputation quality threshold
#'
#' @return
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

  Xfile %>% arrange(CHR, POS) -> X_data
  Yfile %>% arrange(CHR, POS) -> Y_data
  Data = inner_join(X_data, Y_data,
                    by = c("CHR", "POS", "RSID"))

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
  # Test swapping
  all(Data1$A1.x==Data1$A1.y)
  all(Data1$A2.x==Data1$A2.y)
  Data = Data1
  Data$A1 = Data$A1.x
  Data$A2 = Data$A2.x
  #nrow(Data)


  Data = inner_join(Data, RHOfile) #based on rsid / chr / pos
  #nrow(Data)

  LDfile$INFO = NULL #lots of SNPs have different info (UKBB vs UK10K?), needed otherwise only 207,914 SNPs left
  Data = inner_join(Data, LDfile) #based on rsid / chr / pos
  #nrow(Data)

  # Filter out any NA in the effects
  X.na = which(is.na(Data$BETA.x)==T)
  Y.na = which(is.na(Data$BETA.y)==T)
  Data = Data[-unique(c(X.na,Y.na)),]

  return(Data)

}
