#' auprc, fr and ttif calculator
#' 
#' Copyright 2019 PARKER INSTITUTE FOR CANCER IMMUNOTHERAPY
#' Permission is hereby granted, free of charge, to any person obtaining a 
#' copy of this software and associated documentation files (the "Software"), 
#' to deal in the Software without restriction, including without limitation 
#' the rights to use, copy, modify, merge, publish, distribute, sublicense, 
#' and/or sell copies of the Software, and to permit persons to whom the Software 
#' is furnished to do so, subject to the following conditions:
#'   
#' The above copyright notice and this permission notice shall be included in all copies or 
#' substantial portions of the Software.
#' 
#' THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
#' INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR 
#' PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY 
#' CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, 
#' OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE
#' 
#' Questions should be directed to Danny Wells: dwells < at > parkerici.org
#' 



library(PRROC)
library(dplyr)


#' all functions below take in two tables: 
#' @ranked.list a ranked list of pmhc. Must include columns (with the exact name matching)
#' hla: the allele of the pmhc, in format A0201 (i.e, <HLA Gene><Two digit serotype><Two digit allele ID>)
#' rank: the rank of the pmhc. 1 should denote the top ranked peptide and should increase from there. 
#' sequence: the protein sequence of the putative neoepitope
#' 
#' @tested.list an unranked list of tested pmhc. Must include columns with the exact name: 
#' hla: the allele of the pmhc, in format A0201 (i.e, <HLA Gene><Two digit serotype><Two digit allele ID>)
#' sequence: the protein sequence of the tested peptide
#' validated: a boolean variable, if the peptide validated or not. TRUE: validated; FALSE: not validated
#'
#' In both tables, every included pMHC must appear only once.

#' auprc, calculated using PRROC
auprc.calculation <- function(ranked.list, tested.list){
  #merge the lists
  merged.list <- ranked.list %>% 
    inner_join(tested.list, 
               by=c('hla'='hla','sequence'='sequence'))
  if (dim(merged.list)[[1]][1]>0){
    #Calculate AUPRC
    validated.ranks<-merged.list %>% filter(validated==TRUE) %>% .$rank
    not.validated.ranks<-merged.list %>% filter(validated==FALSE) %>% .$rank
    auprc <- PRROC::pr.curve(-validated.ranks, -not.validated.ranks)$auc.integral
    return(auprc)
  } else {
    return(NA)
  }
}

#' fraction ranked (fr)
fr.calculation <- function(ranked.list, tested.list){
  #merge the lists
  merged.list <- ranked.list %>% 
    filter(rank<=100) %>% 
    inner_join(tested.list, 
               by=c('hla'='hla','sequence'='sequence'))
  if (dim(merged.list)[[1]][1]>0){
    #Calculate AUPRC
    fr<-dim(merged.list %>% filter(validated==TRUE))[[1]][1]/dim(tested.list %>% filter(validated==TRUE))[[1]][1]
    return(fr)
  } else {
    return(NA)
  }
}

#' top-twenty immunogenic fraction (ttif)
ttif.calculation <- function(ranked.list, tested.list){
  #merge the lists
  merged.list <- ranked.list %>% 
    filter(rank<20) %>% 
    inner_join(tested.list, 
               by=c('hla'='hla','sequence'='sequence'))
  merged.list %>% filter(validated==TRUE)
  if (dim(merged.list)[[1]][1]>0){
    #Calculate AUPRC
    ttif <- (dim(merged.list %>% filter(validated==TRUE))[[1]][1])/(dim(merged.list)[[1]][1])
    return(ttif)
  } else {
    return(NA)
  }
}

