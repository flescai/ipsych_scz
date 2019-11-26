library(tidyverse)
library(survcomp)
load("qvalueAndCardiff.RData")

res_durv <- read.table("../Result_extTADA_PosteriorAndFDRWed_Sep_27_23_25_38_2017.txt",
                       header = T, stringsAsFactors = F)
res_disruptive <- read.table("../Result_extTADA_PosteriorAndFDRTue_DisruptiveOneCat_Oct_03_10_45_12_2017.txt",
                             header = T, stringsAsFactors = F)
res_cardiff <- read.table("../extada_iontorrent_UK10K_DURVs.tsv",
                          header = T, stringsAsFactors = F)


convertQVal <- function(qvalVector){
  p_values <- qvalVector * rank(qvalVector) / (max(qvalVector) * length(qvalVector))
  return(p_values)
}


res_durv$pval <- convertQVal(res_durv$qvalue)
res_disruptive$pval <- converQVal(res_disruptive$qvalue)

res_durv <- res_durv %>%
  left_join(res_cardiff[c("Gene","Ion_case_N","Ion_con_N","Ion_OR","Ion_P", 
                          "UK10K_Case_N","UK10K_Con_N","UK10K_OR","UK10K_pval")], by = "Gene")

res_disruptive <- res_disruptive %>%
  left_join(res_cardiff[c("Gene","Ion_case_N","Ion_con_N","Ion_OR","Ion_P", 
                          "UK10K_Case_N","UK10K_Con_N","UK10K_OR","UK10K_pval")], by = "Gene")



###############################################
# new weighting according to effective count
###############################################

effectiveCount <- function(Ncases, Ncontrols){
  Neff = 2 * (Ncases * Ncontrols) / (Ncases + Ncontrols)
  return(Neff)
}



# size of datasets used in extTADA (i.e. the clusters)
# taken directly from extTADA model formula
Ncase = sum(1712,1122,458,253,567,1960,1468,203)
Ncontrol = sum(3375,1682,367,196,350,2699,1121,232)

res_durv$ipsychEffectiveCount <- effectiveCount(Ncase, Ncontrol)
res_disruptive$ipsychEffectiveCount <- effectiveCount(Ncase, Ncontrol)


# Iontorrent targeted sequencing study has 5,207 cases and 4,991 controls

res_durv$IonEffectiveCount <- effectiveCount(5207, 4991)
res_disruptive$IonEffectiveCount <- effectiveCount(5207, 4991)



# UK10K exome sequencing dataset consists of 1,352 cases and 4,769 controls.

res_durv$ukEffectiveCount <- effectiveCount(1352, 4769)
res_disruptive$ukEffectiveCount <- effectiveCount(1352, 4769)



#### function to combine p values depending on existence or not on the studies
#### modified to use the effective counts
#### NB - using Z-TRANSFORM METHOD

combineNonMissing <- function(resultsRowData){
  Ion_P <- as.numeric(unname(unlist(resultsRowData[["Ion_P"]])))
  UK10K_pval <- as.numeric(unname(unlist(resultsRowData[["UK10K_pval"]])))
  iPyschPval <- as.numeric(unname(unlist(resultsRowData[["pval"]])))
  
  ipsychWeight <- as.numeric(unname(unlist(resultsRowData[["ipsychEffectiveCount"]])))
  ionWeight <- as.numeric(unname(unlist(resultsRowData[["IonEffectiveCount"]])))
  ukWeight <- as.numeric(unname(unlist(resultsRowData[["ukEffectiveCount"]])))
  
  combinedP <- NA
  
  if(!is.na(Ion_P) & !is.na(UK10K_pval)){
    combinedP <- combine.test(
      p = c(iPyschPval, Ion_P, UK10K_pval),
      weight = c(ipsychWeight, ionWeight, ukWeight),
      method = "z.transform"
    )
  }
  else if(!is.na(Ion_P) & is.na(UK10K_pval)){
    combinedP <- combine.test(
      p = c(iPyschPval, Ion_P),
      weight = c(ipsychWeight, ionWeight),
      method = "z.transform"
    )
  }
  else if(is.na(Ion_P) & !is.na(UK10K_pval)){
    combinedP <- combine.test(
      p = c(iPyschPval, UK10K_pval),
      weight = c(ipsychWeight, ukWeight),
      method = "z.transform"
    )
  }
  return(combinedP)
}


res_durv$combinedP <- unlist(apply(res_durv, 1, combineNonMissing))
res_durv <- res_durv[order(res_durv$qvalue),]

res_disruptive$combinedP <- unlist(apply(res_disruptive, 1, combineNonMissing))
res_disruptive <- res_disruptive[order(res_disruptive$qvalue),]




###################################
# check with another package
###################################

library(metap)


combineNonMissingSumZ <- function(resultsRowData){
  Ion_P <- as.numeric(unname(unlist(resultsRowData[["Ion_P"]])))
  UK10K_pval <- as.numeric(unname(unlist(resultsRowData[["UK10K_pval"]])))
  iPyschPval <- as.numeric(unname(unlist(resultsRowData[["pval"]])))
  
  ipsychWeight <- as.numeric(unname(unlist(resultsRowData[["ipsychEffectiveCount"]])))
  ionWeight <- as.numeric(unname(unlist(resultsRowData[["IonEffectiveCount"]])))
  ukWeight <- as.numeric(unname(unlist(resultsRowData[["ukEffectiveCount"]])))
  
  combinedP <- NA
  result <- NULL
  
  if(!is.na(Ion_P) & !is.na(UK10K_pval)){
    result <- sumz(
      p = c(iPyschPval, Ion_P, UK10K_pval),
      weights = c(ipsychWeight, ionWeight, ukWeight)
    )
  }
  else if(!is.na(Ion_P) & is.na(UK10K_pval)){
    result <- sumz(
      p = c(iPyschPval, Ion_P),
      weights = c(ipsychWeight, ionWeight)
    )
  }
  else if(is.na(Ion_P) & !is.na(UK10K_pval)){
    result <- sumz(
      p = c(iPyschPval, UK10K_pval),
      weight = c(ipsychWeight, ukWeight)
    )
  }
  if (is.numeric(result$p[1,1])){
    combinedP <- result$p[1,1]
  }
  return(combinedP)
}


res_durv$combinedPSumz <- unlist(apply(res_durv, 1, combineNonMissingSumZ))
res_durv <- res_durv[order(res_durv$qvalue),]

res_disruptive$combinedPSumz <- unlist(apply(res_disruptive, 1, combineNonMissingSumZ))
res_disruptive <- res_disruptive[order(res_disruptive$qvalue),]




##################################################
## redo hypergeometric with ASC gene set
##################################################


ASCgenes <- read_tsv("asd_grove_genes.txt")$GENE

ASDgenes <- X102ASC_GeneDiscovery_v2$gene[X102ASC_GeneDiscovery_v2$ASD_vs_DDID=="ASD"]
DDIDgenes <- X102ASC_GeneDiscovery_v2$gene[X102ASC_GeneDiscovery_v2$ASD_vs_DDID=="DDID"]


NewGeneSets <- rbind(
  GenoveseAllSetsUnique,
  data.frame(
    gene = ASCgenes,
    set = "Grove_ASC_genes"
  ),
  data.frame(
    gene = ASDgenes,
    set = "Grove_ASD_prevalent"
  ),
  data.frame(
    gene = DDIDgenes,
    set = "Grove_DDID_prevalent"
  )
)


sigTadaRes = res_durv[res_durv$PP>0.5,]

NewEnrichmentSetsResults <- hypergeoEnrichment(sigTadaRes$Gene, NewGeneSets, genomeMart)
NewEnrichmentSetsResults = numerize(NewEnrichmentSetsResults, "pvalue")
NewEnrichmentSetsResults$padjust = p.adjust(NewEnrichmentSetsResults$pvalue, method = "BH")

NewEnrichmentPlot <- dotplotEnrichment(NewEnrichmentSetsResults, totTestedGenes=23, maxToDisplay = 20, pvalTreshold = 0.05)

pdf("extTADA_significant_enrichmentPlot_20190513.pdf", paper = "a4", height = 0, width = 0)
print(NewEnrichmentPlot)
dev.off()

write_csv(NewEnrichmentSetsResults[order(NewEnrichmentSetsResults$padjust),], "extTADA_sig_newEnrichmentResults_20190513.csv")


save.image("qvalueAndCardiff_20190813.RData")

