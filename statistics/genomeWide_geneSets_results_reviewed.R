

library(tidyverse)
library(survcomp)
library(meta)
load("genomeWide_geneSets_results-reviewed_20190513.RData")


#######################################
# read in all results new and previous
#######################################

allASCResultsCleanWide <- read_tsv("ASCGeneSets_fisherCombined.tab")
genoveseNormalSetsWide <- read_tsv("../allGeneSets_fisherCombined.tab")
genoveseExpressionSetsWide <-read_tsv("../expressionGeneSets_fisherCombined.tab")


allSetsCombinedWide <- rbind(
  genoveseNormalSetsWide,
  genoveseExpressionSetsWide,
  allASCResultsCleanWide
)

allSetsCombinedWide$weightedCombinedP <- unlist(apply(allSetsCombinedWide, 1, fisherPValWeighted))




###################################################
## META-ANALYSIS
## recombine again previous results as well as new
## in order to replot
###################################################

## see file 


### read the sets results
groupSetsRes <- readRDS("../SZMeta_GeneSets_results_groupSets_clean-cov.RData")
expressSetsRes <- readRDS("../SZMeta_GeneSets_results_expressionSets_clean-cov.RData")
ASDgroupsSet <- readRDS("SZMeta_GeneSets_results_ASCSets_clean-cov.RData")
ASCallSet <- readRDS("SZMeta_GeneSets_results_ASCSets-allonly_clean-cov.RData")
geneSetsRes <- rbind(
  groupSetsRes, 
  expressSetsRes,
  ASCallSet,
  ASDgroupsSet
  )

geneSetsRes = numerize(geneSetsRes, "stderror")

levels(geneSetsRes$set)[levels(geneSetsRes$set)=="ASDgenes"] <- "Grove_ASD_prevalent"
levels(geneSetsRes$set)[levels(geneSetsRes$set)=="DDIDgenes"] <- "Grove_DDID_prevalent"


############################################################
############ gene sets functions for meta ##################
############################################################
### tailor calc meta to the sets dataset

calcMetagenSet <- function(dataset, varCat, setName){
  dataSubset = dataset %>%
    filter(type == varCat & analysis == "Odds ratios" & set == setName)
  metaRes = metagen(log(dataSubset$estimate), dataSubset$stderror, sm="OR")
  resData = data.frame(
    dataset = c("Meta - Fixed effects", "Meta - Random effects", "Heterogeneity test"),
    set = setName,
    type = varCat,
    analysis = c("Odds ratios"),
    estimate = c(exp(metaRes$TE.fixed), exp(metaRes$TE.random), metaRes$I2),
    stderror = c(metaRes$seTE.fixed, metaRes$seTE.random, NA),
    lowerConf = c(exp(metaRes$lower.fixed), exp(metaRes$lower.random), metaRes$lower.I2),
    higherConf = c(exp(metaRes$upper.fixed), exp(metaRes$upper.random), metaRes$upper.I2),
    pvalue = c(metaRes$pval.fixed, metaRes$pval.random, pchisq(metaRes$Q, metaRes$df.Q, lower.tail = F)),
    padjusted = c(NA,NA,NA)
  )
  dataSubset = rbind(dataSubset, resData)
  tempRes = list(dataSubset, metaRes)
  return(tempRes)
}




varCatSets = unique(geneSetsRes$type)
varCatSets = as.character(varCatSets)
setNames = as.character(unique(geneSetsRes$set))


metaStatsSetsOR = list()
resultsStatsSetsMeta = data.frame()

totRuns = length(varCatSets) * length(setNames)
runIndex = 1
pb <- txtProgressBar(min = 0, max = totRuns, style = 3)


for (setName in setNames) {
  for (varCat in varCatSets){
    #writeLines(paste0("calc set", setName, " and var cat ", varCat))
    resSet = calcMetagenSet(geneSetsRes, varCat, setName)
    metaStatsSetsOR = c(metaStatsSetsOR, resSet[[2]])
    resultsStatsSetsMeta = rbind(resultsStatsSetsMeta, resSet[[1]])
    setTxtProgressBar(pb, runIndex)
    writeLines("")
    runIndex = runIndex + 1
  }
}


head(resultsStatsSetsMeta)
write_tsv(resultsStatsSetsMeta, "genomewide_gene-sets_meta-analysis_results_20190513.txt")


## order the classification
classifySets$GroupName = factor(classifySets$GroupName, levels = c("Genome-wide_sets", "Organs", "Diseases", "De-novo", "Brain_specific", "Neuron_types", "Synapses"))

## order all sets according to this
classifySets = classifySets[order(classifySets$GroupName),]
resultsStatsSetsMeta$set = factor(resultsStatsSetsMeta$set, levels = c(classifySets$SetName))

resultsStatsSetsMeta = resultsStatsSetsMeta %>%
  left_join(classifySets, by = c("set" = "SetName"))

############################################################
############ selection of significant sets #################
############################################################

signifSets = as.character(unique(resultsStatsSetsMeta[resultsStatsSetsMeta$dataset =="Meta - Fixed effects" & 
                                                        resultsStatsSetsMeta$pvalue < 0.001 &
                                                        resultsStatsSetsMeta$type %in% c("GnURVdistruptive", "GnURVdamaging", "GndURV", "GnURVmissense", "GnURVsynonymous"),][["set"]]))



############################################################
######## PLOT SIGNIFICANT SETS and ALL CATEGORIES ##########
############################################################

resMetaForPlotSets = resultsStatsSetsMeta %>%
  filter(dataset %in% c("DBS", "Clin", "SWE", "Meta - Fixed effects")) %>%
  filter(type %in% c("GnURVdistruptive", "GnURVdamaging", "GndURV", "GnURVmissense", "GnURVsynonymous")) %>%
  filter(set %in% signifSets) %>%
  collect

resMetaForPlotSets$type = factor(resMetaForPlotSets$type, levels = c("GnURVdistruptive", "GnURVdamaging", "GndURV", "GnURVmissense", "GnURVsynonymous"))



sigSetsAllCategoriesPlot = ggplot(resMetaForPlotSets, aes(x=dataset, y=estimate, label=signif(pvalue, digits = 4)))+
  geom_pointrange(stat = "identity", aes(ymin=lowerConf, ymax=higherConf, colour = dataset))+
  geom_text(vjust=-2, size=2)+
  geom_hline(yintercept=1, colour="lightgreen") +
  facet_grid(type~set)

pdf("genomewide_URV_geneSets20190513_significant-sets-all-categories.pdf", paper = "a4r", width = 0, height = 0)
print(sigSetsAllCategoriesPlot)
dev.off()



############################################################
######## PLOT SIGNIFICANT SETS dURV/Syn ONLY ###############
############################################################


resMetaForPlotSetsLimit = resultsStatsSetsMeta %>%
  filter(dataset %in% c("DBS", "Clin", "SWE", "Meta - Fixed effects")) %>%
  filter(type %in% c("GndURV", "GnURVsynonymous")) %>%
  filter(set %in% signifSets[-c(6)]) %>%
  collect


sigSets_dURVonlyPlot = ggplot(resMetaForPlotSetsLimit, aes(x=dataset, y=estimate, label=signif(pvalue, digits = 4)))+
  geom_pointrange(stat = "identity", aes(ymin=lowerConf, ymax=higherConf, colour = dataset))+
  geom_text(vjust=-2, size=2)+
  geom_hline(yintercept=1, colour="lightgreen") +
  facet_grid(type~set) +
  labs(title = "Odds Ratio in most significant gene-sets (dURV)", 
       x="Datasets and Gene sets", 
       y="OR estimate and 95% CI and Variant categories") + 
  theme(plot.title=element_text(hjust = 0.5, 
                                size =14, 
                                face = "bold"),
        axis.text.x=element_text(angle = 30)
  )

pdf("genomewide_URV_geneSets20190513_significant-sets_dURV-only.pdf", paper = "a4r", width = 0, height = 0)
print(sigSets_dURVonlyPlot)
dev.off()



############################################################
###### PLOT SIGNIFICANT SETS Disruptive/Syn ONLY ###########
############################################################


resMetaForPlotSetsDisruptive = resultsStatsSetsMeta %>%
  filter(dataset %in% c("DBS", "Clin", "SWE", "Meta - Fixed effects")) %>%
  filter(type %in% c("GnURVdistruptive", "GnURVsynonymous")) %>%
  filter(set %in% signifSets) %>%
  collect


sigSets_disruptiveURVonlyPlot = ggplot(resMetaForPlotSetsDisruptive, aes(x=dataset, y=estimate, label=signif(pvalue, digits = 4)))+
  geom_pointrange(stat = "identity", aes(ymin=lowerConf, ymax=higherConf, colour = dataset))+
  geom_text(vjust=-2, size=2)+
  geom_hline(yintercept=1, colour="lightgreen") +
  geom_hline(yintercept=1.0340282, colour="lightblue", linetype = "dotdash") +
  facet_grid(type~set) +
  labs(title = "Odds Ratio in most significant gene-sets (disruptiveURV)", 
       x="Datasets and Gene sets", 
       y="OR estimate and 95% CI and Variant categories") + 
  theme(plot.title=element_text(hjust = 0.5, 
                                size =14, 
                                face = "bold"),
        axis.text.x=element_text(angle = 30)
  )



pdf("genomewide_URV_geneSets20190513_significant-sets_disruptiveURVs-only.pdf", paper = "a4r", width = 0, height = 0)
print(sigSets_disruptiveURVonlyPlot)
dev.off()


##################################################################
## GENOME-WIDE PLOT SIGNIFICANT ALL CATEGORIES MISS CONSTRAINED ##
##################################################################

constraineURVRes = resultsStatsSetsMeta %>%
  filter(dataset %in% c("DBS", "Clin", "SWE", "Meta - Fixed effects")) %>%
  filter(type %in% c("GnURVdistruptive", "GnURVdamaging", "GndURV", "GnURVmissense", "GnURVsynonymous")) %>%
  filter(set %in% "constrained") %>%
  mutate(type = factor(type, levels = c("GnURVdistruptive", "GnURVdamaging", "GndURV", "GnURVmissense", "GnURVsynonymous"))) %>%
  mutate(type = plyr::revalue(type, c("GnURVdistruptive" = "Disruptive", "GnURVdamaging" = "Missense damaging", 
                                      "GndURV" = "Damaging and Disruptive", "GnURVmissense" = "Missense non damaging", 
                                      "GnURVsynonymous" = "Synonymous"))) %>%
 ggplot(aes(x=dataset, y=estimate, label=signif(pvalue, digits = 4)))+
  geom_pointrange(stat = "identity", aes(ymin=lowerConf, ymax=higherConf, colour = dataset))+
  geom_text(vjust=-2, size=2)+
  geom_hline(yintercept=1, colour="lightgreen") +
  facet_grid(set~type)+
  labs(title = "Odds Ratio of different URV categories in missense constrained genes", 
       x="Datasets and Variant categories", 
       y="OR estimate and 95% CI") + 
  theme(plot.title=element_text(hjust = 0.5, 
                                size =14, 
                                face = "bold"),
        axis.text.x=element_text(angle = 30)
  )

pdf("genomewide_URV_geneSets20190513_missense-constrained-set_all-URV-categories.pdf", paper = "a4r", width = 0, height = 0)
print(constraineURVRes)
dev.off()


############################################################
######## PLOT ALL TESTED SETS dURV/Syn ONLY ###############
############################################################


dURVonlyRes = resultsStatsSetsMeta %>%
  filter(dataset %in% c("DBS", "Clin", "SWE", "Meta - Fixed effects")) %>%
  filter(type %in% c("GndURV")) %>%
  mutate(dataset = factor(dataset, levels = c("Meta - Fixed effects", "DBS", "Clin", "SWE"))) %>%
  ggplot(aes(x=dataset, y=estimate, label=signif(pvalue, digits = 4)))+
  geom_pointrange(stat = "identity", aes(ymin=lowerConf, ymax=higherConf, colour = dataset))+
  geom_text(vjust=-2, size=2)+
  geom_hline(yintercept=1, colour="lightgreen") +
  geom_hline(yintercept=1.0156060, colour="lightblue", linetype = "dotdash") +
  coord_flip() +
  facet_wrap(~set, ncol = 4)


dURVonlyResOrgans = resultsStatsSetsMeta %>%
  filter(dataset %in% c("DBS", "Clin", "SWE", "Meta - Fixed effects")) %>%
  filter(type %in% c("GndURV")) %>%
  filter(GroupName == "Organs") %>%
  mutate(dataset = factor(dataset, levels = c("Meta - Fixed effects", "DBS", "Clin", "SWE"))) %>%
  ggplot(aes(x=dataset, y=estimate, label=signif(pvalue, digits = 4)))+
  geom_pointrange(stat = "identity", aes(ymin=lowerConf, ymax=higherConf, colour = dataset))+
  geom_text(vjust=-2, size=2)+
  geom_hline(yintercept=1, colour="lightgreen") +
  coord_flip() +
  facet_wrap(~set, ncol = 1)

dURVonlyResOrgansSimplified = resultsStatsSetsMeta %>%
  filter(dataset %in% c("DBS", "Clin", "SWE", "Meta - Fixed effects")) %>%
  filter(type %in% c("GndURV")) %>%
  filter(GroupName == "Organs") %>%
  filter(set %in% c("brain_specific", "kidney_specific", "heart_specific", "colon_specific")) %>%
  mutate(set = factor(set, levels = c("brain_specific", "kidney_specific", "heart_specific", "colon_specific"))) %>%
  mutate(dataset = factor(dataset, levels = c("Meta - Fixed effects", "DBS", "Clin", "SWE"))) %>%
  ggplot(aes(x=dataset, y=estimate, label=signif(pvalue, digits = 4)))+
  geom_pointrange(stat = "identity", aes(ymin=lowerConf, ymax=higherConf, colour = dataset))+
  geom_text(vjust=-2, size=2)+
  geom_hline(yintercept=1, colour="lightgreen") +
  geom_hline(yintercept=1.0156060, colour="lightblue", linetype = "dotdash") +
  coord_flip() +
  facet_wrap(~set, ncol = 1)

dURVonlyResDiseases = resultsStatsSetsMeta %>%
  filter(dataset %in% c("DBS", "Clin", "SWE", "Meta - Fixed effects")) %>%
  filter(type %in% c("GndURV")) %>%
  filter(GroupName == "Diseases") %>%
  mutate(dataset = factor(dataset, levels = c("Meta - Fixed effects", "DBS", "Clin", "SWE"))) %>%
  ggplot(aes(x=dataset, y=estimate, label=signif(pvalue, digits = 4)))+
  geom_pointrange(stat = "identity", aes(ymin=lowerConf, ymax=higherConf, colour = dataset))+
  geom_text(vjust=-2, size=2)+
  geom_hline(yintercept=1, colour="lightgreen") +
  coord_flip() +
  facet_wrap(~set, ncol = 1)

dURVonlyResDenovo = resultsStatsSetsMeta %>%
  filter(dataset %in% c("DBS", "Clin", "SWE", "Meta - Fixed effects")) %>%
  filter(type %in% c("GndURV")) %>%
  filter(GroupName == "De-novo") %>%
  mutate(dataset = factor(dataset, levels = c("Meta - Fixed effects", "DBS", "Clin", "SWE"))) %>%
  ggplot(aes(x=dataset, y=estimate, label=signif(pvalue, digits = 4)))+
  geom_pointrange(stat = "identity", aes(ymin=lowerConf, ymax=higherConf, colour = dataset))+
  geom_text(vjust=-2, size=2)+
  geom_hline(yintercept=1, colour="lightgreen") +
  coord_flip() +
  facet_wrap(~set, ncol = 1)


dURVonlyResBrain = resultsStatsSetsMeta %>%
  filter(dataset %in% c("DBS", "Clin", "SWE", "Meta - Fixed effects")) %>%
  filter(type %in% c("GndURV")) %>%
  filter(GroupName == "Brain_specific") %>%
  mutate(dataset = factor(dataset, levels = c("Meta - Fixed effects", "DBS", "Clin", "SWE"))) %>%
  ggplot(aes(x=dataset, y=estimate, label=signif(pvalue, digits = 4)))+
  geom_pointrange(stat = "identity", aes(ymin=lowerConf, ymax=higherConf, colour = dataset))+
  geom_text(vjust=-2, size=2)+
  geom_hline(yintercept=1, colour="lightgreen") +
  coord_flip() +
  facet_wrap(~set, ncol = 1)


dURVonlyResNeuron = resultsStatsSetsMeta %>%
  filter(dataset %in% c("DBS", "Clin", "SWE", "Meta - Fixed effects")) %>%
  filter(type %in% c("GndURV")) %>%
  filter(GroupName == "Neuron_types") %>%
  mutate(dataset = factor(dataset, levels = c("Meta - Fixed effects", "DBS", "Clin", "SWE"))) %>%
  ggplot(aes(x=dataset, y=estimate, label=signif(pvalue, digits = 4)))+
  geom_pointrange(stat = "identity", aes(ymin=lowerConf, ymax=higherConf, colour = dataset))+
  geom_text(vjust=-2, size=2)+
  geom_hline(yintercept=1, colour="lightgreen") +
  coord_flip() +
  facet_wrap(~set, ncol = 1)


dURVonlyResSynapses = resultsStatsSetsMeta %>%
  filter(dataset %in% c("DBS", "Clin", "SWE", "Meta - Fixed effects")) %>%
  filter(type %in% c("GndURV")) %>%
  filter(GroupName == "Synapses") %>%
  mutate(dataset = factor(dataset, levels = c("Meta - Fixed effects", "DBS", "Clin", "SWE"))) %>%
  ggplot(aes(x=dataset, y=estimate, label=signif(pvalue, digits = 4)))+
  geom_pointrange(stat = "identity", aes(ymin=lowerConf, ymax=higherConf, colour = dataset))+
  geom_text(vjust=-2, size=2)+
  geom_hline(yintercept=1, colour="lightgreen") +
  coord_flip() +
  facet_wrap(~set, ncol = 1)


############################################################
###### PLOT ALL TESTED SETS Disruptive/Syn ONLY ############
############################################################

disruptiveURVonlyRes = resultsStatsSetsMeta %>%
  filter(dataset %in% c("DBS", "Clin", "SWE", "Meta - Fixed effects")) %>%
  filter(type %in% c("GnURVdistruptive")) %>%
  ggplot(aes(x=dataset, y=estimate, label=signif(pvalue, digits = 4)))+
  geom_pointrange(stat = "identity", aes(ymin=lowerConf, ymax=higherConf, colour = dataset))+
  geom_text(vjust=-2, size=2)+
  geom_hline(yintercept=1, colour="lightgreen") +
  coord_flip() +
  facet_wrap(~set, ncol = 6)










save.image("genomeWide_geneSets_results-reviewed_20190513.RData")
