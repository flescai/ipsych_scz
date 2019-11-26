import hail
import sys

### parsing arguments
### NB: IN ORDER TO MAKE THE SCRIPT GENERAL, this needs to pass the VEP config as argument
### arguments passed to the script:
# (1) input file (either a VCF or VDS)
# (2) base name of output
# (3) folder where things will be placed
# (4) samples to include in the analysis
# (5) vep config file
# (6) optional - repartition number integer


input = str(sys.argv[1])
outVCF = str(sys.argv[3]) + '/' + str(sys.argv[2]) + '.vcf.bgz'
outVDS = str(sys.argv[3]) + '/' + str(sys.argv[2]) + '.vds'
tmpDir = str(sys.argv[3]) + '/tmp'
includeSamples = str(sys.argv[4])
vepConfig = str(sys.argv[5])
varTab = str(sys.argv[3]) + '/' + str(sys.argv[2]) + '.filt.vars.tab'
genoTab = str(sys.argv[3]) + '/' + str(sys.argv[2]) + '.filt.geno.tab'

############# partitioning if passed
if (len(sys.argv) >6):
 sys.stderr.write('\n----------- Received argument for repartitioning --------\n')
 partitions = int(sys.argv[6])
 print(partitions)
else:
 sys.stderr.write('\n----------- No repartitioning requested --------\n')

###########################################################################
### main script
hc = hail.HailContext(tmp_dir=tmpDir)
### this part has been changed so the input can be either a VDS or a VCF file
### and the script becomes flexible to handle this
if (input.find('vcf')!=-1):
 sys.stderr.write('\n----------- Reading input VCF --------\n')
 vds = hc.import_vcf(input, skip_bad_ad=True)
else:
 sys.stderr.write('\n----------- Reading input VDS --------\n')
 vds = hc.read(input)


###########################################################################
######## setting filters for variants and genotypes
## this assumes the variants are not split-multi 
## and therefore it is not possible to run variantqc
qc_expressions = [
 'va.callRate = gs.fraction(g => g.isCalled)',
 'va.callStats = gs.callStats(g => v)' ]
 
###########################################################################
## the annotations are prepared in order to extract the information
## on URV as well
## in order to run these annotations, the VDS needs to be annotated with VEP

anno_expressions = [
 'va.info.AN = va.callStats.AN, va.info.AC = va.callStats.AC[1:], va.info.AF = va.callStats.AF[1:]',
 'va.info.nAlleles = va.callStats.AF.length',
 'va.info.MAF = va.callStats.AF.min',
 'va.info.VEPGENE = va.vep.transcript_consequences.find(c => c.consequence_terms.toSet.contains(va.vep.most_severe_consequence)).gene_symbol',
 'va.info.VEPmostSevere = va.vep.most_severe_consequence',
 'va.info.LOF = va.vep.transcript_consequences.find(c => c.consequence_terms.toSet.contains(va.vep.most_severe_consequence)).lof',
 'va.info.LOFflags = va.vep.transcript_consequences.find(c => c.consequence_terms.toSet.contains(va.vep.most_severe_consequence)).lof_flags',
 'va.info.LOFfilter = va.vep.transcript_consequences.find(c => c.consequence_terms.toSet.contains(va.vep.most_severe_consequence)).lof_filter',
 'va.info.SIFTpred = va.vep.transcript_consequences.find(c => c.consequence_terms.toSet.contains(va.vep.most_severe_consequence)).sift_prediction',
 'va.info.POLYpred = va.vep.transcript_consequences.find(c => c.consequence_terms.toSet.contains(va.vep.most_severe_consequence)).polyphen_prediction',
 'va.info.EXACmaf = va.vep.colocated_variants.find(c => c.exac_allele == va.vep.transcript_consequences.find(c => c.consequence_terms.toSet.contains(va.vep.most_severe_consequence)).variant_allele).exac_maf',
 'va.info.EXACdisease = va.vep.colocated_variants.find(c => c.exac_allele == va.vep.transcript_consequences.find(c => c.consequence_terms.toSet.contains(va.vep.most_severe_consequence)).variant_allele).phenotype_or_disease',
 'va.info.IMPACT = va.vep.transcript_consequences.find(c => c.consequence_terms.toSet.contains(va.vep.most_severe_consequence)).impact' ]


###########################################################################
######## filter samples, then filter genotypes, then filter pass variants
###########################################################################

sys.stderr.write("\n------------Filter and QC--------------\n")
vds_filtered = (vds
 .filter_samples_list(includeSamples, keep=True)
 .filter_genotypes(
  '''
   g.dp < 10 ||
   (g.ad[0] + g.ad[1]) / g.dp < 0.9 ||
   (g.isHomRef && (g.ad[0] / g.dp < 0.9 || g.gq < 25)) ||
   (g.isHet && (g.ad[1] / g.dp < 0.25 || g.pl[0] < 25)) ||
   (g.isHomVar && (g.ad[1] / g.dp < 0.9 || g.pl[0] < 25))
  ''', keep=False)
 .annotate_variants_expr(qc_expressions)
 .filter_variants_expr('va.callRate >= 0.9 && va.pass', keep=True))

###########################################################################
######## annotate the file with VEP
###########################################################################

sys.stderr.write("\n------------Annotate the VCF file-------------\n")
sys.stderr.write("Using vep config\n")
print(vepConfig)
vds_annotated = (vds_filtered
 .vep(vepConfig, block_size=1000)
 .annotate_variants_expr(anno_expressions))
sys.stderr.write("\n------------Annotate the Singletons-------------\n")
vds_annotated = vds_annotated.annotate_variants_expr('va.info.isSingleton = gs.filter(g => g.isCalledNonRef).count() == 1')
sys.stderr.write("\n------------Annotate the URVs-------------\n")
vds_annotated = vds_annotated.annotate_variants_expr('va.info.isURV = va.info.isSingleton && !isDefined(va.info.EXACmaf)')
vds_annotated = vds_annotated.annotate_variants_expr('va.info.URVmissense = va.info.isURV && va.info.VEPmostSevere == "missense_variant"')
vds_annotated = vds_annotated.annotate_variants_expr('va.info.URVlofHC = va.info.isURV && va.info.LOF == "HC"')
vds_annotated = vds_annotated.annotate_variants_expr('va.info.URVcategory = if (va.info.isURV && va.info.VEPmostSevere == "synonymous_variant") "synURV" else (if (va.info.isURV && va.info.IMPACT == "HIGH") "disruptive" else (if ((va.info.isURV && va.info.IMPACT == "MODERATE" && va.info.SIFTpred == "deleterious" && "damaging" ~ va.info.POLYpred) || (va.info.isURV && va.info.VEPmostSevere == "inframe_insertion") || (va.info.isURV && va.info.VEPmostSevere == "inframe_deletion")) "damaging" else (if (va.info.isURV && va.info.VEPmostSevere == "missense_variant") "missense" else "other")))')

###########################################################################
########## repartition if needed
###########################################################################

if (len(sys.argv) >6):
 sys.stderr.write('\n----------- Repartitioning Data --------\n')
 vds_annotated = vds_annotated.repartition(partitions, shuffle=True)
else:
 sys.stderr.write('\n----------- No repartitioning requested --------\n')
 
 
###########################################################################
######### Exporting URVs data
###########################################################################

try:
 sys.stderr.write('\n----------- Exporting Variants Tab --------\n')
 vds_annotated.export_variants(varTab, 'VAR = v, CONS = va.info.VEPmostSevere, SIFT = va.info.SIFTpred, POLYPHEN = va.info.POLYpred, IMPACT = va.info.IMPACT, LOF = va.info.LOF, SING = va.info.isSingleton, ExaC = va.info.EXACmaf, URV = va.info.isURV, URVmissense = va.info.URVmissense, URVlof = va.info.URVlofHC, URVcat = va.info.URVcategory')
except:
 sys.stderr.write('\n----------- EXCEPTION caught while exporting variants Tab --------\n')

try:
 sys.stderr.write('\n----------- Exporting Genotypes Tab --------\n')
 vds_annotated.export_genotypes(genoTab, 'SAMPLE = s.id, VAR = v, CONS = va.info.VEPmostSevere, GENE = va.info.VEPGENE, SIFT = va.info.SIFTpred, POLYPHEN = va.info.POLYpred, IMPACT = va.info.IMPACT, LOF = va.info.LOF, SING = va.info.isSingleton, ExaC = va.info.EXACmaf, URV = va.info.isURV, URVmissense = va.info.URVmissense, URVlof = va.info.URVlofHC, URVcat = va.info.URVcategory')
except:
 sys.stderr.write('\n----------- EXCEPTION caught while exporting Genotypes Tab --------\n')


###########################################################################
######## save as VCF file
###########################################################################

try:
 sys.stderr.write('\n---------- Save as VCF ---------\n')
 vds_annotated.export_vcf(outVCF)
except:
 sys.stderr.write('\n----------- EXCEPTION caught while saving VCF --------\n')


###########################################################################
######## save as VDS file
###########################################################################

try:
 sys.stderr.write('\n----------- Save as VDS --------\n')
 vds_annotated.write(outVDS, overwrite=True)
except:
 sys.stderr.write('\n----------- EXCEPTION caught while saving VDS --------\n')

sys.stderr.write("\n------------Completed--------------\n")


###########################################################################