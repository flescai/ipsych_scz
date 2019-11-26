import hail
import sys
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import seaborn
from math import log, isnan

### parsing arguments
### NB: IN ORDER TO MAKE THE SCRIPT GENERAL
### arguments passed to the script:
# (1) input file (either a VCF or VDS)
# (2) folder where things will be placed
# (3) base name of output

input = str(sys.argv[1])
baseDir = str(sys.argv[2])
tmpDir = str(sys.argv[2]) + '/tmp'
baseName = str(sys.argv[3])
urvCounts = baseDir + '/' + baseName + '_samples_URVcounts_byGene.tab'
outVDS = str(sys.argv[3]) + '/' + str(sys.argv[2]) + '.urvCat_byGene.vds'
logFile = baseDir + '/' + baseName + '_URVcounts_byGene.log'

### main script
hc = hail.HailContext(tmp_dir=tmpDir, log=logFile)
### this part has been changed so the input can be either a VDS or a VCF file
### and the script becomes flexible to handle this
if (input.find('vcf')!=-1):
 sys.stderr.write('\n----------- Reading input VCF --------\n')
 vds = hc.import_vcf(input, skip_bad_ad=True)
else:
 sys.stderr.write('\n----------- Reading input VDS --------\n')
 vds = hc.read(input)

sampleAnnotations = [
 'sa.nCalledVariants = gs.filter(g => g.isCalled).count()',
 'sa.nNonRefVariants = gs.filter(g => g.isCalledNonRef()).count()',
 'sa.totLOF = gs.filter(g => va.info.LOF == "HC" && g.isCalledNonRef()).count()',
 'sa.totMissense = gs.filter(g => va.vep.most_severe_consequence == "missense_variant" && g.isCalledNonRef()).count()',
 'sa.nTotalURVs = gs.filter(g => va.info.isURV && g.isCalledNonRef()).count()',
 'sa.nURVdistruptive = gs.filter(g => va.info.URVcategory == "disruptive" && g.isCalledNonRef()).count()',
 'sa.nURVdamaging = gs.filter(g => va.info.URVcategory == "damaging" && g.isCalledNonRef()).count()',
 'sa.nURVmissense = gs.filter(g => va.info.URVcategory == "missense" && g.isCalledNonRef()).count()',
 'sa.nURVsynonymous = gs.filter(g => va.info.URVcategory == "synURV" && g.isCalledNonRef()).count()']

sys.stderr.write('\n----------- Annotating Samples --------\n')
vds_annotated = vds.annotate_samples_expr(sampleAnnotations)

tableHeader = 'SAMPLE = s, nCalledVariants = sa.nCalledVariants, nNonRefVariants = sa.nNonRefVariants, totMissense = sa.totMissense, totLOF = sa.totLOF, totURV = sa.nTotalURVs, nURVdistruptive = sa.nURVdistruptive, nURVdamaging = sa.nURVdamaging, nURVmissense = sa.nURVmissense, nURVsynonymous = sa.nURVsynonymous'

sys.stderr.write('\n----------- Getting List of Genes --------\n')
genes = vds_annotated.query_variants('variants.map(v => va.info.VEPGENE).collect()')
genesUnique = list(set(genes))
genesUnique = filter(None, genesUnique)

sys.stderr.write('\n----------- Annotating the genes ONE at a time to avoid errors --------\n')
for i in genesUnique:
 geneAnnotations = []
 geneLabel = i.replace(".", "_")
 geneLabel = geneLabel.replace("-", "_")
 syn = 'sa.' + geneLabel + '_syn = gs.filter(g => va.info.URVcategory == "synURV" && va.info.VEPGENE == "' + i + '" && g.isCalledNonRef()).count()'
 miss = 'sa.' + geneLabel + '_miss = gs.filter(g => va.info.URVcategory == "missense" && va.info.VEPGENE == "' + i + '" && g.isCalledNonRef()).count()'
 dam = 'sa.' + geneLabel + '_dam = gs.filter(g => va.info.URVcategory == "damaging" && va.info.VEPGENE == "' + i + '" && g.isCalledNonRef()).count()'
 dis = 'sa.' + geneLabel + '_dis = gs.filter(g => va.info.URVcategory == "disruptive" && va.info.VEPGENE == "' + i + '" && g.isCalledNonRef()).count()'
 lofCount = 'sa.' + geneLabel + '_lofCount = gs.filter(g => va.info.LOF == "HC" && va.info.VEPGENE == "' + i + '" && g.isCalledNonRef()).count()'
 URVtot = 'sa.' + geneLabel + '_URVtot = gs.filter(g => va.info.isURV && va.info.VEPGENE == "' + i + '" && g.isCalledNonRef()).count()'
 missTot = 'sa.' + geneLabel + '_missTot = gs.filter(g => va.vep.most_severe_consequence == "missense_variant" && va.info.VEPGENE == "' + i + '" && g.isCalledNonRef()).count()'
 geneAnnotations.append(syn) 
 geneAnnotations.append(miss)
 geneAnnotations.append(dam)
 geneAnnotations.append(dis)
 geneAnnotations.append(lofCount)
 geneAnnotations.append(URVtot)
 geneAnnotations.append(missTot)
 geneAnnotations = [str(x) for x in geneAnnotations]
 try:
  sys.stderr.write('\n----------- Annotating Gene = ' + i + ' --------\n')
  vds_annotated = vds_annotated.annotate_samples_expr(geneAnnotations)
  tableHeader = tableHeader + ', ' + geneLabel + '_syn = sa.' + geneLabel + '_syn' + ', ' + geneLabel + '_miss = sa.' + geneLabel + '_miss' + ', ' + geneLabel + '_dam = sa.' + geneLabel + '_dam' + ', ' + geneLabel + '_dis = sa.' + geneLabel + '_dis'
  tableHeader = tableHeader + ', ' + geneLabel + '_lofCount = sa.' + geneLabel + '_lofCount' + ', ' + geneLabel + '_URVtot = sa.' + geneLabel + '_URVtot' + ', ' + geneLabel + '_missTot = sa.' + geneLabel + '_missTot'
 except:
  sys.stderr.write('\n----------- Gene = ' + i + ' gave an error --------\n')


sys.stderr.write('\n----------- Preview of results --------\n')
samples_counts = vds_annotated.samples_keytable().to_pandas()
print(samples_counts.head(n=10))

sys.stderr.write('\n----------- Writing URV counts --------\n')
vds_annotated.export_samples(urvCounts, tableHeader)


###########################################################################
######## save as VDS file
###########################################################################

sys.stderr.write('\n----------- Save as VDS --------\n')
vds_annotated.write(outVDS, overwrite=True)
sys.stderr.write("\n------------Completed--------------\n")


###########################################################################
