import hail
import sys

### parsing arguments
### NB: IN ORDER TO MAKE THE SCRIPT GENERAL
### arguments passed to the script:
# (1) input file (either a VCF or VDS)
# (2) folder where things will be placed
# (3) base name of output
# (4) variants to be used for PCA

input = str(sys.argv[1])
baseDir = str(sys.argv[2])
tmpDir = str(sys.argv[2]) + '/tmp'
baseName = str(sys.argv[3])
varLDprunedForPCA = str(sys.argv[4])
scoreFile = baseDir + '/' + baseName + '_pca_samples_scores.tab'
loadingsFile = baseDir + '/' + baseName + '_pca_variants_loadings.tab'

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

varTable = hail.KeyTable.import_bed(varLDprunedForPCA)

vds_forPCA = (vds
 .split_multi()
 .filter_variants_table(varTable, keep=True)
 .variant_qc()
 .sample_qc()
 .pca(scores='sa.scores', loadings='va.loadings', as_array=False, k=20))

vds_forPCA.export_samples(scoreFile, 'SAMPLE = s, PC1 = sa.scores.PC1, PC2 = sa.scores.PC2, PC3 = sa.scores.PC3, PC4 = sa.scores.PC4, PC5 = sa.scores.PC5, PC6 = sa.scores.PC6, PC7 = sa.scores.PC7, PC8 = sa.scores.PC8, PC9 = sa.scores.PC9, PC10 = sa.scores.PC10, PC11 = sa.scores.PC11, PC12 = sa.scores.PC12, PC13 = sa.scores.PC13, PC14 = sa.scores.PC14, PC15 = sa.scores.PC15, PC16 = sa.scores.PC16, PC17 = sa.scores.PC17, PC18 = sa.scores.PC18, PC19 = sa.scores.PC19, PC20 = sa.scores.PC20')
vds_forPCA.export_variants(loadingsFile, 'VAR = v, PC1 = va.loadings.PC1, PC2 = va.loadings.PC2, PC3 = va.loadings.PC3, PC4 = va.loadings.PC4, PC5 = va.loadings.PC5, PC6 = va.loadings.PC6, PC7 = va.loadings.PC7, PC8 = va.loadings.PC8, PC9 = va.loadings.PC9, PC10 = va.loadings.PC10, PC11 = va.loadings.PC11, PC12 = va.loadings.PC12, PC13 = va.loadings.PC13, PC14 = va.loadings.PC14, PC15 = va.loadings.PC15, PC16 = va.loadings.PC16, PC17 = va.loadings.PC17, PC18 = va.loadings.PC18, PC19 = va.loadings.PC19, PC20 = va.loadings.PC20')
