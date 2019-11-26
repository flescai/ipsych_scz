#!/bin/sh
#SBATCH --nodes 1
#SBATCH -c 16
#SBATCH --time 1:00:00
#SBATCH --mem 32g
#SBATCH --job-name hail-HWE

## this script is made to export variants in PLINK format
## intended to run SNP pruning for selecting variants to be used
## in PCA analysis

# without this, srun will do the batch script once per node
[ $SLURM_PROCID -ne 0 ] && exit 0

DIR=$1
INPUTFILE=$2
OUTFILE=${INPUTFILE%vds}hwe.tsv

###########################################

export _JAVA_OPTIONS="${_JAVA_OPTIONS} -Djava.io.tmpdir=${DIR}/tmp"
export SPARK_DAEMON_JAVA_OPTS="${_JAVA_OPTIONS} -Djava.io.tmpdir=${DIR}/tmp"


############## Prepare software #########################################
source /com/extra/java/8/load.sh
source /com/extra/shtools/1.0/load.sh
source /com/extra/python/2.7/load.sh
source /com/extra/spark-cluster/1.6.1/load.sh
source /com/extra/hail/2016-09-20/load.sh

############## Compute ##################################################

source spark-cluster-start.sh

## this is filtering samples (if passed) and variants based on QC first
## then based on allele number, call rate and MAF
## in order to perform pruning in plink

spark-submit \
--total-executor-cores 16 \
--executor-memory 24G \
/com/extra/hail/2016-09-20/jar-bin/hail-all-spark.jar \
--tmpdir ${DIR}/tmp \
read -i ${DIR}/${INPUTFILE} \
count \
repartition -n 25 \
exportvariants -o ${DIR}/${OUTFILE} -c 'v, va.qc.AF, va.qc.nHomRef, va.qc.nHet, va.qc.nHomVar, va.qc.rHeterozygosity, va.qc.rHetHomVar, va.qc.rExpectedHetFrequency, va.qc.pHWE'

spark-cluster-stop.sh

### adding the header to the file, which misses from the export
sed -i '1i\variant\tMAF\tnHomRef\tnHet\tnHomVar\trHeterozygosity\trHetHomVar\trExpectedHetFrequency\tpHWE' ${DIR}/${OUTFILE}

echo "ALL DONE"
