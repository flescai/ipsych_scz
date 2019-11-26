#!/bin/sh
#SBATCH --nodes 30
#SBATCH -c 16
#SBATCH --time 11:59:00
#SBATCH --mem 128g
#SBATCH --job-name hailMulti


# without this, srun will do the batch script once per node
[ $SLURM_PROCID -ne 0 ] && exit 0


############## Prepare software #########################################
source /com/extra/java/8/load.sh
source /com/extra/Anaconda-Python/2.2.0-2.7/load.sh
source /com/extra/spark-cluster/2.1.0/load.sh  
source /com/extra/gcc/5.2.0/load.sh

export SPARK_HOME=/com/extra/spark-cluster/2.1.0
export HAIL_HOME=/com/extra/hail/2017-06-30/


############## Compute ##################################################

source /com/extra/spark-cluster/2.1.0/bin/spark-cluster-start.sh
set -x
# alias hail="PYTHONPATH=$SPARK_HOME/python:$SPARK_HOME/python/lib/py4j-0.10.4-src.zip:$HAIL_HOME/python SPARK_CLASSPATH=$HAIL_HOME/build/libs/hail-all-spark.jar python"

PYTHONPATH="$HAIL_HOME/python:$HAIL_HOME/python/lib/pyhail.zip:$SPARK_HOME/python:`ls $SPARK_HOME/python/lib/py4j-*-src.zip`" \
   SPARK_CLASSPATH=$HAIL_HOME/build/libs/hail-all-spark.jar \
   spark-submit \
   --conf spark.shuffle.compress=true \
   --conf spark.kryoserializer.buffer.max=512m \
  --conf spark.driver.memory=128G \
   --conf spark.akka.frameSize=1024 \
   --conf spark.driver.maxResultSize=20G \
   --conf spark.rdd.compress=true \
   --conf spark.ui.port=54054 \
   --conf spark.cores.max=480 \
   --executor-memory 120G \
   --total-executor-cores 480 \
   "$@"
    
source /com/extra/spark-cluster/2.1.0/bin/spark-cluster-stop.sh

