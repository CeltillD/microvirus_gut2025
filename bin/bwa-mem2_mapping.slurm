#!/bin/bash
#SBATCH --job-name=mapp_micro
#SBATCH --output=log/mapp_%j.log
#SBATCH --error=log/mapp_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=5G
#SBATCH --partition=normal
#SBATCH --time=24:00:00

module purge
module load rclone
echo "START job --" $(date)
echo "HOSTNAME --" $(hostname)
WORKDIR=/storage/scratch/$USER    # travail dans le scratch
mkdir -p $WORKDIR

GLOBID=$1                                           # ex : PRJ123_VIR_SRR123
PRJ=$(echo ${GLOBID} |cut -d"_" -f1)                # ex : PRJ123                     
TYPE=$(echo ${GLOBID} |cut -d"_" -f2)               # ex : VIR
ID=$(echo ${GLOBID} |cut -d"_" -f3)                 # ex : SRR123
mkdir -p $WORKDIR/${GLOBID}/stdin
mkdir -p $WORKDIR/${GLOBID}/stdout

rclone copy --include "*${ID}*" \
    microv:microv-gut-data/TRIM/${PRJ}/${TYPE}/${ID}/ $WORKDIR/${GLOBID}/stdin --progress && \
singularity exec --bind /storage:/storage $HOME/singularity_imgs/bwa-mem2_2.2.1--he513fc3_0.sif bwa-mem2 mem \
    -t 4 \
    $HOME/Genomovar_GUT8k.fna \
    <(zcat $WORKDIR/${GLOBID}/stdin/TrimmedReads_Fastp_l80W5M30w4_${GLOBID}_1.fastq.gz) \
    <(zcat $WORKDIR/${GLOBID}/stdin/TrimmedReads_Fastp_l80W5M30w4_${GLOBID}_2.fastq.gz) \
    |gzip - > $WORKDIR/${GLOBID}/stdout/mapp_${GLOBID}.sam.gz && \
rclone copy --include "*" $WORKDIR/${GLOBID}/stdout/ microv:microv-gut-data/MAPPING/${PRJ}/${TYPE}/${ID}/. --progress && \
sleep 10 && \
rm -r $WORKDIR/${GLOBID}
echo "END job --" $(date)
