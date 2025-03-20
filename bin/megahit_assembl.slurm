#!/bin/bash
#SBATCH --job-name=megahit
#SBATCH --output=log/MEGAHIT_%j.log
#SBATCH --error=log/MEGAHIT_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=10G
#SBATCH --partition=normal,lmge
#SBATCH --time=24:00:00

#######################################################################################

#   AUTHORS :   Celtill DUMONT (BIOADAPT,LMGE)
#   DATE    :   14 - 02 - 2025
#   SCRIPT  :   megahit_assembl.slurm
#   TOOLS   :   rclone (v1.55.1-DEV), megahit (v1.2.9) 

#######################################################################################

# sbatch M2_celtill/bin/megahit_assembl.slurm

####################################

module purge 
module load gcc/8.1.0
module load megahit/1.2.9 
module load rclone
echo "START job --" $(date)
echo "HOSTNAME --" $(hostname)
WORKDIR=/storage/scratch/$USER    # travail dans le scratch
mkdir -p $WORKDIR

#NOMBRE ELEMENTS A TRAITER
N=$(rclone ls microv: |grep -P "PRJNA588514|PRJNA804384" |grep "TrimmedReads_Fastp_l80W5M30w4_" |cut -d'_' -f4,5,6 |sort -u |wc -l)

#DEFINITION DE L'ARRAY
if [ -z "$SLURM_ARRAY_TASK_COUNT" ]; then
    sbatch --array=0-$(($N-1)) $0 $@
    exit 0
fi

#VARIABLES
ALLTRIM=($(rclone ls microv:  |grep -P "PRJNA588514|PRJNA804384" |grep "TrimmedReads_Fastp_l80W5M30w4_" |cut -d'_' -f4,5,6 |sort -u))
GLOBID=${ALLTRIM[$SLURM_ARRAY_TASK_ID]}             # ex : PRJ123_VIR_SRR123
PRJ=$(echo ${GLOBID} |cut -d"_" -f1)                # ex : PRJ123                     
TYPE=$(echo ${GLOBID} |cut -d"_" -f2)               # ex : VIR
ID=$(echo ${GLOBID} |cut -d"_" -f3)                 # ex : SRR123
rm -rf $WORKDIR/${GLOBID}

#COPIE DES READS TRIMMED + EXECUTION DU MAPPING
mkdir -p $WORKDIR/${GLOBID}/stdin
rclone copy --include "*${GLOBID}*" \
    microv:microv-gut-data/TRIM/${PRJ}/${TYPE}/${ID}/ $WORKDIR/${GLOBID}/stdin && \
megahit --presets meta-large \
    -1 $WORKDIR/${GLOBID}/stdin/TrimmedReads_Fastp_l80W5M30w4_${GLOBID}_1.fastq.gz \
    -2 $WORKDIR/${GLOBID}/stdin/TrimmedReads_Fastp_l80W5M30w4_${GLOBID}_2.fastq.gz \
    -r <(zcat $WORKDIR/${GLOBID}/stdin/TrimmedReads_Fastp_l80W5M30w4_${GLOBID}_unpaired*) \
    -o $WORKDIR/${GLOBID}/stdout && \
mv $WORKDIR/${GLOBID}/stdout/final.contigs.fa \
    $WORKDIR/${GLOBID}/stdout/Contigs_MEGAHIT_metalarge_${GLOBID}.fna && \
gzip $WORKDIR/${GLOBID}/stdout/*${GLOBID}.fna && \
rclone copy --include "*${GLOBID}*" $WORKDIR/${GLOBID}/stdout/ \
    microv:microv-gut-data/CONTIGS/${PRJ}/${TYPE}/${ID}/. --progress --ignore-checksum && \
sleep 30 && \
rm -rf $WORKDIR/${GLOBID}

echo "END job --" $(date)
