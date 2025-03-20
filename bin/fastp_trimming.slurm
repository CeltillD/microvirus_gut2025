#!/bin/bash
#SBATCH --job-name=fastp_trim
#SBATCH --output=log/trim_%j.log
#SBATCH --error=log/trim_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=5G
#SBATCH --partition=normal
#SBATCH --time=4:00:00

#######################################################################################

#   AUTHOR :   Celtill DUMONT (BIOADAPT,LMGE)
#   DATE    :   14 - 02 - 2025
#   SCRIPT  :   fastp_trimming.slurm
#   TOOLS   :   rclone (v1.55.1-DEV), fastp (v0.20.1)

#######################################################################################

# sbatch M2_celtill/bin/01_Fastp_trimming_CD_14-02-25.sh

#######################################################################################

module purge
module load rclone
echo "START job --" $(date)
echo "HOSTNAME --" $(hostname)
WORKDIR=/storage/scratch/$USER    # travail dans le scratch
mkdir -p $WORKDIR

#NOMBRE ELEMENTS A TRAITER
N=$(rclone ls microv: |grep "RawReads_Curl_x_" |cut -d'_' -f4,5,6 |sort -u |wc -l)

#DEFINITION DE L'ARRAY
if [ -z "$SLURM_ARRAY_TASK_COUNT" ]; then
    sbatch --array=0-$(($N-1)) $0 $@
    exit 0
fi

#VARIABLES
ALLRAW=($(rclone ls microv: |grep "RawReads_Curl_x_" |cut -d'_' -f4,5,6 |sort -u))
GLOBID=${ALLRAW[$SLURM_ARRAY_TASK_ID]}              # ex : PRJ123_VIR_SRR123
PRJ=$(echo ${GLOBID} |cut -d"_" -f1)                # ex : PRJ123                     
TYPE=$(echo ${GLOBID} |cut -d"_" -f2)               # ex : VIR
ID=$(echo ${GLOBID} |cut -d"_" -f3)                 # ex : SRR123

#COPIE DES READS BRUTS + EXECUTION DU TRIMMING
rm -r $WORKDIR/${GLOBID}
mkdir -p $WORKDIR/${GLOBID}/stdin
mkdir -p $WORKDIR/${GLOBID}/stdout
rclone copy --include "*${GLOBID}*" \
    microv:microv-gut-data/RAW/${PRJ}/${TYPE}/${ID}/ $WORKDIR/${GLOBID}/stdin && \
    ls -l $WORKDIR/${GLOBID}/stdin/ && \
singularity exec --bind /storage:/storage $HOME/singularity_imgs/fastp_v0.20.1.sif fastp \
    -i $WORKDIR/${GLOBID}/stdin/RawReads_Curl_x_${GLOBID}_1.fastq.gz \
    -o $WORKDIR/${GLOBID}/stdout/TrimmedReads_Fastp_l80W5M30w4_${GLOBID}_1.fastq.gz \
    -I $WORKDIR/${GLOBID}/stdin/RawReads_Curl_x_${GLOBID}_2.fastq.gz \
    -O $WORKDIR/${GLOBID}/stdout/TrimmedReads_Fastp_l80W5M30w4_${GLOBID}_2.fastq.gz \
    -l 80 \        # taille minimum fixée à 80nt
    -W 5 \         # taille de la fenetre à 5nt
    -M 30 \        # scrore PHRED de 30 exigé 
    -w 4 \         # décalage par pas de 4nt
    --cut_right \  # décalage de la fenetre glissante de gauche à droite (deb-fin)
    --unpaired1 $WORKDIR/${GLOBID}/stdout/TrimmedReads_Fastp_l80W5M30w4_${GLOBID}_unpaired1.fastq.gz \
    --unpaired2 $WORKDIR/${GLOBID}/stdout/TrimmedReads_Fastp_l80W5M30w4_${GLOBID}_unpaired2.fastq.gz && \
rclone copy --include "*" $WORKDIR/${GLOBID}/stdout/ microv:microv-gut-data/TRIM/${PRJ}/${TYPE}/${ID}/. --progress --ignore-checksum && \
sleep 30 && \
rm -r $WORKDIR/${GLOBID}
echo "END job --" $(date)
