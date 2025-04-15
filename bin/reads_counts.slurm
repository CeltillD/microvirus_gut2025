#!/bin/bash
#SBATCH --job-name=counts
#SBATCH --output=log/rcounts_%j.log
#SBATCH --error=log/rcounts_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=100M
#SBATCH --partition=lmge
#SBATCH --time=1:00:00

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
    microv:microv-gut-data/MAPPING/${PRJ}/${TYPE}/${ID}/ $WORKDIR/${GLOBID}/stdin --progress && \
LIB_SIZE=$(grep -vc "^@" <(zcat $WORKDIR/${GLOBID}/stdin/mapp_${GLOBID}.sam.gz)) && \
zcat $WORKDIR/${GLOBID}/stdin/mapp_${GLOBID}.sam.gz \
    |grep -v ^@ \                                                                       # suppression de l'entete
    |awk '$3!="*"&&$5>20' \                                                             # suppression des non-matchs et mauvais mapping (QMAP > 20 == proba erreur < 1%)
    |sort -k1,1 -k5,5nr \                                                               # tri par nom de read et qualité de mapping
    |awk '!seen[$1]++{print $3}' \                                                      # préserve uniquement la première occurence (meilleur QMAP)
    |sort \                                                                             # tri (dispensable)
    |uniq -c \                                                                          # comptage des occurences de chaque génomes
    |awk -v tot=$LIB_SIZE -v name=$ID 'BEGIN{FS=" ";OFS="\t"}{print name,tot,$2,$1}' \  # création TSV tel que : VIROME_NAME / LIB_SIZE / GENOME_TARGET / NB_READS 
    |gzip - > $WORKDIR/${GLOBID}/stdout/reads_counts_$ID.tsv.gz && \
rclone copy --include "*" $WORKDIR/${GLOBID}/stdout/ microv:microv-gut-data/COUNTS/${PRJ}/${TYPE}/${ID}/. --progress && \
cp $WORKDIR/${GLOBID}/stdout/reads_counts_$ID.tsv.gz $HOME/COUNTS/. && \
sleep 5 && \
rm -r $WORKDIR/${GLOBID}
echo "END job --" $(date)
