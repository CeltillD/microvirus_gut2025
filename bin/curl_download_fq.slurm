#!/bin/bash
#SBATCH --job-name=download_FTP
#SBATCH --output=log/download_FTP_%j.log
#SBATCH --error=log/download_FTP_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=1G
#SBATCH --partition=normal
#SBATCH --time=1:00:00


#######################################################################################

#   AUTHORS :   Celtill DUMONT (BIOADAPT,LMGE)
#   DATE    :   13 - 02 - 2025
#   SCRIPT  :   curl_download_fq.slurm
#   TOOLS   :   rclone (v1.55.1-DEV), curl (7.29.0)

#######################################################################################

# sbatch $HOME/M2_celtill/bin/Curl_Download_CD_13-02-25.slurm "$HOME/metadata.csv"

# EX : "$HOME/metadata.csv"
#------------------------------------
#VIR;PRJNA564995;SRR10099602;ftp.sra.ebi.ac.uk/vol1/fastq/SRR100/002/SRR10099602/SRR10099602_1.fastq.gz;ftp.sra.ebi.ac.uk/vol1/fastq/SRR100/002/SRR10099602/SRR10099602_2.fastq.gz
#VIR;PRJNA564995;SRR10099679;ftp.sra.ebi.ac.uk/vol1/fastq/SRR100/079/SRR10099679/SRR10099679_1.fastq.gz;ftp.sra.ebi.ac.uk/vol1/fastq/SRR100/079/SRR10099679/SRR10099679_2.fastq.gz

####################################

module purge
module load rclone
echo "START job --" $(date)
echo "HOSTNAME --" $(hostname)
WORKDIR=/storage/scratch/$USER         # travail dans le scratch
mkdir -p $WORKDIR

####################################

#BOUCLE CONTROLE
if [ -z "$1" ]; then
    echo "Erreur : Argument non-valide" >&2
    exit 1
fi

#NOMBRE ELEMENTS A TRAITER
N=$(cat "$1" |wc -l)                                # 2 liens par lignes ; ex : 6 lignes

#DEFINITION DE L'ARRAY
if [ -z "$SLURM_ARRAY_TASK_COUNT" ]; then
    sbatch --array=0-$(($N+$N-1)) $0 $@             # 0 à (N*2-1) ; ex 0-11 (pour 12 liens)
    exit 0
fi

#DEPLACEMENT DANS LE SCRATCH
echo "cd to $WORKDIR"
cd $WORKDIR

#ATTRIBUTION DES VARIABLES

LINKS=($(cat "$1"|cut -d";" -f4,5 |grep .|tr '\n' ';' |sed 's/;/ /g'))   # array des liens FTP
link=${LINKS[$SLURM_ARRAY_TASK_ID]}                                      # lien à télécharger
ID=$(cat "$1" |grep -w "${link}" |cut -d";" -f3 |sort -u)                # ex : SRR12345 / ERR1345
fwdrev=$(echo ${link}  |cut -d '_' -f2)                                  # ex : _1.fastq.gz / _2.fastq.gz
PRJ=$(cat "$1" |grep -w "${link}" |cut -d";" -f2 |sort -u)               # ex : PRJNA12345
TYPE=$(cat "$1" |grep -w "${link}" |cut -d";" -f1 | sort -u)             # ex : VIR / 16S

mkdir -p $WORKDIR/${ID}/${fwdrev}

echo ${link} ${PRJ} ${TYPE} $WORKDIR/${ID}/${fwdrev}/RawReads_Curl_x_${PRJ}\_${TYPE}\_${ID}_${fwdrev}

#TELECHARGEMENT DES LIENS ET COPIE SUR BUCKET
curl -o $WORKDIR/${ID}/${fwdrev}/RawReads_Curl_x_${PRJ}\_${TYPE}\_${ID}_${fwdrev} ftp://${link} && \
rclone copy --include "*${ID}*fastq.gz" $WORKDIR/${ID}/${fwdrev} microv:microv-gut-data/RAW/${PRJ}/${TYPE}/${ID}/. --progress --ignore-checksum && \
sleep 30 && \
rm -r $WORKDIR/${ID}/${fwdrev}

echo "END job --" $(date)
