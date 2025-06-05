# microvirus_gut2025

Ce répertoire git réunit les principaux scripts que j'ai utilisé pour l'élaboration de mon analyse sur les microvirus intestinaux chez l'Homme. Les 1292 viromes assemblés sont issus de prélèvements médicaux de plusieus centaines de personnes dans le monde. Les objectifs de cette méta-analyse sont dans un premier temps de transmettre un vision large de la diversité des microvirus existante dans cet écosystème mais également de décrire les hôtes bactériens qui se trouvent etre infectés par ces bactériophages. Afin de clarifier un maximum le propos, voici les principales étapes réalisées au cours de cette étude :

<p align="center">
  <img src="https://github.com/CeltillD/microvirus_gut2025/blob/main/WFgut25.png" alt="Aperçu graphique" width="600">
</p>

Afin de faciliter la lecture de mes scripts, voici un bref descriptif :
 - **curl_download_fq.slurm** : utilisé pour l'aquisition des FASTQ de reads des différents viromes par FTP (HPC2)
 - **fastp_trimming.slurm** : controle qualité des reads (HPC2)
 - **megahit_assembl.slurm** : assemblages des reads pour former les contigs (HPC2)
 - **identify_microv.py** : retenir les contigs identifiés comme des génomes complets de microvirus
 - **compute_Viridic.py** : calcul des scores de VIRIDIC pour les couples de microvirus
 - **dereplication_microv.py** : déréplicartion des génomes identiques et formation des génomovars
 - **find_conservative_regions.py** : calcul de mesures comparatives sur la base de Cap et Rep
 - **bwa-mem2_mapping.slurm** : execution du mapping des reads contre les contigs/génomovars (HPC2)
 - **reads_counts.slurm** : comptages des reads préférentiellement mappés (HPC2)
 - **norm_RPKM-CLR.py** : normalisation des données de comptages brutes pour les abondances relatives

*Bonus* *:* **regression_logistique.py** *:* *analyse des genres de micros préférentiellement associé à une condition*
*(en cours d'intégration...)*

L'ensemble des principaux outils utilisés sont décris ci-dessous :
| OUTIL        | VERSION         | PARAMÈTRES                                                                 | FONCTION                                                                 |
|--------------|-----------------|----------------------------------------------------------------------------|--------------------------------------------------------------------------|
| singularity  | V4.1.1          | exec --bind                                                               | Exécuter des outils par le biais de containers/images                   |
| rclone       | V1.55.1-DEV     | -                                                                          | Accès au stockage distant depuis HPC                                    |
| curl         | V7.29           | -                                                                          | Téléchargement des FASTQ bruts en lien FTP                              |
| fastp        | V0.20.1         | -l 80 -W 5 -M 30 -w 4 –-cut_right                                           | Contrôle Qualité – nettoyage des reads bruts                            |
| MvSPAdes     | V4.0.0          | --metaviral --only-assembler                                              | Assemblage adapté aux métagénomes                                       |
| penguiN      | v5.cf8933       | guided_nuclassemble                                                       | Assemblage guidé par la prédiction protéique                            |
| MEGAHIT      | v1.2.9          | --presets meta-large                                                      | Assemblage adapté aux métagénomes                                       |
| prodigal     | v2.6.3          | -p meta                                                                   | Prédiction des protéines                                                 |
| mmseqs       | v15-6f452+ds-2  | easy-search                                                               | Recherche de similarité entre génomes / protéines (~blastx)             |
| blastn       | v2.12.0+        | -outfmt 6 -evalue 1 -word_size 7 -reward 2 -penalty -3 -gapopen 5 -gapextend 2 | Alignement local entre génomes (similaire VIRIDIC)                     |
| clustalo     | V1.2.4          | --iterations=5 --full                                                     | Alignement global de séquences protéiques                               |
| MCL          | V14-137         | -I 2.0                                                                    | Markov Clustering Algorithm – clustering avec scores de similarité      |
| genomad      | V1.11.0         | end-to-end --cleanup                                                      | Prédiction des contigs viraux                                           |
| Bwa-mem2     | V2.2.1          | mem                                                                       | Mapping des lectures contre des contigs / génomes de référence          |
| trimal       | v1.5.rev0       | -gapthreshold 0.2                                                         | Nettoyage des alignements multiples avant phylogénie                    |
| iqtree2      | V2.4.0          | -m VT+F+R7 -alrt 1000 -bb 1000                                            | Phylogénie des genres de micro (Cap+Rep) avec modèle ‘Variable Time’    |

Remarques particulières :
- Les calculs ont été effectués sur les ressources du Mésocentre Clermont Auvergne (HPC2)
