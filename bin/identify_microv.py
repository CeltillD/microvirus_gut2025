#!/usr/bin/env python3

import sys
import os
import subprocess
import warnings
warnings.filterwarnings("ignore", message="Bio.pairwise2 has been deprecated")
from Bio import pairwise2

WORKDIR,CONTIG_FILE,VIROME_NAME=sys.argv[1],sys.argv[2],sys.argv[3]             # repertoire de travail, fichier fasta et nom du virome
                                   
os.system(f'mkdir -p {WORKDIR}STATS')                                           # répertoire pour enregistrer les comptages pour les viromes

FILES = {'Fasta_S1' : f'{WORKDIR}tmp_{VIROME_NAME}.fna',                        # fna contenant les contigs 2 à 10 kb
         'VP1_S2'   : f'{WORKDIR}tmp_VP1_{VIROME_NAME}.tsv',                    # resultats du blastx contre VP1
         'VP4_S2'   : f'{WORKDIR}tmp_VP4_{VIROME_NAME}.tsv',                    # resultats du blastx contre VP4 
         'Fasta_S2' : f'{WORKDIR}tmp_S2_{VIROME_NAME}.fna',                     # fna avec les contigs circulaires
         'FNA_S3'   : f'{WORKDIR}tmp_genes_{VIROME_NAME}.fna',                  # fna des genes prédits (prodigal meta)
         'FAA_S3'   : f'{WORKDIR}tmp_genes_{VIROME_NAME}.faa',                  # faa des genes prédits (prodigal meta)
         'VP1_S3'   : f'{WORKDIR}tmp_proteins_VS_Capsids_{VIROME_NAME}.tsv',    # resultats du premier blastp contre VP1
         'FINAL_FNA': f'{WORKDIR}complete_{VIROME_NAME}.fna',                   # fna des contigs complets décircularisés
         'FAA_S4'   : f'{WORKDIR}tmp_genes-2_{VIROME_NAME}.faa',                # faa des genes prédits (post réarrangement)
         'FNA_S4'   : f'{WORKDIR}tmp_genes-2_{VIROME_NAME}.fna',                # fna des genes prédits (post réarrangement)
         'VP1_S4'   : f'{WORKDIR}tmp_proteins_VS_Capsids_{VIROME_NAME}_2.tsv',  # resulats du second blastp contre VP1
         'FINAL_CAP': f'{WORKDIR}capsides_{VIROME_NAME}.faa',                   # faa des prots de capsides "_1" = meilleur hit 
         '13k_VP1'  : f'{WORKDIR}CLEAN_13k_VP1.faa',                            # faa de référence VP1 (d'après E. Olo Ndela)
         '13k_VP4'  : f'{WORKDIR}CLEAN_13k_VP4.faa'                             # faa de référence VP4 (d'après E. Olo Ndela)
         }

# Le dico SHELL répertorie les commandes bash sous formes de listes executables avec subprocess.run()
# ex : 'Blastx_VP1' -> mmseqs easy-search ~/tmp_virome01.fna ~/CLEAN_13k_VP1.faa ~/tmp_vp1viome01

SHELL = {
    'Blastx_VP1'    : ['mmseqs',    'easy-search', FILES['Fasta_S1'], FILES['13k_VP1'], FILES['VP1_S2'], f'{WORKDIR}tmp_vp1{VIROME_NAME}'],
    'Blastx_VP4'    : ['mmseqs',    'easy-search', FILES['Fasta_S1'], FILES['13k_VP4'], FILES['VP4_S2'], f'{WORKDIR}tmp_vp4{VIROME_NAME}'],
    'Prodigal_1'    : ['prodigal',  '-i',FILES['Fasta_S2'], '-a', FILES['FAA_S3'] , '-d', FILES['FNA_S3'], '-p' ,'meta'],
    'Blastp_VP1_1'  : ['mmseqs',    'easy-search', FILES['FAA_S3'], FILES['13k_VP1'], FILES['VP1_S3'], f'{WORKDIR}tmp_cap_{VIROME_NAME}', '--format-output', 'query,target,evalue,raw,pident'],
    'Prodigal_2'    : ['prodigal',  '-i', FILES['FINAL_FNA'], '-a', FILES['FAA_S4'], '-d', FILES['FNA_S4'], '-p' ,'meta'],
    'Blastp_VP1_2'  : ['mmseqs',    'easy-search', FILES['FAA_S4'], FILES['13k_VP1'], FILES['VP1_S4'], f'{WORKDIR}tmp_cap2_{VIROME_NAME}', '--format-output', 'query,target,evalue,raw,pident']
}

### Définition des dictionnaires :
Dfna = {}                       # Dico reliant identifiant de contig et sa séquence (2 à 10kb)
Dnew = {}                       # // avec des contigs décircularisés
Dmic = {}                       # // avec les contigs contenant VP1 et VP4
Dscore = {}                     # relie chaque contig à son meilleur score de blastp VP1 
DCapside = {}                   # relie une protéine de capside à son contig (premier blastp VP1)
DCapside_2 = {}                 # // pour le sencond blastp VP1
DBestCap = {}                   # relie un contig avec sa protéine ayant le plus au score contre VP1 

### Définition des scores pour un alignement global entre deux séquences :
match_score = 1                 # Score pour une correspondance
mismatch_score = -2             # Pénalité pour une non-correspondance
gap_open_penalty = -5           # Pénalité pour l'ouverture d'un gap
gap_extend_penalty = -2         # Pénalité pour l'extension d'un gap

### Définition des sets :
cvp1 = set()                    # Contigs possédant la capside (VP1)
cvp4 = set()                    # Contigs possédant la protéine de réplication (VP4)
Capsides = set()                # Protéines blastées contre VP1 (1ère itération)
Capsides2 = set()               # Protéines blastées contre VP1 (2e itération)

### Définition des variables de comptages :
NTOT = 0                        # Nombre de contigs en entrée
NLEN = 0                        # nb de contigs de taille 2 à 10 kb 
NCIR = 0                        # nb de contigs circulaire + 2 à 10kb
NCVP = 0                        # nb de contigs avec VP1 et VP4 + circulaire + 2 à 10 kb
NCAP = 0                        # nb de Capsides prédites "en trop"

### Définition des fonctions :

def ExecShellCommand(file, command):
    """
    Fonction pour controler et executer les commandes shell avec subprocess.run()
    Prends en arguments le fichier input (file) et une liste str (cf SHELL)
    Si le fichier n'est pas présent, nettoyage des fichiers temporaires.
    """
    if os.path.isfile(file):
        try:
            subprocess.run(command, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        except subprocess.CalledProcessError:
            os.system(f'rm -rf {WORKDIR}tmp_*{VIROME_NAME}*')
    else:
        print(f'Le fichier {file} est introuvable')

def ReverseComplement(seq):
    """
    Fonction qui retourne la complémentarité d'une séquence en entrée (str)
    """
    complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    Revseq=""
    for nuc in seq:
        Revseq += complement_dict[nuc]
        NewSeq = Revseq[::-1]
    return NewSeq

def step_1() :
    """
    Fonction qui a pour but de créer un FNA avec des séquences décircularisées :
    -> contrôle de la taille [ 2kb - 10kb ]
    -> contrôle de la circularité par alignement global [ debut - fin ]
    -> suppression de la partie redondante ( décircularisation )
    Retourne les comptages pour les contigs : totaux, de bonne taille et circulaires  
    """

    with open(CONTIG_FILE) as f1s1 :
        global NTOT
        global NLEN 
        global NCIR
        id = None
        seq = []
        seq_complete = ""
        for ligne in f1s1:
            ligne = ligne.strip()
            if ligne.startswith(">"):                       # >k127_0 flag=0 multi=9.0 len=225 
                NTOT += 1
                l = ligne.split()
                if id is not None:  
                    seq_complete = "".join(seq)
                    if 2000 <= len(seq_complete) <= 10000 : # critère de taille
                        Dfna[id] = seq_complete
                        NLEN += 1
                id = VIROME_NAME + "_" + str(l[0].replace(">","").split("_")[1]) + "_" + str(len(seq_complete))
                seq = []  
            else:
                seq.append(ligne)  
        if id is not None:                                  # pour la dernière séquence du fasta
            seq_complete = "".join(seq)
            if 2000 <= len(seq_complete) <= 10000:
                Dfna[id] = seq_complete

    with open(FILES['Fasta_S1'],'w') as f2s1 :
        for Contig,Seq in Dfna.items() :
            Motif = Seq[0:10]  		        
            Index = 0
            Index = Seq.rfind(Motif)		            # dernière position du motif		
            if Index != 0:				    
                finContig = Seq[Index:] 				# fragement 2	
                if len(finContig) > 10 and len(finContig) < 1000 and len(Seq) > 2000 :  
                    debContig = Seq[0:len(finContig)]   # fragment 1
                    aln = pairwise2.align.globalms(     # alignement global 
                        debContig,
                        finContig,
                        match_score,
                        mismatch_score,
                        gap_open_penalty,
                        gap_extend_penalty
                        )
                    # on compte ici le nombre de matchs "|" sur l'alignement global :
                    #   ATCGCATCAATCGCC
                    #   ||||||    || ||
                    #   ATCGCAAACGTCTCC
                    ratio = pairwise2.format_alignment(*aln[0]).count('|') / len(debContig)
                    # le ratio est égal au nombre de bases alignées sur la longueur des fragments 
                    if ratio >= 0.95 :
                        NCIR += 1 
                        decircular=Seq[:Index]                                      # séquence décircularisée
                        Dnew[Contig] = decircular
                        print(">"+Contig,decircular,sep="\n", file = f2s1)

    return [str(NTOT),str(NLEN),str(NCIR)]

def step_2() :
    """
    Fonction qui compose un FNA avec les contigs possèdant VP1 et VP4 :
    -> recherche des contigs possedant VP1 (blastx score > 50)
    -> recherche des contigs possedant VP4 (blastx score > 50)
    -> composition d'un dico (Dmic) avec les contigs complets (avec VP1 et VP4)
    -> FNA avec ces séquences + ajout d'une extension (2000nt au debut du contig) à la fin pour ne pas manquer de protéines lors de la prédiction avec Prodigal
    Renvoie le nombre de contigs complets 
    """

    global Dmic
    global NCVP

    if os.path.isfile(FILES['VP1_S2']) :
        with open (FILES['VP1_S2']) as vp1 :
            for l in vp1 :
                l=l.rstrip().split("\t")
                if int(l[11]) >= 50 :
                    cvp1.add(str(l[0]))             # id des contigs possedant VP1
    else :
        print(f'Le fichier {FILES['VP1_S2']} est absent (pas de match contre VP1)')
        os.system(f'rm -rf {WORKDIR}tmp_*{VIROME_NAME}*')

    if os.path.isfile(FILES['VP4_S2']) :
            with open (WORKDIR+'tmp_VP4_'+VIROME_NAME+'.tsv') as vp4 :
                for l in vp4 :
                    l=l.rstrip().split("\t")
                    if int(l[11]) >= 50 :
                        cvp4.add(str(l[0]))         # id des contigs ayant VP4
    else :
        print(f'Le fichier {FILES['VP4_S2']} est absent (pas de match contre VP4)')
        os.system(f'rm -rf {WORKDIR}tmp_*{VIROME_NAME}*')

    ContigsComplets = cvp4.intersection(cvp1)       # id des contigs complets
    NCVP = len(ContigsComplets)
    Dmic = {key: value for key, value in Dnew.items() if key in ContigsComplets}

    with open(FILES['Fasta_S2'],'w') as f1s2 :      # FNA avec extension de 2000nt à la fin (str(s[0:2000]))
        for k,s in Dmic.items() :           
            print(">"+k,s+str(s[0:2000]),sep="\n", file = f1s2) 

    return str(NCVP)

def step_3():
    """
    Fonction qui a pour but de créer le FNA des contigs complets finalisé :
    -> identification des potentielles protéines de capside (VP1) pour chaque contig 
    -> association de chaque contig avec sa proteine de capside associée (meilleur score par contig)
    -> récupération des séquences nt des capsides pour chaque contig
    -> réarrangement du contig en fonction du brin sens de la protéine de capside (reverse complément si nécéssaire)
    -> réécriture de la séquence du contig en démarrant par la protéine de capside 
    -> suppression de l'extention de 2000nt afin d'avoir un contig décircularisé
    -> ajouts d'éléments descriptifs du contigs dans l'entete (taille contig, taille capside, GC%...)
    Chaque contig enregistré est composé d'une capside en début de sa séquence
    """
    DSeqCap = {}        # associe id de la prot à une liste [séquence, posdébut, posfin, sens]
    Dsens = {}          # associe un id de contig au sens de sa capside
    id = 0
    if os.path.isfile(FILES['VP1_S3']) :
        with open(FILES['VP1_S3'],'r') as f1s3 :
            for l in f1s3 :
                l=l.rstrip().split("\t")
                if int(l[3]) > 50 :
                    Contig = "_".join(l[0].split("_")[:-1]) # identifiant du contig
                    IdProt = l[0]                           # identifiant de la protéine 
                    Dscore[IdProt] = int(l[3])
                    if Contig not in Dscore.keys() :
                        Dscore[Contig] = 0                  # initialisation du score max
                    if int(l[3]) > Dscore[Contig] :         # si la proteine courante dépasse le score max du contig
                        Dscore[Contig] = int(l[3])          # cette protéine devient la protéine de capside ayant le meilleur score
                        DBestCap[Contig] = IdProt
                        DCapside[IdProt] = Contig
                    Capsides.add(IdProt)                    
    else :
        print(f'Le fichier {FILES["VP1_S3"]} est absent')
        os.system(f'rm -rf {WORKDIR}tmp_*{VIROME_NAME}*')

    if os.path.isfile(FILES['FNA_S3']) :
        with open(FILES['FNA_S3'],'r') as f2s3 :
            for l in f2s3 :
                l=l.rstrip()
                if l.startswith(">") : # >CONTIG_NOR_43957_CD_2763_5513_17 # 6191 # 6289 # 1 # ;gc_cont=0.475
                    if id in Capsides:
                        DSeqCap[id]=[seq,deb,fin,sens]
                    l=l.replace(" ","").split("#")
                    id = l[0].replace(">","")
                    deb = int(l[1])
                    fin = int(l[2])
                    sens = l[3]
                    seq = ""
                    Dsens['_'.join(id.split('_')[:-1])] = sens   
                else :
                    seq += l
    else :
        print(f'Le fichier {FILES['FNA_S3']} est absent')
        os.system(f'rm -rf {WORKDIR}tmp_*{VIROME_NAME}*')

    with open(FILES['FINAL_FNA'],'w') as f3s3 :
        for cap,info in DSeqCap.items() :
            if cap in set(DBestCap.values()) :
                Contig = DCapside[cap]                      # on retient le contig lié à la protéine courante
                SeqBase = Dmic[Contig]                      # séquence décircularisée du contig
                Ext = str(SeqBase[0:2000])                  # extension de 2000nt 
                SeqLong = SeqBase + Ext                     # ajout de cette extension à la fin de la séquence

                if info[3] != "1" :                         # reverse si la capside est sur le brin negatif
                    FwdSeqLong = ReverseComplement(SeqLong)
                    FwdExt = ReverseComplement(Ext)
                else :
                    FwdSeqLong = SeqLong
                    FwdExt = Ext
                
                # recherche de la position de la capside dans la séquence
                IndCap = FwdSeqLong.find(info[0])           
                # rearrangement : capside au debut + suppresion de l'extension
                FinalContig = str(FwdSeqLong[IndCap:] + FwdSeqLong[:IndCap]).replace(FwdExt,"",1)
                # Comptage du GC % 
                GCprc = (FinalContig.count("C") + FinalContig.count("G")) / len(FinalContig)
                GC = str(round(GCprc,3)).replace("0.","")
                # Enregistrement de la taille du contig
                LenCap = str(info[2] - info[1] + 1 )
                # Composition de la nouvelle entete
                IdContig = "_".join(Contig.split("_")[:-2])
                NumContig = str(Contig.split("_")[-2])
                LenContig = str(len(FinalContig))
                Chevron = f'>{IdContig}_c{NumContig}_VP1-4_Cir{LenContig}_GC{GC}_Cap{LenCap}'
                print(Chevron, FinalContig, sep="\n", file = f3s3)

def step_4():
    """
    Fonction servant à identifier toutes les protéines appariées à VP1 par blastx :
    Les protéines ayant pour suffixes "_1" sont les protéines de capsides (meilleurs scores)
    Les éventuelles protéines supplémentaires sont enregistrés dans un faa avec les capsides principales
    Renvoie le nombre de proteines liées à VP1 par blastp pour l'ensemble des contigs
    """

    global NCAP

    if os.path.isfile(FILES['VP1_S4']) :
        with open(FILES['VP1_S4'],'r') as f1s4 :
            for l in f1s4 :
                l=l.rstrip().split("\t")
                if int(l[3]) > 50 :
                    Capsides2.add(l[0])
        NCAP = len(Capsides2)
    else :
        print(f'Le fichier {FILES['VP1_S4']} est absent')
        os.system(f'rm -rf {WORKDIR}tmp_*{VIROME_NAME}*')
    
    if os.path.isfile(FILES['FAA_S4']) :
        with open(FILES['FINAL_CAP'], 'w') as f3s4 :
            with open(FILES['FAA_S4'],'r') as f2s4 :
                for l in f2s4 :
                    l = l.rstrip()
                    if l.startswith(">") :
                        clause = 0
                        Id = l.split()[0].replace(">","")
                        if Id in Capsides2 :
                            clause = 1
                            print(l, file = f3s4)
                    else :
                        if clause == 1 :
                            print(l, file = f3s4)
    else :
        print(f'Le fichier {FILES['FAA_S4']} est absent')
        os.system(f'rm -rf {WORKDIR}tmp_*{VIROME_NAME}*')
    
    return str(NCAP)

def run_all() :

    print(f'START - {VIROME_NAME}')
    
    counts = step_1()                                           # création FNA contigs 2-10kb et decircularises
    ExecShellCommand(FILES['Fasta_S1'], SHELL['Blastx_VP1'])    # Blastx de ces contigs contre VP1.faa
    ExecShellCommand(FILES['Fasta_S1'], SHELL['Blastx_VP4'])    # Blastx de ces contigs contre VP4.faa
    counts.append(step_2())                                     # création FNA des contigs complets
    ExecShellCommand(FILES['Fasta_S2'], SHELL['Prodigal_1'])    # prédiction des protéines du FNA précédent
    ExecShellCommand(FILES['FAA_S3'], SHELL['Blastp_VP1_1'])    # Blastp pour identifier les capsides potentielles (VP1)
    step_3()                                                    # création du FNA des contigs finalisés
    ExecShellCommand(FILES['FINAL_FNA'], SHELL['Prodigal_2'])   # prédiction des protéines après réarrangement
    ExecShellCommand(FILES['FAA_S4'], SHELL['Blastp_VP1_2'])    # second blastp pour capter les protéines liées à VP1
    counts.append(step_4())                                     # création FAA pour les protéines de capside

    with open(f'{WORKDIR}STATS/{VIROME_NAME}.tsv','w') as f :   # Sauvegarde des comptages
        print(VIROME_NAME,"\t".join(counts),sep = "\t",file=f)  
        # virome01  nb_total    nb_2-10kb   nb_cir  nb_complets nb_capsides

    # NETTOYAGE FICHIERS ET REPERTOIRES

    os.system(f'rm -rf {WORKDIR}tmp_*{VIROME_NAME}*')
    
    print(f'END - {VIROME_NAME}')

run_all()
