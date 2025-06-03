#!/usr/bin/env python3

############

# NOM      : find_conservative_regions.py 
# DATE     : 03-06-2025
# FONCTION : Détécter les régions conservées dans les alignements de Cap et Rep (AA et NT) et calculer les scores d'identités associés 
# AUTEUR   : Dumont C. (M2BI, BIOADAPT, LMGE)
# VERSION  : Python 3.12.2 [GCC 12.3.0]
# PACAKGES : Bio = 1.85
# NB       : Ce script requiert 3 fichiers au préalable (ex : capsidesAA.msa repAA.msa all_genes.fna)

############

from collections import defaultdict
from Bio import AlignIO
import sys

def define_mask(fileMSA):
    """
    Input : MSA_PROTS(cap ou rep)
    Objectifs :
    -> intégrer les alignements multiples des séquences protéiques
    -> définir un "masque" en fonction du seuil (positions conservées)
    Sorties : nombre de colonnes; MSA des Cap/ou/Rep; masque (vecteur booléen)
    """
    proteinsFMA={}
    alignement       = AlignIO.read(fileMSA, "fasta")			    # lecture du fichier MSA
    nombre_sequences = len(alignement) 
    nbCol_MSA        = alignement.get_alignment_length()            # nombre de positions dans l'alignement global
    mask = []                                                       # vecteur booléen indiquant les positions conservées 

    for i in range(nbCol_MSA):										# Parcourir chaque position de l'alignement
        colonne = alignement[:, i]  								# Récupérer tous les caractères de la colonne
        gaps = colonne.count("-")	                                # Récupération du nombre de gaps				
        freqSSgap = (nombre_sequences - gaps) / nombre_sequences    # ratio de fréquence
        if freqSSgap > SEUIL:										# on repère les colonnes du MSA qui ont moins de x% de gaps
            mask.append(True)
        else:
            mask.append(False)

    with open (fileMSA, "r") as faa :                               # enregistrement dans un dico pour l'alignemennt 
        for line in faa :              
            li = line.strip()
            if li.startswith(">"):
                nom=li.split(" ")[0]
                proteinsFMA[nom]=""
            else:
                proteinsFMA[nom]+=li

    return [nbCol_MSA, proteinsFMA, mask]                           

def conservative_regions(Linfo) :                           # Linfo correspond à la liste d'outpout de la fonction define_mask()
    """
    Input : [nbCol_MSA, proteinsFMA, mask] 
    Objectifs : 
    -> pour chaque protéine détecter les régions conservées sur la base du masque
    -> pour chaque position conservée, récupération du codon NT
    -> constitution des séquences conservées de chaque protéines en AA et NT
    Sorties : 
    -> DconsFaa    : [id] = séquences avec régions conservés en AA
    -> DconsFna    : [id] = séquences avec régions conservés en NT
    -> DglobFna    : [id] = séquences globales en NT
    -> proteinsFMA : [id] = séquences globales en AA
    """
    nbCol_MSA   = Linfo[0]                                  # nombre de positions
    proteinsFMA = Linfo[1]                                  # Dico d'alignement multiple global AA
    mask        = Linfo[2]                                  # vecteur bool des positions conservées
    DconsFna    = {}                                        # Dico d'alignement multiple desn régions conservées NUCL
    DconsFaa    = {}                                        # Dico d'alignement multiple desn régions conservées AA
    DglobFna    = {}                                        # Dico d'alignement multiple global NUCL

    for p in proteinsFMA :                                  # application du masque pour chaque prot
        ConsFaa = ""                                        # séquence avec uniquement les régions conservées AA                    
        ConsFna = ""                                        # séquence avec uniquement les régions conservées NUCL
        GlobFna = ""                                        # séquence d'alignement global basée sur la séquence AA convertie en nucl
        pos = -1                                            # compteur des positions non-vides de la séquence en AA

        for i in range(nbCol_MSA) :                         # parcours des colonnes / positions
            if proteinsFMA[p][i] != "-" :                   # si la position courante (i) n'est pas vide
                pos += 1
                GlobFna += DFNA[p][pos]                     # récupération du triplet nucléotidique ayant formé l'AA courant                    
                if mask[i]:                                 # si la position fait partie des positions conservées
                    ConsFaa += proteinsFMA[p][i]            # ajout de l'acide aminé en cette position      
                    ConsFna += DFNA[p][pos]                 # ajout du triplet nt conservé aynt formé l'AA à la position correspondante
                else:
                    ConsFaa += "-"                          # si la position n'est pas conservée 
                    ConsFna += "---"
            else :
                ConsFaa += "-"                              # si la position est vide (pas d'AA dans la position i de cette séquence)
                ConsFna += "---"
                GlobFna += "---"

        DconsFna[p] = ConsFna                               # sauvegarde des séquences crées dans les dicos correspndants
        DconsFaa[p] = ConsFaa
        DglobFna[p] = GlobFna

    return DconsFaa, DconsFna, DglobFna, proteinsFMA        # retourne les 4 dicos (conservés / global avec AA / NT)

def calc_identity(Dfma) :
    """
    Input    : Dico de séquences
    Objectif : Caculer pour chaque couple de séquences un score d'identité
    Sortie   : Dico avec un score d'identité pour chaque couple
    """
    Dscore = defaultdict(float)                                                 # associe un couple de génome à son score d'identité 
    Lseen  = set()                                                              # set de couples déja parcourus 
	
    for x in Dfma :                                                             # parcours des couples
        for y in Dfma :                                                         
            clause  = False                                                     # si il y a un alignement 
            nbAln   = 0                                                         # nombre de bases alignées
            nbId    = 0                                                         # nombre de bases identiques
            ident   = 0                                                         # score d'identité
            gA      = str(x[:x.rfind("_")] if "_" in x else x).replace(">","")  # récupère uniquement le nom du génome (sans extension _x de la prot)
            gB      = str(y[:y.rfind("_")] if "_" in y else y).replace(">","")  
            cpl     = str(gA+"\t"+gB)                                           # nom du couple de génome
            cpli    = str(gB+"\t"+gA)                                           # nom inversé du couple
			
            if x != y and cpli not in Lseen and cpl not in Lseen:               # si nouveau couple de génomes différents
                seqX = Dfma[x]                                                  # récupération des séquences alignées
                seqY = Dfma[y]
                for i in range(len(Dfma[x])) :                                  # parcours des positions de l'alignement (len(Dfma[x] == len(Dfma[y])
                    if seqX[i] != "-" and seqY[i] != "-" :
                        nbAln += 1                                              # si les deux position ne sont pas vides
                        clause = True                                           
                        if seqX[i] == seqY[i] :                                 # si les deux bases à la meme position sont identiques
                            nbId += 1                                           
							
            if clause == True :                                                 # évite une éventuelle division par 0
                ident = round((nbId * 100) / nbAln, 2)                          # calcul de l'identité
                if ident >= 50 :                                                # limite minimale à 50 pour réduire la taille des comparaisons 
                    Dscore[cpl] = ident                                             
			
            Lseen.add(cpli)                                                     # si cpl est vu alors on ajoute cpli pour ne pas refaire la meme comparaison
				
    return Dscore                                                               # renvoie le dico de score d'identité pour ces couples

if __name__ == "__main__":							# lancement du script
    if len(sys.argv) < 4 :
        print("Ustilisation : ./find_conservative_regions.py <int[0-100]> <capsidesAA.msa> <repAA.msa> <all_genes.fna>")
	
    SEUIL       = int(sys.argv[1])/100                                          # seuil minimal pour déterminer les régions conservées
    CAP_GLOB_FAA= sys.argv[2]                                                   # fichier d'alignement multiple des prot de capside
    REP_GLOB_FAA= sys.argv[3]                                                   # fichier d'alignement multiple des prot de replication
    FILE_FNA    = sys.argv[4]                                                   # fichier fasta des gènes (nucléotides)                                                
                                                                                
    COUPLES = set()                                                             # set de couples de génomes 
    DFNA = {}                                                                   # dico de séquences tq : gène01 -> ["ATG", "CCG", ... , "TTC"]
    ALL_DICOS = []                                                              # liste des dicos de résultats

    with open (FILE_FNA, "r") as fna :
        for line in fna :              
            li = line.strip()
            if li.startswith(">"):
                nom=li.split(" ")[0]
            else:
                DFNA[nom] = [li[i:i+3] for i in range(0, len(li), 3)]           # on sépare la séquence nucl en codons 

    CAP_CONS_FAA, CAP_CONS_FNA, CAP_GLOB_FNA, CAP_GLOB_FAA = conservative_regions( define_mask( CAP_GLOB_FAA ) ) 
    REP_CONS_FAA, REP_CONS_FNA, REP_GLOB_FNA, REP_GLOB_FAA = conservative_regions( define_mask( REP_GLOB_FAA ) )

    CGN = calc_identity(CAP_GLOB_FNA)                                           # calcul des distances pour chaque couple génomes
    del CAP_GLOB_FNA                                                            # suppression de l'ancien dico
    CGA = calc_identity(CAP_GLOB_FAA)
    del CAP_GLOB_FAA
    CCN = calc_identity(CAP_CONS_FNA)
    del CAP_CONS_FNA
    CCA = calc_identity(CAP_CONS_FAA)
    del CAP_CONS_FAA
    RGN = calc_identity(REP_GLOB_FNA)
    del REP_GLOB_FNA
    RGA = calc_identity(REP_GLOB_FAA)
    del REP_GLOB_FAA
    RCN = calc_identity(REP_CONS_FNA)
    del REP_CONS_FNA
    RCA = calc_identity(REP_CONS_FAA)
    del REP_CONS_FAA

    ALL_DICOS = [CGN,CGA,CCN,CCA,RGN,RGA,RCN,RCA]

    for d in ALL_DICOS :
        COUPLES.update(d.keys())                                                # set des couples de génomes

    print(f'GenomeA\tGenomeB\tCGN\tCGA\tCCN\tCCA\tRGN\tRGA\tRCN\tRCA')          # entete du TSV
    for cpl in COUPLES :
        print(cpl, end="")                                                      # nom du couple (genomeA = col1 et genomeB = col2)
        for d in ALL_DICOS :
            print(str("\t"+str(d[cpl])), end="")                                # ajout des scores d'identité pour le couple
        print()
