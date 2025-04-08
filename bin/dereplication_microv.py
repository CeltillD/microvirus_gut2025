#!/usr/bin/env python3

import os
import sys
import pandas as pd
from collections import defaultdict

### Initialisation des variables ###

OUTDIR,VIRIDIC_TSV,FNA,PRO = sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4]

### Définition des Fonctions ###

def read_fasta(file_fasta) :
    """
    INPUT : fichier au format fasta
    |
    -> simplification des identifiants de génomes (DnewId {former_id : new_id})
    -> lecture et enregistrement des séquences dans un dictionnaire (Dfna {new_id : seq})
    |
    OUTPUT : dico des identifiants & dico des séquences
    """
    Dfna    = defaultdict(str)
    DnewId  = {}
    ngut    = 0 
    navrc   = 0
    with open(file_fasta) as f1 :                                                      
        for l in f1 : 
            l = l.rstrip()
            if l.startswith(">") : 
                g = l[1:]
                if g.startswith("GUT") :                                            # priorité aux génomes du gut (2025)
                    ngut  += 1                                                      # GUT25_ZUB_82326_HC_c158713_VP1-4_Cir6327_GC366_Cap1884
                    newId  = "G_"+"_".join(g.split("_")[1:4])+"_c"+str(ngut)        # new ID : G_ZUB_82326_HC_c1                                                                  
                    DnewId[g] = newId 

                if g.startswith("AVrC"):                                            # idem pour les génomes de AVrC 
                    navrc += 1
                    newId  = "A_"+"_".join(g.split("_")[1:2])+"_c"+str(navrc)      # AVrC_GutCatV1_XX_c3442_VP1-4_Cir5528_GC418_Cap2046
                    DnewId[g] = newId                                                 # new ID : A_GutCatV1_c3442

                if g.startswith("REF"):                                             # REF13k_EAF_DB_IMGmg_3300010885-c3
                    newId  = g.replace("REF13k","R")                                # new ID : R_EAF_DB_IMGmg_3300010885-c3
                    DnewId[g] = newId
            else :
                Dfna[newId] += l
    return Dfna, DnewId

def prob_genomes(df,DnewId):
    """
    INPUT : dataframe de resultats VIRIDIC (df) et dico des nouveaux identifiants de génome (DnewId)
    |
    -> identifie les prophages et lie les génomes apparentés dans un dico (Dpro)
    -> relève les couples de génomes ayant des critères anormaux (forte identité mais taille très différente)
    -> pour chaque génome on relève son nombre d'occurence dans les couples problématiques (Docc)
    -> pour chaque couple problématique, j'ajoute le génome ayant le plus d'occurence des deux à un set de génomes anormaux (probGenomes)
    |
    OUTPUT : renvoie le set d'identifiants de génomes anormaux (probGenomes)
    """
    Dpro = {}
    Docc = {}
    probGenomes = set()

    with open(PRO,'r') as pro :
        for l in pro :
            l = l.rstrip()
            Dpro["REF13k_"+l] = []
            probGenomes.add("REF13k_"+l)
    
    for _, row in df.iterrows():
        if row['gA'] in Dpro.keys():
            Dpro[row['gA']].append(row['gB'])
        elif row['gB'] in Dpro.keys():
            Dpro[row['gB']].append(row['gA'])

    df_prob = df[
        (df["ANI_AB"] > 70) &                                               # Si ANI entre A -> B élevé :
        ((df["ratioLen"] < 0.70))
    ]
    probPotentiels = pd.concat([df_prob["gA"], df_prob["gB"]]).unique()         # liste des génomes potentiellement problématiques
    for pp in probPotentiels :                                                  # comptage de leur nombre d'occurences problématiques
        Docc[pp] = df_prob['gA'].str.count(pp).sum() + df_prob['gB'].str.count(pp).sum()
    Dprob = df_prob.set_index('gA')['gB'].to_dict()
    for gA, gB in Dprob.items():                                                # pour chaque couple de génomes problématiques, je m'affranchis de celui qui a le plus d'occurences
        if Docc[gA] > Docc[gB] :
            probGenomes.add(DnewId[gA])
        elif  Docc[gA] < Docc[gB] :
            probGenomes.add(DnewId[gB])
        else :                                                                    # si meme nombre d'occurences
            probGenomes.add(DnewId[gA])
            probGenomes.add(DnewId[gB])
    return probGenomes,Dpro

def single_linkage(clusters):                                               # Clustering Single Linkage        
    """
    INPUT : liste de couples, dont chaque couple est une liste de 2 identifiants (clusters) : [ [gA,gB], [gA,gC], [gB,gC] ... ]
    |
    -> algorithme agglomératif qui unie les couples si au moins un des deux identifiants est retrouvé dans un autre couple
    -> tant que la taille des clusters n'est pas fixe on ré-execute le processus
    -> si la taille des clusters est identique à l'itération précédente, on arrete (tous les id sont présent de manière unique dans un cluster) 
    |
    OUTPUT : renvoie une liste de sets (new_clusters), dont chaque set est un cluster de génomes similaires
    """                                   
    while True:                                                             # agglomérer les couples de génomes identiques
        new_clusters = []                           
        merged = set()                                                      
        for i, cluster in enumerate(clusters):
            if i in merged:
                continue
            new_cluster = set(cluster)
            for j, other_cluster in enumerate(clusters[i+1:], start=i+1):
                if j in merged:
                    continue
                if new_cluster & set(other_cluster):
                    new_cluster |= set(other_cluster)
                    merged.add(j)
            new_clusters.append(new_cluster)
        if len(new_clusters) == len(clusters):                              # taille finale des clusters
            return new_clusters
        clusters = new_clusters

def dereplicate(df, DnewId):
    """
    INPUT : dataframe filtrée (df) et dico des identifiants (DnewId)
    |
    -> établie une liste de couples de génomes [ [gA,gB], [gA,gC], [gB,gC] ... ]
    -> éxecution du single linkage sur cette liste (couples)
    -> tri de chaque cluster de façon établir un génome représentatif (issu du GUT de préférence) pour chaque cluster (Dclust { Gref : liste_de_genomes_similaires })
    |
    OUTPUT : renvoie le Dclust qui lie un génome représentatif à son cluster
    """
    couples = []
    Dclust  = {}
    for _, row in df.iterrows():                                       # chaque couple étant identique
        coupleAB = [ DnewId[str(row["gA"])], DnewId[str(row["gB"])] ]
        couples.append(coupleAB)

    clusters = single_linkage(couples)                                          # application du single linkage sur ces couples

    for c in clusters :                                                   # séléction des génomes représentatifs
            c_sorted = sorted(c, key=lambda x: (not x.startswith("G"), x)) # tri du cluster : génome "GUT" en premier   
            Gref  = next(iter(c_sorted))
            Dclust[Gref] = c_sorted
    return Dclust

def crea_files(type, Dclust, Dfna, DnewId, Dpro, df_filtred) :

    os.system(f'mkdir -p {OUTDIR}/{type}')

    pub     = ["CLO", "FER", "JAN", "JUN", "LIA", "NOR", "STO", "TIA", "ZUO", "ZUB"]
    Dfiles  = {
        "name_f1" : f'{OUTDIR}/{type}/{type}_genomes_micro.fna' ,
        "name_f2" : f'{OUTDIR}/{type}/{type}_genomes_ids.tsv',
        "name_f3" : f'{OUTDIR}/{type}/{type}_clusters.tsv',
        "name_f4" : f'{OUTDIR}/{type}/{type}_counts_gut25.tsv',
        "name_f5" : f'{OUTDIR}/{type}/{type}_VIRIDIC.tsv',
        "name_f6" : f'{OUTDIR}/{type}/{type}_prophages.tsv'
    }
    
    with open(Dfiles["name_f1"],'w') as f1 :
        for g in Dclust.keys() :
            print(f'>{g}\n{Dfna[g]}', file = f1)                             # création d'un nouveau fasta avec génomes dérépliqués

    with open(Dfiles["name_f2"],'w') as f2 :              # fichier liant ancien ID et nouvel ID
        for i in DnewId :
            print(i, DnewId[i], sep="\t", file = f2)
    
    with open(Dfiles["name_f3"],'w') as f3 :                   # fichier avec génome de référence -> liste de génomes identiques (dérépliqués)

        with open(Dfiles["name_f4"],'w') as f4 :

            print("GENOME_rep", "CLO", "FER", "JAN", "JUN", "LIA", "NOR", "STO", "TIA", "ZUO", "ZUB", sep = "\t", file = f4)
            for rep in Dclust :
                cluster = ";".join(Dclust[rep])
                print(rep, end = "", file = f4)
                for p in pub :
                    if f'G_{p}' in cluster :
                        print(f'\t1', end = "", file = f4)
                    else : 
                        print(f'\t0', end = "", file = f4)
                print(file = f4)
                print(cluster,sep="\t", file = f3)
    
    with open(Dfiles["name_f5"],'w') as f5 :
        print('gA', 'gB', 'lenA', 'lenB', 'sim%', 'dist%', 'nb_idAB', 'nb_idBA', 'nb_aliA', 'nb_aliB', 'fracAliA', 'fracAliB', 'ratioLen', 'ANI_AB', sep = "\t", file = f5)
        for row in df_filtred.itertuples(index=False):
            ga = DnewId[row[0]]
            gb = DnewId[row[1]]
            if ga in Dclust.keys() and gb in Dclust.keys() and ga != gb:
                filtered_row = row._asdict()  
                del filtered_row['gA']
                del filtered_row['gB']
                print(ga, gb, *filtered_row.values(), sep="\t", file=f5)

    with open(Dfiles['name_f6'],'w') as f6 :
        for prophage in Dpro.keys():
            if len(Dpro[prophage]) > 1 :
                print(prophage, ";".join(Dpro[prophage][1:]), sep="\t", file=f6)
            else :
                print(prophage, file = f6)
                
### EXECUTION DU SCRIPT ###

DFNA,DNEW = read_fasta(FNA)
            
DF_ORI = pd.read_csv(VIRIDIC_TSV, sep='\t')                                     # chargement du TSV en df pandas
DF_ORI['ANI_AB'] = (DF_ORI['nb_idAB'] * 100) / DF_ORI['nb_aliA']                        # ajout d'une colonne pour calculer l'ANI de A sur B (très proche de ANI de B sur A)

PROB,DPRO = prob_genomes(DF_ORI,DNEW)

DF1  = DF_ORI[
    ((~DF_ORI['gA'].isin(PROB))|(~DF_ORI['gB'].isin(PROB))) &                   # on garde seulement les génomes représentatifs                                                                  
    (DF_ORI["ANI_AB"] == 100) &                                                # couples ANI > ou = à 99.5                                            # fraction alignée AsurB >=0.99
    (DF_ORI["ratioLen"] >= 0.999) ]                

DF2  = DF_ORI[
    ((~DF_ORI['gA'].isin(PROB))|(~DF_ORI['gB'].isin(PROB))) &                   # on garde seulement les génomes représentatifs                                                                  
    (DF_ORI["ANI_AB"] >= 99.5) &                                                # couples ANI > ou = à 99.5                                            # fraction alignée AsurB >=0.99
    (DF_ORI["fracAliA"] >= 0.99) &
    (DF_ORI["fracAliB"] >= 0.99) ]                                              # ratio des longueurs égal

crea_files("IDENTICAL", dereplicate(DF1,DNEW), DFNA, DNEW, DPRO, DF1)
crea_files("GENOMOVAR", dereplicate(DF2,DNEW), DFNA ,DNEW, DPRO, DF2)
