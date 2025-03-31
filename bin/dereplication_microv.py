#!/usr/bin/env python3

import os
import sys
import pandas as pd

### Initialisation des Variables ###

OUTDIR,VIRIDIC_TSV,FNA = sys.argv[1],sys.argv[2],sys.argv[3]  

COUPLES  = []                       # liste de couples de génomes identiques (cov>99 & ANI>99)
CLUSTERS = []                       # liste des clusters de génomes identiques (//)

Dfna  = {}                          # dico qui enregistre les séquences du fasta
Docc  = {}                          # dico qui compte le nombre d'occurences d'un génome potentiellement "problématique"
Dprob = {}                          # dico qui lie un couple de génome "problématique"
Didg  = {}                          # dico qui lie l'ancien ID du génome à son nouveau ID
Drep  = {}                          # dico qui lie un génome représentatif avec ses partenaires dans le cluster après filtration

probGenomes = set()                 # set des génomes considérés comme "problématiques" 
reprGenomes = set()                 # set des génomes représentatifs pour chaque cluster

ngut = 0                            # compteur des génomes du gut 2025

os.system(f'mkdir -p {OUTDIR}/INFO')

### Définition des Fonctions ###

def single_linkage(clusters):                                               # Clustering Single Linkage                                           
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
        if len(new_clusters) == len(clusters):                              # taille finale des clusters
            return new_clusters
        clusters = new_clusters

### Lecture du fichiers des génomes de Microvirus ###

with open(FNA) as f1 :                                                      # enregistrement des 49k séquences fasta dans un dico "Dfna"
	g = ""
	for l in f1 : 
		l=l.rstrip()
		if l.startswith(">") : 
			g = l[1:]
			Dfna[g] = ''
		else :
			Dfna[g] += l

### Calcul de l'Average Nucleotide Identity (ANI) des couples de génomes ###

df = pd.read_csv(VIRIDIC_TSV, sep='\t')                                     # chargement du TSV en df pandas
df['ANI_AB'] = (df['nb_idAB'] * 100) / df['nb_aliA']                        # ajout d'une colonne pour calculer l'ANI de A sur B (très proche de ANI de B sur A)

### Gestion des Génomes "problématiques" ###

df_prob = df[
     (df["ANI_AB"] > 80) & (                                                # Si ANI entre A -> B élevé :
     ((df["fracAliA"] > 0.8) & (df["fracAliB"] < 0.6)) |                    # un couple de génome est problématique lorsqu'il y a une différence significative des fractions alignées
     ((df["fracAliA"] < 0.6) & (df["fracAliB"] > 0.8)))
]

probPotentiels = pd.concat([df_prob["gA"], df_prob["gB"]]).unique()         # liste des génomes potentiellement problématiques

for pp in probPotentiels :                                                  # comptage de leur nombre d'occurences problématiques
	Docc[pp] = df_prob['gA'].str.count(pp).sum() + df_prob['gB'].str.count(pp).sum()

Dprob = df_prob.set_index('gA')['gB'].to_dict()

for gA, gB in Dprob.items():                                                # pour chaque couple de génomes problématiques, je m'affranchis de celui qui a le plus d'occurences
	if Docc[gA] > Docc[gB] :
		probGenomes.add(gA)
	elif  Docc[gA] < Docc[gB] :
		probGenomes.add(gB)
	else :                                                                    # si meme nombre d'occurences
		probGenomes.add(gA)
		probGenomes.add(gB)

### Gestion de la Déréplication des Génomes par Clustering Single Linkage ###

df_filtered  = df[
    ((~df['gA'].isin(probGenomes))|(~df['gB'].isin(probGenomes))) &         # on ignore les lignes des génomes problématiques                                                                  
    (df["ANI_AB"] == 100) &                                                 # couples ANI == 100%                                              # fraction alignée AsurB >=0.99
    (df["ratioLen"] >= 0.999)]                                              # ratio des longueurs égal

genomes = pd.concat([df_filtered["gA"], df_filtered["gB"]]).unique()        # noms des génomes valides

for _, row in df_filtered.iterrows():                                       # chaque couple étant identique
    coupleAB = [ str(row["gA"]), str(row["gB"]) ]
    COUPLES.append(coupleAB)

CLUSTERS = single_linkage(COUPLES)                                          # application du single linkage sur ces couples

for cluster in CLUSTERS :                                                   # séléction des génomes représentatifs
        clust = sorted(cluster, key=lambda x: (x.startswith("GUT"), x)) # tri du cluster : génome "GUT" en premier   
        Gref  = next(iter(clust))
        reprGenomes.add(Gref)                                               # on prend le premier génome
        Drep[Gref] = clust            

### Création d'un fichier fasta nettoyé ###

with open(f'{OUTDIR}/genomes_micro.fna','w') as o :
    for g in reprGenomes:

        if g.startswith("GUT") :                                            # priorité aux génomes du gut (2025)
                ngut  += 1                                                  # GUT25_ZUB_82326_HC_c158713_VP1-4_Cir6327_GC366_Cap1884
                newId  = "G_"+"_".join(g.split("_")[1:4])+"_c"+str(ngut)    # new ID : G_ZUB_82326_HC_c1                                                                  
                Didg[g] = newId 

        if g.startswith("AVrC"):                                            # idem pour les génomes de AVrC 
            newId  = "A_"+"_".join(g.split("_")[1:4]).replace("_XX","")     # AVrC_GutCatV1_XX_c3442_VP1-4_Cir5528_GC418_Cap2046
            Didg[g] = newId                                                 # new ID : A_GutCatV1_c3442

        if g.startswith("REF"):                                             # REF13k_EAF_DB_IMGmg_3300010885-c3
            newId  = g.replace("REF13k","R")                                # new ID : R_EAF_DB_IMGmg_3300010885-c3
            Didg[g] = newId

        print(f'>{newId}\n{Dfna[g]}', file = o)                             # création d'un nouveau fasta avec génomes dérépliqués

### Enregistrement des infos supplémentaires ###

with open(f'{OUTDIR}/INFO/genomes_ids.tsv','w') as ids :                    # fichier liant ancien ID et nouvel ID
    for i in Didg :
        print(i, Didg[i], sep="\t", file = ids)

with open(f'{OUTDIR}/INFO/rep_clusters.tsv','w') as r :                     # fichier avec génome de référence -> liste de génomes identiques (dérépliqués)
    for rep in Drep :
        print(rep, ";".join(Drep[rep]), sep="\t", file = r)

with open(f'{OUTDIR}/INFO/derep_VIRIDIC.tsv','w') as v :
    for _, row in df.iterrows():
        if row["gA"] in reprGenomes and row["gB"] in reprGenomes :
            filtered_row = row.drop(['gA', 'gB'])
            print(Didg[row["gA"]],Didg[row["gB"]],*filtered_row, sep="\t", file = v)

#############################################################################
