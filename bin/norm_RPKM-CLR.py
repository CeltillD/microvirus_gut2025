#!/usr/bin/env python3

import sys
import math
from collections import defaultdict
from scipy.stats import gmean

# USAGE : ./norm_RPKM-CLR.py  FINAL_INFO_MAPPING.tsv  clusters_Rep60_Cap50.txt

def calc_rpkm(raw_count,lenMicro,libsize) :             # calcule la valeur du RPKM pour un comptage brut d'un micro dans un virome   calc_rpkm( 12,  5751,   2334115) 
    rpkm = (raw_count * 1e+09) / (lenMicro * libsize)
    return rpkm

def calc_clr(val_rpkm,moy_geo) :                        # calcule la valeur CLR à partir du RPKM et de la moyenne géométrique
    val_clr = round(math.log((val_rpkm / moy_geo)),3)
    return val_clr

INFO_MAPPING = sys.argv[1]                                          # fichier avec tailles de librairies, reads mappés par génomes...
INFO_CLUSTER = sys.argv[2]                                          # fichier indiquant quels génomes sont dans quels genres

DviromeToGenomovarRawCounts = defaultdict(lambda : defaultdict(int)) # comptages bruts par génome dans un virome
DlengthGenomovars = {}                                              # taille des génomes
DlibSize = {}                                                       # taille de la librairie par virome 
DrepGenreToClusters = {}                                            # attribue un génome à son genre (cluster)                    
DrpkmGenomovars = defaultdict(lambda : defaultdict(float))          # valeur de rpkm par génome dans un virome
DgmeanVirome    = {}                                                # lie un virome à sa moyenne géométrique
DViromeToGenusRawCounts = defaultdict(lambda : defaultdict(int))    # comptages bruts par genre dans un virome
DrpkmGenus = defaultdict(lambda : defaultdict(float))               # valeur de rpkm par genre dans un virome

# Exemple INFO_MAPPING
"""
VIROME		LIBSIZE	MICRO			LENGTH	COUNT
JUN_D1Apa_HC	2334115	A_GutCatV1_c1		5029	20
JUN_D1Apa_HC	2334115	A_GutCatV1_c1003	4686	180
"""

with open(INFO_MAPPING,'r') as f1 :
    premiere_ligne = f1.readline() 
    for l in f1 :
        l = l.rstrip().split('\t')
        virome  = l[0]                                      # nom du virome
        libsize = int(l[1])                                 # taille de la librairie (=nb de reads uniques mappés sur TOUS les virus)
        micro   = l[2]                                      # nom du génome
        length  = int(l[3])                                 # taille du génome
        count   = int(l[4])                                 # comptage brut des reads mappés contre ce génome dans ce virome
        DviromeToGenomovarRawCounts[virome][micro] = count  # stockage des infos
        DlibSize[virome] = libsize
        DlengthGenomovars[micro] = length


moy_lg_genome_per_genus = defaultdict(int)
with open(INFO_CLUSTER,'r') as f2 :
    for l in f2 :
        l = l.rstrip().split('\t')
        repGenre        = l[0]
        listeGenomovars = []
        
        for g in l :
            if g in DlengthGenomovars :						            
                moy_lg_genome_per_genus[ repGenre ] +=  DlengthGenomovars[g]  	# pour chaque genre, on fait la moyenne de la taille des génomes
                listeGenomovars.append(g)

        moy_lg_genome_per_genus[ repGenre ] = moy_lg_genome_per_genus[ repGenre ] / len(listeGenomovars)

        DrepGenreToClusters[ repGenre ]     = listeGenomovars                    


for genus in DrepGenreToClusters : 
    for virome in DviromeToGenomovarRawCounts :
        for genomovar in DrepGenreToClusters[genus]:
            if genomovar in DviromeToGenomovarRawCounts[virome] :
                # somme des comptages bruts de génomovars d'un meme genre
                DViromeToGenusRawCounts[virome][genus] += DviromeToGenomovarRawCounts[virome][genomovar] 
        if genus not in DViromeToGenusRawCounts[virome] :
            DViromeToGenusRawCounts[virome][genus] =1                           #         on met le pseudocount si pas ce genre dans ce virome


# calcul des RPKM pour tous les genres dans les viromes
for virome in DviromeToGenomovarRawCounts :                       			# pour chaque virome                 
    for genre in DrepGenreToClusters :							#     pour chaque genre
        lenMicro    = moy_lg_genome_per_genus[genre]					#         on récupère la taille moyenne des génomes de ce genre
        libsize     = DlibSize[virome]
        DrpkmGenus[virome][genre] = calc_rpkm(  DViromeToGenusRawCounts[virome][genre],  lenMicro,   libsize)     # calcul du rpkm
    DgmeanVirome[virome] = gmean(list(DrpkmGenus[virome].values()))            		# calcul du la moyenne geo du virome


list_Virome = sorted(DViromeToGenusRawCounts.keys())

with open("rpkm-clr_matrix_GENUS.tsv",'w') as o2 :

    print("GENUS",end="",file=o2)
    
    for v in list_Virome :
        print("\t",v,sep="",end="",file=o2)
    print(file=o2)
    
    for genuss in DrepGenreToClusters :
        print(genuss,end="",file=o2)
        for virome in list_Virome :
            #print(virome, genuss, DrpkmGenus[virome][genuss],  DgmeanVirome[virome])
            print("\t", calc_clr(DrpkmGenus[virome][genuss],  DgmeanVirome[virome])  ,  sep="",  end="", file=o2)
        print(file=o2)





























