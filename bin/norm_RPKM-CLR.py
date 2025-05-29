#!/usr/bin/env python3

import sys
import math
from collections import defaultdict
from scipy.stats import gmean

def calc_rpkm(raw_count,lenMicro,libsize) :             # calcule la valeur du RPKM pour un comptage brut d'un micro dans un virome
    rpkm = (raw_count * 1e+09) / (lenMicro * libsize)
    return rpkm

def calc_clr(val_rpkm,moy_geo) :                        # calcule la valeur CLR à partir du RPKM et de la moyenne géométrique
    val_clr = round(math.log((val_rpkm / moy_geo)),3)
    return val_clr

INFO_MAPPING = sys.argv[1]                                          # fichier avec tailles de librairies, reads mappés par génomes...
INFO_CLUSTER = sys.argv[2]                                          # fichier indiquant quels génomes sont dans quels genres

DviomeToGenomovarRawCounts = defaultdict(lambda : defaultdict(int)) # comptages bruts par génome dans un virome
DlengthGenomovars = {}                                              # taille des génomes
DlibSize = {}                                                       # taille de la librairie par virome 
DGenomovarsToClusters = {}                                          # attribue un génome à son genre (cluster)                    
DrpkmGenomovars = defaultdict(lambda : defaultdict(float))          # valeur de rpkm par génome dans un virome
DgmeanVirome    = {}                                                # lie un virome à sa moyenne géométrique
DViromeToGenusRawCounts = defaultdict(lambda : defaultdict(int))    # comptages bruts par genre dans un virome
DrpkmGenus = defaultdict(lambda : defaultdict(float))               # valeur de rpkm par genre dans un virome


with open(INFO_MAPPING,'r') as f1 :
    for l in f1 :
        l = l.rstrip().split('\t')
        virome  = l[0]                                      # nom du virome
        libsize = int(l[1])                                 # taille de la librairie (=nb de reads uniques mappés sur TOUS les virus)
        micro   = l[2]                                      # nom du génome
        length  = int(l[3])                                 # taille du génome
        count   = int(l[4])                                 # comptage brut des reads mappés contre ce génome dans ce virome
        DviomeToGenomovarRawCounts[virome][micro] = count   # stockage des infos
        DlibSize[virome] = libsize
        DlengthGenomovars[micro] = length

with open(INFO_CLUSTER,'r') as f2 :
    for l in f2 :
        l = l.rstrip().split('\t')
        DGenomovarsToClusters[l[0]] = l[1].split(';')       # le génome en l[0] a été désigné en amont comme le représentant du cluster (genre)

for virome in DviomeToGenomovarRawCounts :                                          # calcul des RPKM pour tous les génomes dans les viromes
    for genomovar,raw_count in DviomeToGenomovarRawCounts[virome].items():
        lenMicro    = DlengthGenomovars[genomovar]
        libsize     = DlibSize[virome]
        DrpkmGenomovars[virome][genomovar] = calc_rpkm(raw_count,lenMicro,libsize)  # calcul du rpkm
    DgmeanVirome[virome] = gmean(list(DrpkmGenomovars[virome].values()))            # calcul du la moyenne geo du virome

with open("rpkm-clr_matrix_GENOMOVARS.tsv",'w') as o1 :                             # dataframe ts des valeurs normalisées
    print("GENOMOVAR",end="",file=o1)
    for v in DviomeToGenomovarRawCounts :
        print("\t",v,sep="",end="",file=o1)
    print(file=o1)
    for genomovar in DlengthGenomovars :
        print(genomovar,end="",file=o1)
        for virome in DviomeToGenomovarRawCounts :
            if genomovar in DviomeToGenomovarRawCounts[virome] :
                # je calcule directement la valeur normalisée CLR ici à partir du RPKM et GM du virome
                print("\t",calc_clr(DrpkmGenomovars[virome][genomovar],DgmeanVirome[virome]),sep="",end="",file=o1)
            else : 
                # si comptage nul, on met une valeur faible, non-nulle pour execcuter le CLR
                print("\t",calc_clr(1e-03,DgmeanVirome[virome]),sep="",end="",file=o1)
        print(file=o1)

### Meme prcédé ici pour les genres à la différence qu'une somme des comptages brutes par genre est réalisée en amont

for genus in DGenomovarsToClusters : 
    for genomovar in DGenomovarsToClusters[genus]:
        for virome in DviomeToGenomovarRawCounts :
            if genomovar in DviomeToGenomovarRawCounts[virome] :
                # somme des comptages bruts de génomovars d'un meme genre
                DViromeToGenusRawCounts[virome][genus] += DviomeToGenomovarRawCounts[virome][genomovar] 

### ---

for virome in DviomeToGenomovarRawCounts :
    for genomovar,raw_count in DViromeToGenusRawCounts[virome].items():
        lenMicro    = DlengthGenomovars[genomovar]
        libsize     = DlibSize[virome]
        DrpkmGenus[virome][genomovar] = calc_rpkm(raw_count,lenMicro,libsize)
    DgmeanVirome[virome] = gmean(list(DrpkmGenus[virome].values()))

with open("rpkm-clr_matrix_GENUS.tsv",'w') as o2 :
    print("GENUS",end="",file=o2)
    for v in DViromeToGenusRawCounts :
        print("\t",v,sep="",end="",file=o2)
    print(file=o2)
    for genomovar in DGenomovarsToClusters :
        print(genomovar,end="",file=o2)
        for virome in DViromeToGenusRawCounts :
            if genomovar in DViromeToGenusRawCounts[virome] :
                print("\t",calc_clr(DrpkmGenus[virome][genomovar],DgmeanVirome[virome]),sep="",end="",file=o2)
            else : 
                print("\t",calc_clr(1e-03,DgmeanVirome[virome]),sep="",end="",file=o2)
        print(file=o2)
