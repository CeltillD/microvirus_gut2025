#!/usr/bin/python3

import sys
import numpy as np

# Commande bash réalisée en amont du script pour un blastn 50k x 50k
"""
nohup blastn 
	-query				Microv50k.fna 
	-db 				dbmicro50k 
	-out 				m50k_vs_m50k.txt 
	-outfmt 			6 
	-evalue 			1 
	-max_target_seqs 	50000 
	-num_threads 		46 
	-word_size 			7 		
	-reward 			2 
	-penalty 			-3 
	-gapopen 			5 
	-gapextend 			2
	&
"""
ficFna	   = sys.argv[1]													# ex : micro.fasta
fileBLAST  = ficFna[0:-6]+'_BLASTn.tab'						# ex : micro_BLASTn.tab
fileSortie = ficFna[0:-6]+'_Viridic.tsv'					# ex : micro_Viridic.tsv

dseq = {}																          # dico de toutes les séquences fasta
didentAB = {} 																    # dico qui associe chaque couple à son identité
couplec = "debut" 															  # ititialisation du couple courant

with open(ficFna) as f1 : 												# on enregistre les séquences
	nomc =""
	for lig in f1 : 
		lig=lig.rstrip()
		if lig.startswith(">") : 
			nomc = lig[1:]
			dseq[nomc] = ''
		else :
			dseq[nomc] += lig

# création de l'entete du fichier VIRIDIC
fe1 = open(fileSortie, "w")
print( 'gA', 'gB', 'lenA', 'lenB', 'sim%', 'dist%', 'nb_idAB', 'nb_idBA', 'nb_aliA', 'nb_aliB', 'fracAliA', 'fracAliB', 'ratioLen', sep = "\t", file=fe1)
#xxx_c8	Gap1x	4229	5973	54   	45.1	2799    	2793	   3405	      3397	    0.80	    0.56		0.70      

with open(fileBLAST) as f1 : 
	for lig in f1 : 
		lig = lig.rstrip() 
		
		li = lig.split()  
		newCouple = li[0]+" "+li[1]

		if couplec == "debut" : 											                                    # initialisation
			A = np.zeros(len(dseq[li[0]])) 
			couplec = newCouple

		if newCouple != couplec : 				         					                              # si le couple est nouveau en traite couplec
			na      = couplec.split()[0] 									
			nb      = couplec.split()[1]
			idAB    = np.sum(A) 											                                      # addition des %id de chaque position pour avoir le nb de bases identiques  
			fracali = np.count_nonzero(A) 									                                # on récupère le nombre de bases alignées 
			lA      = len(dseq[na])
			lB      = len(dseq[nb])
			sim2    = ((idAB)*100)/(lA) 		                                                # on calcule la similarité en % et on enregistre seulement si > 25
			
			if sim2 > 25 :													                                        # recup le l'id nuc de A sur B et la fraction alignée
				didentAB[couplec] = [idAB, fracali] 						
				
				coupleinv = couplec.split()[1]+" "+couplec.split()[0]
				
				if coupleinv in didentAB : 	
					sim      = ((didentAB[couplec][0]+didentAB[coupleinv][0])*100)/(lA+lB) 		  # on calcule la similarité en %
					dist     = 100 - sim                                                       	# on calcule la distance en %	
					FracAliA = didentAB[couplec][1]   / lA 										                  # on calcule la fraction alignée de A
					FracAliB = didentAB[coupleinv][1] / lB 										                  # on calcule la fraction alignée de B	
					Ratio    = min(lA, lB) / max(lA, lB) 										                    # on calcule le ratio des longueurs de A et B
				
					print( na, nb, lA, lB, sim, dist, didentAB[couplec][0], didentAB[coupleinv][0], didentAB[couplec][1], didentAB[coupleinv][1], FracAliA, FracAliB, Ratio, sep = "\t", file=fe1)
					if na==nb :
						del didentAB[couplec]
					else :
						del didentAB[couplec]
						del didentAB[coupleinv]

			couplec = newCouple
			na      = couplec.split()[0]		
			A       = np.zeros(len(dseq[na])) 									                           # reinit
		 
		if newCouple == couplec : 									                                     # tant que l'on est dans le meme couple courant
			for i in range(min(int(li[6]), int(li[7])+1),max(int(li[6]), int(li[7])+1)) : 
				if A[i-1] == 0 : 											                                       # si la position n'a pas été lue
					A[i-1] = float(li[2])/100

# pour la dernière ligne :

na      = couplec.split()[0] 									
nb      = couplec.split()[1]
idAB    = np.sum(A) 											                                          # addition des %id de chaque position pour avoir le nb de bases identiques  
fracali = np.count_nonzero(A) 									                                    # on récupère le nombre de bases alignées 
lA      = len(dseq[na])
lB      = len(dseq[nb])
sim2    = ((idAB)*100)/(lA) 		                                                    # on calcule la similarité en % et on enregistre seulement si > 25
			
if sim2 > 25 :
	didentAB[couplec] = [idAB, fracali] 						
	coupleinv = couplec.split()[1]+" "+couplec.split()[0]
				
	if coupleinv in didentAB : 	
		sim      = ((didentAB[couplec][0]+didentAB[coupleinv][0])*100)/(lA+lB) 		     # on calcule la similarité en %
		dist     = 100 - sim                                                       	   # on calcule la distance en %	
		FracAliA = didentAB[couplec][1]   / lA 										                     # on calcule la fraction alignée de A
		FracAliB = didentAB[coupleinv][1] / lB 										                     # on calcule la fraction alignée de B	
		Ratio    = min(lA, lB) / max(lA, lB) 										                       # on calcule le ratio des longueurs de A et B
		print( na, nb, lA, lB, sim, dist, didentAB[couplec][0], didentAB[coupleinv][0], didentAB[couplec][1], didentAB[coupleinv][1], FracAliA, FracAliB, Ratio, sep = "\t", file=fe1)				
