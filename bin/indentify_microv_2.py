#!/usr/bin/env python3

import sys
import os
import subprocess
import pandas as pd
from collections import defaultdict
from Bio import pairwise2

ContigFile,ViromeName=sys.argv[2],sys.argv[3]  # fichier fasta et nom du virome
WORKDIR = sys.argv[1]
os.system(f'mkdir -p {WORKDIR}')


Dfna = {}                       # Dico reliant identifiant de contig et sa séquence
Dnew = {}
Dmic = {}
Capsides = set()
DMicroComplets = {}
match_score = 1                 # Score pour une correspondance
mismatch_score = -2             # Pénalité pour une non-correspondance
gap_open_penalty = -5           # Pénalité pour l'ouverture d'un gap
gap_extend_penalty = -2         # Pénalité pour l'extension d'un gap
cvp1 = set()                    # Contigs possédant la capside (VP1)
cvp4 = set()                    # Contigs possédant la protéine de réplication (VP4)
ntot = 0                        # Nombre de contigs en entrée
nlen = 0                        # nb de contigs de taille 2 à 10 kb 
ncir = 0                        # nb de contigs circulaire + 2 à 10kb
ncvp = 0                        # nb de contigs avec VP1 et VP4 + circulaire + 2 à 10 kb
ncv1 = 0                        # nb de contigs avec uniquement VP1 + circulaire + 2 à 10kb
ncv4 = 0                        # nb de contigs avec uniquement VP4 + circulaire + 2 à 10kb
Micro = {}

def ReverseComplement(seq):
	""" Objectif : 
	Reverse complémenter des séquences ADN
	"""
	complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
	Revseq=""
	for nuc in seq:
		Revseq += complement_dict[nuc]
	NewSeq = Revseq[::-1]
	return NewSeq

# Enregistrement des séquences correspondant au critère de taille (2-10kb)
with open (ContigFile) as cf :
        id = None
        seq = []
        seq_complete = ""
        for ligne in cf:
            ligne = ligne.strip()
            if ligne.startswith(">"):                       # >k127_0 flag=0 multi=9.0 len=225 
                l = ligne.split()
                if id is not None:  
                    seq_complete = "".join(seq)
                    if 2000 <= len(seq_complete) <= 10000 : # critère de taille
                        Dfna[id] = seq_complete
                id = ViromeName + "_" + str(l[0].replace(">","").split("_")[1]) + "_" + str(len(seq_complete))
                seq = []  
            else:
                seq.append(ligne)  
        if id is not None:                                  # pour la dernière séquence du fasta
            seq_complete = "".join(seq)
            if 2000 <= len(seq_complete) <= 10000:
                Dfna[id] = seq_complete

# Test de circularité par alignement global de 2 fragments d'un meme contig
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
                ncir += 1 
                decircular=Seq[:Index]                                      # séquence décircularisée
                Dnew[Contig] = decircular
                with open(WORKDIR+"tmp_"+ViromeName+".fna",'a') as stdin :
                    print(">"+Contig,decircular,sep="\n", file = stdin)

#os.system('mmseqs easy-search '+WORKDIR+'tmp_'+ViromeName+'.fna /home/cedumont/Micro/prots/13390_ultimate_VP1s.faa '+WORKDIR+'tmp_VP1_'+ViromeName+'.tsv '+WORKDIR+'tmp_vp1'+ViromeName+' ')
#os.system('mmseqs easy-search '+WORKDIR+'tmp_'+ViromeName+'.fna /home/cedumont/Micro/prots/12933_ultimate_VP4s.faa '+WORKDIR+'tmp_VP4_'+ViromeName+'.tsv '+WORKDIR+'tmp_vp4'+ViromeName+' ')

# Si le fasta des séquences circulaires est présent
if os.path.isfile(WORKDIR+'tmp_'+ViromeName+'.fna') :

    # On Blast les séquences circulaires obtenues contre les fasta des protéines VP1 et VP4
    subprocess.run([
        'mmseqs', 'easy-search',
        f'{WORKDIR}tmp_{ViromeName}.fna',
        '/home/cedumont/Micro/prots/13390_ultimate_VP1s.faa',
        f'{WORKDIR}tmp_VP1_{ViromeName}.tsv',
        f'{WORKDIR}tmp_vp1{ViromeName}'
    ], check=True)

    subprocess.run([
        'mmseqs', 'easy-search',
        f'{WORKDIR}tmp_{ViromeName}.fna',
        '/home/cedumont/Micro/prots/12933_ultimate_VP4s.faa',
        f'{WORKDIR}tmp_VP4_{ViromeName}.tsv',
        f'{WORKDIR}tmp_vp4{ViromeName}'
    ], check=True)

# Si MMSEQS donne un résultat 
if os.path.isfile(WORKDIR+'tmp_VP1_'+ViromeName+'.tsv') :

    # Test si le contig contient la capside (VP1)
    with open (WORKDIR+'tmp_VP1_'+ViromeName+'.tsv') as vp1 :
        for l in vp1 :
            l=l.rstrip().split("\t")
            if int(l[11]) >= 50 :
                cvp1.add(str(l[0]))

if os.path.isfile(WORKDIR+'tmp_VP4_'+ViromeName+'.tsv') :
    # Test si le contig contient la protéine de réplication (VP4)
    with open (WORKDIR+'tmp_VP4_'+ViromeName+'.tsv') as vp4 :
        for l in vp4 :
            l=l.rstrip().split("\t")
            if int(l[11]) >= 50 :
                cvp4.add(str(l[0]))

# Identifiants de contigs contenant les 2 gènes marqueurs
ContigsComplets = cvp4.intersection(cvp1)
ncvp = len(ContigsComplets)

# Composer un nouveau dico avec les potentiels micro
Dmic = {key: value for key, value in Dnew.items() if key in ContigsComplets}

with open(WORKDIR+"tmp_mic_"+ViromeName+".fna",'w') as stdin :
    for k,s in Dmic.items() :
        print(">"+k,s+str(s[0:2000]),sep="\n", file = stdin)

# prodigal -i allcontigs.fna -a allproteins.faa > /dev/null

if os.path.isfile(WORKDIR+'tmp_mic_'+ViromeName+'.fna') :
    subprocess.run([
            'prodigal', '-i',
            f'{WORKDIR}tmp_mic_{ViromeName}.fna',
            '-a',
            f'{WORKDIR}tmp_genes_{ViromeName}.faa',
            '-d',
            f'{WORKDIR}tmp_genes_{ViromeName}.fna',
            '-p' ,'meta'
        ], check=True)

# mmseqs easy-search allproteins.faa {ficherCapside} tmp_proteins_VS_Capsids.tsv tmp --format-output \"query,target,evalue,raw,pident\"

if os.path.isfile(WORKDIR+'tmp_genes_'+ViromeName+'.faa') :
    subprocess.run([
            'mmseqs', 'easy-search',
            f'{WORKDIR}tmp_genes_{ViromeName}.faa',
            '/home/cedumont/Micro/prots/13390_ultimate_VP1s.faa',
            f'{WORKDIR}tmp_proteins_VS_Capsids_{ViromeName}.tsv',
            f'{WORKDIR}tmp_cap_{ViromeName}',
            '--format-output',
            'query,target,evalue,raw,pident'
        ], check=True)

if os.path.isfile(f'{WORKDIR}tmp_proteins_VS_Capsids_{ViromeName}.tsv') :
    Dscore = defaultdict(int)
    DCapside = {}
    DBestCap = {}
    with open(f'{WORKDIR}tmp_proteins_VS_Capsids_{ViromeName}.tsv','r') as cap :
        for l in cap :
            l=l.rstrip().split("\t")
            if int(l[3]) > 50 :
                Contig = "_".join(l[0].split("_")[:-1])
                IdProt = l[0]
                Dscore[IdProt] = int(l[3])
                if int(l[3]) > Dscore[Contig] :
                    Dscore[Contig] = int(l[3])
                    DBestCap[Contig] = IdProt
                    DCapside[IdProt] = Contig
                Capsides.add(IdProt)

if os.path.isfile(f'{WORKDIR}tmp_genes_{ViromeName}.fna') :
    DSeqCap = {}
    Dsens = {}
    id = 0
    with open(f'{WORKDIR}tmp_genes_{ViromeName}.fna','r') as genes :
        for l in genes :
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

with open(f'{WORKDIR}complete_{ViromeName}.fna','w') as fna_f :
    for cap,info in DSeqCap.items() :
        if cap in set(DBestCap.values()) :
            Contig = DCapside[cap]
            SeqBase = Dmic[Contig]
            Ext = str(SeqBase[0:2000])
            SeqLong = SeqBase + Ext
            if info[3] != "1" :
                FwdSeqLong = ReverseComplement(SeqLong)
                FwdExt = ReverseComplement(Ext)
            else :
                FwdSeqLong = SeqLong
                FwdExt = Ext
            IndCap = FwdSeqLong.find(info[0])
            FinalContig = str(FwdSeqLong[IndCap:] + FwdSeqLong[:IndCap]).replace(FwdExt,"",1)
            GCprc = (FinalContig.count("C") + FinalContig.count("G")) / len(Seq)
            GC = str(round(GCprc,3)).replace("0.","")
            LenCap = str(info[2] - info[1])
            IdContig = "_".join(Contig.split("_")[:-2])
            NumContig = str(Contig.split("_")[-2])
            LenContig = str(len(FinalContig))
            Chevron = f'>{IdContig}_c{NumContig}_VP1-4_Cir{LenContig}_GC{GC}_Cap{LenCap}'
            print(Chevron,FinalContig,sep="\n",file=fna_f)

# Nettoyage des éléments temporaires
os.system(f'rm -rf '+WORKDIR+'tmp_*'+ViromeName+'*')
