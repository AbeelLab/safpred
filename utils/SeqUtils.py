#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import numpy as np

def RevComp(s):
    c = s.replace('A', 'X').replace('C','Y').replace('T', 'A').replace('G','C').replace('Y','G').replace('X','T')
    return c[::-1]

def RandSeq(l, A=('A','C','G','T'), p=(.25,.25,.25,.25)):
    """
    Generate a random sequence of length l from an alphabet A with prob distribution
    Input: sequence length, alphabet (default: (A,C,G,T)) and probability dist. (uniform)
    Returns: random sequence
    -- credit: stackoverflow answer 49620846
    """
    return ''.join(np.random.choice(A, p=p) for _ in range(l))

def HamDist(seq1, seq2, toRef=True):
    """
    Compute Hamming distance between 2 sequences of varying lengths, trim if larger
    Input: 2 sequences, flag to use the first sequence as reference
    Returns: Hamming distance
    """
    d = sum(s1 != s2 for (s1,s2) in zip(seq1,seq2))
    if toRef: d = d/len(seq1)
    return d

def GeneDistbak(pos1, pos2, contiglen=None):
    """
    Calculate intergenic distance between two genes on the same contig
    Input: positions of 2 genes (as intervals), circular flag (default: True)
    and contig length (necessary only for circular contigs)
    Returns the intergenic distance
    """
    if pos1[0] > pos2[0]:
        corpos1 = pos2
        corpos2 = pos1
    else:
        corpos1 = pos1
        corpos2 = pos2
    if contiglen:
        d = min(corpos2[0]-corpos1[1], contiglen-corpos2[1]+corpos1[0])
    else:
        d = corpos2[0]-corpos1[1]
    return d

def GeneDist(pos1, pos2, contiglen=None):
    """
    Calculate intergenic distance between two genes on the same contig
    Input: positions of 2 genes (as intervals), circular flag (default: True)
    and contig length (necessary only for circular contigs)
    Returns the intergenic distance
    ref: Mihelčić, M., Šmuc, T. & Supek, F. Patterns of diverse gene functions 
    in genomic neighborhoods predict gene function and phenotype. 
    Sci Rep 9, 19537 (2019). https://doi.org/10.1038/s41598-019-55984-0
    """
    if pos1[0] > pos2[1]:
        corpos1 = pos2
        corpos2 = pos1
    elif pos2[0] > pos1[1]:
        corpos1 = pos1
        corpos2 = pos2
    else:
        return 0
    if contiglen:
        d = min(abs(corpos2[0]-corpos1[1]), abs(corpos2[1]-contiglen-corpos1[0]))
    else:
        d = abs(corpos2[0]-corpos1[1])
    return d
    
def GeneDistStranded(pos1, pos2, contiglen=None):
    """
    Calculate intergenic distance between two genes on the same contig, also 
    takes strandedness into account
    Input: positions of 2 genes (as intervals), circular flag (default: True)
    and contig length (necessary only for circular contigs)
    Returns the intergenic distance
    ref: Mihelčić, M., Šmuc, T. & Supek, F. Patterns of diverse gene functions 
    in genomic neighborhoods predict gene function and phenotype. 
    Sci Rep 9, 19537 (2019). https://doi.org/10.1038/s41598-019-55984-0
    """
    if pos1[-1] != pos2[-1]:
        return 1e6
    if pos1[0] > pos2[1]:
        corpos1 = pos2
        corpos2 = pos1
    elif pos2[0] > pos1[1]:
        corpos1 = pos1
        corpos2 = pos2
    else:
        return 0
    if contiglen:
        d = min(abs(corpos2[0]-corpos1[1]), abs(corpos2[1]-contiglen-corpos1[0]))
    else:
        d = abs(corpos2[0]-corpos1[1])
    return d

def GetSeqID(fastaList, genomeList=None):
    """
    Parse fasta files to collect all sequence entries
    Input: fasta file(s) and a list of corresponding genomes (optional)
    Returns: a dict mapping genome IDs -> list of sequence entries
    """
    from Bio import SeqIO

    if isinstance(fastaList,str): fastaList = [fastaList]
    if not genomeList: genomeList = ['_'.join(line.split('/')[-1].split('_')[:2]) for line in fastaList]
    genome2Seq = {genome: [] for genome in genomeList}
    for (i,fasta) in enumerate(fastaList):
        print("Parsing genome # {}".format(i+1))
        genome2Seq[genomeList[i]] = [rec.name for rec in SeqIO.parse(fasta, 'fasta')]
    seq2Genome = {seqID: k for k,v in genome2Seq.items() for seqID in v}
    return genome2Seq, seq2Genome

def ExtractFeatures(fastaList, annotList, genomeList, extractFeat='cds', keyword=None, outFile=None):
    """
    Extract gene features, their DNA sequence and genome information -- gff only
    Input: list of relative path to fasta files, annotation (gff) files and another
    list of genome IDs, gene keyword, output file to write save the results
    Returns: a dict mapping genome IDs -> genes, dataframe for gene recs
    """
    import pandas as pd
    from Bio import SeqIO

    if keyword: keyword = keyword.lower()
    genome2Gene = {genome: {} for genome in genomeList}

    for (i,gff) in enumerate(annotList):
        genome = genomeList[i]
        print("Genome # {} -- {}".format(i+1,genome))
        with open(gff, 'r') as f:
            seqList = list(SeqIO.parse(fastaList[i], 'fasta'))
            for line in f:
                if line.startswith('#'): continue
                if (keyword) and (keyword not in line.lower()): continue
                recFlds = line.strip().split('\t')
                if recFlds[2].lower() == extractFeat:
                    protID,gene,prod = [],[],[]
                    pos = [*map(int,recFlds[3:5])]
                    attrFlds = recFlds[8].split(';')
                    protID,gene,prod,locustag = '','','',''
                    for attr in attrFlds:
                        x = attr.split('=')
                        if (x[0].lower()=='name') or (x[0].lower()=='genbank'):
                            protID = x[1]
                        elif x[0].lower()=='gene':
                            gene = x[1]
                        elif x[0].lower()=='product':
                            prod = x[1]
                        elif x[0].lower()=='locus_tag':
                            locustag = x[1]
                    if keyword:
                        if bool(not gene and prod):
                            gene = [x for x in prod.lower().split(' ') if 'stx' in x]
                            if gene: gene = gene[0]+'_subunit'
                    if recFlds[0] in genome2Gene[genome].keys():
                        genome2Gene[genome][recFlds[0]][gene] = {'id':protID, 'product':prod, \
                                                                 'locustag':locustag, 'pos':pos}
                    else:
                        genome2Gene[genome][recFlds[0]] = {gene: {'id':protID, 'product':prod, \
                                                                  'locustag':locustag, 'pos':pos}}

            for (k,v) in genome2Gene[genome].items():
                seqRec = [rec for rec in seqList if rec.id == k][0]
                for (gene,attr) in v.items():
                    genome2Gene[genome][k][gene]['seq'] = seqRec.seq[attr['pos'][0]:attr['pos'][1]]
    gene2Genome = {}
    for (k,v) in genome2Gene.items():
        for (acc,genes) in v.items():
            for (gene,attr) in genes.items():
                if attr['id'] not in gene2Genome.keys():
                    gene2Genome[attr['id']] = {'gene': gene, 'product': attr['product'], \
                                               'genomes': [acc], 'seq': str(attr['seq'])}
                else:
                    gene2Genome[attr['id']]['genomes'].append(acc)
    genedf = pd.DataFrame(gene2Genome.values(), index=gene2Genome.keys())
    genedf.loc[:,'gene'] = genedf.apply(lambda x: '_'.join([x.gene.split('_')[0], x['product'].split(' ')[-1]]) \
                                        if 'subunit' in x.gene else x.gene, axis=1)
    genedf.loc[:,'genomes'] = genedf.genomes.apply(lambda x: ';'.join(x))
    genedf.sort_values(by='gene', inplace=True)
    if outFile: genedf.to_csv(outFile, index=True, header=True, sep='\t')

    return genome2Gene, genedf

def ExtractSeqFromFasta(infile, startpos, endpos, strand=None, contig=None):
    """
    Extract a region from fasta file, takes the strand into account as well
    Input: input fasta file path, start and end position (integer), 
    strand ('+' or '-' or None if unknown) and the contig name (if there's more than one)
    Returns: a SeqRecord object with the extracted region
    """
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord
    
    if not os.path.isfile(infile): 
        print("Ooops! File not found")
        return
    fullseq = []
    if contig:
        for rec in SeqIO.parse(infile,'fasta'):
            if rec.id==contig:
                fullseq = rec
        if not fullseq:
            print("This contig does not exist in this fasta file!!")
            return
    else:
        fullseq = SeqIO.read(infile,'fasta')
    
    seq = fullseq.seq[startpos-1:endpos]
    if strand == '-': # Take the reverse complement
        seq = seq.reverse_complement()
    recid = '{} | {}-{} | {}'.format(contig,startpos,endpos,strand)
    extractedseq = SeqRecord(seq, id=recid, name='', description='')
    
    return extractedseq
