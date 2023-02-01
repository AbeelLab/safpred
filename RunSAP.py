#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os,sys,re
import numpy as np
import pandas as pd
import pickle
from copy import deepcopy
from tqdm import tqdm
from Bio import SeqIO

from utils import SeqUtils, SAPUtils

tqdm.pandas()


embmodel = 'esm'

with open('data/goa/go_classes.pkl', 'rb') as f:
    goclasses = pickle.load(f)
with open('data/goa/go_ics.pkl', 'rb') as f:
    goics = pickle.load(f) 
tsvfile = 'data/uniprot/uniprot2go_exp_set.tsv'
truelabels = SAPUtils.ParseTSV(tsvfile)

# Load cluster annotations

print("Loading cluster metadata")

cdf = pd.read_pickle('gtdbrep_allgenes_aa_clustered_addmembers_min10_sorted.pkl.gz')
cdf.update(cdf.members.apply(lambda x: x.split(';')))
cdf.update(cdf.members_contig.apply(lambda x: x.split(';')))
cdf.update(cdf.members_genbank.apply(lambda x: x.split(';')))

print("Loading cluster annotations")

print("Loading embeddings")

if embmodel=='t5':
    with open('data/embeddings/t5_sprot.pkl', 'rb') as f:
        sprotembeddings = pickle.load(f)
elif embmodel=='esm':
    with open('data/embeddings/esm1b_sprot.pkl', 'rb') as f:
        sprotembeddingse = pickle.load(f)
idmap = {k.split('|')[1]: k for k in sprotembeddingst5.keys()}

# Model parameters
goweight = 'distance'
usemaxoperon = True
normsimilarity = True
keepsingletons = True
thset = 'all'
nnset = 'all'
percentile = [99.999]
embmodel = 'esm1b'




# Choose an experiment
bacteria = 'ecoli'
dataset = 'full'
tvalues = np.linspace(0.0, 1.0, num=50)

ogopfile = 'data/operondb/operons.pkl.gz'

with open('data/testprot_{}.txt'.format(bacteria),'r') as f:
    testproteins = [line.strip() for line in f]
with open('datasets/trainprot_{}.txt'.format(bacteria),'r') as f:
    trainproteins = [line.strip() for line in f]

# testembeddings = {idmap[prot]: sprotembeddings[idmap[prot]] for prot in testproteins}
# trainembeddings = {idmap[prot]: sprotembeddings[idmap[prot]] for prot in trainproteins}
testfile = 'data/embeddings/esm1b_{}_test.pkl'.format(bacteria)
trainfile = 'data/embeddings/esm1b_{}_train.pkl'.format(bacteria)
testembeddings, trainembeddings = EmbUtils.GetEmbeddings(testfile, trainfile)

if dataset == 'full':
    keepclusters = cdf.index
else:
    with open('data/trainclusters_{}_{}.txt'.format(bacteria,dataset), 'r') as f:
        keepclusters = [int(line.strip()) for line in f]

print("Edit operon database to remove clusters -- if necessary")

opfile = 'data/operondb/operons_filterintdist_{}{}.pkl.gz'.format(bacteria,dataset)

print("Assigning operons to query points")

nndict, newopdf = SAPUtils.AssignOperons(opfile, embfile, cdf, keepclusters, \
                  testembeddings, normsimilarity=True, keepsingles=True, thset='all', nnset='all')   

print("Predict GO terms")
operonpredictions = SAPUtils.PredictFromAverageOperons(cluster2go, newopdf, nndict)
normoperonpredictions = PredUtils.NormalizePredictions(deepcopy(operonpredictions), goclasses)

