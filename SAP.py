#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from numpy.linalg import norm
from copy import deepcopy
import pandas as pd
import pickle
def ParseOperonFile(opfile, embfile, keepclusters, cdf=None, keepsingles=False, avgoperon=False):
    '''
    Parse the operon file, noticed that I need this more often and I thought I would >.>
    Input: operon file path, list of clusters allowed, parameter to keep singleton operons 
    in the database (default: False)
    Returns the operon database
    '''
    filtopdf = pd.read_pickle(opfile)
    nrfiltopdf = filtopdf.copy()
    nrfiltopdf.update(nrfiltopdf.operon.apply(lambda x: tuple(x)))
    nrfiltopdf.drop_duplicates('operon', inplace=True)
    
    # Remove singleton operons (?)
    if not keepsingles: nrfiltopdf = nrfiltopdf[nrfiltopdf.oplen>1]
    
    # Collect all the clusters that ended up in an operon
    opclusters = set()
    for i,(opnum,row) in enumerate(nrfiltopdf.iterrows()):
        for cnum in row.operon: opclusters.add(cnum)
    # Add the clusters that didn't make it into an operon    
    addclusters = set(keepclusters).difference(opclusters) # add all other clusters 
    addop = [{'operon': tuple([cnum]), 'intergenicdist': [0], 'oplen': 1} for cnum in addclusters]
    addop = pd.DataFrame(pd.DataFrame(addop, index=range(len(addop))))
    newopdf = pd.concat([nrfiltopdf, addop], ignore_index=True, axis=0)
    
    if avgoperon:
        avgemb = []
        keepcdf = cdf.copy().loc[keepclusters]
        with open(embfile, 'rb') as f: allemb = pickle.load(f)
        
        clusterembmatrix = np.array([allemb[row.repid] for cnum,row in keepcdf.iterrows()])
        clusterembmatrix = clusterembmatrix/np.linalg.norm(clusterembmatrix, axis=1, keepdims=True)
        cluster2row = {cnum: i for i,cnum in enumerate(keepcdf.index)}  
        for opnum,row in newopdf.iterrows():
            if row.oplen == 1: 
                avgemb.append(clusterembmatrix[cluster2row[row.operon[0]],:])
            else:
                embvec = np.array([clusterembmatrix[cluster2row[cnum],:] for cnum in row.operon])
                avgemb.append(embvec.mean(axis=0))
        newopdf.loc[:,'avgemb'] = avgemb
            
    return newopdf

def FindNNClusters(opfile, embfile, cdf, keepclusters, testembeddings, \
                   normsimilarity=True, keepsingles=False, thset='operon', nnset='operon'):
    '''
    Run the entire operon based pipeline, used primarily for testing/developing
    Input: operon file path, cluster embeddings file path, cluster metadata, list of 
    clusters allowed, test protein embeddings, parameter to keep singleton operons 
    in the database (default: False), set of clusters to use to calculate NN similarity
    threshold, set of clusters from which we find the NN Clusters
    Returns NN clusters dictionary, and the new operon database
    '''
    percentile = [99.999]
    print("Parsing the operon file")
    filtopdf = pd.read_pickle(opfile)
    nrfiltopdf = filtopdf.copy()
    nrfiltopdf.update(nrfiltopdf.operon.apply(lambda x: tuple(x)))
    nrfiltopdf.drop_duplicates('operon', inplace=True)
    
    # Remove singleton operons (?)
    if not keepsingles: nrfiltopdf = nrfiltopdf[nrfiltopdf.oplen>1]
    
    # Collect all the clusters that ended up in an operon
    opclusters = set()
    for i,(opnum,row) in enumerate(nrfiltopdf.iterrows()):
        for cnum in row.operon: opclusters.add(cnum)
    # Add the clusters that didn't make it into an operon    
    addclusters = set(keepclusters).difference(opclusters) # add all other clusters 
    addop = [{'operon': tuple([cnum]), 'intergenicdist': [0], 'oplen': 1} for cnum in addclusters]
    addop = pd.DataFrame(pd.DataFrame(addop, index=range(len(addop))))
    newopdf = pd.concat([nrfiltopdf, addop], ignore_index=True, axis=0)
    
    # Create the mapping: cluster -> operon
    print("Creating the cluster to operon mapping")
    cnum2op = {}
    for i,(opnum,row) in enumerate(newopdf.iterrows()):
        for cnum in row.operon:
            if cnum in cnum2op.keys():
                cnum2op[cnum].add(opnum)
            else:
                cnum2op[cnum] = {opnum}
                
    # Create the matrix to calculate the threshold
    print("Loading the cluster embeddings matrix")
    with open(embfile, 'rb') as f:
        allemb = pickle.load(f)
    if thset == 'operon':
        thclusters = deepcopy(opclusters)
    elif thset == 'all':
        thclusters = deepcopy(keepclusters)
    keepcdf = cdf.copy().loc[thclusters]
    clusterembmatrix = np.array([allemb[row.repid] for cnum,row in keepcdf.iterrows()])
    clusterembmatrix = clusterembmatrix/np.linalg.norm(clusterembmatrix, axis=1, keepdims=True)
    cluster2row = {cnum: i for i,cnum in enumerate(keepcdf.index)}            
    
    print("Finding NN clusters")
    if nnset == 'operon':
        nnclusters = deepcopy(opclusters)
    elif nnset == 'all':
        nnclusters = deepcopy(keepclusters)
    
    totproteins = len(testembeddings)
    nndict = {queryid: [] for queryid in testembeddings.keys()}
    for i,(queryid,queryemb) in enumerate(testembeddings.items()):
        if i % 5*1e3 == 0: 
            print("Finished {:.2f}% of the test proteins".format(i/totproteins*100))
        dot_products = np.matmul(clusterembmatrix, queryemb)
        similarities = dot_products / np.linalg.norm(queryemb)  
        t = np.percentile(similarities, percentile)
        for cnum in nnclusters:
            similarity = similarities[cluster2row[cnum]]
            if similarity >= t:
                if (normsimilarity) and (min(similarities) != max(similarities)):
                    similarity = float(similarity - min(similarities)) / (max(similarities) - min(similarities))
                alloperons = list(cnum2op[cnum])
                nndict[queryid].append((similarity, cnum, alloperons))
    nonn = [k for k,v in nndict.items() if len(v)==0]
    print("Couldn't find NN clusters for {} test proteins".format(len(nonn)))
    return nndict, newopdf

def AssignOperons(opfile, embfile, cdf, keepclusters, testembeddings, \
                  normsimilarity=True, keepsingles=False, thset='operon', nnset='operon'):
    '''
    Assign query points to the most similar operon in the database
    Input: operon file path, cluster embeddings file path, cluster metadata, list of 
    clusters allowed, test protein embeddings, parameter to keep singleton operons 
    in the database (default: False), flag to calculate similarity threshold based on 
    operons or clusters (default: operons), flag to define unit of NNs (default: operons, which
    will exclude singletons)
    Returns assigned operons, and the new operon database
    '''
    percentile = [99.999]
    print("Parsing the operon file")
    newopdf = ParseOperonFile(opfile, embfile, keepclusters, cdf, keepsingles, avgoperon=True)
    
    # Create the mapping: cluster -> operon
    print("Creating the cluster to operon mapping")
    cnum2op = {}
    for i,(opnum,row) in enumerate(newopdf.iterrows()):
        for cnum in row.operon:
            if cnum in cnum2op.keys():
                cnum2op[cnum].add(opnum)
            else:
                cnum2op[cnum] = {opnum}
                
    # Create the matrix to calculate the threshold
    print("Creating the embeddings matrix")
    if thset == 'operon': # use only the nonsingleton operons to calculate the threshold
        thoperons = deepcopy(newopdf[newopdf.oplen>1].index)
    elif thset == 'all':
        thoperons = deepcopy(newopdf.index)
    embmatrix = np.array(newopdf.loc[thoperons,'avgemb'].to_list())
    embmatrix = embmatrix/np.linalg.norm(embmatrix, axis=1, keepdims=True)
    operon2row = {opnum: i for i,opnum in enumerate(thoperons)}
    
    print("Finding NN clusters")
    if nnset == 'operon': # assign onlythe nonsingleton operons
        nnoperons = deepcopy(newopdf[newopdf.oplen>1].index)
    elif nnset == 'all':
        nnoperons = deepcopy(newopdf.index)
    
    totproteins = len(testembeddings)
    nndict = {queryid: [] for queryid in testembeddings.keys()}
    for i,(queryid,queryemb) in enumerate(testembeddings.items()):
        if i % 5*1e3 == 0: 
            print("Finished {:.2f}% of the test proteins".format(i/totproteins*100))
        dotproducts = np.matmul(embmatrix, queryemb)
        similarities = dotproducts / np.linalg.norm(queryemb)  
        t = np.percentile(similarities, percentile)
        for opnum in nnoperons:
            similarity = similarities[operon2row[opnum]]
            if similarity >= t:
                if (normsimilarity) and (min(similarities) != max(similarities)):
                    similarity = float(similarity - min(similarities)) / (max(similarities) - min(similarities))
                nndict[queryid].append((similarity, opnum))
    nonn = [k for k,v in nndict.items() if len(v)==0]
    print("Couldn't assign operons for {} test proteins".format(len(nonn)))
    return nndict, newopdf

def PredictFromAverageOperons(cluster2go, operondf, nndict):
    """
    Predict GO terms from predicted operons
    Input: dictionary mapping protein clusters in operons to go terms, operon database
    as a DataFrame (operons are zero-based indexed), dictionary of protein clusters 
    from the operons database that are most similar to the query protein (key: protein ID)
    Returns a dictionary of predictions from the operon database
    """
    testsize = len(nndict)
    operon2len, operon2intergenicdist, operon2cluster = ExtractOperonStats(operondf)
    # predictions = {queryid.split('|')[1]: {} for queryid in nndict.keys()}
    predictions = dict()
    for i,(queryid,nn) in enumerate(nndict.items()):
        if i % 1e3 == 0:
            print("Predicted {:.2f}% of the test set".format(i/testsize*100))
        gopred = dict()
        operon2gofreq = dict()
        operon2sim = dict()
        for opsim,opnum in nn:
            operon2sim[opnum] = opsim
            operon2gofreq[opnum] = dict()
            clusters = operon2cluster[opnum]
            for cnum in clusters:
                for go in cluster2go.get(cnum,[]):
                    if go not in operon2gofreq[opnum].keys():
                        operon2gofreq[opnum][go] = float(1 / len(clusters))
                    else:
                        operon2gofreq[opnum][go] = operon2gofreq[opnum][go] + float(1 / len(clusters))
        for opsim,opnum in nn:
            for go,gofreq in operon2gofreq[opnum].items():
                if go not in gopred.keys():
                    gopred[go] = float(gofreq*opsim)
                else:
                    gopred[go] = max(gopred[go], float(gofreq*opsim))
        predictions[queryid.split('|')[1]] = gopred
    
    return predictions
