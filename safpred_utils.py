#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from numpy.linalg import norm
from copy import deepcopy
import pandas as pd
import pickle

def parse_database_file(db_path, emb_path, keep_clusters=None, keep_singletons=False, avg_synteny=False):
    '''
    Parse the database file, noticed that I need this more often and I thought I would >.>
    Input: database file path, embeddings file path, list of clusters allowed, parameter to keep singleton regions 
    in the database (default: False)
    Returns the synteny database
    '''
    # Remove duplicate regions
    db_df = pd.read_pickle(db_path)
    nr_db_df = db_df.copy()
    nr_db_df.update(nr_db_df.region.apply(lambda x: tuple(x)))
    nr_db_df.drop_duplicates('region', inplace=True)
    
    # Remove singleton regions
    if not keep_singletons: 
        nr_db_df = nr_db_df[nr_db_df.region_len > 1]
    
    # Collect all the clusters that ended up in a region
    syteny_clusters = set()
    for i, (rnum, row) in enumerate(nr_db_df.iterrows()):
        for cnum in row.region: 
            synteny_clusters.add(cnum)
    # Add the clusters that didn't make it into a region
    if not keep_clusters:
        keep_clusters = deepcopy(synteny_clusters)
    add_clusters = set(keep_clusters).difference(synteny_clusters) # add all other clusters 
    add_region = [{'region': tuple([cnum]), 'intergenic_dist': [0], 'region_len': 1} for cnum in add_clusters]
    add_region = pd.DataFrame(pd.DataFrame(add_region, index=range(len(add_region))))
    new_db_df = pd.concat([nr_db_df, add_region], ignore_index=True, axis=0)
    
    if avg_synteny:
        avg_emb = []
        # keep_cluster_df = cluster_df.copy().loc[keep_clusters]
        with open(emb_path, 'rb') as f: 
            all_emb = pickle.load(f)

        # Make the cluster embedding matrix from the clusters allowed
        cluster_emb_matrix = np.array([all_emb[cnum] for cnum in keep_clusters])
        cluster_emb_matrix = cluster_emb_matrix/np.linalg.norm(cluster_emb_matrix, axis=1, keepdims=True)
        cluster2row = {cnum: i for i, cnum in enumerate(keep_clusters)}  
        for rnum, row in new_db_df.iterrows():
            if row.region_len == 1: # only 1 cluster in the region
                avg_emb.append(cluster_emb_matrix[cluster2row[row.region[0]], :])
            else: # multiple clusters: take the average
                emb_vec = np.array([cluster_emb_matrix[cluster2row[cnum], :] for cnum in row.region])
                avg_emb.append(emb_vec.mean(axis=0))
        new_db_df.loc[:, 'avg_emb'] = avg_emb
            
    return new_db_df

def extract_db_stats(db_df):
    """
    Helper function to extract basic stats from synteny database dataframe
    Input: synteny database as a DataFrame 
    Returns two dictionaries mapping region IDs to region length and intergenic distance
    """
    df = db_df.copy()
    region2len = df.region_len.to_dict()
    region2intergenicdist = df.intergenic_dist.to_dict()
    region2cluster = df.region.to_dict()
    
    return region2len, region2intergenicdist, region2cluster    

def find_nn_clusters(db_path, emb_path, keep_clusters, test_embeddings, norm_sim=True, 
                     keep_singletons=False, th_set='synteny', nn_set='synteny'):
    '''
    Run the entire synteny based pipeline, used primarily for testing/developing
    Input: database file path, cluster embeddings file path, list of clusters allowed, 
    test protein embeddings, parameter to normalize cos similarity (default: True), 
    parameter to keep singleton regions in the database (default: False), 
    set of clusters to use to calculate NN similarity threshold ('synteny' or 'all'),
    set of clusters from which we find the NN Clusters ('synteny' or 'all')
    Returns NN clusters dictionary, and the new synteny database
    '''
    percentile = [99.999]
    print("Parsing the database file")
    db_df = pd.read_pickle(db_path)
    nr_db_df = db_df.copy()
    nr_db_df.update(nr_db_df.region.apply(lambda x: tuple(x)))
    nr_db_df.drop_duplicates('region', inplace=True)
    
    # Remove singleton regions
    if not keep_singletons: 
        nr_db_df = nr_db_df[nr_db_df.oplen > 1]
    
    # Collect all the clusters that ended up in a region
    syteny_clusters = set()
    for i, (rnum, row) in enumerate(nr_db_df.iterrows()):
        for cnum in row.region: 
            synteny_clusters.add(cnum)
    # Add the clusters that didn't make it into a region    
    add_clusters = set(keep_clusters).difference(synteny_clusters) # add all other clusters 
    add_region = [{'region': tuple([cnum]), 'intergenic_dist': [0], 'region_len': 1} for cnum in add_clusters]
    add_region = pd.DataFrame(pd.DataFrame(add_region, index=range(len(add_region))))
    new_db_df = pd.concat([nr_db_df, add_region], ignore_index=True, axis=0)
    
    # Create the mapping: cluster -> syntenic region
    print("Creating the cluster to syntenic region mapping")
    cnum2region = {}
    for i, (rnum, row) in enumerate(new_db_df.iterrows()):
        for cnum in row.region:
            if cnum in cnum2region.keys():
                cnum2region[cnum].add(rnum)
            else:
                cnum2region[cnum] = {rnum}
                
    # Create the matrix to calculate the NN threshold
    print("Loading the cluster embeddings matrix")
    with open(emb_path, 'rb') as f:
        all_emb = pickle.load(f)
    if th_set == 'synteny':
        th_clusters = deepcopy(synteny_clusters)
    elif th_set == 'all':
        th_clusters = deepcopy(keep_clusters)
    cluster_emb_matrix = np.array([all_emb[cnum] for cnum in keep_clusters])
    cluster_emb_matrix = cluster_emb_matrix/np.linalg.norm(cluster_emb_matrix, axis=1, keepdims=True)
    cluster2row = {cnum: i for i, cnum in enumerate(keep_clusters)}  
                         
    print("Finding NN clusters")
    if nn_set == 'synteny':
        nn_clusters = deepcopy(synteny_clusters)
    elif nn_set == 'all':
        nn_clusters = deepcopy(keep_clusters)
    
    num_test = len(test_embeddings)
    nn_dict = {query_id: [] for query_id in test_embeddings.keys()}
    for i, (query_id, query_emb) in enumerate(test_embeddings.items()):
        if i % 5*1e3 == 0: 
            print("Finished {:.2f}% of the test proteins".format(i/num_test*100))
        dot_products = np.matmul(cluster_emb_matrix, query_emb)
        similarities = dot_products / np.linalg.norm(query_emb)  
        t = np.percentile(similarities, percentile)
        for cnum in nn_clusters:
            similarity = similarities[cluster2row[cnum]]
            if similarity >= t:
                if (norm_sim) and (min(similarities) != max(similarities)):
                    similarity = float(similarity - min(similarities)) / (max(similarities) - min(similarities))
                all_regions = list(cnum2region[cnum])
                nn_dict[query_id].append((similarity, cnum, all_regions))
    
    no_nn = [k for k, v in nn_dict.items() if len(v) == 0]
    print("Couldn't find NN clusters for {} test proteins".format(len(no_nn)))
                         
    return nn_dict, new_db_df

def assign_regions(db_path, emb_path, keep_clusters, test_embeddings, norm_sim=True, 
                   keep_singletons=False, ths_et='synteny', nn_set='synteny'):
    '''
    Assign query points to the most similar region in the database
    Input: database file path, cluster embeddings file path, list of clusters allowed, 
    test protein embeddings, parameter to normalize cos similarity (default: True), 
    parameter to keep singleton regions in the database (default: False), 
    set of clusters to use to calculate NN similarity threshold ('synteny' or 'all'),
    set of clusters from which we find the NN Clusters ('synteny' or 'all')
    Returns assigned regions, and the new synteny database
    '''
    percentile = [99.999]
    print("Parsing the database file")
    new_db_df = parse_database_file(db_path, emb_path, keep_clusters, keep_singletons, avg_synteny=True)
    
    # Create the mapping: cluster -> syntenic region
    print("Creating the cluster to syntenic region mapping")
    cnum2region = {}
    for i, (rnum, row) in enumerate(new_db_df.iterrows()):
        for cnum in row.region:
            if cnum in cnum2region.keys():
                cnum2region[cnum].add(rnum)
            else:
                cnum2region[cnum] = {rnum}

    # Create the matrix to calculate the threshold
    print("Creating the embeddings matrix")
    if th_set == 'synteny': # use only the nonsingleton operons to calculate the threshold
        th_regions = deepcopy(new_db_df[new_db_df.region_len > 1].index)
    elif th_set == 'all':
        th_regions = deepcopy(new_db_df.index)
    emb_matrix = np.array(new_db_df.loc[th_regions, 'avg_emb'].to_list())
    emb_matrix = cluster_emb_matrix/np.linalg.norm(emb_matrix, axis=1, keepdims=True)
    region2row = {rnum: i for i, rnum in enumerate(th_regions)}                         
                       
    
    print("Finding NN clusters")
    if nn_set == 'synteny': # assign only the nonsingleton regions
        nn_regions = deepcopy(new_db_df[new_db_df.region_len > 1].index)
    elif nn_set == 'all':
        nn_regions = deepcopy(new_db_df.index)

    num_test = len(test_embeddings)
    nn_dict = {query_id: [] for query_id in test_embeddings.keys()}
    for i, (query_id, query_emb) in enumerate(test_embeddings.items()):
        if i % 5*1e3 == 0: 
            print("Finished {:.2f}% of the test proteins".format(i/num_test*100))
        dot_products = np.matmul(emb_matrix, query_emb)
        similarities = dot_products / np.linalg.norm(query_emb)  
        t = np.percentile(similarities, percentile)
        for rnum in nn_regions:
            similarity = similarities[region2row[rnum]]
            if similarity >= t:
                if (norm_sim) and (min(similarities) != max(similarities)):
                    similarity = float(similarity - min(similarities)) / (max(similarities) - min(similarities))
                nn_dict[query_id].append((similarity, rnum))
                
    no_nn = [k for k, v in nn_dict.items() if len(v) == 0]
    print("Couldn't assign regions for {} test proteins".format(len(no_nn)))
                       
    return nn_dict, new_db_df

def predict_from_avg_synteny(cluster2go, db_df, nn_dict):
    """
    Predict GO terms from the assigned regions
    Input: dictionary mapping protein clusters in regions to GO terms, synteny database
    as a DataFrame (regions are zero-based indexed), dictionary of protein clusters 
    from the database that are most similar to the query protein (key: protein ID)
    Returns a dictionary of predictions from the synteny database
    """
    num_test = len(nn_dict)
    region2len, _region2intergenic_dist, region2cluster = extract_db_stats(db_df)

    predictions = dict()
    for i, (query_id, nn) in enumerate(nn_dict.items()):
        if i % 1e3 == 0:
            print("Predicted {:.2f}% of the test set".format(i/num_test*100))
        go_pred = dict()
        region2gofreq = dict()
        region2sim = dict()
        for rsim, rnum in nn:
            region2sim[rnum] = rsim
            region2gofreq[rnum] = dict()
            clusters = region2cluster[rnum]
            for cnum in clusters:
                for go in cluster2go.get(cnum, []): # collect go terms for each cluster
                    if go not in region2gofreq[rnum].keys():
                        region2gofreq[rnum][go] = float(1 / len(clusters))
                    else:
                        region2gofreq[rnum][go] = region2gofreq[rnum][go] + float(1 / len(clusters))
        for rsim, rnum in nn:
            for go, gofreq in region2gofreq[rnum].items():
                if go not in gopred.keys():
                    gopred[go] = float(gofreq*rsim) # scale go term freq with the similarity
                else:
                    gopred[go] = max(gopred[go], float(gofreq*rsim))
        predictions[query_id] = gopred
    
    return predictions
