import numpy as np
from numpy.linalg import norm
from copy import deepcopy
import pandas as pd
import pickle

def NaivePredictor(tsvfile, trainproteins, testproteins):
    """
    Run the baseline naive predictor. It's simply the frequency of each GO term
    in the training set
    Input: proteins in the training set and proteins in the test set
    Returns naive predictions
    """
    from utils import SeqUtils
    annotmap = SeqUtils.ParseTSV(tsvfile)
    numtrainproteins = float(len(trainproteins))
    # Collect go term frequencies
    gopred = dict()
    for prot in trainproteins:
        for goterm in annotmap.get(prot,[]):
            if goterm in gopred.keys():
                gopred[goterm] = gopred[goterm] + float(1.0 / numtrainproteins)
            else:
                gopred[goterm] = float(1.0 / numtrainproteins)    
    # Create the predictions dictionary
    predictions = dict()
    for queryid in testproteins:
        predictions[queryid] = gopred
    
    return predictions

def PredictFromBlast(blastdf, annotmap, alltestIDs):
    """
    Run baseline BLAST predictor based on the max % identity approach
    Input: tabular output file from running blast against a database, dictionary
    mapping sequence IDs to go terms, list of sequence IDs of proteins in the test set
    Returns baseline BLAST predictions 
    """
    predictions = dict()
    for _,row in blastdf.iterrows():
        queryID = row['queryid']
        trainID = row['targetid']
        pident = float(row['pident']) / 100.0 
        if queryID not in predictions:
            predictions[queryID] = dict()       
        traingoterms = annotmap[trainID]
        for goterm in traingoterms:
            if goterm not in predictions[queryID]:
                predictions[queryID][goterm] = pident
            else:
                predictions[queryID][goterm] = max(pident, predictions[queryID][goterm])

    for queryID in alltestIDs:
        if queryID not in predictions:
            predictions[queryID] = dict()
                
    return predictions

def NormalizePredictions(predictions, goclasses):
    normpredictions = deepcopy(predictions)
    min_probability = dict()
    max_probability = dict()
    for category in ["MF", "CC", "BP"]:
        min_probability[category] = 2
        max_probability[category] = -2

    for test_protein_id, test_protein_predictions in normpredictions.items():
        # Determine range of prediction probabilities
        for go_term, probability in test_protein_predictions.items():
            category = go_classes[go_term]
            min_probability[category] = min(min_probability[category], probability)
            max_probability[category] = max(max_probability[category], probability)

    # Normalize per class
    for test_protein_id, test_protein_predictions in normpredictions.items():
        for go_term, probability in test_protein_predictions.items():
            category = go_classes[go_term]
            if (min_probability[category] < 2) and (abs(max_probability[category] - min_probability[category]) > 0.0000000001):
                test_protein_predictions[go_term] = (probability - min_probability[category]) / (max_probability[category] - min_probability[category])
    return normpredictions
