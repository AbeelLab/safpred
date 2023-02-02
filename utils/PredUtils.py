import numpy as np
from copy import deepcopy
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
    minprob = dict()
    maxprob = dict()
    for category in ["MF", "CC", "BP"]:
        minprob[category] = 2
        maxprob[category] = -2

    for queryprotID, querypred in normpredictions.items():
        # Determine range of prediction probabilities
        for goterm, probability in querypred.items():
            category = goclasses[goterm]
            minprob[category] = min(minprob[category], probability)
            maxprob[category] = max(maxprob[category], probability)

    # Normalize per class
    for queryprotID, querypred in normpredictions.items():
        for goterm, probability in querypred.items():
            category = goclasses[goterm]
            if (minprob[category] < 2) and (abs(maxprob[category] - minprob[category]) > 0.0000000001):
                normpredictions[goterm] = (probability - minprob[category]) / (maxprob[category] - minprob[category])
    return normpredictions

def GetEmbeddings(testfile, trainfile, percentile=99.999):
    """
    Extract embeddings for train and test proteins from pickle files
    Input: test and train embeddings pickle file paths
    Returns train and test embeddings dicitonary mapping protein IDs to embedding vectors
    """
    with open(testfile, "rb") as f:
        testembeddings = pickle.load(f)
    with open(trainfile, "rb") as f:
        trainembeddings = pickle.load(f)
    return testembeddings, trainembeddings

def CreateTrainMatrix(trainembeddings):
    numemb = len(trainembeddings)
    embsize = len(list(trainembeddings.values())[0])

    embmatrix = np.zeros((numemb, embsize))

    ID2row = dict()

    for i, (trainprotID,trainemb) in enumerate(trainembeddings.items()):
        embmatrix[i] = trainemb
        ID2row[trainprotID] = i

    return embmatrix/np.linalg.norm(embmatrix, axis=1, keepdims=True), ID2row  

def RunKNN(tsvfile, testfile, trainfile, percentile=99.0):
    from utils import SeqUtils
    annotmap = SeqUtils.ParseTSV(tsvfile)
    testembeddings, trainembeddings = GetEmbeddings(testfile, trainfile)

    embeddingmatrix, ID2row = CreateTrainMatrix(trainembeddings)
    
    # Iterate over test set and make predictions
    testsize = len(testembeddings)
    predictions = dict()

    for i,(queryprotID, queryemb) in enumerate(testembeddings.items()):
        if i % 1e3 == 0:
            print("Predicted {} out of {} -- {} to go".format(i,testsize,testsize-i))
   
        # To each test protein, associate a dictionary of GO terms and their associated probabilities
        gopred = dict()

        # Determine similarities
        dotproducts = np.matmul(embeddingmatrix, queryemb)
        similarities = dotproducts / np.linalg.norm(queryemb)

        t = np.percentile(similarities, percentile)

        for trainprotID, trainemb in trainembeddings.items():
            similarity = similarities[ID2row[trainprotID]]

            # Only consider neighbors with similarity over the threshold
            if similarity >= t:
                # Fetch GO terms associated with train protein
                trainprotgoterms = annotmap[trainprotID]

                # Iterate over the train protein go terms
                for goterm in trainprotgoterms:
                    if goterm in gopred:
                        # If this GO term has been predicted before, store the maximum similarity
                        gopred[goterm] = max(gopred[goterm], similarity)
                    else:
                        gopred[goterm] = similarity

        # Make prediction for this test protein
        predictions[queryprotID] = gopred
    return predictions
