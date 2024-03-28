import numpy as np
from copy import deepcopy
import pickle

def naive_predictor(annot_file_path, train_proteins, test_proteins):
    """
    Run the baseline naive predictor. It's simply the frequency of each GO term
    in the training set
    
    Parameters
    ----------
    annot_file_path : str
        Path to tabular annotation file, maps protein IDs to GO term annotations
    train_proteins : list
        List of proteins in the training set
    test_proteins : list
        List of proteins in the test set

    Returns
    -------
    predictions : dict
        A dictionary of predictions, maps protein IDs to predicted GO terms and their probability
    """
    
    from utils import seq_utils
    annot_map = seq_utils.load_annot_file(annot_file_path)
    num_trainproteins = float(len(train_proteins))
    # Collect go term frequencies
    go_pred = dict()
    for prot in train_proteins:
        for go_term in annot_map.get(prot, []):
            if go_term in go_pred.keys():
                go_pred[go_term] = go_pred[go_term] + float(1.0 / num_trainproteins)
            else:
                go_pred[go_term] = float(1.0 / num_trainproteins)    
    # Create the predictions dictionary
    predictions = dict()
    for query_id in test_proteins:
        predictions[query_id] = go_pred
    
    return predictions

def PredictFromBlast(blast_df, annot_file_path, test_proteins):
    """
    Run the baseline BLAST predictor based on the max % identity approach    
    Parameters
    ----------
    blast_df : pandas Dataframe
        Dataframe of running blastp against the training set
    annot_file_path : str    
        Path to tabular annotation file, maps protein IDs to GO term annotations
    test_proteins : list
        List of proteins in the test set

    Returns
    -------
    predictions : dict
        A dictionary of predictions, maps protein IDs to predicted GO terms and their probability
    """
    from utils import seq_utils
    annot_map = seq_utils.load_annot_file(annot_file_path)
    
    predictions = dict()
    for _, row in blast_df.iterrows():
        query_id = row['queryid']
        train_id = row['targetid']
        pident = float(row['pident']) / 100.0 
        if query_id not in predictions:
            predictions[query_id] = dict()       
        train_go_terms = annot_map[train_id]
        for go_term in train_go_terms:
            if go_term not in predictions[query_id]:
                predictions[query_id][go_term] = pident
            else:
                predictions[query_id][go_term] = max(pident, predictions[query_id][go_term])

    for query_id in test_proteins:
        if query_id not in predictions:
            predictions[query_id] = dict()
                
    return predictions

def normalize_prediction(predictions, go_classes):
    """ Helper function to normalize predictions within a GO ontology class
    """
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
