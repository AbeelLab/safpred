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

def blast_baseline(blast_df, annot_file_path, test_proteins):
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

def normalize_predictions(predictions, go_classes):
    """ Helper function to normalize predictions within a GO ontology class
    """
    norm_predictions = deepcopy(predictions)
    min_prob = dict()
    max_prob = dict()
    for category in ["MF", "CC", "BP"]:
        min_prob[category] = 2
        max_prob[category] = -2

    for query_id, query_pred in norm_predictions.items():
        # Determine range of prediction probabilities
        for go_term, go_prob in query_pred.items():
            category = go_classes[go_term]
            min_prob[category] = min(min_prob[category], go_prob)
            max_prob[category] = max(max_prob[category], go_prob)

    # Normalize per class
    for query_id, query_pred in norm_predictions.items():
        for go_term, go_prob in query_pred.items():
            category = go_classes[go_term]
            if (min_prob[category] < 2) and (abs(max_prob[category] - min_prob[category]) > 0.0000000001):
                norm_predictions[go_term] = (go_prob - min_prob[category]) / (max_prob[category] - min_prob[category])
    return norm_predictions

def load_embeddings(train_file_path, test_file_path, percentile=99.999):
    """
    Extract embeddings for train and test proteins from pickle files
    Input: test and train embeddings pickle file paths
    Returns train and test embeddings dictionary mapping protein IDs to embedding vectors
    """
    with open(train_file_path, "rb") as f:
        train_embeddings = pickle.load(f)
    with open(test_file_path, "rb") as f:
        test_embeddings = pickle.load(f)        
    return train_embeddings, test_embeddings

def create_train_matrix(train_embeddings):
    num_emb = len(train_embeddings)
    emb_size = len(list(train_embeddings.values())[0])

    emb_matrix = np.zeros((num_emb, emb_size))

    id2row = dict()

    for i, (train_id, train_emb) in enumerate(train_embeddings.items()):
        emb_matrix[i] = train_emb
        id2row[train_id] = i

    return emb_matrix/np.linalg.norm(emb_matrix, axis=1, keepdims=True), id2row  

def SAFPrednn(annot_file_path, train_file_path, test_file_path, percentile=99.0):
    from utils import seq_utils
    
    annot_map = seq_utils.load_annot_file(annot_file_path)
    train_embeddings, test_embeddings = load_embeddings(train_file_path, test_file_path)

    train_emb_matrix, id2row = create_train_matrix(train_embeddings)
    
    # Iterate over test set and make predictions
    num_test = len(test_embeddings)
    predictions = dict()

    for i, (query_id, query_emb) in enumerate(test_embeddings.items()):
        if i % 1e3 == 0:
            print("Predicted {} out of {} -- {} to go".format(i, num_test, num_test - i))
   
        # To each test protein, associate a dictionary of GO terms and their associated probabilities
        go_pred = dict()

        # Determine similarities
        dot_products = np.matmul(train_emb_matrix, query_emb)
        similarities = dot_products / np.linalg.norm(query_emb)

        t = np.percentile(similarities, percentile)

        for train_id, train_emb in train_embeddings.items():
            similarity = similarities[id2row[train_id]]

            # Only consider neighbors with similarity over the threshold
            if similarity >= t:
                # Fetch GO terms associated with train protein
                train_prot_goterms = annot_map[train_id]

                # Iterate over the train protein go terms
                for go_term in train_prot_goterms:
                    if go_term in go_pred:
                        # If this GO term has been predicted before, store the maximum similarity
                        go_pred[go_term] = max(go_pred[go_term], similarity)
                    else:
                        go_pred[go_term] = similarity

        # Make prediction for this test protein
        predictions[query_id] = go_pred
        
    return predictions
