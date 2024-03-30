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
                query_pred[go_term] = (go_prob - min_prob[category]) / (max_prob[category] - min_prob[category])
    
    return norm_predictions

def propagate_predictions(predictions, remove_root=False, go_path='data/go/go.obo', name_space=False):
    """
    Propagate GO predictions, retains the max confidence   
    Parameters
    ----------
    predictions : dict
        A dictionary of predictions, maps protein IDs to predicted GO terms and their probability
    remove_root : bool
        Parameter to remove the root term. Default is False
    go_path : str
        Path to the GO ontology obo file. Default is 'data/go/go.obo'      
    name_space : str    
        If you want to want to propagate within this namespace only. Default is False, i.e.
        propagates for all namespaces

    Returns
    -------
    predictions : dict
        A dictionary of propagated predictions, maps protein IDs to predicted GO terms and their probability
    """
    from utils.misc_utils import Ontology
    
    # Remove obsolete terms
    obsolete_terms = {'GO:0052312', 'GO:1902586', 'GO:2000775'}
    # Remove the root term for evaluations
    if remove_root:
        root_terms = ['GO:0008150','GO:0005575','GO:0003674']
    else:
        root_terms = []
    root_terms = obsolete_terms.union(root_terms)    
                  
    prop_predictions = dict()
    go = Ontology(go_path, with_rels=True)
    for test_prot, test_predictions in predictions.items():
        new_predictions = deepcopy(test_predictions)
        for go_term, go_prob in test_predictions.items():
            anc_go_set = set([go_term])
            anc_go_set |= go.get_ancestors(go_term)
            if name_space: # propagate within the ontology of this NS
                cur_ns = go.get_namespace(go_term)
                new_go_set = set([anc_go for anc_go in anc_go_set if go.get_namespace(anc_go) == cur_ns])
            else:
                new_go_set = anc_go_set
            if remove_root: # remove root terms
                new_go_set.difference(root_terms)
            for new_go in new_go_set:
                if new_go in new_predictions.keys():
                    new_predictions[new_go] = max(new_predictions[new_go], go_prob)
                else:
                    new_predictions[new_go] = go_prob
        prop_predictions[test_prot] = new_predictions
    
    return prop_predictions 

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

