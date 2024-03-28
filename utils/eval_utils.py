#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import numpy as np
from copy import deepcopy

def calc_fmax(real_predictions, ytrue, thresholds, test_proteins, go_classes, 
              eval_go_class, remove_root=False, eval_mode='partial'):
    """
    Calculate Fmax, similar to CAFA challenges
    
    Parameters
    ----------
    real_predictions : dict
        Dict of predictions, maps protein IDs to GO terms assigned and their probability
    ytrue : dict
        Dict of real GO terms, maps protein IDs to GO terms
    thresholds : list
        List of prob thresholds to calculate the F1 score for
    test_proteins : list
        List of test protein IDs
    go_classes : list
        List of all GO terms in the ontology
    eval_go_class : str
        GO class that the Fmax will be calculated for 
    remove_root : bool
        Parameter to remove the root term. Default is False
    eval_mode : str
        Parameter to calculate recall over all proteins ('full') or only those that
        have a GO term ('partial') as it is done in CAFA. Default is 'partial'

    Returns
    -------
    f_max : float
        The maximum F1 score achieved over all thresholds
    pr : list
        List of all precision values calculated over all thresholds
    rc : list
        List of all recall values calculated over all thresholds
    coverage : float
        Prediction coverage at the threshold value that maximizes F1-score
    th_max : float
        The threshold value that maximizes F1-score
    """                
    num_test = len(test_proteins)
    predictions = deepcopy(real_predictions)
    # Remove obsolete terms
    obsolete_terms = {'GO:0052312', 'GO:1902586', 'GO:2000775'}
    # Remove the root term for evaluations
    if remove_root:
        root_terms = ['GO:0008150','GO:0005575','GO:0003674']
    else:
        root_terms = []
    root_terms = obsolete_terms.union(root_terms)  
    pr = []
    rc = []
    ff = []
    cov_list = []

    fmax = -np.inf
    num_prot = 0

    for t in thresholds:
        precision = 0.0
        recall = 0.0

        # Total number of test proteins for which any prediction was made and 
        # those with a prediction probability >= t
        num_predany = 0
        num_predthreshold = 0

        for test_prot in test_proteins:
            pred_flag = 0 # variable to check if there were any predictions for this query
            
            go_pred = deepcopy(predictions[test_prot]) # predictions
            go_true = set(deepcopy(ytrue[test_prot])) # real values
            go_true_class = set(filter(lambda x: (go_classes[x] == eval_go_class) and 
                                       (x not in root_terms), go_true))
            if len(go_true_class) == 0:
                continue
            intersection_size = 0
            num_prot += 1

            # Keep track of the number of terms that belong to the desired go class
            num_classterms = 0
            for go_term, go_prob in go_pred.items():
                if (go_classes[go_term] == eval_go_class) and (go_term not in root_terms):
                    num_classterms = num_classterms + 1

                    # Check threshold condition
                    if go_prob >= t:
                        pred_flag = 1
                        if go_term in go_true_class:
                            intersection_size = intersection_size + 1

            # If t = 0, override precision -- all terms in the ontology are predicted with probability 0
            if t == 0:
                precision = precision + (float(len(go_true_class)) / len(go_classes))
            elif (num_classterms > 0):
                precision = precision + (float(intersection_size) / num_classterms)
                
            if len(go_true_class) > 0:
                recall = recall + (float(intersection_size) / len(go_true_class))
            
            num_predthreshold = num_predthreshold + pred_flag
            if num_classterms > 0:
                num_predany = num_predany + 1

        if num_predthreshold > 0:
            precision = float(precision) / num_predthreshold
        else:
            precision = 0.0
        if eval_mode == 'full': # evaluate over all proteins that have a GO term
            recall = float(recall) / num_prot
        elif eval_mode == 'partial': # evaluate over proteins that we predicted only
            # Saw this option in the CAFA assessment tool
            recall = float(recall) / num_predany

        # If t is 0, override recall as 1
        if t == 0:
            recall = 1.0    
        
        pr.append(precision)
        rc.append(recall)
        cov_list.append(float(num_predany) / float(num_prot))

        if precision + recall > 0:
            ffcurr = (2 * precision * recall) / (precision + recall)
            ff.append((2 * precision * recall) / (precision + recall))
            fmax = max(fmax, ffcurr)
        else:
            ffcurr = 0.0
            ff.append(ffcurr)
            fmax = max(fmax, ffcurr)
    fmax_idx = np.argmax(ff)
    th_max = thresholds[fmax_idx]
    coverage = cov_list[fmax_idx]
    
    return fmax, pr, rc, coverage, th_max


def calc_smin(real_predictions, ytrue, thresholds, test_proteins, go_classes, 
              go_ics, eval_go_class, remove_root=False):
    """
    Calculate Smin, similar to CAFA challenges
    
    Parameters
    ----------
    real_predictions : dict
        Dict of predictions, maps protein IDs to GO terms assigned and their probability
    ytrue : dict
        Dict of real GO terms, maps protein IDs to GO terms
    thresholds : list
        List of prob thresholds to calculate the F1 score for
    test_proteins : list
        List of test protein IDs
    go_classes : list
        List of all GO terms in the ontology
    go_ics : dict
        Dict mapping GO terms to their information content
    eval_go_class : str
        GO class that the Fmax will be calculated for 
    remove_root : bool
        Parameter to remove the root term. Default is False

    Returns
    -------
    s_min : float
        The minimum semantic distance achieved over all thresholds
    ru_list : list
        List of all remaining uncertainty values calculated over all thresholds
    mi_list : list
        List of all mutual information values calculated over all thresholds
    """                
    import math
    num_test = len(test_proteins)
    predictions = deepcopy(real_predictions)
    # Remove obsolete terms
    obsolete_terms = {'GO:0052312', 'GO:1902586', 'GO:2000775'}
    # Remove the root term for evaluations
    if remove_root:
        root_terms = ['GO:0008150','GO:0005575','GO:0003674']
    else:
        root_terms = []
    root_terms = obsolete_terms.union(root_terms)  
    num_prot = 0 # number of proteins with a GO term in this ontology    
    ru_list = []
    mi_list = []
    smin = np.inf

    for t in thresholds:
        ru = 0.0
        mi = 0.0

        for test_prot in test_proteins:            
            go_true = set(deepcopy(ytrue[test_prot])) # real values
            go_true_class = set(filter(lambda x: (go_classes[x] == eval_go_class) and 
                                       (x not in root_terms), go_true))
            if len(go_true_class) == 0:
                continue
            num_prot = num_prot + 1
            # Fetch predictions for this test protein
            go_pred = deepcopy(predictions[test_prot]) # predictions
            go_pred_threshold = set()
            num_classterms = 0
            for go_term, go_prob in go_pred.items():
                if (go_classes[go_term] == eval_go_class) and (go_term not in root_terms):
                    num_classterms = num_classterms + 1
                    # Check threshold condition
                    if go_prob >= t:
                        go_pred_threshold.add(go_term)
            for true_term in go_true_class:
                if true_term not in go_pred_threshold:
                    ru = ru + go_ics[true_term]
            for pred_term in go_pred_threshold:
                if pred_term not in go_true_class:
                    mi = mi + go_ics[pred_term]
        ru = float(ru) / num_prot
        mi = float(ru) / num_prot
        
        if (ru > 0) or (mi > 0):
            scur = math.sqrt((ru ** 2) + (mi ** 2))
            smin = min(smin, scur)

        ru_list.append(ru)
        mi_list.append(mi)

    return smin, ru_list, mi_list

def calc_aupr(pr, rc):
    """
    Calculate aupr, area under the precision recall curve
    
    Parameters
    ----------
    pr : list
        List of all precision values calculated over all thresholds
    rc : list
        List of all recall values calculated over all thresholds

    Returns
    -------
    aupr : float
        Area under the precision recall curve
    """  
    pr_arr = np.array(pr)
    rc_arr = np.array(rc)

    sorted_idx = np.argsort(rc_arr)
    rc_arr_sorted = rc_arr[sorted_idx]
    pr_arr_sorted = pr_arr[sorted_idx]
    aupr = np.trapz(pr_arr_sorted, rc_arr_sorted)
    return aupr
