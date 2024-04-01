#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  1 14:07:03 2024

@author: aysun
"""

import pickle
import numpy as np
from utils import eval_utils, seq_utils     

predictions_file = 'out/safpred_predictions_normalized.pkl'
annot_file = 'data/sprot_db_current_metadata.tsv'
test_seq_file = 'data/input/test_seq_ecoli_small.fasta'
go_class_file = 'data/go/go_classes.pkl'
go_ics_file = 'data/go/go_ics.pkl'
eval_go_class = 'bp'
remove_root = True
thresholds = np.linspace(0.0, 1.0, num=50)

test_proteins = seq_utils.get_seq_id(test_seq_file)

with open(go_class_file, 'rb') as f:
    go_classes = pickle.load(f)
with open(go_ics_file, 'rb') as f:
    go_ics= pickle.load(f)    
annot_map = seq_utils.load_annot_file(annot_file)
true_labels = {test_prot: annot_map[test_prot] for test_prot in test_proteins}
with open(predictions_file, 'rb') as f:
    predictions = pickle.load(f)


print("Calculating Fmax")
fmax, pr, rc, coverage, th_max = eval_utils.calc_fmax(predictions, true_labels, 
                                                      thresholds, test_proteins, 
                                                      go_classes, eval_go_class, 
                                                      remove_root=remove_root)

print("Fmax is {:.2f} at threshold {:.2f}".format(fmax, th_max))

print("Calculating Smin")
smin, ru_list, mi_list = eval_utils.calc_smin(predictions, true_labels, thresholds, 
                                              test_proteins, go_classes, go_ics, 
                                              eval_go_class, remove_root=remove_root)

print("Smin is {:.2f}".format(smin))
