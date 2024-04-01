#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  1 17:03:26 2024

@author: aysun
"""

import os
import pickle
from utils import pred_utils, seq_utils

train_seq_file = 'data/input/train_seq_file.fasta'
test_seq_file = 'data/input/test_seq_file.fasta'
annot_file_path = 'data/uniprot/sprot_db_current_metadata.tsv'
predictions_path = 'out/blast_predictions_normalized.pkl'

print("Create blast database from the training fasta sequences")
cmd = 'makeblastdb -in {} -parse_seqids -dbtype prot'.format(train_seq_file)
os.system(cmd)

print("Run blastp")
blast_db = 'data/input/train_seq_file'
test_proteins = seq_utils.get_seq_id(test_seq_file)
with open('data/go/go_classes.pkl', 'rb') as f:
    go_classes = pickle.load(f)

blast_df = seq_utils.blast_wrapper(test_seq_file, blast_db)
predictions = pred_utils.blast_baseline(blast_df, annot_file_path, test_proteins)

prop_predictions = pred_utils.propagate_predictions(predictions, remove_root=True)    
norm_predictions = pred_utils.normalize_predictions(prop_predictions, go_classes)

print("Saving the normalized predictions to {}".format(predictions_path))
with open(predictions_path, 'wb') as f:
    pickle.dump(norm_predictions, f)
