#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import pickle
from copy import deepcopy
from tqdm import tqdm
from Bio import SeqIO
from utils import embed_utils, seq_utils, safpred_utils, pred_utils

tqdm.pandas()

def main():
    parser = argparse.ArgumentParser(description = 'run SAFPred')
    parser.add_argument('--train_seq', '-i', type=str, required=True, 
                        help = 'Path to the training set sequences file in fasta format')
    parser.add_argument('--test_seq', '-t', type=str, required=True, 
                        help = 'Path to the test set sequences file in fasta format')   
    parser.add_argument('--annot_file', '-a', type=str, required=True, 
                        help = 'Path to the annotation file, mapping all proteins (train and test) to GO terms')
    parser.add_argument('--output', '-o', type=str, required=True, 
                        help='Output directory where the prediction results will be stored')    
    parser.add_argument('--db_path', '-d', type=str, default='data/safpreddb/safpreddb.pkl.gz', 
                        help = 'Path to the synteny database file')
    parser.add_argument('--db_emb_path', '-b', type=str, default='data/safpreddb/safpreddb_emb.pkl', 
                        help = 'Path to the synteny database embeddings file')
    parser.add_argument('--emb_model', '-e', type=str, default='esm', 
                        help = 'Embedding model: "esm", "t5" or "none". which requires train and test embedding files as input')
    parser.add_argument('--train_emb_file', '-f', type=str, default=None, 
                        help = 'Path to the training embeddings, required only if you want to use your own')
    parser.add_argument('--test_emb_file', '-g', type=str, default=None, 
                        help = 'Path to the test embeddings, required only if you want to use your own')
    parser.add_argument('--nn_percentile', '-p', type=float, default=99.999, 
                        help = 'Percentile for the NN thresholds. Default is 99.999')
    parser.add_argument('--norm_sim', '-n', action='store_false', 
                        help = 'Normalize NN similarities. Default is True')
    parser.add_argument('--keep_singletons', '-k', action='store_false', 
                        help = 'Keep singleton entries in the synteny database. Default is True')    

    args = parser.parse_args()
    train_seq_file = args.train_seq
    test_seq_file = args.test_seq
    annot_file_path = args.annot_file
    db_path = args.db_path
    db_emb_path = args.db_emb_path
    emb_model = args.emb_model
    output_path = args.output

    predictions_path = os.path.join(output_path, 'safpred_predictions_normalized.pkl')

    if emb_model == 'none':
        if args.train_emb_file:
            print("Using pre-calculated embeddings")
            if args.test_emb_file:
                train_embeddings, test_embeddings = pred_utils.load_embeddings(args.train_emb_file, args.test_emb_file)
            else:
                print("You have to input test embeddings too! --- Mission aborted")
    else:
        if emb_model == 'esm':
            print("Extracting ESM-1b embeddings for the training sequences")
            train_embeddings = embed_utils.esm1b_embed(train_seq_file, 'data/models/esm1b/')
            print("Extracting ESM-1b embeddings for the test sequences")
            test_embeddings = embed_utils.esm1b_embed(test_seq_file, 'data/models/esm1b/')
        elif emb_model == 't5':
            print("Extracting T5-XL-U50 embeddings for the training sequences")
            train_embeddings = embed_utils.t5xlu50_embed(train_seq_file, 'data/models/t5/')
            print("Extracting T5-XL-U50 embeddings for the test sequences")
            test_embeddings = embed_utils.t5xlu50_embed(test_seq_file, 'data/models/t5/')
        else:
            print("Invalid embeddings model name --- Mission aborted")  

    train_proteins = seq_utils.get_seq_id(train_seq_file)
    test_proteins = seq_utils.get_seq_id(test_seq_file)

    # Load cluster to GO term mapping dictionary
    with open('data/safpreddb/cluster2go.pkl', 'rb') as f:
        cluster2go = pickle.load(f)
    with open('data/go/go_classes.pkl', 'rb') as f:
        go_classes = pickle.load(f)
    
    nn_dict, db_df = safpred_utils.assign_regions(db_path, db_emb_path, test_embeddings=test_embeddings, 
                                                  keep_clusters=None, norm_sim=args.norm_sim, 
                                                  keep_singletons=args.keep_singletons, th_set='synteny', nn_set='synteny')
    safprednn_predictions = safpred_utils.safprednn(annot_file_path, train_embeddings, test_embeddings)
    safpredsynteny_predictions = safpred_utils.predict_from_avg_synteny(cluster2go, db_df, nn_dict)
    predictions = safpred_utils.combine_predictors(test_proteins, safprednn_predictions, safpredsynteny_predictions)
    prop_predictions = pred_utils.propagate_predictions(predictions, remove_root=True)    
    norm_predictions = pred_utils.normalize_predictions(prop_predictions, go_classes)

    print("Saving the normalized predictions to directory {}".format(output_path))
    with open(predictions_path, 'wb') as f:
        pickle.dump(norm_predictions, f)

if __name__ == '__main__':
    main()
