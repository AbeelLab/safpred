#!/usr/bin/env python3
# -*- coding: utf-8 -*-

def calc_ham_dist(seq1, seq2, to_ref=True):
    """
    Calculate Hamming distance between 2 sequences of varying lengths, trim if larger
    Input: 2 sequences, flag to use the first sequence as reference
    Returns: Hamming distance
    """
    d = sum(s1 != s2 for (s1, s2) in zip(seq1, seq2))
    if to_ref: d = d/len(seq1)
    return d
    
def calc_intergenic_dist(pos1, pos2, contig_len=None):
    """
    Calculate intergenic distance between two genes on the same contig, also 
    takes strandedness into account
    Input: positions of 2 genes (as intervals), circular flag (default: True)
    and contig length (necessary only for circular contigs)
    Returns the intergenic distance - or None if the genes are on different strands
    ref: Mihelčić, M., Šmuc, T. & Supek, F. Patterns of diverse gene functions 
    in genomic neighborhoods predict gene function and phenotype. 
    Sci Rep 9, 19537 (2019). https://doi.org/10.1038/s41598-019-55984-0
    """
    if pos1[-1] != pos2[-1]:
        return None
    if pos1[0] > pos2[1]:
        corpos1 = pos2
        corpos2 = pos1
    elif pos2[0] > pos1[1]:
        corpos1 = pos1
        corpos2 = pos2
    else:
        return 0
    if contig_len:
        d = min(abs(corpos2[0] - corpos1[1]), abs(corpos2[1] - contig_len - corpos1[0]))
    else:
        d = abs(corpos2[0] - corpos1[1])

    return d

def load_annot_file(annot_file_path):
    """
    Parse tabular annotation file mapping protein IDs to go terms
    Input: tsv file path
    Returns a dictionary mapping protein sequence IDs to go terms
    """
    import pandas as pd
    import numpy as np
    # Obsolete go terms are removed
    obsolete_terms = {'GO:0052312', 'GO:1902586', 'GO:2000775'}
    annot_df = pd.read_csv(annot_file_path, sep="\t", na_filter=False)
    annot_map = dict()

    for row in annot_df.itertuples():
        seq_id = row.id
        go_terms = np.array([str(i) for i in row.go_exp.split(';') if str(i) != obsolete_terms])
        # goterms = np.array([str(i) for i in row.go_exp.split(';')])
        annot_map[seq_id] = go_terms
    
    return annot_map
    
def get_seq_id(fasta_file_path):
    """
    Get sequence IDs from a fasta file
    Input: fasta file path
    Returns a list sequence IDs
    """
    from Bio import SeqIO

    seq_ids = []
    for rec in SeqIO.parse(fasta_file_path, 'fasta'):
        seq_id = rec.id
        seq_ids.append(seq_id)

    return seq_ids

def run_wget(cmd, verbose=False, *args, **kwargs):
    import subprocess

    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                               text=True, shell=True)
    std_out, std_err = process.communicate()
    if verbose:
        print(std_out.strip(), std_err)
    pass

def download_sprot(output_dir, release='current'):
    """
    Download a specific release of the SwissProt database, it will fetch both the protein sequences
    and their metadata

    Parameters
    ----------
    output_dir : str
        Directory to write the protein sequences in a fasta file and their metadata in a tabular file
        using the prefix "sprot_db_<release>"
    release : str
        SwissProt database release to download. Default is "current"

    Returns
    -------
    None.    
    """
    import gzip
    import os
    import pandas as pd
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq

    if release == 'current':
        sprot_url = 'ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz'
        go_url = 'http://purl.obolibrary.org/obo/go.obo'  
        
        cmd = 'wget -P {} {}'.format(os.path.join(output_dir, 'uniprot'), sprot_url)
        run_wget(cmd)
        cmd = 'wget -P {} {}'.format(os.path.join(output_dir, 'go'), go_url)
        run_wget(cmd)
    else:
        print("This script doesn't work for any releases other than current yet --- TODO")
        return
    
    
    exp_codes = ['EXP','IDA','IPI','IMP','IGI','IEP','HTP','HDA','HMP','HGI',
                 'HEP','IBA','IBD','IKR','IRD','IC','TAS']
    prot_list = []
    acc_list = []
    seq_list = []
    annot_list = []
    with gzip.open(os.path.join(output_dir, 'uniprot/uniprot_sprot.dat.gz'), 'rt') as f:
        prot = ''
        acc = ''
        seq = ''
        annots = []
        for line in f:
            flds = line.strip().split('   ')
            if (flds[0].lower() == 'id') and (len(flds) > 1): # found a new protein entry
                if len(prot) > 0: # but first, save the entry from previous iteration
                    prot_list.append(prot)
                    acc_list.append(acc)
                    seq_list.append(seq)
                    annot_list.append(annots)
                prot = flds[1]
                annots = []
                seq = ''
            elif (flds[0].lower() == 'ac') and (len(flds) > 1):
                acc = [f.strip(';') for f in flds[1].split()]
            elif (flds[0].lower() == 'dr') and (len(flds) > 1):
                flds = flds[1].split('; ')
                if flds[0].lower() == 'go':
                    go_id = flds[1]
                    go_code = flds[3].split(':')[0]
                    annots.append((go_id, go_code))
            elif flds[0].lower() == 'sq':
                seq = next(f).strip().replace(' ', '')
                while True:
                    sq = next(f).strip().replace(' ', '')
                    if sq == '//':
                        break
                    else:
                        seq += sq
        prot_list.append(prot)
        acc_list.append(acc)
        seq_list.append(seq)
        annot_list.append(annots)
    
    # Create a pandas dataframe to store the SwissProt metadata
    sprot_df = pd.DataFrame({'protein': prot_list, 'acc': acc_list, 'seq': seq_list, 'go': annot_list})
    keep_exp_prots = []
    keep_exp_annots = []
    for idx, row in sprot_df.iterrows():
        exp_annots = []
        for go_id, go_code in row.go:
            if go_code in exp_codes:
                exp_annots.append(go_id)
        if len(exp_annots) > 0:
            keep_exp_prots.append(idx)
            keep_exp_annots.append(';'.join(exp_annots))
    sprot_exp_df = sprot_df.loc[keep_exp_prots]
    sprot_exp_df.loc[:,'go_exp'] = keep_exp_annots
    sprot_exp_df.to_csv(os.path.join(output_dir, 'uniprot/sprot_db_current_metadata.tsv'), sep='\t', index=True, header=True)

    print("Writing the SwissProt database outputs")
    rec_list = []
    with open(os.path.join(output_dir, 'uniprot/sprot_db_current_metadata.tsv'), 'w') as f:
        f.write('id' + '\t' + 'go_exp' + '\n')
        for idx, row in sprot_exp_df.iterrows():
            for acc in row.acc:
                f.write(acc + '\t' + row.go_exp + '\n')
                rec = SeqRecord(Seq(row.seq), id=acc, name='', description='')
                rec_list.append(rec)
    SeqIO.write(rec_list, os.path.join(output_dir, 'uniprot/sprot_db_current.fasta'), 'fasta')

def blast_wrapper(test_seq_file, blast_db, e_value=0.001, word_size=11, 
                  out_fmt=['6', 'qaccver', 'saccver', 'qstart','qend', 'pident',
                           'length','bitscore','qcovs', 'evalue']):
    """    
    Python wrapper script to run blast against any database
    
    Parameters
    ----------
    test_seq_file : str
        Path to test sequences in fasta format
    blast_db : str
        Path to the blast database, created from the training set sequences
    e_value : float, optional
        E-value cutoff for the blast hits. The default is 0.001.
    word_size : int, optional
        Word size parameter in blastp. The default is 11.
    out_fmt : str, optional
        The output format in tabular form and the columns we want in the blast output. 
        The default is ['6', 'qaccver', 'saccver', 'qstart','qend', 'pident', 'length', 
                        'bitscore','qcovs', 'evalue'].

    Returns
    -------
    df : pandas.DataFrame
        Blast hits stored in a neat pandas dataframe, if the blast ran successfully
    """
    import subprocess, re
    import pandas as pd

    out_fmt = ' '.join(out_fmt)
    cmd = ('blastp', '-task', 'blastp-short', '-db', blast_db, '-query', test_seq_file, 
           '-outfmt', out_fmt, '-evalue', str(e_value), '-word_size', str(word_size))
    output = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, \
                            universal_newlines=True)
    if output.returncode == 0:
        col_names = re.sub('\d', '', out_fmt).strip().split(' ')
        result = [x.split('\t') for x in output.stdout.split('\n') if x]
        df = pd.DataFrame(result, columns=col_names)
        float_cols = {'pident', 'bitscore', 'qcovs', 'evalue'}.intersection(col_names)
        int_cols = {'length', 'qstart', 'qend'}.intersection(col_names)
        df.loc[:, float_cols] = df[float_cols].applymap(float)
        df.loc[:, int_cols] = df[int_cols].applymap(int)
        return df
    else:
        print("Ooops something went wrong")
        return 
    