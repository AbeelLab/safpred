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
    # Obsolete go terms are removed
    obsolete_terms = {'GO:0052312', 'GO:1902586', 'GO:2000775'}
    annot_df = pd.read_csv(annot_file_path, sep="\t", na_filter=False)
    annot_map = dict()

    for row in annot_df.itertuples():
        seq_id = row.id
        go_terms = np.array([str(i) for i in row.go_exp.split(';') if str(i) != obsolete_terms])
        # goterms = np.array([str(i) for i in row.go_exp.split(';')])
        annot_map[seq_id] = go_terms
    
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
