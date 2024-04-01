import pickle
import numpy as np    
from Bio import SeqIO

def esm1b_embed(in_path, model_dir, output_path=None, max_len=1000):
    """ 
    Extract ESM-1b embeddings for input aminoacid sequences
    Parameters
    ----------
    in_path : str
        path to fasta file storing input sequences
    model_dir : str
        Directory where pretrained model weights are stored
    output_path : str, optional
        path to write the extracted embedding vectors. The default is None.
    max_len : int, optional
        maximum sequence length. The default is 1000.

    Returns
    -------
    Embedding vectors
    """
    from bio_embeddings.embed import ESM1bEmbedder
    print("Loading the ESM-1b model")
    # Load the ESM-1b model from model_dir so we avoid downloading the weights
    # from scratch each time
    esm1b = ESM1bEmbedder(model_directory=model_dir, device='cpu')

    enc_seqs = []
    seq_ids = []
    
    print("Extracting embeddings")
    for i, rec in enumerate(SeqIO.parse(in_path, 'fasta')):
        seq = str(rec.seq)
        seq = seq.replace('*','') # remove any missing aminoacids
        seq_id = rec.id
        if i % 1e2 == 0:
            print("Processed {} sequences so far".format(i + 1))
        if len(seq) <= max_len:
            enc_seq = esm1b.embed(seq)
        else:
            # ESM-1b input is limited to 1024 aa, so we need to break long 
            # protein sequences into smaller chunks
            chunks = [seq[i:i + max_len] for i in range(0, len(seq), max_len)]
            enc_seq = np.vstack([esm1b.embed(chunk) for chunk in chunks])
            enc_seq = esm1b.reduce_per_protein(enc_seq)
        enc_seq = esm1b.reduce_per_protein(enc_seq)
        enc_seqs.append(enc_seq)
        seq_ids.append(seq_id)

    enc_dict = {'seq_ids': seq_ids, 'enc_seqs': enc_seqs}
    if output_path:
        with open(output_path, 'wb') as f:
            pickle.dump(enc_dict, f)
        print("Pickled embedding vectors in " + output_path)
        
    return enc_dict

def t5xlu50_embed(in_path, model_dir, output_path=None):
    """ 
    Extract T5-XL-U50 embeddings for input aminoacid sequences
    Parameters
    ----------
    in_path : str
        path to fasta file storing input sequences
    model_dir : str
        Directory where pretrained model weights are stored
    output_path : str, optional
        path to write the extracted embedding vectors. The default is None.

    Returns
    -------
    Embedding vectors
    """
    from bio_embeddings.embed import ProtTransT5XLU50Embedder
    
    print("Loading the ESM-1b model")
    # Load the T5-XL-U50 model from model_dir so we avoid downloading the weights
    # from scratch each time
    t5xlu50 = ProtTransT5XLU50Embedder(model_directory=model_dir, 
                                       half_model=True, device='cpu')

    enc_seqs = []
    seq_ids = []
    
    print("Extracting embeddings")
    for i, rec in enumerate(SeqIO.parse(in_path, 'fasta')):
        seq = str(rec.seq)
        seq = seq.replace('*','') # remove any missing aminoacids
        seq_id = rec.id
        if i % 1e2 == 0:
            print("Processed {} sequences so far".format(i + 1))
        enc_seq = t5xlu50.reduce_per_protein(t5xlu50.embed(seq))
        enc_seqs.append(enc_seq)
        seq_ids.append(seq_id)

    enc_dict = {'seq_ids': seq_ids, 'enc_seqs': enc_seqs}
    if output_path:
        with open(output_path, 'wb') as f:
            pickle.dump(enc_dict, f)
        print("Pickled embedding vectors in " + output_path)
        
    return enc_dict
