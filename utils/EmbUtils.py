
def ExtractT5Embeddings(infile, outfile, modeldir='../data/models/prottrans_t5_xl_u50/''):
    import os,sys,re
    import numpy as np
    from Bio import SeqIO
    from bio_embeddings.embed import ProtTransT5XLU50Embedder
    t5xlu50 = ProtTransT5XLU50Embedder(model_directory=modeldir, \
                                       half_model=True, device='cpu')
    fout = open(outfile, 'w+')
    for i,rec in enumerate(SeqIO.parse(inFile, 'fasta')):
        seq = str(rec.seq)
        seq = seq.replace('*','')
        recid = rec.id
        emb = t5xlu50.reduce_per_protein(t5xlu50.embed(seq))
        f.write("{}\t{}\n".format(recid,' '.join(map(str,emb)))

def ExtractESMEmbeddings(infile, outfile, maxlen=1000, modeldir='../data/models/esm1b'):

    import os,sys,re
    import numpy as np
    from Bio import SeqIO
    from bio_embeddings.embed import ESM1bEmbedder
    esm1b = ESM1bEmbedder(model_directory=modeldir, device='cpu')

    fout = open(outfile, 'w+')
    for i,rec in enumerate(SeqIO.parse(infile, 'fasta')):
        seq = str(rec.seq)
        seq = seq.replace('*','')
        recid = rec.id
        if len(seq) <= maxlen:
            emb = esm1b.embed(seq)
            emb = esm1b.reduce_per_protein(emb)
        else:
            chunks = [seq[i:i + maxlen] for i in range(0, len(seq), maxlen)]
            emb = np.vstack([esm1b.embed(chunk) for chunk in chunks])
            emb = esm1b.reduce_per_protein(emb)
        f.write("{}\t{}\n".format(recid,' '.join(map(str,emb)))

