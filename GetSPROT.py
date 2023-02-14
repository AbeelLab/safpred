def GetSPROT(release='current'):
    import gzip, wget
    import pandas as pd
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord
    from utils.MiscUtils import Ontology
        expcodes = ['EXP','IDA','IPI','IMP','IGI','IEP','HTP','HDA','HMP','HGI',\
            'HEP','IBA','IBD','IKR','IRD','IC','TAS']
        protlist, acclist, seqlist, annotlist = [], [], [], []
        with gzip.open(infile, 'rt') as f:
            prot, acc, seq, annots = '', '', '', []
            for line in f:
                flds = line.strip().split('   ')
                if flds[0].lower()=='id' and len(flds)>1: # found a new entry
                    if len(prot)>0: # save the entry from previous iteration
                        protlist.append(prot)
                        acclist.append(acc)
                        seqlist.append(seq)
                        annotlist.append(annots)
                    prot = flds[1]
                    annots = []
                    seq = ''
                elif flds[0].lower()=='ac' and len(flds)>1:
                    acc = [f.strip(';') for f in flds[1].split()]
                elif flds[0].lower()=='dr' and len(flds)>1:
                    flds = flds[1].split('; ')
                    if flds[0].lower()=='go':
                        goid = flds[1]
                        gocode = flds[3].split(':')[0]
                        annots.append((goid,gocode))
                elif flds[0].lower()=='sq':
                    seq = next(f).strip().replace(' ', '')
                    while True:
                        sq = next(f).strip().replace(' ', '')
                        if sq == '//':
                            break
                        else:
                            seq += sq
            protlist.append(prot)
            acclist.append(acc)
            seqlist.append(seq)
            annotlist.append(annots)
        
        # Create a pandas dataframe
        sprotdf = pd.DataFrame({'protein': protlist, 'acc': acclist, 'seq': seqlist, 'go': annotlist})
        keepexp = []
        keepexpannots = []
        for idx, row in sprotdf.iterrows():
            expannots = []
            for goid,gocode in row.go:
                if gocode in expcodes:
                    expannots.append(goid)
            if len(expannots)>0:
                keepexp.append(idx)
                keepexpannots.append(';'.join(expannots))
        expsprotdf = sprotdf.loc[keepexp]
        expsprotdf.loc[:,'go_exp'] = keepexpannots
        return expsprotdf
    
    if release=='current':
        sproturl = 'ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz'
        GOurl = 'http://purl.obolibrary.org/obo/go.obo'       
        wget.download(sproturl, out='uniprot_sprot.dat.gz')
        wget.download(GOurl, out='go.obo')
        
    expsprotdf = LoadSPROT('uniprot_sprot.dat.gz')
    reclist = []
    with open('uniprot2go_exp.tsv', 'w') as f:
        f.write('id' + '\t' + 'go_exp' + '\n')
        for idx,row in expsprotdf.iterrows():
            for acc in row.acc:
                f.write(acc + '\t' + row.go_exp + '\n')
                rec = SeqRecord(row.seq, id=acc, name='', description='')
                reclist.append(rec)
    SeqIO.write(reclist, 'uniprot_sprot.fasta', 'fasta')
