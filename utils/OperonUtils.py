def EditOperonDB(outfile, clusterdbpath, operondbpath, keepclusters, maxintergenicdist=300, addsingletons=True):

    import numpy as np
    import pandas as pd
    import pickle
    from copy import deepcopy
    from tqdm import tqdm
    from Bio import SeqIO

    from utils import SeqUtils

    tqdm.pandas()


    cdf = pd.read_pickle('clusterdbpath')
    cdf.update(cdf.members.apply(lambda x: x.split(';')))
    cdf.update(cdf.members_contig.apply(lambda x: x.split(';')))
    cdf.update(cdf.members_genbank.apply(lambda x: x.split(';')))

    # Load cluster metadata
    print("Loading base stats")
    with open('../data/operondb/basestats.pkl.gz', 'rb') as f:
        basestats = pickle.load(f)
    gene2pos = basestats['gene2pos']
    gene2cluster = basestats['gene2cluster']
    gene2contig = basestats['gene2contig']
    contig2len = basestats['contig2len']
    cluster2genes = basestats['cluster2genes']
    allclusters = set(cluster2genes.keys())

    print("Loading the original operon database")
    opdf = pd.read_pickle(operondbpath)
    nropdf = opdf.copy()
    nropdf.update(nropdf.operon.apply(lambda x: tuple(x)))
    nropdf.drop_duplicates('operon', inplace=True)

    print("Removing the clusters (and their operons) from the original database")
    opdict = {}
    totdone = len(nropdf)
    totremove = 0
    for i,(opnum,row) in enumerate(nropdf.iterrows()):
        if i % 1e4 == 0:
            # print("Finished {:.2f}% of the potential operons".format(i/totdone*100))
        keepidx = []
        keepcnum = []
        for idx,cnum in enumerate(row.operon):
            if cnum in keepclusters:
                keepidx.append(idx)
                keepcnum.append(cnum)
        d = deepcopy(row.minintergenicdist)
        if len(keepidx)==1:
            newd = [0]
        elif len(keepidx)==2:
            newd = [sum(d[keepidx[0]:keepidx[1]])]
        else:
            newd = []
            for idx1,idx2 in zip(keepidx,keepidx[1:]):
                newd.append(sum(d[idx1:idx2]))
        oplen = len(keepidx)
        if oplen > 0:
            opdict[opnum] = {'operon': keepcnum, 'oplen': oplen, 'minintergenicdist': newd}
        else:
            totremove = totremove + 1
            if totremove % 1e4 == 0:
                print("Removing operon {} --- removed {} operons so far".format(opnum,totremove))
    smallopdf = pd.DataFrame(opdict.values(), index=opdict.keys())


    # Break up operons with large intergenic distance
    print("Breaking up operons with large intergenic distance")
    totop = len(smallopdf)
    opcount = 0
    newopdict = {}
    for i,(idx,row) in enumerate(smallopdf.iterrows()):
        if i % 1e4==0:
            print("Finished {:.2f}% of the potential operons".format(i/totop*100))
        addop = []
        addintd = []
        for j,(c1,c2) in enumerate(zip(row.operon,row.operon[1:])):
            intd = row.minintergenicdist[j]
            if intd <= maxdist:
                if len(addop)==0:
                    addop.append(c1)
                addop.append(c2)
                addintd.append(intd)
            else:
                if len(addop)>0:
                    newopdict[opcount] = {'operon': addop, 'intergenicdist': addintd}
                    addop = []
                    addintd = []
                    opcount = opcount+1
        if len(addop)>0:
            newopdict[opcount] = {'operon': addop, 'intergenicdist': addintd}
            opcount = opcount+1
    newopdf = pd.DataFrame(newopdict.values(), index=newopdict.keys())
    if addsingletons:
    # Need to add the singleton operons
    singletons = smallopdf[smallopdf.apply(lambda x: (x.oplen==1) and (x.operon[0] in keepclusters), axis=1)].copy()
    ss = singletons.copy()
    ss.update(ss.operon.apply(lambda x: tuple(x)[0]))
    ss.drop_duplicates('operon', inplace=True)
    ss.update(ss.operon.apply(lambda x: [x]))
    ss.rename(columns={'minintergenicdist': 'intergenicdist'}, inplace=True)
    newopdf.loc[:,'oplen'] = newopdf.operon.apply(lambda x: len(x))
    newopdf = newopdf.append(ss, ignore_index=True)
    
    print("Saving the new operon database")
    newopdf.to_pickle(outfile)
