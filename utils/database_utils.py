import pandas as pd
from copy import deepcopy


def remove_redundant_entries(db_df):
    """ Helper function to remove redundant entries in a synteny database
    """
    nr_db_df = db_df.copy()
    nr_db_df.update(nr_db_df.region.apply(lambda x: tuple(x)))
    nr_db_df.drop_duplicates('region', inplace=True)
    return nr_db_df

def edit_database(indb_path, outdb_path, keep_clusters, max_intergenic_dist=300, add_singletons=True):
    """
    Edit an existing synteny database to remove gene clusters and the syntenic regions
    associated with them
    
    Parameters
    ----------
    indb_path : str
        Path to the existing synteny database, should be a pickle file
    outdb_path : str
        Path to write the edited synteny database out in a pickle file
    keep_clusters : list
        List of clusters (or genes) to keep in the database
    max_intergenic_dist : int, optional
        Maximum intergenic distance allowed within a syntenic region. Default is 300
    add_singletons : bool, optional
        Flag to include singleton regions in the database. The default is True.

    Returns
    -------
    out_db_df : pandas DataFrame
        The new synteny database described in a pandas dataframe
    """

    print("Loading the input database and removing redundant entries")
    db_df = pd.read_pickle(indb_path)
    nr_db_df = db_df.copy()
    nr_db_df.update(nr_db_df.region.apply(lambda x: tuple(x)))
    nr_db_df.drop_duplicates('region', inplace=True)

    print("Removing the gene clusters (and their regions) from the input database")
    syn_dict = {}
    num_regions = len(nr_db_df)
    num_removed = 0
    for i, (rnum, row) in enumerate(nr_db_df.iterrows()):
        if i % 1e4 == 0:
            print("Finished {:.2f}% of the potential syntenic regions".format(i/num_regions*100))
        new_region = []
        new_clusters = []
        for idx, cluster_num in enumerate(row.region):
            if cluster_num in keep_clusters:
                new_region.append(idx)
                new_clusters.append(cluster_num)
        d = deepcopy(row.intergenic_dist)
        if len(new_region) == 1: # Only 1 gene in the region: set the distance to 0 
            new_dist = [0]
        elif len(new_region) == 2:
            new_dist = [sum(d[new_region[0]:new_region[1]])]
        else: # More than 2 entries: calculate pairwise distance for all
            new_dist = []
            for idx1, idx2 in zip(new_region, new_region[1:]):
                new_dist.append(sum(d[idx1:idx2]))
        region_len = len(new_region)
        if region_len > 0: # Found a region we can keep: save it to the synteny dict
            syn_dict[rnum] = {'region': new_clusters, 'region_len': region_len, 
                              'intergenic_dist': new_dist}
        else: # No members left in the region: removed from the database
            num_removed = num_removed + 1
            if num_removed % 1e4 == 0:
                print("Removing region {} --- removed {} regions so far".format(rnum, num_removed))
    reduced_db_df = pd.DataFrame(syn_dict.values(), index=syn_dict.keys())
    
    # Break up regions with large intergenic distance
    print("Breaking up regions with large intergenic distance")
    num_regions = len(reduced_db_df)
    region_count = 0
    out_syn_dict = {}
    for i, (rnum, row) in enumerate(reduced_db_df.iterrows()):
        if i % 1e4==0:
            print("Finished {:.2f}% of the potential syntenic regions".format(i/num_regions*100))
        add_region = []
        add_intergenic_dist = []
        for j, (c1, c2) in enumerate(zip(row.region, row.region[1:])):
            intergenic_dist = row.intergenic_dist[j]
            if intergenic_dist <= max_intergenic_dist:
                if len(add_region) == 0:
                    add_region.append(c1)
                add_region.append(c2)
                add_intergenic_dist.append(intergenic_dist)
            else:
                if len(add_region) > 0:
                    out_syn_dict[region_count] = {'region': add_region, 
                                                  'intergenic_dist': add_intergenic_dist}
                    add_region = []
                    add_intergenic_dist = []
                    region_count = region_count + 1
        if len(add_region) > 0:
            out_syn_dict[region_count] = {'region': add_region, 
                                          'intergenic_dist': add_intergenic_dist}
            region_count = region_count + 1
    out_db_df = pd.DataFrame(out_syn_dict.values(), index=out_syn_dict.keys())
    
    if add_singletons:
    # Need to add the singleton operons
        singletons = reduced_db_df[reduced_db_df.apply(lambda x: (x.region_len == 1) 
                                                       and (x.region[0] in keep_clusters), 
                                                       axis=1)].copy()
        singletons.update(singletons.region.apply(lambda x: tuple(x)[0]))
        singletons.drop_duplicates('region', inplace=True)
        singletons.update(singletons.region.apply(lambda x: [x]))
        # Recalculate region stats after adding singletons
        out_db_df = pd.concat([out_db_df, singletons], ignore_index=True)
        out_db_df.loc[:,'region_len'] = out_db_df.region.apply(lambda x: len(x))
    
    print("Saving the new synteny database")
    out_db_df.to_pickle(outdb_path)
    
    return out_db_df
