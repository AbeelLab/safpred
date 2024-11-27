# SAFPred: Synteny-aware gene function prediction for bacteria using protein embeddings

SAFPred is a novel synteny-aware gene function prediction tool based on protein embeddings, to annotate bacterial species. SAFPred distinguishes itself from existing methods for annotating bacteria in two ways: (i) it uses embedding vectors extracted from state-of-the-art protein language models and (ii) it incorporates conserved synteny across the entire bacterial kingdom using a novel synteny-based approach proposed in our work.

In this repository the scripts we used to run and evaluate SAFPred.

# Table of Contents

- [Dependencies](#dependencies)
- [Running SAFPred](#safpred)
  - [Usage](#usage)
  - [Miscellaneous scripts](#misc-scripts)
- [Data](#data)
  - [SwissProt database and GO ontology](#sprot-go)
  - [SAFPredDB](#safpreddb)
  - [Train and test sequences](#train-test)
  - [Protein language models](#plm)
- [Citation](#citation)

# Dependencies <a name="dependencies"></a>

- pandas

- biopython

- bio_embeddings

- tqdm

You can create a new conda environment and install all dependencies with the environment file in this repository.

```bash
conda env create -f environment.yml
```

# Data <a name="data"></a>

### SwissProt database and GO ontology <a name="sprot-go"></a>

- The latest release of GO ontology can be downloaded from [the Gene Ontology website](http://geneontology.org/).
  - To reproduce our results, use release 2021-11-16 from the GO archives.
- You can download the SwissProt data used to evaluate SAFPred [from the UniProt website](https://www.uniprot.org/help/downloads). 
  - To reproduce the results, use release 2021-04 from the [previous releases page](https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/).
  - You can run the script `download_sprot_wrapper.py` to download and process the data to apply the filtering procedure as we describe in the manuscript.
    1. Sequence length is kept within [40,1000]
    2. Proteins with at least one experimental GO annotation are retained (evidence codes: EXP, IDA, IPI, IMP, IGI, IEP, HTP, HDA, HMP, HGI, HEP, IBA, IBD, IKR, IRD, IC, TAS).

```bash
./scripts/download_sprot_wrapper.py -h

usage: download_sprot_wrapper.py [-h] --output_dir OUTPUT_DIR --release RELEASE

Download SwissProt database and preprocess it

optional arguments:
  -h, --help            show this help message and exit
  --output_dir OUTPUT_DIR, -o OUTPUT_DIR
                        Directory to write the output files
  --release RELEASE, -r RELEASE
                        Database release date to download -- only the "current" option is
                        supported at the moment
```

   

- Here is an example run:
  
  ```bash
  ./download_sprot_wrapper.py -o data -r current
  ```

    which will create the following directory

```bash
data
├── go
│   └── go.obo
└── uniprot
    └── sprot_db_current.fasta
    └── sprot_db_current_metadata.tsv
    └── uniprot_sprot.dat.gz
```

1. `sprot_db_current_metadata.tsv` is a tab-separated file, mapping uniprot entries to experimental GO annotations
2. `sprot_db_current.fasta` is a fasta file that contains the corresponding proteins sequences for the proteins in `sprot_db_current_metadata.tsv`.

### SAFPredDB <a name="safpreddb"></a>

The synteny database, SAFPredDB, is based on the entire Genome Taxonomy Database (GTDB Release 202, retrieved on 31/03/2022). You can find the current release, as well as the previous versions of GTDB [here](https://gtdb.ecogenomic.org/). Our synteny database will be uploaded to the TU Delft repositories soon. 

- You can download the full SAFPredDB and the corresponding embeddings [from 4TU.ResearchData](https://doi.org/10.4121/ac84802e-853f-46f1-9786-b9d29c0f7557.v1). You can find the accompanying scripts to load and edit an existing database, but also create your own in the SAFPredDB Github repository [here](https://github.com/AbeelLab/SAFPredDB)). On the 4TU.ResearchData you will find 3 files:
  
  - `safpreddb_full_nr.pkl.tar.gz`
    
    The database itself, stored as a pickled pandas DataFrame. The syntenic regions are numbered (zero-based index), each one has the columns:
    
    - `region`: the gene clusters within the syntenic region, you can find which genes were clustered in which cluster IDs in the `safpreddb_cluster_dict.pkl` file
    
    - `intergenic_dist`: intergenic distance between the clusters in a syntenic region. Since there are multiple genes in each cluster, and thus multiple intergenic distance values between any two cluster in a region, we report the minimum of these values
    
    - `region_len`: number of clusters in a region
  
  - `safpreddb_full_emb.pkl`
      A dictionary mapping each region to the average synteny vector, which was calculated by taking the average of embedding vectors (extracted using the ESM-1b model) of each cluster in the region
  
  - `safpreddb_cluster_dict.pkl`
    
    A dictionary mapping each cluster to the genes it caontains, and the genomes these genes are located on as `{cluster_id: {'genes': genes, 'genomes': genomes}}`

- To test our code, we provide a small portion of the full database which you can download along with the corresponding embeddings from [here](https://surfdrive.surf.nl/files/index.php/s/vNizLfgLL4gqUeZ).

### Train and test sequences <a name="train-test"></a>

To test our code, we created a small subset of our full *E. coli* experiment. In this repository, you can find both the fasta sequences and the embeddings:

```bash
data/input
├── test_embeddings_ecoli_small.pkl
├── test_proteins_ecoli_small.txt
├── test_seq_ecoli_small.fasta
├── train_embeddings_ecoli_small.pkl
├── train_proteins_ecoli_small.txt
└── train_seq_ecoli_small.fasta
```

### Protein Language Models <a name="plm"></a>

If you want to extract protein embeddings using the language models we used in our work (instead of providing your own embeddings) you need to download the model weights and files from their own repositories, [Evolutionary Scale Modeling (esm): Pretrained language models for proteins](https://github.com/facebookresearch/esm)) and [T5-XL-U50: Get protein embeddings from protein sequences](https://github.com/sacdallago/bio_embeddings). 

Once you have the model data, make sure you place them under the `models` directory

```bash
data
├── models
│   ├── esm1b
│   └── t5
```

# Running SAFPred <a name="safpred"></a>

To run SAFPred successfully, a few annotation and database files are needed, and you should preserve the directory structure as we provide in this repository:

```bash
data
├── go
│   ├── go_classes.pkl
│   ├── go_ics.pkl
│   └── go.obo
├── safpreddb
│   ├── cluster2go_full.pkl
│   ├── cluster2go.pkl
```

Note that the `cluster2go.pkl` file is a subset of the full mapping in `cluster2go_full.pkl`, it was created to accompany the toy example in this repository. If you want to use the full SAFPredDB database from 4TU.ResearchData, you need the full mapping, `cluster2go_full.pkl`. 

### Usage <a name="usage"></a>

```bash
./safpred.py -h
usage: safpred.py [-h] --train_seq TRAIN_SEQ --test_seq TEST_SEQ --annot_file ANNOT_FILE --output
                  OUTPUT [--db_path DB_PATH] [--db_emb_path DB_EMB_PATH] [--emb_model EMB_MODEL]
                  [--train_emb_file TRAIN_EMB_FILE] [--test_emb_file TEST_EMB_FILE]
                  [--nn_percentile NN_PERCENTILE] [--norm_sim] [--keep_singletons]

run SAFPred

optional arguments:
  -h, --help            show this help message and exit
  --train_seq TRAIN_SEQ, -i TRAIN_SEQ
                        Path to the training set sequences file in fasta format
  --test_seq TEST_SEQ, -t TEST_SEQ
                        Path to the test set sequences file in fasta format
  --annot_file ANNOT_FILE, -a ANNOT_FILE
                        Path to the annotation file, mapping all proteins (train and test) to GO
                        terms
  --output OUTPUT, -o OUTPUT
                        Output directory where the prediction results will be stored
  --db_path DB_PATH, -d DB_PATH
                        Path to the synteny database file
  --db_emb_path DB_EMB_PATH, -b DB_EMB_PATH
                        Path to the synteny database embeddings file
  --emb_model EMB_MODEL, -e EMB_MODEL
                        Embedding model: "esm", "t5" or "none". which requires train and test
                        embedding files as input
  --train_emb_file TRAIN_EMB_FILE, -f TRAIN_EMB_FILE
                        Path to the training embeddings, required only if you want to use your own
  --test_emb_file TEST_EMB_FILE, -g TEST_EMB_FILE
                        Path to the test embeddings, required only if you want to use your own
  --nn_percentile NN_PERCENTILE, -p NN_PERCENTILE
                        Percentile for the NN thresholds. Default is 99.999
  --norm_sim, -n        Normalize NN similarities. Default is True
  --keep_singletons, -k
                        Keep singleton entries in the synteny database. Default is True
```

The following example run will make predictions for 50 *E. coli* proteins, using embedding vectors we already created before, and write the final prediction output into the directory `out`:

```bash
./safpred.py -i data/input/train_seq_ecoli_small.fasta -t data/input/test_seq_ecoli_small.fasta \
    -o out -a data/uniprot/sprot_db_current_metadata.tsv -d data/safpreddb/safpreddb.pkl.gz \
    -b data/safpreddb/safpreddb_emb.pkl -e none -f data/input/train_embeddings_ecoli_small.pkl \
    -g data/input/test_embeddings_ecoli_small.pkl 
```

### Miscellaneous scripts <a name="misc-scripts"></a>

Under the directory `scripts` you will find miscellaneous scripts we used when running the experiments for our manuscript. 

- `download_sprot_wrapper.py`: python wrapper to download the current release of the SwissProt database and GO ontology, and prepare the databases to run SAFPred.

- `evaluate_predictions.py`_: python wrapper to evaluate SAFPred prediction outputs

- `run_blast_baseline.py`: python wrapper to (i) create a blast db from the training set sequences, (ii) predict using the [frequency blast](https://doi.org/10.1186/s13059-019-1835-8) approach and (iii) propagate and normalize predictions.

# Citation <a name="citation"></a>

If you find our method or any of the original script in this repository useful, please cite the manuscript:

```python
@article {urhan2024safpred,
    author = {Aysun Urhan and Bianca-Maria Cosma and Ashlee M. Earl and Abigail L. Manson and Thomas Abeel},
    title = {SAFPred: Synteny-aware gene function prediction for bacteria using protein embeddings},
    year = {2024},
    volume = {40},
    number = {6},
    pages = {btae328},
    month = {05},
    doi = {10.1093/bioinformatics/btae328},
    journal = {Bioinformatics}
}
```

If you use the full SAFPredDB database from the [4TU.ResearchData](https://doi.org/10.4121/ac84802e-853f-46f1-9786-b9d29c0f7557.v1) please cite 

```python
@misc{urhan2024safpreddb,
  doi = {10.4121/AC84802E-853F-46F1-9786-B9D29C0F7557},
  url = {https://data.4tu.nl/datasets/ac84802e-853f-46f1-9786-b9d29c0f7557},
  author = {Urhan, Aysun and Cosma, Bianca-Maria and Earl, Ashlee M. and Manson, Abigail L. and Abeel, Thomas},
  keywords = {Microbiology, FOS: Biological sciences, Genetics, Biological Sciences, bionformatics, microbial genomics, genomics, protein language model, bacterial genomics, comparative genomics, protein embeddings, sequence analysis, bacterial synteny},
  language = {en},
  title = {SAFPredDB: Bacterial synteny database},
  publisher = {4TU.ResearchData},
  year = {2024},
  copyright = {Creative Commons Attribution Non Commercial 4.0 International}
}
```
