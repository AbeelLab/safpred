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

- You can download the full SAFPredDB [here](https://surfdrive.surf.nl/files/index.php/s/KrWoxvEYdP38Xid) and the corresponding embeddings [here]. You can find the accompanying scripts to load and edit an existing database, but also create your own in the SAFPredDB Github repository [here](https://github.com/AbeelLab/SAFPredDB)).

- To test our code, we provide a small portion of the full database which you can download along with the corresponding embeddings from [here](https://surfdrive.surf.nl/files/index.php/s/vNizLfgLL4gqUeZ) 

### Train and test sequences <a name="train-test"></a>

To test our code, we created a small subset of our full *E. coli* experiment. In this repository, you can find both the fasta sequences and the embeddings:

```bash
data/input
├── test_embeddings_ecoli_small.pkl
├── test_seq_ecoli_small.fasta
├── train_embeddings_ecoli_small.pkl
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
│   ├── cluster2go.pkl
```

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

If you find our method or any of the original script in this repository useful, please cite

```python
@article {urhan2023safpred,
    author = {Aysun Urhan and Bianca-Maria Cosma and Ashlee M. Earl and Abigail L. Manson and Thomas Abeel},
    title = {SAP: Synteny-aware gene function prediction for bacteria using protein embeddings},
    elocation-id = {2023.05.02.539034},
    year = {2023},
    doi = {10.1101/2023.05.02.539034},
    publisher = {Cold Spring Harbor Laboratory},
    URL = {https://www.biorxiv.org/content/early/2023/11/21/2023.05.02.539034},
    eprint = {https://www.biorxiv.org/content/early/2023/11/21/2023.05.02.539034.full.pdf},
    journal = {bioRxiv}
}
```
