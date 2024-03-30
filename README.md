# SAFPred: Synteny-aware gene function prediction for bacteria using protein embeddings

SAFPred is a novel synteny-aware gene function prediction tool based on protein embeddings, to annotate bacterial species. SAFPred distinguishes itself from existing methods for annotating bacteria in two ways: (i) it uses embedding vectors extracted from state-of-the-art protein language models and (ii) it incorporates conserved synteny across the entire bacterial kingdom using a novel synteny-based approach proposed in our work.

In this repository the scripts we used to run and evaluate SAFPred.

### Dependencies

- numpy

- pandas

- biopython

- bio_embeddings_

### Data

- The latest release of GO ontology can be downloaded from [the Gene Ontology website](http://geneontology.org/).
  - To reproduce our results, use release 2021-11-16 from the GO archives.
- You can download the SwissProt data used to evaluate SAFPred [from the UniProt website](https://www.uniprot.org/help/downloads). 
  - To reproduce the results, use release 2021-04 from the [previous releases page](https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/).
  - You can run the script `download_sprot_wrapper.py` to download and process the data to apply the filtering procedure as we describe in the manuscript.
    1. Sequence length is kept within [40,1000]
    2. Proteins with at least one experimental GO annotation are retained (evidence codes: EXP, IDA, IPI, IMP, IGI, IEP, HTP, HDA, HMP, HGI, HEP, IBA, IBD, IKR, IRD, IC, TAS).

```bash
./download_sprot_wrapper.py -h

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
- The synteny database, SAFPredDB, is based on the entire Genome Taxonomy Database (GTDB Release 202, retrieved on 31/03/2022). You can find the current release, as well as the previous versions of GTDB [here](https://gtdb.ecogenomic.org/). Our synteny database will be uploaded to TU Delft repositories soon. To test our code, we provide a small portion of the full database, you can find the accompanying scripts load and edit an existing database, but also create your own in the SAFPredDB Github repository [here](https://github.com/AbeelLab/SAFPredDB)).
