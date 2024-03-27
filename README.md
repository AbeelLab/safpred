# SAFPred: Synteny-aware gene function prediction for bacteria using protein embeddings
SAP is a novel synteny-aware gene function prediction tool based on protein embeddings, to annotate bacterial species. SAFPred distinguishes itself from existing methods for annotating bacteria in two ways: (i) it uses embedding vectors extracted from state-of-the-art protein language models and (ii) it incorporates conserved synteny across the entire bacterial kingdom using a novel synteny-based approach proposed in our work.

In this repository the scripts we used to run and evaluate SAFPred.

### Dependencies

### Data
- The latest release of GO ontology can be downloaded from [the Gene Ontology website](http://geneontology.org/).
  - To reproduce our results, use release 2021-11-16 from the GO archives.
- You can download the SWISSPROT data used to evaluate SAFPred [from the UniProt website](https://www.uniprot.org/help/downloads). 
  - To reproduce the results, use release 2021-04 from the [previous releases page](https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/).
  - You can run the script `GetSPROT.py` to download and process the data to apply the filtering procedure as we describe in the manuscript.
    1. Sequence length is kept within [40,1000]
    2. Proteins with at least one experimental GO annotation are retained (evidence codes: EXP, IDA, IPI, IMP, IGI, IEP, HTP, HDA, HMP, HGI, HEP, IBA, IBD, IKR, IRD, IC, TAS).
  - `GetSPROT.py` will create 2 text files:
    1. `uniprot2go_exp.tsv` is a tab-separated file, mapping uniprot entries to experimental GO annotations
    2. `uniprot2fasta_exp.fasta` is a fasta file that contains the corresponding proteins sequences for the proteins in `uniprot2go_exp.tsv`.
- The operon database is based on the entire Genome Taxonomy Database (GTDB Release 202, retrieved on 31/03/2022). You can find the current release, as well as the previous versions of GTDB [here](https://gtdb.ecogenomic.org/). Our operon database will be uploaded to TU Delft repositories soon.
