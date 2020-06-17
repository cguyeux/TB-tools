# TB-tools
------------------

This **README.md** gives you the gist of the TB-tools package.

## Main content
--------------------------

- TB-annotator: Collect and annotate Mycobacterium tuberculosis WGS data.    
- CRISPRbuilder-TB: de novo reconstruction of the CRISPR locus of tuberculosis.


## Purpose of the package
--------------------------

The main objective of these tools is to provide various genotyping tools for the
*Mycobacterium tuberculosis* complex. These tools work only for Whole Genome
Sequencing (WGS) data, more specifically for Illumina paired end sequences.
The choice to consider only WGS data comes from the fact that the assembled 
genomes of *M. tuberculosis* often have many errors, either from areas difficult 
to reconstruct (e.g., CRISPR locus), multiple repeated sequences, or the use 
of an unrepresentative strain as a reference (h37Rv).
Both programs in this package receive either a Sequence Read Archive (SRA) 
accession number, or a list of such SRAs in a text file (one number per line), 
or the ID of an NCBI BioProject.

**tb-annotator** is used to collect information from a given WGS *tuberculosis*
genome. SRA files are downloaded if needed, and name of the strain, location, 
bioproject, biosample, etc. are first collected from the NCBI server, when 
available:  Read files are then systematically studied, to provide the following 
information:
 - number of reads, their length, and the mean coverage;
 - 43-spacers based and 98-spacers based *in silico*, and *in vitro* infered, 
 spoligotypes, plus the SIT number;
 - lineages based on various state-of-the-art SNP collections: PGG, Coll, 
 Palittapongarnpim (L1), Shitikov (L2), and Stucki (L4);
 - the number and position (in h37Rv) of numerous Insertion Sequences, including
 IS6110.


## Requirements
---------------

TB-tools needs the following dependencies to work:

* python >= "3.4"
* fastq-dump from sra-tools (NCBI) : https://github.com/ncbi/sra-tools
* blast+ (NCBI) : https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download

To install the required python libraries: `pip3 install -r requirements.txt`

## How to use TB-tools
----------------------

To launch a first analysis: `python tbannotator.py -sra ERR037527`

Help about this script: `python tbannotator.py -h`

By default, only basic information (e.g., Coll et al. lineage) is printed in
the terminal (to save execution time), and results are stored in a pickled 
python dictionary in the 'sequence' directory. To read such a data:

`from pickle import load`
`with open('sequences/ERR037527/ERR037527.pkl', 'rb') as f:`
`    dico = load(f.read())`


## Citation

>Christophe Guyeux, Christophe Sola, Guislaine Refrégier. CRISPRbuilder-TB: “CRISPR-Builder for tuberculosis”. Exhaustive reconstruction of the CRISPR locus in Mycobacterium tuberculosis complex using SRA. Submitted article (2020).