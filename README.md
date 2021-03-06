# TB-tools

## Purpose of the package

- TB-annotator: Collect and annotate *Mycobacterium tuberculosis* WGS data.    
- CRISPRbuilder-TB: de novo reconstruction of the CRISPR locus of *tuberculosis*
(coming soon).

The main objective of these tools is to provide various genotyping tools for the
*Mycobacterium tuberculosis* complex. These tools work only for Whole Genome
Sequencing (WGS) data in fasta files.
The choice to consider only WGS data comes from the fact that the assembled 
genomes of *M. tuberculosis* often have many errors, either from areas difficult 
to reconstruct (e.g., CRISPR locus), multiple repeated sequences, or the use 
of an unrepresentative strain as a reference (h37Rv).
Both programs in this package receive either a Sequence Read Archive (SRA) 
accession number, or a list of such SRAs in a text file (one number per line), 
or the ID of an NCBI BioProject.

**TB-annotator** is used to collect information from a given WGS *tuberculosis*
genome. SRA files are downloaded if needed, and name of the strain, location, 
bioproject, biosample, etc. are first collected from the NCBI server, when 
available:  Read files are then systematically studied, to provide the following 
information:
 - number of reads, their length, and the mean coverage;
 - 43-spacers based and 98-spacers based *in silico*, and *in vitro* inferred, 
 spoligotypes, plus the SIT number;
 - lineages based on various state-of-the-art SNP collections: PGG, Coll, 
 Palittapongarnpim (L1), Shitikov (L2), and Stucki (L4);
 - the number and position (in h37Rv) of numerous Insertion Sequences, including
 IS6110.

**CRISPRbuilder-TB** reconstructs the whole CRISPR locus starting from WGS data.
This allows to deduce the true spoligotype, to detect mutants of spacers and
direct repeats, insertion of mobile elements, and duplications. The Cas locus 
is reconstructed too. This is a semi automatic approach that leads to a set of 
contigs to assemble manually. Depending on the number, length, and quality of
SRAs, the number of contigs can range from 1-2 patterns, in the best case scenario
where the good quality of sequences allows an automatic reconstruction of the 
CRISPR cut in mobile element positions, to several contigs difficult to process,
for too short or polluted reads.


## Requirements

TB-tools needs Python 3, and the following dependencies to work:

* *fastq-dump* from sra-tools (NCBI) : https://github.com/ncbi/sra-tools
* *blastn* and *makeblastdb* from blast+ (NCBI) : https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download

A 'bin' directory as been added to the project with the 3 required executables, 
and for the GNU/Linux, Mac and Windows platforms, but the preferred method is
a clean installation from the source.


## Installation

You can either download a zip file by clicking on the 'Clone or download' green
button, or clone the repository by writing in a terminal:
<pre>
git clone https://github.com/cguyeux/TB-tools.git
</pre>

To install the required python libraries: 
<pre>
pip3 install -r requirements.txt
</pre>


## How to use TB-tools

To launch a first analysis: 
<pre>
python tbannotator.py -sra ERR037527
</pre>
Help about this script: 
<pre>
python tbannotator.py -h
</pre>
By default, only basic information (e.g., Coll et al. lineage) is printed in
the terminal (to save execution time), and results are stored in a pickled 
python dictionary in the 'sequence' directory. To read such a data using 
Python:
<pre>
from pickle import load
with open('sequences/ERR037527/ERR037527.pkl', 'rb') as f:
    dico = load(f.read())
</pre>

## Citation

>Christophe Guyeux, Christophe Sola, Guislaine Refrégier. CRISPRbuilder-TB: “CRISPR-Builder for tuberculosis”. Exhaustive reconstruction of the CRISPR locus in Mycobacterium tuberculosis complex using SRA. Submitted article (2020).