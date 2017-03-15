# HLA*PRG:LA

HLA\*PRG:LA stands for "HLA*PRG, linear approximation". An accompanying blog post will explain what exactly this means, but the basic idea is to seed graph alignments with linear alignments to the sequences that the graph consists of.

HLA\*PRG:LA should be faster and less resource-intensive than HLA*PRG.

## Installing HLA*PRG:LA

### Prerequisites

g++ with support for C++1 (e.g. 4.7.2)  
Boost >= 1.59  
Bamtools (if you run into issues, try this commit: https://github.com/pezmaster31/bamtools/commit/2d7685d2aeedd11c46ad3bd67886d9ed65c30f3e)  
libz

### Compilation

Create a directory structure for HLA*PRG-LA:

`mkdir HLA-PRG-LA HLA-PRG-LA/bin HLA-PRG-LA/src HLA-PRG-LA/obj HLA-PRG-LA/temp HLA-PRG-LA/working HLA-PRG-LA/graphs`

Clone this repository into HLA-PRG-LA/src:

`cd HLA-PRG-LA; git clone https://github.com/AlexanderDilthey/HLA-PRG-LA.git .`

Compile: modify the makefile and then

`make all`

Test that an executable has been created by executing

`../bin/HLA-PRG-LA`

... you should receive an error message like

terminate called after throwing an instance of 'std::runtime_error'  
what():  Please specify arguments --bwa_bin

If you receive errors about shared libraries, modify your LD_LIBRARY_PATH accordingly.

### Download data package

### Index graph and first run


## Running HLA*PRG:LA

`./inferHLATypes.pl --BAM /path/to/indexed.bam --graph PRG_MHC_GRCh38_withIMGT --sampleID $mySampleID --maxThreads 7`

A few notes:
* All output goes into ../working/$mySampleID (where $mySampleID is a variable). Use a unique sample ID for each sample.
* HLA*PRG:LA compares the reference genome underlying your BAM with a database of known references. This database tells HLA*PRG:LA which regions in which reference genome are relevant for HLA typing, and reads from these are extracted and processed. We currently have support for various versions of B37 and for the 1000 Genomes reference file for GRCh38. If the program complains that it cannot find a compatible entry in its internal database, please get in touch - adding more references is easy (see below), and we want to support as wide a range of popular references as possible!

## Interpreting the output from HLA*PRG:LA

For each sample with ID $mySampleID, the main output file is ../working/$mySampleID/hla/R1_bestguess_G.txt. Columns:
* Locus: Should be clear
* Chromosome: 1 or 2. Calls are not phased across loci!
* Allele: Called allele at G group resolution (exons 2/3 for class I genes, exon 2 for class II genes)
* Q1: Quality score (~ probability), between 0 and 1.
* Q2: Ignore
* AverageCoverage: Locus average coverage
* CoverageFirstDecile: First decile coverage
* MinimumCoverage: Minimum column coverage for locus
* proportionkMersCovered: Proportion k-mers belonging to the called allele observed in locus input data
* LocusAvgColumnError: Average alignment error for this locus
* NColumns_UnaccountedAllele_fGT0.2: Number of columns with evidence for the presence of alleles not accounted for by the called alleles
* perfectG: Perfect translation from HLA*PRG:LA call to G grop resolution? (1 = Yes)

Checking AverageCoverage for each sample is an important sanity check - AverageCoverage fluctuates between samples and loci, but its expectation depends on average whole-genome coverage. Also, plotting AverageCoverage for all samples and all loci in your cohort is a good idea!

There are three additional quality indicators: Q1, proportionkMersCovered, NColumns_UnaccountedAllele_fGT0.2. Q1 is usually equal to or very close to 1. proportionkMersCovered should always be 1, at least for high-coverage Illumina samples. NColumns_UnaccountedAllele_fGT0.2 is also usually 0, but NColumns_UnaccountedAllele_fGT0.2 != 0 is not a reliable indicator for false calls. Deviations from the "usual" values sometimes indicate the presence of novel alleles. See the [original HLA*PRG paper](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005151), in which we discussed the predictive value of these indicators in some detail.

If perfectG != 0, you might want to check ../working/$mySampleID/hla/R1_bestguess.txt, which contains the un-translated output from the internal model of the algorithm.

Also, ../working/$mySampleID/reads_per_level.txt gives coverage across the MHC - the coordinate system is that of the PRG itself, i.e. it is not identical to normal genomic coordinates. The file also contains string identifiers for each level, indicating e.g. gene names. We haven't tested the output of this file in any way, but it might be interesting for some applications.

### Citing HLA*PRG:LA

Please cite the [original HLA*PRG paper](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005151).
