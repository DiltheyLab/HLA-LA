# HLA*PRG:LA

HLA\*PRG:LA stands for "HLA*PRG, linear approximation". An accompanying [blog post](https://genomeinformatics.github.io/HLA-PRG-LA/) explains what exactly this means, but the basic idea is to seed graph alignments with linear alignments to the sequences that the graph consists of.

HLA\*PRG:LA is faster and less resource-intensive than HLA*PRG.

See [this blog post](https://genomeinformatics.github.io/HLA-PRG-LA/) for a description of the algorithm.

## Installing HLA*PRG:LA

### Prerequisites

g++ with support for C++11 (e.g. 4.7.2)  
Boost >= 1.59  
Bamtools (if you run into issues, try this commit: https://github.com/pezmaster31/bamtools/commit/2d7685d2aeedd11c46ad3bd67886d9ed65c30f3e)  
libz

bwa >= 0.7.12  
samtools >= 1.3  
picard

### Compilation

Create a directory structure for HLA\*PRG-LA:

`mkdir HLA-PRG-LA HLA-PRG-LA/bin HLA-PRG-LA/src HLA-PRG-LA/obj HLA-PRG-LA/temp HLA-PRG-LA/working HLA-PRG-LA/graphs`

Clone this repository into HLA-PRG-LA/src:

`cd HLA-PRG-LA/src; git clone https://github.com/AlexanderDilthey/HLA-PRG-LA.git .`

Compile: modify the paths to libraries and includes in the makefile and then

`make all`

Test that an executable has been created by executing

`../bin/HLA-PRG-LA --action testBinary`

... and you should the following message:

`HLA*PRG:LA binary functional!` 

Any other message indicates that there is a problem - if you receive errors about shared libraries, modify your `LD_LIBRARY_PATH` accordingly.

### Download the data package

Download the data package (http://www.well.ox.ac.uk/PRG_MHC_GRCh38_withIMGT.tar.gz, 2.3G) and extract it into `HLA-PRG/graphs`, i.e.

~~~~
cd HLA-PRG-LA/graphs
wget http://www.well.ox.ac.uk/PRG_MHC_GRCh38_withIMGT.tar.gz
tar -xvzf PRG_MHC_GRCh38_withIMGT.tar.gz
~~~~

md5sum for PRG_MHC_GRCh38_withIMGT.tar.gz is 525a8aa0c7f357bf29fe2c75ef1d477d.

### Modifying paths.ini

HLA\*PRG:LA makes use of bwa, samtools and picard for various steps of the inference process. It is recommended to manually specify the paths to the right executables in the file `HLA-PRG-LA/src/paths.ini` - the cloned repository will contain an example file, and the format should be self-explanatory.

Note that you can specify multiple alternatives per program, for example for running HLA\*PRG:LA in a heterogeneous environment (HLA\*PRG:LA will always use the first alternative present, and also try a `which` if none of the specified alternatives are present).

### Index graph

Finally, pre-compute the graph index structure - this can take a few hours:

`../bin/HLA-PRG-LA --action prepareGraph --PRG_graph_dir ../graphs/PRG_MHC_GRCh38_withIMGT`

### Test run

Download and index the NA12878 test CRAM file from https://gembox.cbcb.umd.edu/shared/NA12878.cram (63G; md5sum 3dc5d3ace0102e8c1f1ffd2270c1359d), run HLA\*PRG:LA, and compare the output with https://github.com/AlexanderDilthey/HLA-PRG-LA/blob/master/NA12878_example_output_G.txt.

All allele calls should agree, and `Q` should be 1 for all calls.

## Running HLA\*PRG:LA

`./inferHLATypes.pl --BAM /path/to/indexed.bam --graph PRG_MHC_GRCh38_withIMGT --sampleID $mySampleID --maxThreads 7`

A few notes:
* All output goes into `../working/$mySampleID` (where `$mySampleID` is a variable). Use a unique sample ID for each sample.
* You can also specify a CRAM file.
* Both CRAM and BAM files need to be indexed.
* Modify `--maxThreads 7` according to your needs.
* HLA\*PRG:LA tries to automatically figure out the right reference genome for your BAM. It compares the index of your BAM/CRAM with a database of known references, that contains the regions relevant for HLA typing. Reads from these regions are extracted and processed. We currently have support for various versions of B37 and for the 1000 Genomes GRCh38 reference. If the program complains that it cannot find a compatible entry in its internal database, please get in touch - adding more references is easy (see below), and we want to support as wide a range of popular references as possible!
* You do *not* need to modify the utilized graph depending on whether your BAM uses GRCh37 or GRCh38. 

### CRAM files

If you use CRAM input, make sure that your CRAM file contains *all* of the original sample reads, including the unmapped ones (which are typically enriched for HLA-derived reads). We've sometimes come CRAM files for which this hasn't quite been the case; and the resulting HLA calls were not very good (the coverage statistics in the call file are sometimes, but now always, indicative of such problems).

## Interpreting the output from HLA*PRG:LA

For each sample with ID $mySampleID, the main output file is ../working/$mySampleID/hla/R1_bestguess_G.txt. Columns:
* `Locus`: Locus (HLA-A, HLA-B, etc..).
* `Chromosome`: 1 or 2. Calls are not phased across loci!
* `Allele`: Called allele at G group resolution (exons 2/3 for class I genes, exon 2 for class II genes).
* `Q1`: Quality score (~ probability), between 0 and 1.
* `Q2`: Ignore.
* `AverageCoverage`: Locus average coverage.
* `CoverageFirstDecile`: First decile coverage.
* `MinimumCoverage`: Minimum column coverage for locus.
* `proportionkMersCovered`: Proportion k-mers belonging to the called allele observed in locus input data.
* `LocusAvgColumnError`: Average alignment error for this locus.
* `NColumns_UnaccountedAllele_fGT0.2`: Number of columns with evidence for the presence of alleles not accounted for by the called alleles.
* `perfectG`: Perfect translation from HLA*PRG:LA call to G grop resolution? (1 = Yes).

Checking `AverageCoverage` for each sample is an important sanity check - `AverageCoverage` fluctuates between samples and loci, but its expectation depends on average whole-genome coverage. Also, plotting `AverageCoverage` for all samples and all loci in your cohort is a good idea!

There are three additional quality indicators: `Q1`, `proportionkMersCovered`, `NColumns_UnaccountedAllele_fGT0.2`. `Q1` is usually equal to or very close to 1. `proportionkMersCovered` should always be 1, at least for high-coverage Illumina samples. `NColumns_UnaccountedAllele_fGT0.2` is also usually 0, but `NColumns_UnaccountedAllele_fGT0.2` != 0 is not a reliable indicator for false calls. Deviations from the "usual" values sometimes indicate the presence of novel alleles. See the [original HLA\*PRG paper](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005151), in which we discussed the predictive value of these indicators in some detail.

If `perfectG` != 0, you might want to check `../working/$mySampleID/hla/R1_bestguess.txt`, which contains the un-translated output from the internal model of the algorithm.

Also, `../working/$mySampleID/reads_per_level.txt` gives coverage across the MHC - the coordinate system is that of the PRG itself, i.e. it is not identical to normal genomic coordinates. The file also contains string identifiers for each level, indicating e.g. gene names. We haven't tested the output of this file in any way, but it might be interesting for some applications.

## Adding further references

Whenever you apply HLA\*PRG:LA to a BAM file, the algorithm compares the reference genome of the BAM (contig identifiers and contig lengths, using `samtools idxstats`) to a database of known references. This database tells the algorithm which BAM regions are relevant for HLA typing; it will extract reads from these regions and use them for typing.

If, however, your utilized reference genome is not in the internal database (to re-iterate this point: we're looking for exact matches in terms of contig identifiers and contig lengths!), you will receive an error message.

The easiest thing to do is then to send the output from samtools idxstats to us, and we'll add your reference to the database.

You can also, however, modify the database yourself:
* Each graph has its own reference database. Typically your reference directory is `HLA-PRG-LA/graphs/PRG_MHC_GRCh38_withIMGT/knownReferences`.
* Each file in this directory represents one reference. The files are tab-separated with header lines and each line shall have the full number of fields (even if they are empty). The required fields are:
  * contigID: Name of the contig
  * contigLength: Length of the contig
  * ExtractCompleteContig: Requires a value, either 0 (no) or 1 (yes)
  * PartialExtraction_Start: Mututally incompatible with ExtractCompleteContig = 1 - if partial extraction, start of the extraction (1-based coordinates).
  * PartialExtraction_Stop: Stop coordinate for partial extraction (see PartialExtraction_Start)
* For the MHC from chromosome 6, we typically use partial extraction with the MHC coordinates.
* MHC ALTs and HLA genomic sequences (if any) are usually completely extracted.
* Unmapped reads should usually be extracted. Use '\*' (without quotes) as contigID, 0 for contigLength, and 1 for ExtractCompleteContig.

If you create your own extraction files, have a look at the existing files first - they will give you a good example to work from.

Important: if you miss regions, the reads aligned to these are lost for the inference process - make sure to catch all alternative MHC contigs and all HLA references!

## Citing HLA*PRG:LA

Please cite the [original HLA*PRG paper](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005151).
