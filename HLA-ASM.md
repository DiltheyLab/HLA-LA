# HLA*ASM

## Basics

HLA*ASM can find HLA gene exon coordinates in long read-based assemblies and carry out HLA typing at G group resolution.

It can also carry out a comparison to known sample HLA types, which can be a useful assembly quality criterion.

## Installing HLA*ASM

HLA*ASM is part of [HLA*LA](README.md). 

### Prerequisites

[bwa](https://github.com/lh3/bwa) >= 0.7.12
nucmer (part of [MUMmer3.23](http://mummer.sourceforge.net/))
List::MoreUtils (Perl module)
Text::LevenshteinX (Perl module)
Bio::DB::HTS (Perl module)

#### Installing the Perl modules

Installing `Text::LevenshteinXS` and `List::MoreUtils` should be straighforward and be as simple as `perl -MCPAN -e shell`, followed by `install Text::LevenshteinXS` and `install List::MoreUtils`.

Installing [Bio::DB::HTS](http://search.cpan.org/~rishidev/Bio-DB-HTS/lib/Bio/DB/HTS.pm) can be a bit more tricky. You can try the CPAN approach, but might want to refer to the [GitHub repository](https://github.com/Ensembl/Bio-DB-HTS) for more information.

### Download the data package

Download the data package (https://www.dropbox.com/s/mnkig0fhaym43m0/reference_HLA_ASM.tar.gz, md5sum `4e85ddb1ca1130710815f9c0e5e52e22`) and install it into the main directory of HLA*LA:

~~~~
cd HLA-LA/src
wget https://www.dropbox.com/s/mnkig0fhaym43m0/reference_HLA_ASM.tar.gz
tar -xvzf reference_HLA_ASM.tar.gz
~~~~

### Modifying paths.ini

Like HLA*LA, HLA*ASM uses the `path.ini` mechanism to find its dependencies. Refer to [README.md](README.md) for details. HLA*ASM has a separate working directory.

## Running HLA*ASM

`HLA-ASM.pl --assembly_fasta /path/to/assembly.fasta --sampleID MySample --truth /path/to/HLAreference.txt`

Parameters
* `--assembly_fasta`: Path to a single combined FASTA file of the assembly.
* `--sampleID`: A unique sample ID. Calls with different sample IDs can be run in parallel.
* `--truth`: An HLA truth file (optional). For format see the examples in the `reference_HLA_ASM` directory.

All output goes by default into `output_HLA_ASM/$mySampleID` (where `$mySampleID` is the value if `--sampleID`). Use a unique sample ID for each sample.

## Interpreting the output from HLA*LA

The main output file is `output_HLA_ASM/$mySampleID/summary.txt`.

Columns:
* `contigID`: contig ID.
* `locus`: HLA locus.
* `calledGenotypes`: called HLA genotypes, G group resolution.
* `components`: exons used for determination of the HLA genotype.
* `editDistance_calledGenotypes_assembly`: edit distance between the called HLA genotypes and the assembly.
* `minEditDistance_assembly_truth`: minimum edit distance between the assembly and any of the provided true HLA types for the locus.
* `minEditDistance_calledGenotype_truth`: minimum edit distance between the called HLA genotypes and any of the provided true HLA types for the locus.
* `minEditDistance_assembly_truth_whichAlleles`: alleles (from the truth set) underlying the field `minEditDistance_assembly_truth`.
* `minEditDistance_calledGenotype_truth_whichAlleles`: allele pairs (called HLA genotypes / truth type) underlying the field `minEditDistance_assembly_truth`.				

We also generate additional *.tab files in `output_HLA_ASM/$mySampleID/temp_alignments_separate`, the most interesting of which is `genePositions.tab`. It contains the positions of found genes and exons. These can be used to determine gene presence/absence from the assembly or can serve as the basis for higher-resolution typing.

## Citing HLA*ASM

For the time being, please cite the [original HLA*PRG paper](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005151).






