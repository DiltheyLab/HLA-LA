#!/usr/bin/env perl

# ./inferHLATypes.pl --BAM /gpfs1/well/gsk_hla/bam_output/AA02O9Q_Z2.bam --graph PRG_MHC_GRCh38_withIMGT --sampleID NA12878Direct

use warnings;
use strict;
use FindBin;
use File::Spec;
use Getopt::Long;
use Data::Dumper; 
use Sys::Hostname;
use Cwd qw/getcwd abs_path/;
use List::MoreUtils qw/mesh/;
use List::Util qw/all shuffle sum/;


$| = 1;
my $this_bin_dir = $FindBin::RealBin;

my $_BAM;
my $graph;
my $sampleID;
my $sampleID_tumor;
my $qsub;
my $samtools_bin;
my $bwa_bin;
my $java_bin;
my $picard_bin;
my $GATK_bin;
# my $WhatsHap_bin;
# my $bcftools_bin;
# my $spades_bin = '/software/SPAdes/3.13.0/skylake/bin/spades.py';
my $workingDir_param;
my $maxThreads = 1;
my $moreReferencesDir;
my $extractExonkMerCounts;
my $longReads = 0;
my $testing = 0;
my $samtools_T;
my $action = 'call';
my $twoStageReadExtraction = 0;

GetOptions (
	'BAM:s' => \$_BAM,
	'graph:s' => \$graph,
	'sampleID:s' => \$sampleID,
	'sampleID_tumor:s' => \$sampleID_tumor,
	'action:s' => \$action,
	'qsub:s' => \$qsub,
	'workingDir:s' => \$workingDir_param,
	'samtools_bin:s' => \$samtools_bin,
	'bwa_bin:s' => \$bwa_bin,
	'java_bin:s' => \$java_bin,
	'picard_bin:s' => \$picard_bin,
	# 'GATK_bin:s' => \$GATK_bin,
	# 'WhatsHap_bin:s' => \$WhatsHap_bin,
	# 'BCFTools_bin:s' => \$bcftools_bin,
	'maxThreads:s' => \$maxThreads,
	'moreReferencesDir:s' => \$moreReferencesDir,
	'extractExonkMerCounts:s' => \$extractExonkMerCounts,
	'longReads:s' => \$longReads,
	'testing:s' => \$testing,
	'samtools_T:s' => \$samtools_T,
	'twoStageReadExtraction:s' => \$twoStageReadExtraction,
);

if(1 == 0)
{
	my $referenceGenome_for_GATK = '/gpfs/project/dilthey/projects/HLA-LA-devel/working/NA12878_mini/remap/remapped_and_projected.bam.projected.SAM.ref.fa';
	my $refGenome_href = readFASTA($referenceGenome_for_GATK);
	my $phaseSets_href = getPhaseSetsPerGene('/gpfs/project/dilthey/projects/HLA-LA-devel/working/NA12878_mini/remap/HLA.phased.VCF.bck');
	my @genes = sort map {$_ =~ s/\*ref//g; $_} grep {$_ =~ /\*ref/} keys %$refGenome_href;
	print "Detected genes: " , join(", ", @genes), "\n";
	print "Phase set statistics\n";
	foreach my $gene (@genes)
	{
		die unless(exists $refGenome_href->{$gene . '*ref'});
		my $ref_length = length($refGenome_href->{$gene . '*ref'});
		print "\t", $gene, " [length: $ref_length]\n";
		if(exists $phaseSets_href->{$gene})
		{
			my @phaseSets = sort keys %{$phaseSets_href->{$gene}};
			foreach my $PS (@phaseSets)
			{
				
				print "\t\tPhase set ", $PS, ", length ", $phaseSets_href->{$gene}{$PS}{coordinates}[1] - $phaseSets_href->{$gene}{$PS}{coordinates}[0] + 1, " [", $phaseSets_href->{$gene}{$PS}{coordinates}[0], " - ", $phaseSets_href->{$gene}{$PS}{coordinates}[1], "]\n";
			}
		}
	}
	exit;
}

if ($extractExonkMerCounts)
{
	warn "--extractExonkMerCounts 1: This feature is experimental and unlikely to work on your machine";
	die unless(-e '../exonkMerExtraction/GRCh38.forkMers');
	die unless(-e '../exonkMerExtraction/exonCoordinates_manual.txt.forExtraction');	
}

unless((not $longReads) or ($longReads eq 'ont2d') or ($longReads eq 'pacbio'))
{
	die "Please specify --longReads ont2d or --longReads pacbio";
}

my %paths_ini;
my $paths_ini = $this_bin_dir . '/paths.ini';
if(-e $paths_ini)
{
	open(INI, '<', $paths_ini) or die "Cannot open $paths_ini";
	while(<INI>)
	{
		chomp;
		next unless($_);
		$_ =~ s/[\n\r]//g;
		next if($_ =~ /^\s+$/);
		die "Invalid line format in $paths_ini -- expect all lines to be either empty or key=value pairs" unless($_ =~ /^(.+)=(.*)$/);
		my $id = $1;
		my @alts = split(/,/, $2);
		$paths_ini{$id} = \@alts;
	}
	close(INI);
}
$samtools_bin = find_path('samtools_bin', $samtools_bin, 'samtools');
$bwa_bin = find_path('bwa_bin', $bwa_bin, 'bwa');
$java_bin = find_path('java_bin', $java_bin, 'java');
$picard_bin = find_path('picard_bin', $picard_bin, 'picard');
#$GATK_bin = find_path('GATK', $GATK_bin, 'gatk');
#$WhatsHap_bin = find_path('WhatsHap', $WhatsHap_bin, 'whatshap');
#$bcftools_bin = find_path('bcftools_bin', $bcftools_bin, 'bcftools');

# debug commands
# trimContigs('/gpfs/project/dilthey/projects/HLA-LA-devel/working/NA12878_mini/remap/assembled_A_complete.untrimmed.fasta', '/gpfs/project/dilthey/projects/HLA-LA-devel/working/NA12878_mini/remap/assembled_A_complete.untrimmed.fasta.reads.bam.pileup.txt', 'test.fa');
# alignContigs_geneAware('/gpfs/project/dilthey/projects/HLA-LA-devel/working/NA12878_mini/remap/remapped_and_projected.bam.projected.SAM.ref.fa', '/gpfs/project/dilthey/projects/HLA-LA-devel/working/NA12878_mini/remap/assembled_all_contigs.fa', 'testOut.bam', $bwa_bin, $samtools_bin);
# exit;

my $FASTQ_extraction_command_part1;
if($picard_bin =~ /SamToFastq\.jar$/)
{
	die "Please use a recent Picard version and specify the path to picard.jar";
}
elsif($picard_bin =~ /picard-tools$/)
{
	$FASTQ_extraction_command_part1 = qq($picard_bin SamToFastq);	
#	die "Please use a recent Picard version and specify the path to picard.jar";
}
elsif($picard_bin =~ /picard\.jar$/)
{
	$FASTQ_extraction_command_part1 = qq($java_bin -Xmx10g -XX:-UseGCOverheadLimit -jar $picard_bin SamToFastq);
}
elsif($picard_bin =~ /picard$/)
{
	die "Please use a recent Picard version and specify the path to picard.jar";
}
else
{
	die "I can't interpret the specified Picard command: $picard_bin";
}

if($testing)
{
	my $previous_dir = getcwd;
	chdir($this_bin_dir) or die "Cannot cd into $this_bin_dir";
	die "Binary ../bin/HLA-LA (from $this_bin_dir) not there!" unless(-e '../bin/HLA-LA');
	my $cmd_test = qq(../bin/HLA-LA	--action testBinary);
	system($cmd_test) and die "HLA*LA test command $cmd_test failed";
	chdir($previous_dir) or die "Cannot chdir into $previous_dir";
	
	print "HLA-LA.pl test\n\n";
	print "\t", "samtools_bin", ": ", $samtools_bin, "\n";
	print "\t", "bwa_bin", ": ", $bwa_bin, "\n";
	print "\t", "java_bin", ": ", $java_bin, "\n";
	print "\t", "picard_bin", ": ", $picard_bin, "\n";
	# print "\t", "GATK", ": ", $GATK_bin, "\n";
	# print "\t", "WhatsHap", ": ", $WhatsHap_bin, "\n";
	# print "\t", "bcftools_bin", ": ", $bcftools_bin, "\n";
	exit 0;
}


my $working_dir;
if($paths_ini{workingDir}[0] and not defined $workingDir_param)
{
	$working_dir = $paths_ini{workingDir}[0];
	$working_dir =~ s/\$HLA\-LA\-DIR/$this_bin_dir/;
}
else
{
	unless(defined $workingDir_param)
	{
		die "\n\nPlease specify a working directory via --workingDir.\n\nOutput for sample with ID \$sampleID will go a correspondingly named sub-directory of the working directory.\n\nFor example, if --workingDir is /path/working, and --sampleID is mySample, then the output will go into directory /path/working/mySample.\n\n";
	}
	$working_dir = $workingDir_param;
}

unless(-d $working_dir)
{
	die "\n\nSpecified working directory $working_dir is either not present or not a directory.\n\nYou might have specified an invalid path via --workingDir.\n\n";
}

$working_dir = abs_path($working_dir);

unless($sampleID =~ /^\w+$/)
{
	die "Please use only alphanumeric characters - \\w+ - for --sampleID";
}
if($action eq 'somatic')
{
	unless($sampleID_tumor =~ /^\w+$/)
	{
		die "Please use only alphanumeric characters - \\w+ - for --sampleID_tumor";
	}	
	if($sampleID eq $sampleID_tumor)
	{
		die "You provided the same value for --sampleID and --sampleID_tumor - this is not supported (and probably doesn't make sense)";
	}
}

my $working_dir_thisSample = $working_dir . '/' . $sampleID;
my $working_dir_tumorSample;
if($action eq 'somatic')
{
	$working_dir_tumorSample = $working_dir . '/' . $sampleID_tumor;
}
print "HLA-LA.pl\n\n";

print "Identified paths:\n";
print "\t", "samtools_bin", ": ", $samtools_bin, "\n";
print "\t", "bwa_bin", ": ", $bwa_bin, "\n";
print "\t", "java_bin", ": ", $java_bin, "\n";
print "\t", "picard_bin", ": ", $picard_bin, "\n";
# print "\t", "GATK_bin", ": ", $GATK_bin, "\n";
if($action eq 'somatic')
{
	print "\t", "Sample-specific working directory (NORMAL)", ": ", $working_dir_thisSample, "\n";
	print "\t", "Sample-specific working directory (TUMOR)", ": ", $working_dir_tumorSample, "\n";
}
else
{
	print "\t", "General working directory", ": ", $working_dir, "\n";
}

print "\n";


unless(-d $working_dir_thisSample)
{
	mkdir($working_dir_thisSample) or die "Cannot mkdir $working_dir_thisSample";
}
 
my $samtools_version = `$samtools_bin --version` ;
die "Can't parse samtools version output" unless($samtools_version =~ /samtools ([\d\.]+)/);
$samtools_version = $1;
my $samtools_version_numeric = $samtools_version;
if($samtools_version_numeric =~ /^(\d+)\.(\d+)\.(\d+)$/)
{
	$samtools_version_numeric = $1 . '.' . $2 . $3;
}
unless($samtools_version_numeric >= 1.3)
{
	if($samtools_version_numeric =~ /^(\d+)\.(\d+)$/)
	{
		my $samtools_version_secondField = $2;
		unless($samtools_version_secondField >= 3)
		{
			die "I need samtools >=1.3";
		}
	}
	else
	{
		die "I need samtools >=1.3";
	}
}

my $BAM;
my $BAM_tumor;
my $call2_HLAtypes;
my %call2_hla_relevant_readIDs_primaryBAM;
my %somatic_hla_relevant_readIDs_tumorBAM;
my $call2_processing_dir;
my $call2_fn_mapping;
my $call2_fn_mapping_whichHaplotype;
my $somatic_processing_dir_tumor;

if($action eq 'call')
{
	$BAM = File::Spec->abs2rel($_BAM);
	unless(-e $BAM)
	{
		die "BAM (or CRAM) $BAM (inferred from input $_BAM) not existing" unless(-e $BAM);
		unless((-e $BAM . '.bai') or (-e $BAM . '.crai'))
		{
			die "File $BAM does not appear to be indexed";
		}
	}
}
elsif(($action eq 'call2') or ($action eq 'somatic')) 
{
	$call2_HLAtypes = $working_dir_thisSample . '/hla/R1_bestguess.txt';
	$call2_HLAtypes = File::Spec->abs2rel($call2_HLAtypes);	
	unless(-e $call2_HLAtypes)
	{
		die "Expected file $call2_HLAtypes not found - it looks like you have not yet carried out the first step (--action call)";
	}
	
	$BAM = $working_dir_thisSample . '/remapped_with_a.bam';
	$BAM = File::Spec->abs2rel($BAM);
	unless(-e $BAM)
	{
		die "Weird - expected BAM $BAM not found";
	}
	
	$call2_processing_dir = $working_dir_thisSample . '/remap';
	unless(-d $call2_processing_dir)
	{
		mkdir($call2_processing_dir) or die "Cannot mkdir $call2_processing_dir";
	}
	$call2_fn_mapping = $call2_processing_dir . '/ref_for_remap.fa';
	$call2_fn_mapping_whichHaplotype = $call2_processing_dir . '/ref_for_remap.fa.seq2Hap';
	
	if($action eq 'somatic')
	{
		my $call2_tumor_HLAtypes = $working_dir_tumorSample . '/hla/R1_bestguess.txt';
		$call2_tumor_HLAtypes = File::Spec->abs2rel($call2_HLAtypes);	
		unless(-e $call2_tumor_HLAtypes)
		{
			die "Expected file $call2_tumor_HLAtypes not found - it looks like you have not yet carried out the first step (--action call) on the tumor sample";
		}
		
		$BAM_tumor = $working_dir_tumorSample . '/remapped_with_a.bam';
		$BAM_tumor = File::Spec->abs2rel($BAM);
		unless(-e $BAM_tumor)
		{
			die "Weird - expected (tumor) BAM $BAM_tumor not found";
		}
		
		$somatic_processing_dir_tumor = $working_dir_thisSample . '/remap_tumor';
		unless(-d $somatic_processing_dir_tumor)
		{
			mkdir($somatic_processing_dir_tumor) or die "Cannot mkdir $somatic_processing_dir_tumor";
		}
	}
}
else
{
	die "Unknown --action $action - try call or call2";
}

my $full_graph_dir = $FindBin::RealBin . '/../graphs/' . $graph;
my $known_references_dir = $full_graph_dir . '/knownReferences';
unless(-e $full_graph_dir)
{
	die "Graph directory $full_graph_dir not found - valid graph names are subdirectories of the graphs directory in the HLA-LA root";
}
unless((-e $full_graph_dir . '/sequences.txt') and ((-e $full_graph_dir . '/extendedReferenceGenomePath.txt') or (-e $full_graph_dir . '/extendedReferenceGenome/extendedReferenceGenome.fa')) and (-d $known_references_dir))
{
	die "Graph directory $full_graph_dir does not seem to be complete - does this directory specify a valid graph for HLA-LA?";
}

unless((-e $full_graph_dir . '/serializedGRAPH'))
{
	die "It seems that you have not yet indexed graph $graph - suggested command: ../bin/HLA-LA --action prepareGraph --PRG_graph_dir ../graphs/${graph}";
}

my $target_FASTQ_1 = $working_dir_thisSample . '/R_1.fastq';
my $target_FASTQ_2 = $working_dir_thisSample . '/R_2.fastq';
my $target_FASTQ_U = $working_dir_thisSample . '/R_U.fastq';

my $target_FASTQ_1_tumor;
my $target_FASTQ_2_tumor;
my $target_FASTQ_U_tumor;
if($action eq 'somatic')
{
	$target_FASTQ_1_tumor = $working_dir_thisSample . '/R_1_tumor.fastq';
	$target_FASTQ_2_tumor = $working_dir_thisSample . '/R_2_tumor.fastq';
	$target_FASTQ_U_tumor = $working_dir_thisSample . '/R_U_tumor.fastq';
}

if(($action eq 'call2') or ($action eq 'somatic'))
{
	my $hla_working_dir = $working_dir_thisSample . '/hla';
	$hla_working_dir = File::Spec->abs2rel($hla_working_dir);	
	
	my %calledHLA;
	{
		open(CALLS, '<', $call2_HLAtypes) or die "Cannot open $call2_HLAtypes";
		my $headerLine = <CALLS>;
		chomp($headerLine);
		my @header_fields = split(/\t/, $headerLine);
		while(<CALLS>)
		{
			my $line = $_;
			chomp($line);
			next unless($line);
			my @line_fields = split(/\t/, $line, -1);
			die unless(scalar(@line_fields) == scalar(@header_fields));
			my %line = (mesh @header_fields, @line_fields);
			die unless(defined $line{'Locus'});
			die unless(defined $line{'Chromosome'});
			die unless(defined $line{'Allele'});
			my @alleles = split(/;/, $line{'Allele'});
			$calledHLA{$line{'Locus'}}{$line{'Chromosome'}} = \@alleles;
		}
		close(CALLS);
	}
	
	{
	
		my $fn_relevant_reads = $working_dir_thisSample . '/reads_2_gene.txt';

		open(READS, '<', $fn_relevant_reads) or die "Cannot open $fn_relevant_reads";
		while(<READS>)
		{
			my $line = $_;
			chomp($line);
			next unless($line);
			my @line_fields = split(/\t/, $line);
			die unless(scalar(@line_fields) == 2);
			$call2_hla_relevant_readIDs_primaryBAM{$line_fields[0]}++;
		}
		close(READS);
		
		print "Identified ", scalar(keys %call2_hla_relevant_readIDs_primaryBAM), " relevant read IDs.\n";		
			
		if($action eq 'somatic')
		{
			my $fn_relevant_reads_tumor = $working_dir_tumorSample . '/reads_2_gene.txt';
		
			open(TUMORREADS, '<', $fn_relevant_reads_tumor) or die "Cannot open $fn_relevant_reads_tumor";
			while(<TUMORREADS>)
			{
				my $line = $_;
				chomp($line);
				next unless($line);
				my @line_fields = split(/\t/, $line);
				die unless(scalar(@line_fields) == 2);
				$somatic_hla_relevant_readIDs_tumorBAM{$line_fields[0]}++;
			}
			close(TUMORREADS);	

			print "Identified ", scalar(keys %somatic_hla_relevant_readIDs_tumorBAM), " relevant tumor read IDs.\n";
		}
	}
	
	my %missing_gene_calls;
	my %have_raw_sequences_for_genes;
	my @files_raw_sequences = glob($full_graph_dir . '/pseudoGenomic_fullLengthMapping/raw_HLA-*.fa');
	foreach my $file (@files_raw_sequences)
	{
		die unless($file =~ /pseudoGenomic_fullLengthMapping\/raw_HLA-(\w+).fa/);
		my $gene = $1;
		$have_raw_sequences_for_genes{$gene}++;
		unless(exists $calledHLA{$gene})
		{
			$missing_gene_calls{$gene}++;
		}
	}
	
	open(REMAP, '>', $call2_fn_mapping) or die "Cannot open $call2_fn_mapping";
	open(REMAPHAP, '>', $call2_fn_mapping_whichHaplotype) or die "Cannot open $call2_fn_mapping_whichHaplotype";
	
	my %added_alleles;	
	foreach my $locusID (keys %calledHLA)
	{
		my $fn_call_sequences = $full_graph_dir . '/pseudoGenomic_fullLengthMapping/raw_HLA-' . $locusID . '.fa';
		unless(-e $fn_call_sequences)
		{
			die "File $fn_call_sequences missing";
			next;
		}
		my $locus_seq_href = readFASTA($fn_call_sequences);
		foreach my $chromosome (keys %{$calledHLA{$locusID}})
		{
			my @alleles = @{$calledHLA{$locusID}{$chromosome}};
			foreach my $allele (@alleles)
			{
				die "Allele $allele not existing in $fn_call_sequences" unless(exists $locus_seq_href->{$allele});
				print REMAPHAP join("\t", $locusID, $chromosome, $allele), "\n";				
				unless($added_alleles{$allele})
				{
					print REMAP '>', $allele, "\n", $locus_seq_href->{$allele}, "\n";
					$added_alleles{$allele}++;
				}
			}
		}
		
	}
	
	foreach my $locusID (keys %missing_gene_calls)
	{
		my $fn_locus_sequences = $full_graph_dir . '/pseudoGenomic_fullLengthMapping/raw_HLA-' . $locusID . '.fa';
		unless(-e $fn_locus_sequences)
		{
			die "File $fn_locus_sequences missing";
			next;
		}	
		my $locus_seq_href = readFASTA($fn_locus_sequences);
		foreach my $alleleID (keys %$locus_seq_href)
		{
			next unless($alleleID =~ /\*/);
			print REMAPHAP join("\t", $locusID, '?', $alleleID), "\n";						
			unless($added_alleles{$alleleID})
			{
				print REMAP '>', $alleleID, "\n", $locus_seq_href->{$alleleID}, "\n";
				$added_alleles{$alleleID}++;
			}		
		}
	}
	
	my $fn_pgf_n_masked = $full_graph_dir . '/pseudoGenomic_fullLengthMapping/PGF_with_Ns.fa';
	my $pgf_href = readFASTA($fn_pgf_n_masked);
	die unless(scalar(keys %$pgf_href) == 1);
	foreach my $seqID (keys %$pgf_href)
	{
		print REMAP '>', $seqID, "\n", $pgf_href->{$seqID}, "\n";	
	}
	print REMAP "\n";
	close(REMAP);
	close(REMAPHAP);
	
	print "\nGenerated remapping file: $call2_fn_mapping\n";
	
	# next steps
	# - reconstruct eight haplotypes from graph segment sequences
	# - for each called exonic allele, find the closest genomic allele(s)
	# - put together a file that has
	#   - for all called exonic sequences, the corresponding closest genomic allele sequences, and the corresponding PGF alignment, and the graph level / start position in the global MSA (?)
	#   - the eight B38 haplotypes, with the corresponding HLA alleles substituted -
	#              + ideally: with the genomic allele at the locus that is most closely related to the one already in the haplotype
	#              + (do we need a distance matrix)?
	#              + plus the pairwise alignment to PGF
	
}

my $mapAgainstCompleteGenome;
extractRelevantReadsFromBAM(
	$known_references_dir, 
	$BAM,
	$working_dir_thisSample,
	$target_FASTQ_1,
	$target_FASTQ_2,
	$target_FASTQ_U,
	$longReads,
	\$mapAgainstCompleteGenome,
);

my $mapAgainstCompleteGenome_tumor;
if($action eq 'somatic')
{
	extractRelevantReadsFromBAM(
		$known_references_dir, 
		$BAM_tumor,
		$somatic_processing_dir_tumor,
		$target_FASTQ_1_tumor,
		$target_FASTQ_2_tumor,
		$target_FASTQ_U_tumor,
		$longReads,
		\$mapAgainstCompleteGenome_tumor,
	);
}

if(($action eq 'call') and ($twoStageReadExtraction))
{
	print "INFO: Carry out two-stage read extraction\n";
	if(not $mapAgainstCompleteGenome)
	{
		print "INFO: You activated --twoStageReadExtraction, but for this input file no unmapped reads are being extracted -- the two-stage extraction process should not do any harm, but will also not speed things up.\n"; # todo remove
	}	
	
	my $fn_mapping_PRGonly = $full_graph_dir. '/mapping_PRGonly/referenceGenome.fa';
	die "Missing file: $fn_mapping_PRGonly" unless(-e $fn_mapping_PRGonly);
	
	my $BWA_index_file_bwt = $fn_mapping_PRGonly . '.bwt';
	die "File $fn_mapping_PRGonly does not seem to be bwa-indexed - please run 'bwa index $fn_mapping_PRGonly' and try again" unless(-e $BWA_index_file_bwt);
	
	my $bwa_x = $longReads ? ' -x ' . $longReads . ' ' : '';
	my $bwa_cmd;
	if($longReads)
	{
		$bwa_cmd = "$bwa_bin mem -t $maxThreads $bwa_x $fn_mapping_PRGonly $target_FASTQ_U";
	}
	else
	{
		$bwa_cmd = "$bwa_bin mem -t $maxThreads $bwa_x $fn_mapping_PRGonly $target_FASTQ_1 $target_FASTQ_2";
	}
	print $bwa_cmd, "\n";
	my %readIDs_filter;
	open(BWAPIPE, $bwa_cmd . '|') or die "Cannot open pipe to BWA command: $bwa_cmd";
	while(<BWAPIPE>)
	{
		my $line = $_;
		chomp($line);
		next unless($line);
		next if(substr($line, 0, 1) eq '@');
		my @line_fields = split(/\t/, $line);
		my $FLAGS = $line_fields[1];
		die unless(defined $FLAGS);
		my $readID = $line_fields[0];
		my $isUnmapped = (($FLAGS & 4) ? 1 : 0);
		if(! $isUnmapped)
		{
			$readIDs_filter{$readID}++;
		}
	}
	close(BWAPIPE) or die "Could not execute BWA command: $bwa_cmd";
	
	print "INFO: Two-stage extraction targeting " . scalar(keys %readIDs_filter) . " read IDs.\n";
	
	my $target_FASTQ_1_preFilter = $target_FASTQ_1 . '.prefilter';
	my $target_FASTQ_2_preFilter = $target_FASTQ_2 . '.prefilter';
	my $target_FASTQ_U_preFilter = $target_FASTQ_U . '.prefilter';
	my $mv_command = qq(mv $target_FASTQ_1 $target_FASTQ_1_preFilter && mv $target_FASTQ_2 $target_FASTQ_2_preFilter && mv $target_FASTQ_U $target_FASTQ_U_preFilter);
	system($mv_command) and die "Move command $mv_command failed";
	
	filterReadIDs([$target_FASTQ_1_preFilter, $target_FASTQ_2_preFilter, $target_FASTQ_U_preFilter], \%readIDs_filter, '.filtered', 1);
	my $target_FASTQ_1_postFilter = $target_FASTQ_1_preFilter . '.filtered';
	my $target_FASTQ_2_postFilter = $target_FASTQ_2_preFilter . '.filtered';
	my $target_FASTQ_U_postFilter = $target_FASTQ_U_preFilter . '.filtered';	
	die "Missing file $target_FASTQ_1_postFilter" unless(-e $target_FASTQ_1_postFilter);
	die "Missing file $target_FASTQ_2_postFilter" unless(-e $target_FASTQ_2_postFilter);
	die "Missing file $target_FASTQ_U_postFilter" unless(-e $target_FASTQ_U_postFilter);
	
	my $mv_command_2 = qq(mv $target_FASTQ_1_postFilter $target_FASTQ_1 && mv $target_FASTQ_2_postFilter $target_FASTQ_2 && mv $target_FASTQ_U_postFilter $target_FASTQ_U);
	system($mv_command_2) and die "Move command $mv_command_2 (II) failed";	
	
}	

if($action eq 'call')
{
	if($extractExonkMerCounts)
	{
		die if($longReads);
		die unless(-e '../exonkMerExtraction/GRCh38.forkMers');
		die unless(-e '../exonkMerExtraction/exonCoordinates_manual.txt');	
		
		my $command_extraction = qq(perl extractkMerCounts.pl --sampleID $sampleID --outputDirectory $working_dir_thisSample --referenceGenome ../exonkMerExtraction/GRCh38.forkMers --exonCoordinates ../exonkMerExtraction/exonCoordinates_manual.txt --FASTQ1 $target_FASTQ_1 --FASTQ2 $target_FASTQ_2 --bwa_bin $bwa_bin --samtools_bin $samtools_bin --maxThreads $maxThreads);
		print "Now executing: $command_extraction\n";
		system($command_extraction) and die "Command $command_extraction failed\n";
	}
	else
	{
		my $host = hostname();
		# my $MHC_PRG_2_bin = (($host =~ /rescomp/) or ($host =~ /comp[ABC]/)) ? '../bin_cluster3/HLA-LA' : '../bin/HLA-LA';
		my $MHC_PRG_2_bin = '../bin/HLA-LA';

		my $previous_dir = getcwd;
		chdir($this_bin_dir) or die "Cannot cd into $this_bin_dir";

		die "Binary $MHC_PRG_2_bin not there!" unless(-e $MHC_PRG_2_bin);
		my $command_MHC_PRG = qq($MHC_PRG_2_bin --action HLA --maxThreads $maxThreads --sampleID $sampleID --outputDirectory $working_dir_thisSample --PRG_graph_dir $full_graph_dir --FASTQU $target_FASTQ_U --FASTQ1 $target_FASTQ_1 --FASTQ2 $target_FASTQ_2 --bwa_bin $bwa_bin --samtools_bin $samtools_bin --mapAgainstCompleteGenome $mapAgainstCompleteGenome --longReads $longReads);
		
		print "\nNow executing:\n$command_MHC_PRG\n";
		if(system($command_MHC_PRG) != 0) 
		{
			die "HLA-LA execution not successful. Command was $command_MHC_PRG\n";
		}

		chdir($previous_dir) or die "Cannot cd into $previous_dir";
	}
}
elsif(($action eq 'call2') or ($action eq 'somatic'))
{
	die "Long-read mode not supported yet for modes 'call2' or 'somatic'" if($longReads);

	unless(-e 'Perl/callsToVCF.pl')
	{
		die "Script Perl/callsToVCF.pl not found - are you calling me from the right directory (i.e. the 'src' directory of HLA*LA)?";
	}
	my $call2_HLAtypes_VCF = $call2_HLAtypes . '.VCF';
	my $cmd_generate_VCF_normalCalls = qq(perl Perl/callsToVCF.pl --graph $graph --callFile $call2_HLAtypes --VCFout $call2_HLAtypes_VCF --sampleID $sampleID);
	system($cmd_generate_VCF_normalCalls) and die "VCF generation command (from G group calls) $cmd_generate_VCF_normalCalls failed";

	my $cmd_bwa_index = qq($samtools_bin faidx $call2_fn_mapping && $bwa_bin index $call2_fn_mapping);
	system($cmd_bwa_index) and die "Could not execute: $cmd_bwa_index";
	
	my $BAM_PROJECTED_NORMAL = $call2_processing_dir . '/remapped_and_projected.bam';
	my $SAM_UNPROJECTED_NORMAL = $call2_processing_dir . '/remapped_and_projected.bam.unprojected.SAM';
	my $referenceGenome_for_GATK = filterReadsAndGenerateProjectedBAM($target_FASTQ_1, $target_FASTQ_2, $target_FASTQ_U, \%call2_hla_relevant_readIDs_primaryBAM, $BAM_PROJECTED_NORMAL, $call2_fn_mapping, $longReads, $mapAgainstCompleteGenome);
	die "Missing file: $SAM_UNPROJECTED_NORMAL" unless(-e $SAM_UNPROJECTED_NORMAL);
	
	my $call2_outputPrefix = $call2_processing_dir . '/haplotypeHMMInput';

	my $dict = $referenceGenome_for_GATK;
	$dict =~ s/\.fa$/.dict/;
	if(-e $dict)
	{
		unlink($dict) or die "Cannot unlink $dict";
	}
	
	my $picard_firstPart = (-x $picard_bin) ? $picard_bin : qq(java -jar $picard_bin);
	my $cmd_Picard_1 = qq($picard_firstPart CreateSequenceDictionary R=$referenceGenome_for_GATK O=$dict);
	print "Now executing: $cmd_Picard_1 \n";
	system($cmd_Picard_1) and die "Command $cmd_Picard_1 failed";
	
	my $cmd_samtools_index = qq($samtools_bin faidx $referenceGenome_for_GATK);
	print "Now executing: $cmd_samtools_index \n";
	system($cmd_samtools_index) and die "Command $cmd_samtools_index failed";
			
	my $BAM_PROJECTED_NORMAL_RG = $BAM_PROJECTED_NORMAL . '.rg.bam';		
	my $cmd_Picard_2 = qq($picard_firstPart AddOrReplaceReadGroups I=$BAM_PROJECTED_NORMAL O=$BAM_PROJECTED_NORMAL_RG  RGID=4  RGLB=lib1  RGPL=illumina  RGPU=unit1  RGSM=$sampleID; $samtools_bin index $BAM_PROJECTED_NORMAL_RG);
	print "Now executing: $cmd_Picard_2 \n";
	system($cmd_Picard_2) and die "Command $cmd_Picard_2 failed";
	
	my $gene_2_readIDs_href = extract_gene2readID_fromProjectedBAM($BAM_PROJECTED_NORMAL_RG, $samtools_bin);
	
	if($action eq 'call2')
	{
		my $haplotypeHMM_inputDataPreparation_cmd = qq(perl Perl/localReassembly.pl --graph PRG_MHC_GRCh38_withIMGT --inputSAM $SAM_UNPROJECTED_NORMAL --reference $call2_fn_mapping --samtools_bin $samtools_bin --outputPrefix $call2_outputPrefix); 
		print "Now executing: $haplotypeHMM_inputDataPreparation_cmd \n";
		system($haplotypeHMM_inputDataPreparation_cmd) and die "Command $haplotypeHMM_inputDataPreparation_cmd failed";	
				
		my $MHC_PRG_2_bin = '../bin/HLA-LA';

		my $previous_dir = getcwd;
		chdir($this_bin_dir) or die "Cannot cd into $this_bin_dir";

		die "Binary $MHC_PRG_2_bin not there!" unless(-e $MHC_PRG_2_bin);
		my $command_MHC_PRG = qq($MHC_PRG_2_bin --action readHMM --inputPrefix $call2_outputPrefix);
				
		print "\nNow executing:\n$command_MHC_PRG\n";
		if(system($command_MHC_PRG) != 0) 
		{
			die "HLA-LA execution not successful. Command was $command_MHC_PRG\n";
		}
		
		chdir($previous_dir) or die "Cannot cd into $previous_dir";				
		
		if(1 == 0)
		{				
			# my $GATK_bin;
			my $WhatsHap_bin;
			my $bcftools_bin;
			my $spades_bin = '/software/SPAdes/3.13.0/skylake/bin/spades.py';			
					
			if(1 == 0)
			{				
				my $call2_fn_VCF = $call2_processing_dir . '/HLA.VCF';
				my $cmd_GATK = qq($GATK_bin --java-options "-Xmx4g" HaplotypeCaller -R $referenceGenome_for_GATK -I $BAM_PROJECTED_NORMAL_RG -O $call2_fn_VCF);
				print "Now executing: $cmd_GATK \n";
				system($cmd_GATK) and die "Command $cmd_GATK failed";	
				
				print "Generated (unphased) VCF: $call2_fn_VCF\n";

				my $call2_fn_VCF_canonical = $call2_processing_dir . '/HLA.canonical.VCF';
				my $cmd_make_canonical = qq(perl Perl/convertToCanonicalVCF.pl --graph $graph --VCFin $call2_fn_VCF --VCFout $call2_fn_VCF_canonical);
				system($cmd_make_canonical) and die "Projection command $cmd_make_canonical failed";
				
				print "Generated canonical (unphased) VCF: $call2_fn_VCF_canonical\n";
				
				my $call2_fn_phased_VCF = $call2_processing_dir . '/HLA.phased.VCF';
				my $cmd_WhatsHap = qq($WhatsHap_bin phase --indels -o $call2_fn_phased_VCF --reference $referenceGenome_for_GATK $call2_fn_VCF $BAM_PROJECTED_NORMAL_RG);
				print "Now executing: $cmd_WhatsHap \n";
				system($cmd_WhatsHap) and die "Command $cmd_WhatsHap failed";	
				
				print "Generated (phased) VCF: $call2_fn_phased_VCF\n";
				
				my $call2_fn_phased_gzipped_VCF = $call2_fn_phased_VCF . '.gz';			
				my $gzip_VCF_command = qq(bgzip -c $call2_fn_phased_VCF > $call2_fn_phased_gzipped_VCF && tabix -p vcf $call2_fn_phased_gzipped_VCF);
				system($gzip_VCF_command) and die "Gzip command failed: $gzip_VCF_command";
				die "Gzipped VCF $call2_fn_phased_gzipped_VCF missing" unless(-e $call2_fn_phased_gzipped_VCF);
				
				my $call_BAM_readsPartitioned = $call2_processing_dir . '/HLA.phased.VCF.reads.bam';
				my $call_list_readsPartitioned = $call2_processing_dir . '/HLA.phased.VCF.reads.listByHaplotag';
				my $cmd_WhatsHap_partition_1 = qq($WhatsHap_bin haplotag --output $call_BAM_readsPartitioned --reference $referenceGenome_for_GATK --output-haplotag-list $call_list_readsPartitioned $call2_fn_phased_gzipped_VCF $BAM_PROJECTED_NORMAL_RG);
				print "Now executing: $cmd_WhatsHap_partition_1 \n";
				system($cmd_WhatsHap_partition_1) and die "Command $cmd_WhatsHap_partition_1 failed";	
				
				my $forAssembly_FASTQ_1 = $call2_processing_dir . '/readsForAssembly_R1.fastq';
				my $forAssembly_FASTQ_2 = $call2_processing_dir . '/readsForAssembly_R2.fastq';
				my $forAssembly_FASTQ_U = $call2_processing_dir . '/readsForAssembly_RU.fastq';
				
				my $cp_if_e = sub {my $input = shift; my $output = shift; if(-e $input){system("cp $input $output") and die "cp failed";}};
				$cp_if_e->($target_FASTQ_1, $forAssembly_FASTQ_1);
				$cp_if_e->($target_FASTQ_2, $forAssembly_FASTQ_2);
				$cp_if_e->($target_FASTQ_U, $forAssembly_FASTQ_U);
				
				my $refGenome_href = readFASTA($referenceGenome_for_GATK);
				my $phaseSets_href = getPhaseSetsPerGene($call2_fn_phased_VCF);
				my @genes = sort map {$_ =~ s/\*ref//g; $_} grep {$_ =~ /\*ref/} keys %$refGenome_href;
				print "Detected genes: " , join(", ", @genes), "\n";
				print "Phase set statistics\n";
				my %assembled_contigs;
				foreach my $gene (@genes)
				{
					die unless(exists $refGenome_href->{$gene . '*ref'});
					my $ref_length = length($refGenome_href->{$gene . '*ref'});
					print "\t", $gene, " [length: $ref_length]\n";
					
					if(1 == 0)
					{
						my $createTrimmedSpadesAssembly = sub {
							my $spades_outputFile = shift;
							my $FASTQ_1 = shift;
							my $FASTQ_2 = shift;
							my $FASTQ_U = shift;
							my $assembled_contigs_href = shift;
							my $contigID_prefix = shift;
							
							my $spades_outputDir = $spades_outputFile . '_spadesDir';
							if(-e $spades_outputDir)
							{
								my $cmd_del = qq(rm -rf $spades_outputDir);
								system($cmd_del);
							}
							
							mkdir($spades_outputDir) or die "Cannot mkdir $spades_outputDir";
									
							my $spades_outputFile_untrimmed = $spades_outputFile . '.untrimmed.fasta';
							
							my $BAM_reads_mapped_to_gene = $spades_outputFile_untrimmed . '.reads.bam';
							my $pileupFile_BAM_reads_mapped_to_gene = $spades_outputFile_untrimmed . '.reads.bam.pileup.txt'; 
							
							my $careful = ($spades_outputFile !~ /_complete_/) ? ' --careful ' : '';
							my $spades_cmd = qq($spades_bin -o $spades_outputDir $careful -1 $FASTQ_1 -2 $FASTQ_2);
							if(system($spades_cmd) == 0)
							{
								my $contigs_fn = $spades_outputDir . '/contigs.fasta';
								die "Contigs file $contigs_fn missing" unless(-e $contigs_fn);
								my $copy_cmd = qq(cp $contigs_fn $spades_outputFile_untrimmed);
								system($copy_cmd) and die "Copying command $copy_cmd failed";
								
								my $cmd_bwa_map_toGene = qq($bwa_bin index $spades_outputFile_untrimmed && $bwa_bin mem $spades_outputFile_untrimmed $FASTQ_1 $FASTQ_2 | $samtools_bin sort --reference $spades_outputFile_untrimmed --threads 1 -O BAM - > $BAM_reads_mapped_to_gene && $samtools_bin index $BAM_reads_mapped_to_gene);
								system($cmd_bwa_map_toGene) and die "Could not execute: $cmd_bwa_map_toGene";		

								my $cmd_pileup = qq(perl Perl/BAMfrequencies.pl --BAM $BAM_reads_mapped_to_gene --referenceGenome $spades_outputFile_untrimmed --samtools_bin $samtools_bin --outputFile $pileupFile_BAM_reads_mapped_to_gene);
								system($cmd_pileup) and die "Pileup command $cmd_pileup failed";
								
								trimContigs($spades_outputFile_untrimmed, $pileupFile_BAM_reads_mapped_to_gene, $spades_outputFile);
								
								if($assembled_contigs_href)
								{
									my $fasta_href = readFASTA($spades_outputFile);
									foreach my $contigID (keys %$fasta_href)
									{
										my $newContigID = $contigID_prefix . $contigID . '_trimmed';
										$assembled_contigs_href->{$newContigID} = $fasta_href->{$contigID};
									}
								}
								
								my $cmd_del = qq(rm -rf $spades_outputDir);
								# system($cmd_del);						
							}
							else
							{
								warn "Spades cmd failed: $spades_cmd";
							}
							
							
						};
						
						my $spades_outputFile_wholeGene = $call2_processing_dir . '/assembled_' . $gene . '_complete';
						if(exists $gene_2_readIDs_href->{$gene})
						{
							filterReadIDs([$forAssembly_FASTQ_1, $forAssembly_FASTQ_2, $forAssembly_FASTQ_U], $gene_2_readIDs_href->{$gene}, '_filteredFor' . $gene . '.fastq');
								
							my $forAssembly_oneGene_FASTQ_1 = $call2_processing_dir . '/readsForAssembly_R1.fastq' . '_filteredFor' . $gene . '.fastq';
							my $forAssembly_oneFene_FASTQ_2 = $call2_processing_dir . '/readsForAssembly_R2.fastq' . '_filteredFor' . $gene . '.fastq';
							my $forAssembly_oneGene_FASTQ_U = $call2_processing_dir . '/readsForAssembly_RU.fastq' . '_filteredFor' . $gene . '.fastq';
							
							$createTrimmedSpadesAssembly->($spades_outputFile_wholeGene, $forAssembly_oneGene_FASTQ_1, $forAssembly_oneFene_FASTQ_2, $forAssembly_oneGene_FASTQ_U, \%assembled_contigs, 'Gene' . $gene . '_complete_');					
						}

						if(exists $phaseSets_href->{$gene})
						{
							my @phaseSets = sort keys %{$phaseSets_href->{$gene}};
							foreach my $PS (@phaseSets)
							{
								
								print "\t\tPhase set ", $PS, ", length ", $phaseSets_href->{$gene}{$PS}{coordinates}[1] - $phaseSets_href->{$gene}{$PS}{coordinates}[0] + 1, " [", $phaseSets_href->{$gene}{$PS}{coordinates}[0], " - ", $phaseSets_href->{$gene}{$PS}{coordinates}[1], "]\n";
								
								filterReadsForPhaseSet($forAssembly_FASTQ_1, $forAssembly_FASTQ_2, $forAssembly_FASTQ_U, $gene . '_PS' . $PS, $call_list_readsPartitioned, $PS);
								
								foreach my $H (qw/H1 H2/)
								{
									my $fn_filtered_R1 = $forAssembly_FASTQ_1 . $gene . '_PS' . $PS . '_' . $H . '.fq';
									my $fn_filtered_R2 = $forAssembly_FASTQ_2 . $gene . '_PS' . $PS . '_' . $H . '.fq';				
									my $fn_filtered_U = $forAssembly_FASTQ_U . $gene . '_PS' . $PS . '_' . $H . '.fq';

									my $spades_outputFile = $call2_processing_dir . '/assembled_' . $gene . '_PS' . $PS . '_' . $H . '.fasta'; 

									$createTrimmedSpadesAssembly->($spades_outputFile, $fn_filtered_R1, $fn_filtered_R2, $fn_filtered_U, \%assembled_contigs, 'Gene' . $gene . '_PS' . $PS . '_' . $H . '_');
								}						
							}
						}
					}
				}
				
				print "\nReassembly done. Have " . scalar(keys %assembled_contigs) . " contigs.\n";
				
				my $FASTA_all_contigs = $call2_processing_dir . '/assembled_all_contigs.fa';
				writeFASTA($FASTA_all_contigs, \%assembled_contigs);
				my $BAM_all_contigs = $call2_processing_dir . '/assembled_all_contigs.fa.bam';
				
				my $VCF_all_contigs = $call2_processing_dir . '/assembled_all_contigs.fa.vcf';

				alignContigs_geneAware($referenceGenome_for_GATK, $FASTA_all_contigs, $BAM_all_contigs, $bwa_bin, $samtools_bin);
				
				my $cmd_bwa_allContigs_ref = qq($bwa_bin index $referenceGenome_for_GATK && $bwa_bin mem $referenceGenome_for_GATK $FASTA_all_contigs | $samtools_bin sort --reference $referenceGenome_for_GATK --threads 1 -O BAM - > $BAM_all_contigs && $samtools_bin index $BAM_all_contigs);
				system($cmd_bwa_allContigs_ref) and die "Could not execute: $cmd_bwa_allContigs_ref";
				
				my $cmd_variantCalling_BCFTools = qq($bcftools_bin mpileup -f $referenceGenome_for_GATK $BAM_all_contigs | $bcftools_bin call -m -Ov -o $VCF_all_contigs);
				print "BCFTools: $cmd_variantCalling_BCFTools \n\n";
				system($cmd_variantCalling_BCFTools) and die "Could not execute: $cmd_variantCalling_BCFTools";

				foreach my $R (1, 2)
				{
					my $inputFASTQ = ($R == 1) ? $target_FASTQ_1 : $target_FASTQ_2;
					my $call_BAM_readsPartitioned_H1 = $call2_processing_dir . '/HLA.phased.VCF.reads.fastq_H1_R' . $R . '.fastq';
					my $call_BAM_readsPartitioned_H2 = $call2_processing_dir . '/HLA.phased.VCF.reads.fastq_H2_R' . $R . '.fastq';
					my $cmd_WhatsHap_partition_2 = qq($WhatsHap_bin split --output-h1 $call_BAM_readsPartitioned --output-h2 $call_BAM_readsPartitioned_H2 $inputFASTQ $call_list_readsPartitioned);
					print "Now executing: $cmd_WhatsHap_partition_2 \n";
					system($cmd_WhatsHap_partition_2) and die "Command $cmd_WhatsHap_partition_2 failed";			
				}	

				my $call2_fn_phased_VCF_canonical = $call2_processing_dir . '/HLA.phased.canonical.VCF';
				my $cmd_make_canonical_phased = qq(perl Perl/convertToCanonicalVCF.pl --graph $graph --VCFin $call2_fn_phased_VCF --VCFout $call2_fn_phased_VCF_canonical);
				system($cmd_make_canonical_phased) and die "Projection command $cmd_make_canonical_phased failed";
				
				print "Generated canonical (phased) VCF: $call2_fn_phased_VCF_canonical\n";
				
				my $genotypes_originalCalls_href = readVCF($call2_HLAtypes_VCF);
				my $genotypes_GATK_href = readVCF($call2_fn_phased_VCF_canonical);
				
				my $originalCalls_VCF_total = 0;
				my $originalCalls_VCF_missing = 0;
				my $originalCalls_VCF_matching = 0;
				my $originalCalls_VCF_mismatching = 0;
				foreach my $position (keys %$genotypes_originalCalls_href)
				{ 
					$originalCalls_VCF_total++;
					if(exists $genotypes_GATK_href->{$position})
					{
						my @gt_1 = sort @{$genotypes_originalCalls_href->{$position}};
						my @gt_2 = sort @{$genotypes_GATK_href->{$position}};
						die unless(scalar(@gt_1) == 2);
						die unless(scalar(@gt_2) == 2);
						if(($gt_1[0] eq $gt_2[0]) and ($gt_1[1] eq $gt_2[1]))
						{
							$originalCalls_VCF_matching++;
						}
						else
						{
							$originalCalls_VCF_mismatching++;
							warn "VCF mismatch $position: " . join('/', @gt_1) . " in original calls, " . join('/', @gt_2) . " in GATK-based VCF.\n";
						}
					}
					else
					{
						$originalCalls_VCF_missing++;
					}
				}
				
				print "\nVCF comparison $call2_HLAtypes_VCF (original calls) v/s $call2_fn_phased_VCF_canonical (GATK-based calls):\n";
				print "\t", "Total positions in original calls: $originalCalls_VCF_total\n";
				print "\t", "... missing from GATK calls: $originalCalls_VCF_missing\n";
				print "\t", "... consistent with GATK calls: $originalCalls_VCF_matching\n";
				print "\t", "... mismatch against GATK calls: $originalCalls_VCF_mismatching\n";
				print "\n";		
			}
			
			if(1 == 0)
			{
				my $call2_fn_VCF = $call2_processing_dir . '/HLA_FreeBayes.VCF';
				my $cmd_FreeBayes = qq(freebayes -f $referenceGenome_for_GATK $BAM_PROJECTED_NORMAL_RG > $call2_fn_VCF);
				print "Now executing: $cmd_FreeBayes \n";
				system($cmd_FreeBayes) and die "Command $cmd_FreeBayes failed";	
				
				print "Generated (unphased) VCF: $call2_fn_VCF\n";

				my $call2_fn_VCF_canonical = $call2_processing_dir . '/HLA_FreeBayes.canonical.VCF';
				my $cmd_make_canonical = qq(perl Perl/convertToCanonicalVCF.pl --graph $graph --VCFin $call2_fn_VCF --VCFout $call2_fn_VCF_canonical);
				system($cmd_make_canonical) and die "Projection command $cmd_make_canonical failed";
				
				print "Generated canonical (unphased) VCF: $call2_fn_VCF_canonical\n";
				
				my $call2_fn_phased_VCF = $call2_processing_dir . '/HLA_FreeBayes.phased.VCF';
				my $cmd_WhatsHap = qq($WhatsHap_bin phase --indels -o $call2_fn_phased_VCF --reference $referenceGenome_for_GATK $call2_fn_VCF $BAM_PROJECTED_NORMAL_RG);
				print "Now executing: $cmd_WhatsHap \n";
				system($cmd_WhatsHap) and die "Command $cmd_WhatsHap failed";	
				
				print "Generated (phased) VCF: $call2_fn_phased_VCF\n";

				my $call2_fn_phased_VCF_canonical = $call2_processing_dir . '/HLA_FreeBayes.phased.canonical.VCF';
				my $cmd_make_canonical_phased = qq(perl Perl/convertToCanonicalVCF.pl --graph $graph --VCFin $call2_fn_phased_VCF --VCFout $call2_fn_phased_VCF_canonical);
				system($cmd_make_canonical_phased) and die "Projection command $cmd_make_canonical_phased failed";
				
				print "Generated canonical (phased) VCF: $call2_fn_phased_VCF_canonical\n";
				
				my $genotypes_originalCalls_href = readVCF($call2_HLAtypes_VCF);
				my $genotypes_GATK_href = readVCF($call2_fn_phased_VCF_canonical);
				
				my $originalCalls_VCF_total = 0;
				my $originalCalls_VCF_missing = 0;
				my $originalCalls_VCF_matching = 0;
				my $originalCalls_VCF_mismatching = 0;
				foreach my $position (keys %$genotypes_originalCalls_href)
				{ 
					$originalCalls_VCF_total++;
					if(exists $genotypes_GATK_href->{$position})
					{
						my @gt_1 = sort @{$genotypes_originalCalls_href->{$position}};
						my @gt_2 = sort @{$genotypes_GATK_href->{$position}};
						die unless(scalar(@gt_1) == 2);
						die unless(scalar(@gt_2) == 2);
						if(($gt_1[0] eq $gt_2[0]) and ($gt_1[1] eq $gt_2[1]))
						{
							$originalCalls_VCF_matching++;
						}
						else
						{
							$originalCalls_VCF_mismatching++;
							warn "VCF mismatch $position: " . join('/', @gt_1) . " in original calls, " . join('/', @gt_2) . " in FreeBayes-based VCF.\n";
						}
					}
					else
					{
						$originalCalls_VCF_missing++;
					}
				}
				
				print "\nVCF comparison $call2_HLAtypes_VCF (original calls) v/s $call2_fn_phased_VCF_canonical (FreeBayes-based calls):\n";
				print "\t", "Total positions in original calls: $originalCalls_VCF_total\n";
				print "\t", "... missing from FreeBayes calls: $originalCalls_VCF_missing\n";
				print "\t", "... consistent with FreeBayes calls: $originalCalls_VCF_matching\n";
				print "\t", "... mismatch against FreeBayes calls: $originalCalls_VCF_mismatching\n";
				print "\n";		
			}
		}			
	}
	elsif($action eq 'somatic')
	{
		my $BAM_PROJECTED_TUMOR = $somatic_processing_dir_tumor . '/remapped_and_projected_tumor.bam';
		filterReadsAndGenerateProjectedBAM($target_FASTQ_1_tumor, $target_FASTQ_2_tumor, $target_FASTQ_U_tumor, \%somatic_hla_relevant_readIDs_tumorBAM, $BAM_PROJECTED_TUMOR, $call2_fn_mapping, $longReads, $mapAgainstCompleteGenome);

		my $BAM_PROJECTED_TUMOR_RG = $BAM_PROJECTED_TUMOR . '.rg.bam';					
		my $cmd_Picard_3 = qq($picard_firstPart AddOrReplaceReadGroups I=$BAM_PROJECTED_TUMOR O=$BAM_PROJECTED_TUMOR_RG  RGID=5  RGLB=lib1  RGPL=illumina  RGPU=unit1  RGSM=$sampleID_tumor; $samtools_bin index $BAM_PROJECTED_TUMOR_RG);
		print "Now executing: $cmd_Picard_3 \n";
		system($cmd_Picard_3) and die "Command $cmd_Picard_3 failed";

		my $somatic_fn_VCF = $call2_processing_dir . '/HLA_somatic.VCF';
		my $cmd_GATK = qq($GATK_bin --java-options "-Xmx4g" Mutect2 -R $referenceGenome_for_GATK -I $BAM_PROJECTED_TUMOR_RG -tumor $sampleID_tumor -I $BAM_PROJECTED_NORMAL_RG -normal $sampleID -O $somatic_fn_VCF);
		print "Now executing: $cmd_GATK \n";
		system($cmd_GATK) and die "Command $cmd_GATK failed";	
		 
		print "Generated VCF: $somatic_fn_VCF\n";

		my $somatic_fn_VCF_canonical = $call2_processing_dir . '/HLA_somatic.canonical.VCF';
		my $cmd_make_canonical = qq(perl Perl/convertToCanonicalVCF.pl --graph $graph --VCFin $somatic_fn_VCF --VCFout $somatic_fn_VCF_canonical);
		system($cmd_make_canonical) and die "Projection command $cmd_make_canonical failed";
		
		print "Generated canonical VCF: $somatic_fn_VCF_canonical\n";
	}
	else
	{
		die;
	}
}	

sub alignContigs_geneAware
{
	my $referenceGenome_for_GATK = shift;
	my $FASTA_all_contigs = shift;
	my $BAM_all_contigs = shift;
	my $bwa_bin = shift;
	my $samtools_bin = shift;
	
	my $reference_href = readFASTA($referenceGenome_for_GATK);
	my $contigs_href = readFASTA($FASTA_all_contigs);
	
	my %contigsByGene;
	foreach my $contigID (keys %$contigs_href)
	{
		die "alignContigs_geneAware(..): Cannot parse contig ID: '$contigID'" unless($contigID =~ /Gene(\w+?)_/);
		my $geneID = $1;
		$contigsByGene{$geneID}{$contigID}++;
		my $refID = $geneID . '*ref';
		die "Missing reference contig ID: '$refID' for $contigID" unless(exists $reference_href->{$refID});
	}
	
	my @BAMs;
	foreach my $geneID (keys %contigsByGene)
	{
		my $refFile = $BAM_all_contigs . '.byGene.' . $geneID . '.ref.fa';
		my $contigsFile = $BAM_all_contigs . '.byGene.' . $geneID . '.contigs.fa';
		my $contigsFileBAM = $BAM_all_contigs . '.byGene.' . $geneID . '.contigs.fa.bam';
		
		my $refID = $geneID . '*ref';
		
		open(REF, '>', $refFile) or die "Cannot open $refFile";
		print REF '>', $refID, "\n";
		print REF $reference_href->{$refID}, "\n";
		close(REF);
		
		open(TOMAP, '>', $contigsFile) or die "Cannot open $contigsFile";
		foreach my $contigID (keys %{$contigsByGene{$geneID}})
		{
			print TOMAP '>', $contigID, "\n";
			print TOMAP $contigs_href->{$contigID}, "\n";
		}
		close(TOMAP);
		
		my $cmd_bwa_allContigs_ref = qq($bwa_bin index $refFile && $bwa_bin mem $refFile $refFile | $samtools_bin sort --reference $refFile --threads 1 -O BAM - > $contigsFileBAM && $samtools_bin index $contigsFileBAM);
		system($cmd_bwa_allContigs_ref) and die "Command '$cmd_bwa_allContigs_ref' failed";
		
		push(@BAMs, $contigsFileBAM);
	}
	
	my $BAM_all_contigs_unsorted = $BAM_all_contigs . '.unsorted.bam';
	if(-e $BAM_all_contigs_unsorted)
	{
		unlink($BAM_all_contigs_unsorted) or die "Cannot unlink $BAM_all_contigs_unsorted";
	}	
	my $samtools_merge_cmd = qq($samtools_bin merge --reference $referenceGenome_for_GATK $BAM_all_contigs_unsorted ) . join(' ', @BAMs);
	system($samtools_merge_cmd) and die "Samtools command '$samtools_merge_cmd' failed";
	
	if(-e $BAM_all_contigs)
	{
		unlink($BAM_all_contigs) or die "Cannot unlink $BAM_all_contigs";
	}	
	my $samtools_sort_cmd = qq($samtools_bin sort --reference $referenceGenome_for_GATK -o $BAM_all_contigs $BAM_all_contigs_unsorted && $samtools_bin index $BAM_all_contigs);
	system($samtools_sort_cmd) and die "Samtools command '$samtools_sort_cmd' failed";
}
			
sub filterReadsForPhaseSet
{
	my $target_FASTQ_1 = shift;
	my $target_FASTQ_2 = shift;
	my $target_FASTQ_U = shift;
	my $suffixForFilter = shift;
	my $call_list_readsPartitioned = shift;
	my $PS = shift;
	
	die unless(defined $suffixForFilter);
	die unless(defined $call_list_readsPartitioned);
	die unless(defined $PS);
	
	my @targetReadIDs_H1;
	my @targetReadIDs_H2;
	open(READ2PS, '<', $call_list_readsPartitioned) or die "Cannot open $call_list_readsPartitioned";
	my $READ2PS_headerLine = <READ2PS>;
	chomp($READ2PS_headerLine);
	my @header_fields = split(/\t/, $READ2PS_headerLine);
	die unless(scalar(@header_fields) == 4);
	die unless($header_fields[0] eq '#readname');
	die unless($header_fields[1] eq 'haplotype');
	die unless($header_fields[2] eq 'phaseset');
	die unless($header_fields[3] eq 'chromosome');
	while(<READ2PS>)
	{
		my $line = $_;
		chomp($line);
		next unless($line);
		my @line_fields = split(/\t/, $line);
		die unless(scalar(@line_fields) == scalar(@header_fields));
		my $readID = $line_fields[0];
		my $haplotype = $line_fields[1];
		my $phaseSet = $line_fields[2];
		if($phaseSet eq $PS)
		{
			if($haplotype eq 'H1')
			{
				push(@targetReadIDs_H1, $readID);
			}
			else
			{
				die unless($haplotype eq 'H2');
				push(@targetReadIDs_H2, $readID);
			}
		}
	}
	close(READ2PS);
	
	my @filesToFilter = grep {defined $_} ($target_FASTQ_1, $target_FASTQ_2, $target_FASTQ_U);
	filterReadIDs(\@filesToFilter, {map {$_ => 1} @targetReadIDs_H1}, $suffixForFilter . '_H1.fq');
	filterReadIDs(\@filesToFilter, {map {$_ => 1} @targetReadIDs_H2}, $suffixForFilter . '_H2.fq');
}


sub getPhaseSetsPerGene
{
	my $VCF = shift;
	my %phaseSets_byGene;
	
	open(VCF, '<', $VCF) or die "Cannot open $VCF";
	while(<VCF>)
	{
		my $line = $_;
		chomp($line);
		$line =~ s/[\n\r]//g;
		next unless($line);
		next if(substr($line, 0, 2) eq '##');
		if(substr($line, 0, 1) eq '#')
		{
			my @header_fields = split(/\t/, $line);
			die "Unexpected header in $VCF" unless($header_fields[0] eq '#CHROM');  
			die unless($header_fields[1] eq 'POS');
			die unless($header_fields[2] eq 'ID');
			die unless($header_fields[3] eq 'REF');
			die unless($header_fields[4] eq 'ALT');
		}
		else
		{			
			my @line_fields = split(/\t/, $line);
			
			my $CHROM = $line_fields[0];
			next unless($CHROM =~ /^(\w+)\*ref$/);
			my $GENE = $1;
			my $POS = $line_fields[1];
			my $REF = $line_fields[3];
			my $ALT = $line_fields[4];

			my $FORMAT = $line_fields[8];
			my @FORMAT = split(/:/, $FORMAT);
			die unless($FORMAT[0] eq 'GT');
			
			my @idx_PS = grep {$FORMAT[$_] eq 'PS'} (0 .. $#FORMAT);

			my $gt = $line_fields[9];
			my @GT = split(/:/, $gt);
			die unless(scalar(@GT) == scalar(@FORMAT));
			
			die "Invalid genotype: $gt" unless($GT[0] =~ /^(\d+)([\/\|])(\d+)$/);
			my $sep = $2;			
			if($sep eq '|')
			{
				die unless(scalar(@idx_PS));
				my $phaseSet = $GT[$idx_PS[0]];
				
				if(not defined $phaseSets_byGene{$GENE}{$phaseSet}{coordinates}[0])
				{
					$phaseSets_byGene{$GENE}{$phaseSet}{coordinates}[0] = $POS;
					$phaseSets_byGene{$GENE}{$phaseSet}{coordinates}[1] = $POS;
				}
				die unless($phaseSets_byGene{$GENE}{$phaseSet}{coordinates}[0] <= $POS);
				if($POS > $phaseSets_byGene{$GENE}{$phaseSet}{coordinates}[1])
				{
					$phaseSets_byGene{$GENE}{$phaseSet}{coordinates}[1] = $POS;
				}				
			}
		}		
	}
	close(VCF);
	
	return \%phaseSets_byGene;
}

sub readVCF
{
	my $VCF = shift;
	my %gt_forReturn;
	
	my $pos_total = 0;
	my $pos_phased = 0;
	my %phaseSets;
	
	open(VCF, '<', $VCF) or die "Cannot open $VCF";
	while(<VCF>)
	{
		my $line = $_;
		chomp($line);
		$line =~ s/[\n\r]//g;
		next unless($line);
		next if(substr($line, 0, 2) eq '##');
		if(substr($line, 0, 1) eq '#')
		{
			my @header_fields = split(/\t/, $line);
			die unless($header_fields[0] eq '#CHROM');
			die unless($header_fields[1] eq 'POS');
			die unless($header_fields[2] eq 'ID');
			die unless($header_fields[3] eq 'REF');
			die unless($header_fields[4] eq 'ALT');
		}
		else
		{
			$pos_total++;
			
			my @line_fields = split(/\t/, $line);
			
			my $CHROM = $line_fields[0];
			my $POS = $line_fields[1];
			my $REF = $line_fields[3];
			my $ALT = $line_fields[4];
			
			my $POS_stop = $POS + length($REF) - 1;
			my $pos_key = $CHROM . ':' . $POS . '-' . $POS_stop;
			
			my @ALT = split(/,/, $ALT);
			my @alleles = ($REF, @ALT);
						
			my $FORMAT = $line_fields[8];
			my @FORMAT = split(/:/, $FORMAT);
			die unless($FORMAT[0] eq 'GT');
			
			my @idx_PS = grep {$FORMAT[$_] eq 'PS'} (0 .. $#FORMAT);

			my $gt = $line_fields[9];
			my @GT = split(/:/, $gt);
			die unless(scalar(@GT) == scalar(@FORMAT));
			
			die "Invalid genotype: $gt" unless($GT[0] =~ /^(\d+)([\/\|])(\d+)$/);
			my $a1 = $1;
			my $a2 = $3;
			my $sep = $2;
			$pos_phased++ if($sep eq '|');
			
			my @genotypes = map {die unless(defined $alleles[$_]); $alleles[$_]} ($a1, $a2);

			if(scalar(@idx_PS))
			{
				die unless(scalar(@idx_PS) == 1);
				$phaseSets{$GT[$idx_PS[0]]}++;
			}
			
			$gt_forReturn{$pos_key} = \@genotypes;
		}
		
	}
	close(VCF);
	
	print "Statistics $VCF\n";
	print "\t", "Positions: $pos_total \n";
	print "\t", "Positions phased: $pos_phased \n";
	print "\t", "Phase sets: ", scalar(keys %phaseSets), "\n";
	print "\n";
		
	return \%gt_forReturn;
	
}

sub readTextPileup
{
	my $pileupFile = shift;
	my %forReturn;
	open(PILEUP, '<', $pileupFile) or die "Cannot open $pileupFile";
	while(<PILEUP>)
	{
		my $line = $_;
		chomp($line);
		next unless($line);
		my @line_fields = split(/\t/, $line, -1);
		die Dumper("Weird number of fields in pileup file $pileupFile", \@line_fields) unless(scalar(@line_fields) == 4);
		my $contigID = $line_fields[0];
		my $pos = $line_fields[1];
		my $refBase = $line_fields[2];
		my $pileUpString = $line_fields[3];
		
		$forReturn{$contigID}{$pos}{ref} = $refBase;
		
		my @pileUpDetails = split(/;/, $pileUpString);
		foreach my $pileUpDetail (@pileUpDetails)
		{
			die unless($pileUpDetail =~ /^(.+)=(.+)$/); 
			my $char = $1;
			my $count = $2;
			next if($char =~ /\-/);
			die "Weird character in pileup ($pileupFile): $char" unless($char =~ /^[ACGT\(\)]+$/i);
			my $char_noStrand = uc($char);
			$char_noStrand =~ s/[\(\)]/-/g;
			$forReturn{$contigID}{$pos}{counts}{stranded}{$char} = $count;
			$forReturn{$contigID}{$pos}{counts}{nonstranded}{$char_noStrand} += $count;
			$forReturn{$contigID}{$pos}{coverage} += $count;			
		}
		
	}
	close(PILEUP);
	
	return \%forReturn;
}

sub trimContigs
{
	my $fasta = shift;
	my $pileupFile = shift;
	my $outputFile = shift;
	
	my $threshold = 5;
	my $slidingWindowSize = 5;
	
	my $pileup_href = readTextPileup($pileupFile);
	
	my $contigs_href = readFASTA($fasta);
	my %trimmed_contigs;
	foreach my $contigID (keys %$contigs_href)
	{
		my $contigSeq = $contigs_href->{$contigID};
		my $contigOriginalLength = length($contigSeq);
		my @coverages;
		$#coverages = length($contigSeq) - 1;
		if(exists $pileup_href->{$contigID})
		{
			for(my $pos_0Based = 0; $pos_0Based < length($contigSeq); $pos_0Based++)
			{
				my $refC = uc(substr($contigSeq, $pos_0Based, 1));
				my $pos_1Based = $pos_0Based + 1; 
				my $coverage = (exists $pileup_href->{$contigID}{$pos_1Based}{coverage}) ? $pileup_href->{$contigID}{$pos_1Based}{coverage} : 0;
				my $nonMajorityBase = 0;
				if($coverage)
				{
					my $supposedRefC = $pileup_href->{$contigID}{$pos_1Based}{ref};
					die "Reference character mismatch - $refC v/s $supposedRefC - contig $contigID, position $pos_1Based (1-based)" unless($refC eq $supposedRefC);
					warn Dumper("trimContigs surprise: '$refC' not present at position $pos_1Based of $contigID pileup", $pileup_href->{$contigID}{$pos_1Based}) unless($pileup_href->{$contigID}{$pos_1Based}{counts}{nonstranded}{$refC});
					my $refC_coverageProportion = $pileup_href->{$contigID}{$pos_1Based}{counts}{nonstranded}{$refC} / $pileup_href->{$contigID}{$pos_1Based}{coverage};
					if($refC_coverageProportion < 0.75)
					{
						$nonMajorityBase = 1;
						# print "trimContigs warning: $refC at 1-based position $pos_1Based of $contigID has coverage proportion $refC_coverageProportion - " . Dumper($pileup_href->{$contigID}{$pos_1Based});
					}
				}
				if(($coverage < $threshold) or $nonMajorityBase)
				{
					substr($contigSeq, $pos_0Based, 1) = 'N';
				}
				$coverages[$pos_0Based] = $coverage;
			}
			die unless(scalar(@coverages) == length($contigSeq));
			
			my $firstOKIndex = undef;
			for(my $pos_0Based = 0; $pos_0Based <= $#coverages; $pos_0Based++)
			{
				my $minIndex = $pos_0Based - $slidingWindowSize;
				my $maxIndex = $pos_0Based + $slidingWindowSize;
				$minIndex = 0 if($minIndex < 0);
				$maxIndex = $#coverages if($maxIndex > $#coverages);
				my @avgValues = @coverages[$minIndex .. $maxIndex];
				my $avgCoverage = sum(@avgValues)/scalar(@avgValues);
				if($avgCoverage > $threshold)
				{
					$firstOKIndex = $pos_0Based;
					last;
				}
			}
			
			my $lastOKIndex = undef;
			for(my $pos_0Based = $#coverages; $pos_0Based >= 0; $pos_0Based--)
			{
				my $minIndex = $pos_0Based - $slidingWindowSize;
				my $maxIndex = $pos_0Based + $slidingWindowSize;
				$minIndex = 0 if($minIndex < 0);
				$maxIndex = $#coverages if($maxIndex > $#coverages);
				my @avgValues = @coverages[$minIndex .. $maxIndex];
				my $avgCoverage = sum(@avgValues)/scalar(@avgValues);
				if($avgCoverage > $threshold)
				{
					$lastOKIndex = $pos_0Based;
					last;
				}
			}		
			
			if((defined $firstOKIndex) and (defined $lastOKIndex))
			{
				$contigSeq = substr($contigSeq, $firstOKIndex, $lastOKIndex - $firstOKIndex + 1);
				(my $contigSeq_N = $contigSeq) =~ s/[^N]//gi;
				my $number_remaining_Ns = length($contigSeq_N);
				$trimmed_contigs{$contigID} = $contigSeq;
				print "Contig trimming report for $contigID ($fasta) of length $contigOriginalLength: Extract from $firstOKIndex to $lastOKIndex, $number_remaining_Ns Ns remaining.\n";
			}
			else
			{
				$trimmed_contigs{$contigID . '_untrimmed'} = $contigSeq;				
				warn "Contig trimming report for $contigID ($fasta) of length $contigOriginalLength: No valid start and end positions found\n";
			}
		}
		else
		{
			$trimmed_contigs{$contigID . '_untrimmed'} = $contigSeq;							
			warn "Contig ID '$contigID' missing from pileup file $pileupFile" 
		}
	}
	
	writeFASTA($outputFile, \%trimmed_contigs);
}

sub randomizeReadOrder
{
	my $inputFiles_aref = shift;
	die unless(scalar(@$inputFiles_aref) >= 1);
	
	my @readData_per_file;
	foreach my $f (@$inputFiles_aref)
	{
		my @readData_thisFile;
		open(FASTQIN, '<', $f) or die "Cannot open $f";
		while(<FASTQIN>)
		{
			my $line = $_;
			chomp($line);
			next unless($line);
			my $readID_with_at = $line;
			die "Weird FASTQ format" unless(substr($readID_with_at, 0, 1) eq '@');
			my $readID = substr($readID_with_at, 1);
			$readID =~ s/\/[12]$//;

			my $seq = <FASTQIN>; chomp($seq);
			my $plus = <FASTQIN>; die unless(substr($plus, 0, 1) eq '+');
			my $qual = <FASTQIN>; chomp($qual);
			die "Weird FASTQ format (II)" unless(length($qual) == length($seq));

			push(@readData_thisFile, [$readID, $readID_with_at, $seq, $plus, $qual]);
		}
		close(FASTQIN);
		push(@readData_per_file, \@readData_thisFile);
	}
	die unless(scalar(@readData_per_file) == scalar(@$inputFiles_aref));
	
	my $n_reads = scalar(@{$readData_per_file[0]});
	die Dumper("Undefined elements in \@readData_per_file", $inputFiles_aref, [map {defined $_ ? 1 : 0} @readData_per_file]) unless(all {defined($_)} @readData_per_file);
	die unless(all {scalar(@{$_}) == $n_reads} @readData_per_file);

	my @readIDs = map {$_->[0]} @{$readData_per_file[0]};
	foreach my $readData_oneFile (@readData_per_file)
	{
		my @readIDs_thisFile = map {$_->[0]} @{$readData_oneFile};
		die unless(all {$readIDs[$_] eq $readIDs_thisFile[$_]} (0 .. $#readIDs));
	}
	
	my @readIDs_indices_shuffled = shuffle((0 .. $#readIDs));
	die unless(scalar(@readIDs_indices_shuffled) == scalar(@readIDs));
	for(my $fileI = 0; $fileI <= $#readData_per_file; $fileI++)
	{
		my $f = $inputFiles_aref->[$fileI];
		open(FASTQOUT, '>', $f) or die "Cannot open $f";
		foreach my $readIndex (@readIDs_indices_shuffled)
		{
			my $read_aref = $readData_per_file[$fileI][$readIndex];
			print FASTQOUT $read_aref->[1], "\n", $read_aref->[2], "\n", $read_aref->[3], $read_aref->[4], "\n"			
		}
		close(FASTQOUT);
	}
}

sub filterReadIDs
{
	my $inputFiles_aref = shift;
	my $filter_IDs_href = shift;
	my $addSuffix = shift;
	my $failOnError = shift;
	die unless(defined $addSuffix);
	
	my %total_readIDs;
	my %got_IDs;
	foreach my $f (@$inputFiles_aref)
	{
		my $fn_out = $f . $addSuffix;
		open(FASTQIN, '<', $f) or die "Cannot open $f";
		open(FASTQOUT, '>', $fn_out) or die "Cannot open $fn_out";
		while(<FASTQIN>)
		{
			my $line = $_;
			chomp($line);
			next unless($line);
			my $readID_with_at = $line;
			die "Weird FASTQ format" unless(substr($readID_with_at, 0, 1) eq '@');
			my $readID = substr($readID_with_at, 1);
			$readID =~ s/\/[12]$//;
			$total_readIDs{$readID}++;
			my $seq = <FASTQIN>; chomp($seq);
			my $plus = <FASTQIN>; die unless(substr($plus, 0, 1) eq '+');
			my $qual = <FASTQIN>; chomp($qual);
			die "Weird FASTQ format (II)" unless(length($qual) == length($seq));
			if($filter_IDs_href->{$readID})
			{
				print FASTQOUT $readID_with_at, "\n", $seq, "\n", $plus, $qual, "\n";
				$got_IDs{$readID}++;
			}
		}
		close(FASTQIN);
		close(FASTQOUT);
	}
	
	my $irregularities = 0;
	foreach my $readID (keys %$filter_IDs_href)
	{
		unless((exists $got_IDs{$readID} ) and (($got_IDs{$readID} == 1) or ($got_IDs{$readID} == 2)))
		{
			$irregularities++;
			print "Missing read ID: $readID\n";
			if($failOnError)
			{
				die "Fatal error: Missing read ID: $readID\n";
			}
		}
	}
	if($irregularities)
	{
		warn "filterReadIDs(..): observed $irregularities irregularities when filtering for " . scalar(keys %$filter_IDs_href) . " read IDs.";
	}
	
	return scalar(keys %total_readIDs);
}

sub filterReadIDs_nonTarget
{
	my $inputFiles_aref = shift;
	my $filter_IDs_href = shift;
	my $addSuffix = shift;
	my $samplingRate = shift;
	
	die unless(defined $samplingRate);
	die unless($samplingRate >= 0);
	die unless($samplingRate <= 1);
	
	my %included_readIDs;
	my %excluded_readIDs;
	foreach my $f (@$inputFiles_aref)
	{
		my $fn_out = $f . $addSuffix;
		open(FASTQIN, '<', $f) or die "Cannot open $f";
		open(FASTQOUT, '>', $fn_out) or die "Cannot open $fn_out";
		while(<FASTQIN>)
		{
			my $line = $_;
			chomp($line);
			next unless($line);
			my $readID_with_at = $line;
			die "Weird FASTQ format" unless(substr($readID_with_at, 0, 1) eq '@');
			my $readID = substr($readID_with_at, 1);
			$readID =~ s/\/[12]$//;
			my $seq = <FASTQIN>; chomp($seq);
			my $plus = <FASTQIN>; die unless(substr($plus, 0, 1) eq '+');
			my $qual = <FASTQIN>; chomp($qual);
			die "Weird FASTQ format (II)" unless(length($qual) == length($seq));
			
			if($filter_IDs_href->{$readID})
			{
				next;
			}
			
			if($excluded_readIDs{$readID})
			{
				next;
			}

			my $randomFraction =  rand(1);
			die unless(($randomFraction >= 0) and ($randomFraction <= 1));
			
			if($included_readIDs{$readID} or ($randomFraction <= $samplingRate))
			{
				$included_readIDs{$readID}++;
				print FASTQOUT $readID_with_at, "\n", $seq, "\n", $plus, $qual, "\n";
			}
			else
			{
				$excluded_readIDs{$readID}++;
			}
		}
		close(FASTQIN);
		close(FASTQOUT);
	}
	
	my $n_total_reads = scalar(keys %included_readIDs) + scalar(keys %excluded_readIDs);
	my $achieved_sampling_rate = sprintf("%.2f", scalar(keys %included_readIDs)/$n_total_reads);
	print "Info: achieved sampling rate (of non-targeted reads) $achieved_sampling_rate, target " .sprintf("%.2f", $samplingRate) . "\n";
}

sub extractRelevantReadsFromBAM
{
	my $known_references_dir = shift;
	my $BAM = shift;
	my $working_dir_thisSample = shift;
	my $extractInto_FASTQ1 = shift;
	my $extractInto_FASTQ2 = shift;
	my $extractInto_FASTQU = shift;
	my $longReads = shift;
	my $mapAgainstCompleteGenome_forRet_ref = shift;
	
	die unless (-d $working_dir_thisSample);
	
	my $command_idxstats = qq($samtools_bin idxstats $BAM);
	my $command_idxstats_output = `$command_idxstats`;
	unless($command_idxstats_output)
	{
		die "Didn't get a sensible output from command $command_idxstats";
	}
	my %BAM_idx;
	my @BAM_idx_contigOrder;
	my @idxstats_lines = split(/\n/, $command_idxstats_output);
	foreach my $l (@idxstats_lines)
	{
		my @l_fields = split(/\t/, $l);
		die unless(scalar(@l_fields) >= 1);
		die if(exists $BAM_idx{$l_fields[0]});
		$BAM_idx{$l_fields[0]} = $l_fields[1];
		push(@BAM_idx_contigOrder, $l_fields[0]);
	}

	my @files_references = glob($known_references_dir . '/*.txt');
	die "No known reference files in knownReferences ($full_graph_dir)?" unless(@files_references);
	my $additional_references_dir = $this_bin_dir . '/additionalReferences/' . $graph;
	if(-e $additional_references_dir)
	{
		my @additional_files_references = glob($additional_references_dir . '/*.txt');
		push(@files_references, @additional_files_references);
	}
	if($moreReferencesDir and (-e $moreReferencesDir))
	{
		my $mD = $moreReferencesDir . '/'  . $graph;
		if(-d $mD)
		{
			my @even_more_files_references = glob($mD . '/*.txt');
			push(@files_references, @even_more_files_references);			
		}
	}

	my @compatible_files;
	my %extractContigs_complete_byFile;
	my %extractContigs_partial_byFile;
	foreach my $f (@files_references)
	{
		open(F, '<', $f) or die "Cannot open $f";
		my $F_firstLine = <F>;
		chomp($F_firstLine);
		my @firstLine_fields = split(/\t/, $F_firstLine);
		my @expected_firstLine_fields = qw/contigID contigLength ExtractCompleteContig PartialExtraction_Start PartialExtraction_Stop/;
		die "Incorrect header for $f ($#firstLine_fields vs $#expected_firstLine_fields)" unless($#firstLine_fields == $#expected_firstLine_fields);
		for(my $i = 0; $i <= $#firstLine_fields; $i++)
		{
			die "Incorrect header for $f" unless($firstLine_fields[$i] eq $expected_firstLine_fields[$i]);
		}
		
		my $n_contigs = 0;
		my $n_contigs_matching = 0;
		while(<F>)
		{
			my $line = $_;
			chomp($line);
			next unless($line);
			my @line_fields = split(/\t/, $line, -1);
			$n_contigs++;
			if((exists $BAM_idx{$line_fields[0]}) and ($BAM_idx{$line_fields[0]} == $line_fields[1]))
			{
				$n_contigs_matching++;
			}
			else
			{
				# print "Mismatch $line_fields[0] - $line_fields[1] - " . $BAM_idx{$line_fields[0]} . "\n";
			}
			
			die if(($line_fields[2]) and (($line_fields[3]) and ($line_fields[4])));
			
			if($line_fields[2])
			{
				if($line_fields[0] eq '*')
				{
					$extractContigs_complete_byFile{$f}{$line_fields[0]} = 1;
				}
				else
				{
					$extractContigs_partial_byFile{$f}{$line_fields[0]} = [1, $line_fields[1]];						
				}
			}
			
			if(($line_fields[3]) and ($line_fields[4]))
			{	
				die if($line_fields[0] eq '*');
				die "Coordinate field $line_fields[3] has non-numeric characters" unless($line_fields[3] =~ /^\d+$/);
				die "Coordinate field $line_fields[4] has non-numeric characters" unless($line_fields[4] =~ /^\d+$/);
				
				$extractContigs_partial_byFile{$f}{$line_fields[0]} = [$line_fields[3], $line_fields[4]];			
			}
		}	
		close(F);

		if(($n_contigs_matching == $n_contigs) and ($n_contigs == scalar(@BAM_idx_contigOrder)))
		{
			push(@compatible_files, $f);
		}
	}

	if(scalar(@compatible_files) == 0)
	{
		die "Have found no compatible reference specifications in $known_references_dir - create a file for this BAM file ($BAM) and try again.";
	}
	if(scalar(@compatible_files) > 1)
	{
		die "Found more than one compatible reference file in $known_references_dir - a duplicate?\n\n".Dumper(\@compatible_files);
	}

	my $compatible_reference_file = $compatible_files[0];

	my @refIDs_for_extraction;
	foreach my $refID (@BAM_idx_contigOrder)
	{
		next if($refID eq '*');
		die if((exists $extractContigs_complete_byFile{$compatible_reference_file}{$refID}) and (exists $extractContigs_partial_byFile{$compatible_reference_file}{$refID}));
		die unless((not defined $extractContigs_complete_byFile{$compatible_reference_file}{$refID}) or ($extractContigs_complete_byFile{$compatible_reference_file}{$refID} eq '0') or ($extractContigs_complete_byFile{$compatible_reference_file}{$refID} eq '1'));
		if($extractContigs_complete_byFile{$compatible_reference_file}{$refID})
		{
			push(@refIDs_for_extraction, $refID);
		}
		if($extractContigs_partial_byFile{$compatible_reference_file}{$refID})
		{
			push(@refIDs_for_extraction, $refID . ':' . $extractContigs_partial_byFile{$compatible_reference_file}{$refID}[0] . '-' . $extractContigs_partial_byFile{$compatible_reference_file}{$refID}[1]);
		}
	}

	die "No contigs for extraction specified in $compatible_reference_file?" unless(scalar(@refIDs_for_extraction));

	my $target_extraction = $working_dir_thisSample . '/extraction.bam';

	my $threads_minus_1 = $maxThreads - 1;
	die unless($threads_minus_1 >= 0);

	my $target_extraction_mapped = $working_dir_thisSample . '/extraction_mapped.bam';
	if($samtools_T)
	{
		die "File $samtools_T specified via --samtools_T not existing" unless(-e $samtools_T);
	}
	my $view_T_switch = ($samtools_T) ? " -T $samtools_T " : "";
	my $extraction_command = qq($samtools_bin view -\@ $threads_minus_1 $view_T_switch -bo $target_extraction_mapped $BAM ).join(' ', @refIDs_for_extraction);
	print "Extract reads from ", scalar(@refIDs_for_extraction), " regions...\n";
	if(system($extraction_command) != 0)
	{
		die "Extraction command $extraction_command failed";
	}

	if($extractContigs_complete_byFile{$compatible_reference_file}{'*'})
	{
		my $target_extraction_unmapped = $working_dir_thisSample . '/extraction_unmapped.bam';
		
		my $extraction_command_unmapped = qq($samtools_bin view -\@ $threads_minus_1 $view_T_switch $BAM '*' | awk '{if (\$3 == "*") print \$0}' | $samtools_bin view -bo $target_extraction_unmapped -);
		print "Extract unmapped reads...\n";
		
		if(system($extraction_command_unmapped) != 0)
		{
			my $extraction_count_unmapped = qq($samtools_bin view -\@ $threads_minus_1 $view_T_switch $BAM '*' | awk '{if (\$3 == "*") print \$0}' | wc);
			my $count_output = `$extraction_count_unmapped`;
			chomp($count_output);
			if($count_output =~ /^\s*0\s+/)
			{
				unlink($target_extraction) if (-e $target_extraction);		
				my $copy_command = qq(cp $target_extraction_mapped $target_extraction );
				if(system($copy_command) != 0)
				{
					die "Copy command $copy_command failed";
				}
			}
			else
			{
				die "Extraction command $extraction_command_unmapped for unmapped reads failed, and unmapped read count does not seem to be 0 ('$count_output')";
			}
		}
		else
		{
			unlink($target_extraction) if (-e $target_extraction);
			my $extraction_command_merge = qq($samtools_bin merge $target_extraction $target_extraction_mapped $target_extraction_unmapped);
			print "Merging...\n";		
			if(system($extraction_command_merge) != 0)
			{
				die "Merge command $extraction_command_merge failed";
			}
		}		
	}
	else
	{
		my $mv_command = qq(mv $target_extraction_mapped $target_extraction);
		if(system($mv_command) != 0)
		{
			die "Move command $mv_command failed";
		}	
	}

	my $index_command = qq($samtools_bin index $target_extraction);
	print "Indexing...\n";
	if(system($index_command) != 0)
	{
		die "Index command $index_command failed";
	}

	#my $target_FASTQ_1 = $working_dir_thisSample . '/R_1.fastq';
	#my $target_FASTQ_2 = $working_dir_thisSample . '/R_2.fastq';
	my $target_FASTQ_U_nonSplit = $extractInto_FASTQU . '.nonSplit';
	
	my $FASTQ_extraction_command = qq($FASTQ_extraction_command_part1 VALIDATION_STRINGENCY=LENIENT I=$target_extraction F=$extractInto_FASTQ1 F2=$extractInto_FASTQ2 FU=$target_FASTQ_U_nonSplit 2>&1);

	print "Extract FASTQ...\n\t$FASTQ_extraction_command\n";
	my $FASTQ_extraction_output = `$FASTQ_extraction_command`;
	#if(($FASTQ_extraction_output =~ /Exception/) or ($FASTQ_extraction_output !~ /net.sf.picard.sam.SamToFastq done/))
	if(($FASTQ_extraction_output !~ /picard.sam.SamToFastq done/))
	{
		die "Picard output: \n\n" . $FASTQ_extraction_output . "\n\nExtraction command $FASTQ_extraction_command $! \n\nAbort because the Picard FASTQ extraction process might have failed. I think so because I could not find the string 'picard.sam.SamToFastq done' in the Picard output.\n\n";
	}

	if($longReads)
	{
		die "You activated --longReads, but the two files $extractInto_FASTQ1 and $extractInto_FASTQ2 (which store paired-end reads) are not empty - this is weird, and I will abort." unless(((-s $target_FASTQ_1) == 0) && ((-s $target_FASTQ_2) == 0));
		unless((-s $target_FASTQ_U_nonSplit) > 0)
		{
			die "You activated --longReads, but my attempt to extract long unpaired reads from the specified input BAM failed. Abort.";
		}
		
		open(INLONG, '<', $target_FASTQ_U_nonSplit) or die;
		open(OUTLONG, '>', $extractInto_FASTQU) or die;
		while(<INLONG>)
		{
			chomp;
			my $readID = $_; 
			die unless(substr($readID, 0, 1) eq '@');
			my $sequence = <INLONG>; chomp($sequence);
			my $plus = <INLONG>;
			my $qualities = <INLONG>; chomp($qualities);
			die unless(substr($plus, 0, 1) eq '+');
			die unless(length($sequence) == length($qualities));
			
			
			if(length($sequence) < 50000)
			{
				print OUTLONG $readID, "\n", $sequence, "\n", "+\n", $qualities, "\n";
			}
			else
			{
				my $runningI = 0;
				while(length($sequence))
				{
					die unless(length($sequence) == length($qualities));
				
					my $thisPart_readID = '>rP' .  $runningI . substr($readID, 1);
					
					my $extractionLength = (length($sequence) > 50000) ? 50000 : length($sequence);
					print OUTLONG $thisPart_readID, "\n", substr($sequence, 0, $extractionLength), "\n", "+\n", substr($qualities, 0, $extractionLength), "\n";

					substr($sequence, 0, $extractionLength) = '';
					substr($qualities, 0, $extractionLength) = '';
					
					$runningI++;
				}
			}
		}
		close(INLONG);
		close(OUTLONG);

		
	}	
	else
	{
		die "You didn't activate --longReads, but the two files $extractInto_FASTQ1 and $extractInto_FASTQ2 (which store paired-end reads) are empty - this is weird, and I will abort." unless(((-s $extractInto_FASTQ1) > 0) && ((-s $extractInto_FASTQ2) > 0));	
		
		my $mv_cmd = qq(mv $target_FASTQ_U_nonSplit $extractInto_FASTQU);
		system($mv_cmd) and die "Could not execute move command (III): $mv_cmd";
	}

	#if(system($FASTQ_extraction_command) != 0)
	#{
	#}

	$$mapAgainstCompleteGenome_forRet_ref = ($extractContigs_complete_byFile{$compatible_reference_file}{'*'}) ? 1 : 0;
}

sub find_path
{
	my $id = shift;
	my $supplied_value = shift;
	my $forWhich = shift;
	
	if(defined $supplied_value)
	{
		die "Command-line supplied value/file for parameter $id not existing" unless(-e $supplied_value);
		return $supplied_value;
	}
	
	if(exists $paths_ini{$id})
	{
		foreach my $alternative (@{$paths_ini{$id}})
		{
			if(-e $alternative)
			{
				return $alternative;
			}
		}
	}	

	if($forWhich)
	{
		my $which_output = `which $forWhich`;
		$which_output =~ s/[\n\r]//g;
		if($which_output and (-e $which_output))
		{
			return $which_output;
		}
	}
	
	die "I couldn't figure out a path for ${id}. Order of precedence: check for parameter --${id}; check paths.ini in $this_bin_dir for a key named $id; parse, if command string defined, the output of the command 'which ${forWhich}' ('which ' means that the command string is not defined).";
}


sub readFASTA
{
	my $file = shift;	
	my $cut_sequence_ID_after_whitespace = shift;
	
	my %R;
	
	open(F, '<', $file) or die "Cannot open $file";
	my $currentSequence;
	while(<F>)
	{		
		my $line = $_;
		chomp($line);
		$line =~ s/[\n\r]//g;
		if(substr($line, 0, 1) eq '>')
		{
			if($cut_sequence_ID_after_whitespace)
			{
				$line =~ s/\s+.+//;
			}
			$currentSequence = substr($line, 1);
			$R{$currentSequence} = '';
		}
		else
		{
			die "Weird input in $file" unless (defined $currentSequence);
			$R{$currentSequence} .= $line;
		}
	}	
	close(F);
		
	return \%R;
}

sub filterReadsAndGenerateProjectedBAM
{
	my $target_FASTQ_1_thisFunction = shift;
	my $target_FASTQ_2_thisFunction = shift;
	my $target_FASTQ_U_thisFunction = shift;
	my $call2_hla_relevant_readIDs_primaryBAM_href = shift;
	my $outputBAM_thisFunction = shift;
	my $refGenome = shift;
	my $longReads_thisFunction = shift;
	my $mapAgainstCompleteGenome_thisFunction = shift;
	
	die "Long-read mode not supported yet" if($longReads_thisFunction);
	
	my $fn_call2_SAM = $outputBAM_thisFunction . '.unprojected.SAM';
	my $fn_call2_SAM_projected = $outputBAM_thisFunction . '.projected.SAM';

	my $total_readIDs_n = filterReadIDs([$target_FASTQ_1_thisFunction, $target_FASTQ_2_thisFunction, $target_FASTQ_U_thisFunction], $call2_hla_relevant_readIDs_primaryBAM_href, '.filtered_call2', 1);
	my $target_FASTQ_1_postFiltering = $target_FASTQ_1_thisFunction . '.filtered_call2';
	my $target_FASTQ_2_postFiltering = $target_FASTQ_2_thisFunction . '.filtered_call2';
	my $target_FASTQ_U_postFiltering = $target_FASTQ_U_thisFunction . '.filtered_call2';
		

	if(! $longReads_thisFunction)
	{
		die "This is weird" if($mapAgainstCompleteGenome_thisFunction); # we check this because this means that we don't extract unmapped reads
		my $total_nonExtracted_reads = $total_readIDs_n - scalar(keys %$call2_hla_relevant_readIDs_primaryBAM_href);
		if($total_nonExtracted_reads > 0)
		{
			my $targetReads_additionalExtraction;
			if($total_nonExtracted_reads >= (0.3 * scalar(keys %$call2_hla_relevant_readIDs_primaryBAM_href)))
			{
				$targetReads_additionalExtraction = (0.3 * scalar(keys %$call2_hla_relevant_readIDs_primaryBAM_href));
			}
			else
			{
				print STDERR "\n";
				warn "WARNING: There are only $total_nonExtracted_reads non-HLA reads, but I would like to have " . (0.3 * scalar(keys %$call2_hla_relevant_readIDs_primaryBAM_href)) . " -- continue, but downstream mapping quality may suffer due to incorrectly estimated insert sizes.";
				print STDERR "\n";
				
				$targetReads_additionalExtraction = $total_nonExtracted_reads;
			}
			my $samplingRate = $targetReads_additionalExtraction / $total_nonExtracted_reads;
			die unless($samplingRate > 0);
			die unless($samplingRate <= 1);
			filterReadIDs_nonTarget([$target_FASTQ_1_thisFunction, $target_FASTQ_2_thisFunction, $target_FASTQ_U_thisFunction], $call2_hla_relevant_readIDs_primaryBAM_href, '.filtered_call2_additional', $samplingRate);
			
			my $cmd_cat_1 = qq(cat ${target_FASTQ_1_thisFunction}.filtered_call2_additional >> $target_FASTQ_1_postFiltering);
			my $cmd_cat_2 = qq(cat ${target_FASTQ_2_thisFunction}.filtered_call2_additional >> $target_FASTQ_2_postFiltering);
			my $cmd_cat_U = qq(cat ${target_FASTQ_U_thisFunction}.filtered_call2_additional >> $target_FASTQ_U_postFiltering);
			system($cmd_cat_1) and die "Cat coommand (I) failed: $cmd_cat_1";
			system($cmd_cat_2) and die "Cat coommand (II) failed: $cmd_cat_2";
			system($cmd_cat_U) and die "Cat coommand (III) failed: $cmd_cat_U";
			
			randomizeReadOrder([$target_FASTQ_1_postFiltering, $target_FASTQ_2_postFiltering]);
			randomizeReadOrder([$target_FASTQ_U_postFiltering]);
		}
		else
		{
			warn "Apparently the extracted FASTQs don't contain any non-HLA reads";
		}
	}
	
	my $cmd_bwa_map = qq($bwa_bin mem -t $maxThreads $refGenome $target_FASTQ_1_postFiltering $target_FASTQ_2_postFiltering > $fn_call2_SAM);
	if(system($cmd_bwa_map))
	{
		warn "Could not execute: $cmd_bwa_map, try again...";
		system("$bwa_bin mem $refGenome $target_FASTQ_1_postFiltering $target_FASTQ_2_postFiltering > $fn_call2_SAM") and die "Second attempt also failed: $cmd_bwa_map";
	}	
	print "\nGenerated SAM file: $fn_call2_SAM\n\n";
	
	my $cmd_projection = qq(perl Perl/projectSAM.pl --graph PRG_MHC_GRCh38_withIMGT --inputSAM $fn_call2_SAM --reference $refGenome --outputSAM $fn_call2_SAM_projected --samtools_bin $samtools_bin);
	system($cmd_projection) and die "Projection command $cmd_projection failed";
		
	my $initialBAM = $fn_call2_SAM_projected . '.bam';
	my $initialBAM_bai = $fn_call2_SAM_projected . '.bam.bai';
	die "Missing file $initialBAM" unless(-e $initialBAM);
	die "Missing file $initialBAM_bai" unless(-e $initialBAM_bai);
	
	my $mv_command_I = qq(mv $initialBAM $outputBAM_thisFunction);
	if(system($mv_command_I) != 0)
	{
		die "Move command $mv_command_I failed";
	}

	my $mv_command_II = qq(mv $initialBAM_bai ${outputBAM_thisFunction}.bai);
	if(system($mv_command_II) != 0)
	{
		die "Move command $mv_command_II failed";
	}	
	
	my $referenceGenome_for_GATK = $fn_call2_SAM_projected . '.ref.fa';
	die "Missing file: $referenceGenome_for_GATK" unless(-e $referenceGenome_for_GATK);
	
	return $referenceGenome_for_GATK;
}


sub extract_gene2readID_fromProjectedBAM
{
	my $BAM = shift;
	my $samtools_bin = shift;
	die unless(defined $samtools_bin);
	
	my $n_alignments_noGene = 0;	
	my $n_reads_unambiguous = 0;
	my $n_reads_multipleGenes = 0;
	
	my %read_2_gene;	
	my $view_cmd = qq($samtools_bin view $BAM);
	open(SAMTOOLS, "$view_cmd |") or die "Cannot open samtools view command $view_cmd";
	while(<SAMTOOLS>)
	{
		my $line = $_;
		chomp($line);
		my @line_fields = split(/\t/, $line);
		my $readID = $line_fields[0];
		my $referenceID = $line_fields[2];
		if($referenceID =~ /^(\w+)\*.+$/)
		{
			my $geneID = $1;
			$read_2_gene{$readID}{$geneID}++;
		}
		else
		{
			$n_alignments_noGene++;
		}
		
	}
	close(SAMTOOLS);
	
	my %forReturn;	
	foreach my $readID (keys %read_2_gene)
	{
		my @genes = keys %{$read_2_gene{$readID}};
		if(scalar(@genes) == 1)
		{
			$forReturn{$genes[0]}{$readID} = 1;
			$n_reads_unambiguous++;
		}
		else
		{
			$n_reads_unambiguous++;
		}
	}
	
	print "extract_gene2readID_fromProjectedBAM(..) summary:\n";
	print "\tBAM: $BAM\n";
	print "\t", "n_alignments_noGene", ": ", $n_alignments_noGene, "\n";
	print "\t", "n_reads_unambiguous", ": ", $n_reads_unambiguous, "\n";
	print "\t", "n_reads_multipleGenes", ": ", $n_reads_multipleGenes, "\n";
	
	return \%forReturn;
}
	
	
sub writeFASTA
{
	my $file = shift;
	# print "Writing $file\n";
	my $href = shift;
	open(F, '>', $file) or die "Cannot open $file";
	foreach my $key (sort keys %$href)
	{
		my $seq = $href->{$key};
		print F '>', $key, "\n";
		# print "\t", $key, "\t", length($seq), "\n";
		while($seq)
		{
			my $toPrint;
			if(length($seq) > 50)
			{
				$toPrint = substr($seq, 0, 50);
				substr($seq, 0, 50) = '';
			}
			else
			{
				$toPrint = $seq;
				$seq = '';
			}	
			print F $toPrint, "\n";
		}
	}
	close(F);	
}
