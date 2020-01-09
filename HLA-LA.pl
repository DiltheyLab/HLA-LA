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
use List::Util qw/all shuffle/;

$| = 1;
my $this_bin_dir = $FindBin::RealBin;

my $_BAM;
my $graph;
my $sampleID;
my $qsub;
my $samtools_bin;
my $bwa_bin;
my $java_bin;
my $picard_bin;
my $GATK_bin;
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
	'action:s' => \$action,
	'qsub:s' => \$qsub,
	'workingDir:s' => \$workingDir_param,
	'samtools_bin:s' => \$samtools_bin,
	'bwa_bin:s' => \$bwa_bin,
	'java_bin:s' => \$java_bin,
	'picard_bin:s' => \$picard_bin,
	'GATK_bin:s' => \$GATK_bin,
	'maxThreads:s' => \$maxThreads,
	'moreReferencesDir:s' => \$moreReferencesDir,
	'extractExonkMerCounts:s' => \$extractExonkMerCounts,
	'longReads:s' => \$longReads,
	'testing:s' => \$testing,
	'samtools_T:s' => \$samtools_T,
	'twoStageReadExtraction:s' => \$twoStageReadExtraction,
);

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
$GATK_bin = find_path('GATK', $GATK_bin, 'gatk');

my $FASTQ_extraction_command_part1;
if($picard_bin =~ /SamToFastq\.jar$/)
{
	die "Please use a recent Picard version and specify the path to picard.jar";
}
elsif($picard_bin =~ /picard-tools$/)
{
	die "Please use a recent Picard version and specify the path to picard.jar";
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
	print "\t", "GATK_bin", ": ", $GATK_bin, "\n";
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
my $working_dir_thisSample = $working_dir . '/' . $sampleID;

print "HLA-LA.pl\n\n";

print "Identified paths:\n";
print "\t", "samtools_bin", ": ", $samtools_bin, "\n";
print "\t", "bwa_bin", ": ", $bwa_bin, "\n";
print "\t", "java_bin", ": ", $java_bin, "\n";
print "\t", "picard_bin", ": ", $picard_bin, "\n";
print "\t", "GATK_bin", ": ", $GATK_bin, "\n";
print "\t", "General working directory", ": ", $working_dir, "\n";
print "\t", "Sample-specific working directory", ": ", $working_dir_thisSample, "\n";

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
	die "I need samtools >=1.3";
}

my $BAM;
my $call2_HLAtypes;
my %call2_hla_relevant_readIDs_primaryBAM;
my $call2_processing_dir;
my $call2_fn_mapping;
my $call2_fn_VCF;
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
elsif($action eq 'call2')
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

if($action eq 'call2')
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

			
		# my @files_readIDs = glob($hla_working_dir . '/*_readIDs_*');
		# foreach my $readIDs_file (@files_readIDs)
		# {
			# open(READS, '<', $readIDs_file) or die "Cannot open $readIDs_file";
			# while(<READS>)
			# {
				# my $line = $_;
				# chomp($line);
				# next unless($line);
				# $call2_hla_relevant_readIDs_primaryBAM{$line}++;
			# }
			# close(READS);
		# }
		
		print "Identified ", scalar(keys %call2_hla_relevant_readIDs_primaryBAM), " relevant read IDs.\n";
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
	close(REMAP);
	
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
		my $isUnmapped = (($FLAGS | 4) ? 1 : 0);
		if(! $isUnmapped)
		{
			$readIDs_filter{$FLAGS}++;
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
elsif($action eq 'call2')
{
	my $fn_call2_SAM = $call2_processing_dir . '/primaryReads_remapped.sam';
	my $fn_call2_SAM_projected = $call2_processing_dir . '/primaryReads_remapped.projected.sam';
	my $fn_call2_BAM_projected = $fn_call2_SAM_projected . '.bam';
	
	my $cmd_bwa_index = qq($bwa_bin index $call2_fn_mapping &> /dev/null);
	system($cmd_bwa_index) and die "Could not execute: $cmd_bwa_index";
	
	die "Long-read mode not supported yet" if($longReads);
	 
	my $total_readIDs_n = filterReadIDs([$target_FASTQ_1, $target_FASTQ_2, $target_FASTQ_U], \%call2_hla_relevant_readIDs_primaryBAM, '.filtered_call2', 1);
	my $target_FASTQ_1_postFiltering = $target_FASTQ_1 . '.filtered_call2';
	my $target_FASTQ_2_postFiltering = $target_FASTQ_2 . '.filtered_call2';
	my $target_FASTQ_U_postFiltering = $target_FASTQ_U . '.filtered_call2';
	
	if(! $longReads)
	{
		die "This is weird" if($mapAgainstCompleteGenome); # we check this because this means that we don't extract unmapped reads
		my $total_nonExtracted_reads = $total_readIDs_n - scalar(keys \%call2_hla_relevant_readIDs_primaryBAM);
		if($total_nonExtracted_reads > 0)
		{
			my $targetReads_additionalExtraction;
			if($total_nonExtracted_reads >= (0.3 * scalar(keys %call2_hla_relevant_readIDs_primaryBAM)))
			{
				$targetReads_additionalExtraction = (0.3 * scalar(keys %call2_hla_relevant_readIDs_primaryBAM));
			}
			else
			{
				print STDERR "\n";
				warn "WARNING: There are only $total_nonExtracted_reads non-HLA reads, but I would like to have " . (0.3 * scalar(keys %call2_hla_relevant_readIDs_primaryBAM)) . " -- continue, but downstream mapping quality may suffer due to incorrectly estimated insert sizes.";
				print STDERR "\n";
				
				$targetReads_additionalExtraction = $total_nonExtracted_reads;
			}
			my $samplingRate = $targetReads_additionalExtraction / $total_nonExtracted_reads;
			die unless($samplingRate > 0);
			die unless($samplingRate <= 1);
			filterReadIDs_nonTarget([$target_FASTQ_1, $target_FASTQ_2, $target_FASTQ_U], \%call2_hla_relevant_readIDs_primaryBAM, '.filtered_call2_additional', $samplingRate);
			
			my $cmd_cat_1 = qq(cat ${target_FASTQ_1}.filtered_call2_additional >> $target_FASTQ_1_postFiltering);
			my $cmd_cat_2 = qq(cat ${target_FASTQ_2}.filtered_call2_additional >> $target_FASTQ_2_postFiltering);
			my $cmd_cat_U = qq(cat ${target_FASTQ_U}.filtered_call2_additional >> $target_FASTQ_U_postFiltering);
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
	
	my $cmd_bwa_map = qq($bwa_bin mem -t $maxThreads $call2_fn_mapping $target_FASTQ_1_postFiltering $target_FASTQ_2_postFiltering > $fn_call2_SAM);
	system($cmd_bwa_map) and die "Could not execute: $cmd_bwa_map";
	
	print "\nGenerated SAM file: $fn_call2_SAM\n\n";
	
	my $cmd_projection = qq(perl Perl/projectSAM.pl --graph PRG_MHC_GRCh38_withIMGT --inputSAM $fn_call2_SAM --outputSAM $fn_call2_SAM_projected);
	system($cmd_projection) and die "Projection command $cmd_projection failed";
	die unless(-e $fn_call2_BAM_projected);
	
	# GATK
	{
		my $referenceGenome_for_GATK = $fn_call2_SAM_projected . '.ref.fa';
		die "Missing file: $referenceGenome_for_GATK" unless(-e $referenceGenome_for_GATK);
		
		my $BAM = $fn_call2_BAM_projected;
		my $BAM_RG = $BAM . '.rg.bam';		
		my $dict = $referenceGenome_for_GATK;
		$dict =~ s/\.fa$/.dict/;
		if(-e $dict)
		{
			unlink($dict) or die "Cannot unlink $dict";
		}
		
		$call2_fn_VCF = $call2_processing_dir . '/HLA.VCF';
 
		my $cmd_Picard_1 = qq(java -jar $picard_bin CreateSequenceDictionary R=$referenceGenome_for_GATK O=$dict);
		print "Now executing: $cmd_Picard_1 \n";
		system($cmd_Picard_1) and die "Command $cmd_Picard_1 failed";
		
		my $cmd_samtools_index = qq($samtools_bin faidx $referenceGenome_for_GATK);
		print "Now executing: $cmd_samtools_index \n";
		system($cmd_samtools_index) and die "Command $cmd_samtools_index failed";
				
		my $cmd_Picard_2 = qq(module load Picard; java -jar $picard_bin AddOrReplaceReadGroups I=$BAM O=$BAM_RG  RGID=4  RGLB=lib1  RGPL=illumina  RGPU=unit1  RGSM=20; $samtools_bin index $BAM_RG);
		print "Now executing: $cmd_Picard_2 \n";
		system($cmd_Picard_2) and die "Command $cmd_Picard_2 failed";
		
		my $cmd_GATK = qq($GATK_bin --java-options "-Xmx4g" HaplotypeCaller -R $referenceGenome_for_GATK -I $BAM_RG -O $call2_fn_VCF);
		print "Now executing: $cmd_GATK \n";
		system($cmd_GATK) and die "Command $cmd_GATK failed";	
	}
	
	print "Generated VCF: $call2_fn_VCF\n";
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
			die "Extraction command $extraction_command_unmapped failed";
		}
		
		unlink($target_extraction) if (-e $target_extraction);
		my $extraction_command_merge = qq($samtools_bin merge $target_extraction $target_extraction_mapped $target_extraction_unmapped);
		print "Merging...\n";		
		if(system($extraction_command_merge) != 0)
		{
			die "Merge command $extraction_command_merge failed";
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