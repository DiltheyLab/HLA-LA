#!/usr/bin/perl

# ./inferHLATypes.pl --BAM /gpfs1/well/gsk_hla/bam_output/AA02O9Q_Z2.bam --graph PRG_MHC_GRCh38_withIMGT --sampleID NA12878Direct

use warnings;
use strict;
use FindBin;
use File::Spec;
use Getopt::Long;
use Data::Dumper;
use Sys::Hostname;

$| = 1;
my $this_bin_dir = $FindBin::RealBin;

my $_BAM;
my $graph;
my $sampleID;
my $qsub;
my $samtools_bin;
my $bwa_bin;
my $java_bin;
my $picard_sam2fastq_bin;
my $maxThreads = 1;
GetOptions (
	'BAM:s' => \$_BAM,
	'graph:s' => \$graph,
	'sampleID:s' => \$sampleID,
	'qsub:s' => \$qsub,
	
	'samtools_bin:s' => \$samtools_bin,
	'bwa_bin:s' => \$bwa_bin,
	'java_bin:s' => \$java_bin,
	'picard_sam2fastq_bin:s' => \$picard_sam2fastq_bin,
	'maxThreads:s' => \$maxThreads,
);

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
		die unless($_ =~ /^(.+)=(.+)$/);
		my $id = $1;
		my @alts = split(/,/, $2);
		$paths_ini{$id} = \@alts;
	}
	close(INI);
}
$samtools_bin = find_path('samtools_bin', $samtools_bin, 'samtools');
$bwa_bin = find_path('bwa_bin', $bwa_bin, 'bwa');
$java_bin = find_path('java_bin', $java_bin, 'java');
$picard_sam2fastq_bin = find_path('picard_sam2fastq_bin', $picard_sam2fastq_bin, undef);

print "inferHLATypes.pl\n\n";

print "Identified paths:\n";
print "\t", "samtools_bin", ": ", $samtools_bin, "\n";
print "\t", "bwa_bin", ": ", $bwa_bin, "\n";
print "\t", "java_bin", ": ", $java_bin, "\n";
print "\t", "picard_sam2fastq_bin", ": ", $picard_sam2fastq_bin, "\n";
print "\n";

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

unless($sampleID =~ /^\w+$/)
{
	die "Please use only alphanumeric characters - \\w+ - for --sampleID";
}

my $BAM = File::Spec->abs2rel($_BAM);
unless(-e $BAM)
{
	die "BAM (or CRAM) $BAM (inferred from input $_BAM) not existing" unless(-e $BAM);
	unless((-e $BAM . '.bai') or (-e $BAM . '.crai'))
	{
		die "File $BAM does not appear to be indexed";
	}
}

my $full_graph_dir = $FindBin::RealBin . '/../graphs/' . $graph;
my $known_references_dir = $full_graph_dir . '/knownReferences';
unless(-e $full_graph_dir)
{
	die "Graph directory $full_graph_dir not found - valid graph names are subdirectories of the graphs directory in the HLA-PRG-LA root";
}
unless((-e $full_graph_dir . '/sequences.txt') and ((-e $full_graph_dir . '/extendedReferenceGenomePath.txt') or (-e $full_graph_dir . '/extendedReferenceGenome/extendedReferenceGenome.fa')) and (-d $known_references_dir))
{
	die "Graph directory $full_graph_dir does not seem to be complete - does this directory specify a valid graph for HLA-PRG-LA?";
}

# get index for this BAM
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
	die "Incorrect header for $f" unless($#firstLine_fields == $#expected_firstLine_fields);
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
	die "Have found no compatible reference specifications in $known_references_dir - create a file for this BAM file and try again.";
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

my $working_dir = $this_bin_dir . '/../working/' . $sampleID;
unless(-d $working_dir)
{
	mkdir($working_dir) or die "Cannot mkdir $working_dir";
}

my $target_extraction = $working_dir . '/extraction.bam';


my $target_extraction_mapped = $working_dir . '/extraction_mapped.bam';
my $extraction_command = qq($samtools_bin view -bo $target_extraction_mapped $BAM ).join(' ', @refIDs_for_extraction);
print "Extract reads from ", scalar(@refIDs_for_extraction), " regions...\n";
if(system($extraction_command) != 0)
{
	die "Extraction command $extraction_command failed";
}

if($extractContigs_complete_byFile{$compatible_reference_file}{'*'})
{
	my $target_extraction_unmapped = $working_dir . '/extraction_unmapped.bam';
	
	my $extraction_command_unmapped = qq($samtools_bin view $BAM '*' | awk '{if (\$3 == "*") print \$0}' | $samtools_bin view -bo $target_extraction_unmapped -);
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

my $target_FASTQ_1 = $working_dir . '/R_1.fastq';
my $target_FASTQ_2 = $working_dir . '/R_2.fastq';
my $FASTQ_extraction_command = qq($java_bin -Xmx10g -XX:-UseGCOverheadLimit -jar $picard_sam2fastq_bin VALIDATION_STRINGENCY=LENIENT I=$target_extraction F=$target_FASTQ_1 F2=$target_FASTQ_2 2>&1);
print "Extract FASTQ...\n\t$FASTQ_extraction_command\n";
my $FASTQ_extraction_output = `$FASTQ_extraction_command`;
#if(($FASTQ_extraction_output =~ /Exception/) or ($FASTQ_extraction_output !~ /net.sf.picard.sam.SamToFastq done/))
if(($FASTQ_extraction_output !~ /net.sf.picard.sam.SamToFastq done/))
{
	die "Extraction command $FASTQ_extraction_command failed $! -- output\n\n" . $FASTQ_extraction_output;
}

#if(system($FASTQ_extraction_command) != 0)
#{
#}

my $mapAgainstCompleteGenome = ($extractContigs_complete_byFile{$compatible_reference_file}{'*'}) ? 1 : 0;
$mapAgainstCompleteGenome = 0;

my $host = hostname();
my $MHC_PRG_2_bin = (($host =~ /rescomp/) or ($host =~ /comp[ABC]/)) ? '../bin_cluster3/HLA-PRG-LA' : '../bin/HLA-PRG-LA';
die "Binary $MHC_PRG_2_bin not there!" unless(-e $MHC_PRG_2_bin);
my $command_MHC_PRG = qq($MHC_PRG_2_bin --action HLA --maxThreads $maxThreads --sampleID $sampleID --outputDirectory $working_dir --PRG_graph_dir $full_graph_dir --FASTQ1 $target_FASTQ_1 --FASTQ2 $target_FASTQ_2 --bwa_bin $bwa_bin --samtools_bin $samtools_bin --mapAgainstCompleteGenome $mapAgainstCompleteGenome);

print "\nNow executing:\n$command_MHC_PRG\n";

if(system($command_MHC_PRG) != 0)
{
	die "HLA-PRG-LA execution not successful. Command was $command_MHC_PRG\n";
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

