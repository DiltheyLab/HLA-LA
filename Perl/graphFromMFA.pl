#!/usr/bin/perl

use strict;
use Data::Dumper;
use Getopt::Long;
use File::Basename;
use Cwd 'abs_path';

my $MFA;
my $graphName;
my $linearSequenceName;
GetOptions ('MFA:s' => \$MFA,
'graphName:s' => \$graphName,
'linearSequenceName:s' => \$linearSequenceName,
);

unless($MFA)
{
	die "Please specify input MFA in FASTA format with --MFA";
}

unless($graphName)
{
	die "Please specify output graph name --graphName";
}

unless(-e $MFA)
{
	die "Specified MFA $MFA not existing!";
}

my $baseDir_graphs_MHCPRG1 = '/Net/birch/data/dilthey/MHC-PRG/tmp2/GS_nextGen/';
unless(-e $baseDir_graphs_MHCPRG1)
{
	die "MHC-PRG-1 graph dir not existing - $baseDir_graphs_MHCPRG1.";
}
my $baseDir_graphs_MHCPRG2 = '/gpfs1/well/gsk_hla/MHC-PRG-2/graphs/';
my $baseDir_graphs_MHCPRG2_linearReferences = '/gpfs1/well/gsk_hla/MHC-PRG-2/linearReferenceALTs/';
unless(-e $baseDir_graphs_MHCPRG2)
{
	die "MHC-PRG-2 graph dir not existing - $baseDir_graphs_MHCPRG2.";
}
unless(-e $baseDir_graphs_MHCPRG2_linearReferences)
{
	die "baseDir_graphs_MHCPRG2_linearReferences graph dir not existing - $baseDir_graphs_MHCPRG2_linearReferences.";
}


my $seq_href = readFASTA($MFA);

my $MFA_basename = basename($MFA);

unless(scalar(keys %$seq_href))
{
	die "No sequences in $MFA?";
}
my $L;
my %id_output_2_id_input;
my %S_for_output;
foreach my $seqID (keys %$seq_href)
{
	my $S = $seq_href->{$seqID};
	unless(defined $L)
	{
		$L = length($S);
	}
	die "Length problem with $seqID - not equal to $L, which is assumed consensus length" unless($L == length($S));
	
	die "Weird characteters in $seqID" unless($S =~ /^[ACGTN\-]+$/i);
	$S = uc($S);
	$S =~ s/-/_/g;
	my $id_for_output = $seqID;
	$id_for_output =~ s/\s/_/g;	
	$S_for_output{$id_for_output} = $S;
	
	$id_output_2_id_input{$id_for_output} = $seqID;
}

my $outputDir_MHC_PRG_1 = $baseDir_graphs_MHCPRG1 . '/' . $graphName;
mkdir($outputDir_MHC_PRG_1);
die unless(-e $outputDir_MHC_PRG_1);

my @positionIDs;
my $f_positions = $outputDir_MHC_PRG_1 . '/positions.txt';
open(POSITIONS, '>', $f_positions) or die "Cannot open $f_positions";
for(my $pI = 0; $pI < $L; $pI++)
{
	my $pName = 'P_' . $pI . '_' . ($pI+1);
	print POSITIONS join(' ', $pName, ($pI+1)), "\n";
	push(@positionIDs, $pName);
}
close(POSITIONS);

my $f_segments = $outputDir_MHC_PRG_1 . '/segments.txt';
my $fn_alignment = 'PRG_' . $MFA_basename;

open(SEGMENTS, '>', $f_segments) or die "Can't open $f_segments";
print SEGMENTS $fn_alignment, "\n";
close(SEGMENTS);

my $f_alignment = $outputDir_MHC_PRG_1 . '/' . $fn_alignment;
open(ALIGNMENTOUT, '>', $f_alignment) or die "Cannot open $f_alignment";
print ALIGNMENTOUT join(' ', 'IndividualID', @positionIDs), "\n";
foreach my $key (keys %S_for_output)
{
	my @f = split(//, $S_for_output{$key});
	die unless($#f == $#positionIDs);
	print ALIGNMENTOUT join(' ', $key, @f), "\n";
}
close(ALIGNMENTOUT);

my $f_genomicMapping = $outputDir_MHC_PRG_1 . '/genomicMapping.txt';
open(F, '>', $f_genomicMapping) or die "Cannot open $f_genomicMapping";
print F "\n";
close(F);


print "Generated input for graph construction with MHC-PRG-1 $outputDir_MHC_PRG_1\n";

my $cmd_suggestion = qq(cd ${baseDir_graphs_MHCPRG1}../../src/ \n../bin/MHC-PRG domode createConcatenatedVariationGraphs ${baseDir_graphs_MHCPRG1}/${graphName} --noPGFprotection);

print "\nSuggested commands\n";
print $cmd_suggestion, "\n\n";

my $outputDir_MHC_PRG_2;
if($linearSequenceName)
{
	$outputDir_MHC_PRG_2 = $baseDir_graphs_MHCPRG2_linearReferences . '/' . $linearSequenceName . '/genePRG/';
	unless(-e $outputDir_MHC_PRG_2)
	{
		die "Directory expected but not existing: $outputDir_MHC_PRG_2"
	}
}
else
{
	$outputDir_MHC_PRG_2 = $baseDir_graphs_MHCPRG2 . '/' . $graphName;
}
mkdir($outputDir_MHC_PRG_2);
die "Directory $outputDir_MHC_PRG_2 not existing" unless(-e $outputDir_MHC_PRG_2);

my $MHCPRG2_PRG_dir = $outputDir_MHC_PRG_2 . '/PRG';
mkdir($MHCPRG2_PRG_dir);
writeToFile($MHCPRG2_PRG_dir . '/ReadMe.txt', "Please execute graph building with MHC-PRG-1 and copy here");

my $MHCPRG2_PRG_refGenomeDir = $outputDir_MHC_PRG_2 . '/referenceGenome';
mkdir($MHCPRG2_PRG_refGenomeDir);
die unless(-e $MHCPRG2_PRG_refGenomeDir);


my $MHCPRG2_mappingDir = $outputDir_MHC_PRG_2 . '/mapping';
mkdir($MHCPRG2_mappingDir);
die unless(-e $MHCPRG2_mappingDir);

my $MHCPRG2_translationDir = $outputDir_MHC_PRG_2 . '/translation';
mkdir($MHCPRG2_translationDir);
die unless(-e $MHCPRG2_translationDir);

my $MHCPRG2_alignedTestSequences = $outputDir_MHC_PRG_2 . '/alignedTestSequences';
mkdir($MHCPRG2_alignedTestSequences);
die unless(-e $MHCPRG2_alignedTestSequences);

my $MHCPRG2_sequences = $outputDir_MHC_PRG_2 . '/sequences.txt';
open(SEQUENCES, '>', $MHCPRG2_sequences) or die "Cannot open $MHCPRG2_sequences";
print SEQUENCES join("\t", qw/SequenceID Name FASTAID Chr Start_1based Stop_1based/), "\n";
my %refGenome;

my %S_alignedTest;
my $PRG_counter = 0;
foreach my $seqID (keys %S_for_output)
{
	$PRG_counter++;
	my $S = $S_for_output{$seqID};
	die "Weird characteters in $seqID" unless($S =~ /^[ACGTN_]+$/i);
	$S = uc($S);
	$S =~ s/_//g;
	die unless($S =~ /^[ACGTN]+$/);
	
	my $PRG_id = 'PRG_' . $PRG_counter;
	$refGenome{$PRG_id} = $S;
	
	print SEQUENCES join("\t", $PRG_counter, $seqID, $PRG_id, '', '', ''), "\n";
	
	my %oneFA = ($PRG_id => $S);
	my $fn_oneFA = $MHCPRG2_mappingDir . '/' . $PRG_counter . '.fa';
	writeFASTA($fn_oneFA, \%oneFA);
	
	my $id_alignment = $id_output_2_id_input{$seqID};
	die unless($id_alignment);
	my $printed_translations = 0;
	my $PRG_position = -1;
	my $alignedSequence = $seq_href->{$id_alignment};
	die unless($alignedSequence);
	die unless($alignedSequence =~ /^[ACGTN\-]+$/i);
	my $fn_translation = $MHCPRG2_translationDir . '/' . $PRG_counter . '.txt';
	open(TRANSLATION, '>', $fn_translation) or die "Cannot open $fn_translation";
	for(my $i = 0; $i < length($alignedSequence); $i++)
	{		
		my $c = substr($alignedSequence, $i, 1);
		if($c ne '-')
		{
			$printed_translations++;
			print TRANSLATION $i, "\n";
		}
	}
	close(TRANSLATION);
	die Dumper("N translation error", $printed_translations, length($S)) unless(length($S) == $printed_translations);
	
	my $alignedSequence_forTest = $alignedSequence;
	$alignedSequence_forTest = uc($alignedSequence_forTest);
	$alignedSequence_forTest =~ s/\-/_/g;
	die unless($alignedSequence_forTest =~ /^[ACGTN_]+$/);
	$S_alignedTest{$id_alignment} = $alignedSequence_forTest;
	
}
my $alignedTestSequences_fn = $MHCPRG2_alignedTestSequences . '/sequences.mfa';
writeFASTA($alignedTestSequences_fn, \%S_alignedTest);
close(SEQUENCES);

my $referenceGenome_fa = $MHCPRG2_PRG_refGenomeDir . '/ref.fa';
writeFASTA($referenceGenome_fa, \%refGenome);
close(SEQUENCES);

my $extendedReferenceGenome_absolutePath = abs_path($referenceGenome_fa);
writeToFile($outputDir_MHC_PRG_2 . '/extendedReferenceGenomePath.txt', $extendedReferenceGenome_absolutePath);

print "Output for MHC-PRG-2 went to $outputDir_MHC_PRG_2\n";
sub readFASTA
{
	my $file = shift;	
	my %R;
	
	open(F, '<', $file) or die "Cannot open $file";
	my $currentSequence;
	while(<F>)
	{
		if(($. % 1000000) == 0)
		{
		# 	print "\r", $.;
		}
		
		my $line = $_;
		chomp($line);
		$line =~ s/[\n\r]//g;
		if(substr($line, 0, 1) eq '>')
		{
			$currentSequence = substr($line, 1);
		}
		else
		{
			$R{$currentSequence} .= $line;
		}
	}	
	close(F);
		
	return \%R;
}

sub writeFASTA
{
	my $file = shift;
	# print "Writing $file\n";
	my $href = shift;
	open(F, '>', $file) or die "Cannot open $file";
	foreach my $key (keys %$href)
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
	
sub writeToFile
{
	my $f = shift;
	my $content = shift;
	open(F, '>', $f) or die "Cannot open $f";
	print F $content;
	close(F);
}