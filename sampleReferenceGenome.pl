#!/usr/bin/perl

use strict;
use warnings;
use List::MoreUtils qw/mesh/;
use Data::Dumper;

my $n_samples = 8;

my $graphDir = qq(/gpfs1/well/gsk_hla/MHC-PRG-2/graphs/PRG_MHC_GRCh38_withIMGT);
my $outputDir = $graphDir . '/sampledReferenceGenomes';
unless(-e $outputDir)
{
	mkdir($outputDir) or die "Cannot mkdir $outputDir";
}


open(F, '<', $graphDir . '/extendedReferenceGenomePath.txt') or die "Can' t open ".$graphDir . '/extendedReferenceGenomePath.txt';
my $referenceGenomePath = <F>; chomp($referenceGenomePath);
close(F);

my $referenceGenome_href = readFASTA($referenceGenomePath, 1);

my %PRGsequences;
open(SEQUENCES, '<', $graphDir . '/sequences.txt') or die "Can't open ".$graphDir . '/extendedReferenceGenomePath.txt';
my $headerLine_sequences = <SEQUENCES>;
chomp($headerLine_sequences);
my @headerFields_sequences = split(/\t/, $headerLine_sequences);
while(<SEQUENCES>)
{
	my $line = $_;
	chomp($line);
	next unless($line);
	
	my @fields = split(/\t/, $line);
	die unless($#fields == $#headerFields_sequences);
	
	my %line = (mesh @headerFields_sequences, @fields);
	
	my $refID;
	if($line{'Chr'})
	{
		$refID = $line{'Chr'};
	}
	else
	{
		$refID = $line{'FASTAID'};
	}
	
	unless(exists $referenceGenome_href->{$refID})
	{
		die Dumper("Can't find reference genome sequence", \%line);
	}
	
	$PRGsequences{$refID} = $referenceGenome_href->{$refID};
}
close(SEQUENCES);

my $fn_genomes = $graphDir . '/sampledReferenceGenomes.txt';
open(SAMPLEDGENOMES, '>', $fn_genomes) or die "Cannot open $fn_genomes";
my $n_random_samples = $n_samples - 1;
my $sampling_proportion = 2/$n_random_samples;
for(my $sampleI = 1; $sampleI <= $n_samples; $sampleI++)
{
	my %refGenome_for_output;
	my $fn_for_output = $outputDir . '/' . $sampleI . '.fa';
	print SAMPLEDGENOMES $fn_for_output, "\n";
	
	if($sampleI == 1)
	{
		%refGenome_for_output = %PRGsequences;
	}
	else
	{
		foreach my $sequenceID (keys %PRGsequences)
		{
			if(rand(1) <= $sampling_proportion)
			{
				$refGenome_for_output{$sequenceID} = $PRGsequences{$sequenceID};
			}
		}
	}
	
	print "Sample $sampleI, generate reference genome with ", scalar(keys %refGenome_for_output), " sequences.\n";
	writeFASTA($fn_for_output, \%refGenome_for_output);
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
			$R{$currentSequence} .= uc($line);
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