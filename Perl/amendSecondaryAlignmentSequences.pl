use strict;

#!/usr/bin/env perl

## Author: Alexander Dilthey (HHU/UKD, NHGRI-NIH), Evan Biederstedt (NYGC), Nathan Dunn (LBNL), Nancy Hansen (NIH), Aarti Jajoo (Baylor), Jeff Oliver (Arizona), Andrew Olsen (CSHL)
## License: The MIT License, https://github.com/NCBI-Hackathons/Graph_Genomes/blob/master/LICENSE

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;   
use List::Util qw/max all/;
use List::MoreUtils qw/mesh/;
use Bio::DB::HTS;

$| = 1;

## Usage:
## amendSecondaryAlignmentSequences.pl --input <path to SAM>
##             --output <path to SAM> 

my $input;
my $output;

GetOptions (
	'input:s' => \$input, 
	'output:s' => \$output, 
);

unless($input)
{
	die "Please specify --input (SAM file)";
}

unless($output)
{
	die "Please specify --output (SAM file)";
}

my $currentReadID;
my $currentReadSequence_positive;
my %processedReadIDs;
open(OUTPUT, '>', $output) or die "Cannot open $output";
open(INPUT, '<', $input) or die "Cannot open $input";
while(<INPUT>)
{
	my $line = $_;
	next unless($line);
	if(substr($line, 0, 1) eq '@')
	{ 
		print OUTPUT $line;
	}
	else
	{ 
		my @line_fields = split(/\t/, $line);
		die unless(scalar(@line_fields) >= 11);
		my $Name = $line_fields[0];
		my $FLAG = $line_fields[1];
		my $rName = $line_fields[2];
		my $POS = $line_fields[3];
		my $mapQ = $line_fields[4];
		my $CIGAR = $line_fields[5];
		my $SEQ = $line_fields[9];
		my $QUAL = $line_fields[10];

		if($SEQ eq '*')
		{
			die unless($currentReadID eq $Name);
			my $currentReadSequence_rightStrand;
			die "Read $Name line $. in $input, we apparently don't have a primary read sequence: $currentReadSequence_positive" unless(defined $currentReadSequence_positive);
			if(isReverseComplemented($FLAG))
			{
				$currentReadSequence_rightStrand = reverseComplement($currentReadSequence_positive);
			}
			else
			{
				$currentReadSequence_rightStrand = $currentReadSequence_positive;
			}
			
			die unless($CIGAR =~ /^((\d+)H)?(.+?)((\d+)H)?$/);
			my $hardclipping_left = $2;
			my $hardclipping_right = $5;
			my $middleComponent = $3;
			die if($middleComponent =~ /H/i);
			
			$hardclipping_left = 0 unless($hardclipping_left);
			$hardclipping_right = 0 unless($hardclipping_right);
			
			die "Weird hardclipping value: $hardclipping_left $CIGAR" unless($hardclipping_left =~ /^\d+$/);
			die unless($hardclipping_right =~ /^\d+$/);
			
			substr($currentReadSequence_rightStrand, 0, $hardclipping_left) = '';
			substr($currentReadSequence_rightStrand, length($currentReadSequence_rightStrand) - $hardclipping_right, $hardclipping_right) = '';
			
			die unless(length($currentReadSequence_rightStrand) == (length($currentReadSequence_positive) - $hardclipping_left - $hardclipping_right));
						
			$line_fields[9] = $currentReadSequence_rightStrand;
			
			print OUTPUT join("\t", @line_fields);
		}
		else
		{
			print OUTPUT $line;
		}
		
		if((not defined $currentReadID) or ($currentReadID ne $Name))
		{
			die "Reads are not ordered by read ID" if($processedReadIDs{$Name});
			$processedReadIDs{$Name}++;
			die unless(isPrimary($FLAG));
			die if($SEQ eq '*');
			die if($CIGAR =~ /H/);
			if(isReverseComplemented($FLAG))
			{
				die unless(defined $SEQ);
				$currentReadSequence_positive = reverseComplement($SEQ);
			}
			else
			{
				$currentReadSequence_positive = $SEQ;
			}
		}
		$currentReadID = $Name;
	}
	
}
close(INPUT);
close(OUTPUT);

print "\n\nDone. Generated $output.\n\n";


sub isReverseComplemented
{
	my $FLAG = shift;
	die unless(defined $FLAG);
	return ($FLAG & 16);
}

sub isNotPrimary
{
	my $FLAG = shift;
	die unless(defined $FLAG);
	return ($FLAG & 256);
}


sub isPrimary
{
	my $FLAG = shift;
	die unless(defined $FLAG);
	return (! isNotPrimary($FLAG));
}

sub reverseComplement
{
	my $kMer = shift;
	$kMer =~ tr/ACGT/TGCA/;
	return reverse($kMer);
	return $kMer;
}